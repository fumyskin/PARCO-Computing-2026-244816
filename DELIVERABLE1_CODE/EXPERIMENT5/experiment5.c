
/******************************************************************************
 * How to build:
 *
 * GCC
 *      gcc mergebased_spmv.c -lm -O3 -fopenmp -o spmv
 *
 * Intel
 *      icc mergebased_spmv.c -openmp -O3 -lrt -fno-alias -xHost -lnuma -o spmv
 *      export KMP_AFFINITY=granularity=core,scatter
 *   
 ******************************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <stdint.h>

/* Include your sparse_matrix.h and utils.h headers here */
/* These would need to be converted to C as well */

//---------------------------------------------------------------------
// Globals, constants, and type declarations
//---------------------------------------------------------------------

static bool g_quiet = false;
static bool g_verbose = false;
static bool g_verbose2 = false;
static int g_omp_threads = -1;
static int g_expected_calls = 1000000;

//---------------------------------------------------------------------
// Utility types
//---------------------------------------------------------------------


// Check types !!! (is long ok? )
typedef struct {
    int x;
    int y;
} int2_t;

// CountIterator alternative in C (attempt)
typedef struct {
    long value;
} CountingIterator;

// Constructor
CountingIterator make_iter(long start) {
    CountingIterator it;
    it.value = start;
    return it;
}

// dereference 
long iter_get(const CountingIterator *it) {
    return it->value;
}

// increment 
void iter_inc(CountingIterator *it) {
    it->value++;
}

// add n 
CountingIterator iter_add(const CountingIterator *it, long n) {
    return make_iter(it->value + n);
}

// subtract n 
CountingIterator iter_sub(const CountingIterator *it, long n) {
    return make_iter(it->value - n);
}

// distance 
long iter_distance(const CountingIterator *a, const CountingIterator *b) {
    return a->value - b->value;
}

// comparisons 
int iter_equal(const CountingIterator *a, const CountingIterator *b) {
    return a->value == b->value;
}


//---------------------------------------------------------------------
// Helper macros
//---------------------------------------------------------------------

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

//---------------------------------------------------------------------
// MergePath Search
//---------------------------------------------------------------------

/*
Computes the begin offsets into A and B for the specific diagonal
*/
int2_t merge_path_search(
    int diagonal,
    const int* restrict a,
    int a_len,
    int b_len)
{
    // diagonal search range (in x coordinate space)
    int x_min = MAX(diagonal - b_len, 0);
    int x_max = MIN(diagonal, a_len);

    // 2D binary search along the diagonal search range
    while (x_min < x_max) {
        int x_pivot = (x_min + x_max) >> 1;
        /* b[diagonal - x_pivot - 1] is implicitly the counting iterator value */
        if (a[x_pivot] <= (diagonal - x_pivot - 1)) {
            // keep top right half of diagonal range
            x_min = x_pivot + 1;
        } else {
            // keep bottom left half of diagonal range
            x_max = x_pivot;
        }
    }

    return int2_t(
        x = MIN(x_min, a_len);  // x coordinate in A
        y = diagonal - x_min;   // y coordinate in B
    )
 
}

//---------------------------------------------------------------------
// SpMV verification
//---------------------------------------------------------------------
/**
 * Compute reference SpMV y = alpha*Ax + beta*y
 * For single precision
 */
void spmv_gold_float(
    int num_rows,
    const int* restrict row_offsets,
    const int* restrict column_indices,
    const float* restrict values,
    const float* restrict vector_x,
    const float* restrict vector_y_in,
    float* restrict vector_y_out,
    float alpha,
    float beta)
{
    for (int row = 0; row < num_rows; ++row) {
        float partial = beta * vector_y_in[row];
        for (int offset = row_offsets[row]; offset < row_offsets[row + 1]; ++offset) {
            partial += alpha * values[offset] * vector_x[column_indices[offset]];
        }
        vector_y_out[row] = partial;
    }
}

//---------------------------------------------------------------------
// CPU merge-based SpMV (single precision)
//---------------------------------------------------------------------

void omp_merge_csrmv(
    int num_threads,
    Sparse_CSR* a,
    int* restrict row_end_offsets,
    int* restrict column_indices,
    float* restrict values,
    float* restrict vector_x,
    float* restrict vector_y_out)
{
    // Temporary storage for inter-thread fix-up after load-balanced work
    int* row_end_offsets = a->row_ptr + 1;
    CountingIterator nonzero_indices = make_iter(0);  // initialize counter
    int num_merge_items = a->n_rows + a->values;
    int items_per_thread = (num_merge_items + num_threads - 1) / num_threads;
    static int row_carry_out[256];
    static float value_carry_out[256];

    #pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int tid = 0; tid < num_threads; tid++) {

        int start_diagonal = MIN(items_per_thread * tid, num_merge_items);
        int end_diagonal = MIN(start_diagonal + items_per_thread, num_merge_items);

        int2_t thread_coord = merge_path_search(start_diagonal, row_end_offsets, num_rows, num_nonzeros, &thread_coord);
        int2_t thread_coord_end = merge_path_search(end_diagonal, row_end_offsets, num_rows, num_nonzeros, &thread_coord_end);

        /* Consume whole rows */
        for (; thread_coord.x < thread_coord_end.x; ++thread_coord.x) {
            float running_total = 0.0f;
            for (; thread_coord.y < row_end_offsets[thread_coord.x]; ++thread_coord.y) {
                running_total += values[thread_coord.y] * vector_x[column_indices[thread_coord.y]];
            }
            vector_y_out[thread_coord.x] = running_total;
        }

        /* Consume partial portion of thread's last row */
        float running_total = 0.0f;
        for (; thread_coord.y < thread_coord_end.y; ++thread_coord.y) {
            running_total += values[thread_coord.y] * vector_x[column_indices[thread_coord.y]];
        }

        /* Save carry-outs */
        row_carry_out[tid] = thread_coord_end.x;
        value_carry_out[tid] = running_total;
    }

    /* Carry-out fix-up (rows spanning multiple threads) */
    for (int tid = 0; tid < num_threads - 1; ++tid) {
        if (row_carry_out[tid] < num_rows) {
            vector_y_out[row_carry_out[tid]] += value_carry_out[tid];
        }
    }
}

//---------------------------------------------------------------------
// Timing utilities
//---------------------------------------------------------------------

/* Conversion of timing structure needed */
typedef struct {
    struct timespec start;
    struct timespec stop;
} cpu_timer_t;

void timer_start(cpu_timer_t* timer) {
    clock_gettime(CLOCK_MONOTONIC, &timer->start);
}

void timer_stop(cpu_timer_t* timer) {
    clock_gettime(CLOCK_MONOTONIC, &timer->stop);
}

double timer_elapsed_ms(const cpu_timer_t* timer) {
    double start_ms = timer->start.tv_sec * 1000.0 + timer->start.tv_nsec / 1e6;
    double stop_ms = timer->stop.tv_sec * 1000.0 + timer->stop.tv_nsec / 1e6;
    return stop_ms - start_ms;
}

//---------------------------------------------------------------------
// Display utilities
//---------------------------------------------------------------------

void display_perf(double setup_ms, double avg_ms, int num_nonzeros, int num_rows, size_t value_size) {
    double nz_throughput, effective_bandwidth;
    size_t total_bytes = (num_nonzeros * (value_size * 2 + sizeof(int))) +
                         (num_rows) * (sizeof(int) + value_size);

    nz_throughput = (double)num_nonzeros / avg_ms / 1.0e6;
    effective_bandwidth = (double)total_bytes / avg_ms / 1.0e6;

    if (!g_quiet) {
        printf("fp%d: %.4f setup ms, %.4f avg ms, %.5f gflops, %.3lf effective GB/s\n",
               (int)(value_size * 8),
               setup_ms,
               avg_ms,
               2 * nz_throughput,
               effective_bandwidth);
    } else {
        printf("%.5f, %.5f, %.6f, %.3lf, ",
               setup_ms, avg_ms,
               2 * nz_throughput,
               effective_bandwidth);
    }
    fflush(stdout);
}

//---------------------------------------------------------------------
// Main (example usage - would need full integration with matrix loading)
//---------------------------------------------------------------------

int main(int argc, char** argv) {
    printf("Merge-based SpMV implementation in C\n");
    printf("Full integration requires sparse_matrix.h and utils.h to be converted as well\n");
    
    /* Example usage would go here */
    /* Parse command line arguments */
    /* Load matrix */
    /* Run tests */
    
    return 0;
}