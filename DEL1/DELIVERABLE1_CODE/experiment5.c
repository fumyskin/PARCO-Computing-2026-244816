#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <pthread.h>
#include <omp.h>
#include <math.h>
#include <limits.h>
#include "mmio.h"
#include "specifications.h"

//define fma block operation
#if defined(__x86_64__) && defined(__FMA__)
static inline double fma_fallback(double a, double b, double c) {
    __m128d A = _mm_set_sd(a);
    __m128d B = _mm_set_sd(b);
    __m128d C = _mm_set_sd(c);
    __m128d R = _mm_fmadd_sd(A, B, C);
    return _mm_cvtsd_f64(R);
}
#else
static inline double fma_fallback(double a, double b, double c) {
    return a * b + c; // may be fused by compiler with -O3 -ffp-contract=fast
}
#endif

// A small struct to return coordinates (row_index, nz_index)
typedef struct {
    unsigned row;   // x coordinate (number of row-end events consumed)
    unsigned nz;    // y coordinate (number of nonzero events consumed)
} MergeCoord;


// function to initialize a struct COO given the data extracted from .mtx file
Sparse_Coordinate* initialize_COO(
    unsigned n_rows,
    unsigned n_cols,
    unsigned nnz,
    unsigned* row_indices,
    unsigned* col_indices,
    double* values
)
{
    Sparse_Coordinate* struct_COO = surely_malloc(sizeof(Sparse_Coordinate));
    struct_COO->n_rows = n_rows;
    struct_COO->n_cols = n_cols;
    struct_COO->nnz = nnz;
    struct_COO->row_indices = row_indices;
    struct_COO->col_indices = col_indices;
    struct_COO->values = values;

    return struct_COO;
}

// function to perform Spmv on COO
void SpMV_COO(Sparse_Coordinate* COO, double* vec, double* res){
    for(unsigned i = 0; i < COO->n_rows; i++){
        res[i] = 0;
    }

    for(unsigned nnz_id = 0; nnz_id < COO->nnz; nnz_id++){
        unsigned i = COO->row_indices[nnz_id];
        unsigned j = COO->col_indices[nnz_id];
        double val = COO->values[nnz_id];

        res[i] += val * vec[j]; 
    } 

    return;
}


// MergePathSearch:
// diagonal : target diagonal index in the logical merge path (0..num_rows+nnz)
// a_len    : number of row-end entries (num_rows)
// b_len    : number of nonzeros (nnz)
// row_end_offsets : pointer to row end offsets for rows 0..(num_rows-1)
//                   (this is A.row_ptr + 1 in the paper: the end index for each row)
static inline MergeCoord MergePathSearch(
    unsigned diagonal,
    const unsigned *row_end_offsets,
    unsigned a_len,
    unsigned b_len)
{
    int x_min = (int)diagonal - (int)b_len;
    if (x_min < 0) x_min = 0;
    int x_max = (int)diagonal;
    if (x_max > (int)a_len) x_max = (int)a_len;

    while (x_min < x_max) {
        int pivot = (x_min + x_max) >> 1;
        int b_index = (int)diagonal - pivot - 1;

        int cond;
        if (b_index < 0) {
            cond = 0;
        } else {
            cond = (row_end_offsets[pivot] <= (unsigned)b_index);
        }

        if (cond) {
            x_min = pivot + 1;
        } else {
            x_max = pivot;
        }
    }

    MergeCoord c;
    c.row = (unsigned)x_min;
    c.nz  = diagonal - (unsigned)x_min;
    return c;
}



// Merge-path based parallel CSR SpMV
// m     : pointer to your Sparse_CSR
// x     : dense input vector (length m->n_cols)
// y     : output vector (length m->n_rows) - must be allocated by caller
// num_threads: if <=0 uses omp_get_max_threads()
void merge_path_csr_mv(const Sparse_CSR *m, const double *x, double *y, int num_threads)
{
    if (!m || !x || !y) return;

    unsigned num_rows = m->n_rows;
    unsigned nnz = m->row_ptr[num_rows]; // last entry is total number of nonzeros
    const unsigned *row_ptr = m->row_ptr;
    const unsigned *col_ind = m->col_ind;
    const double   *val     = m->values;

    if (num_threads <= 0) num_threads = omp_get_max_threads();

    // Prepare the "row_end_offsets" array: row_ptr + 1 in the paper (end offset for each row)
    // We can just use &row_ptr[1] as pointer, but MergePathSearch expects length = num_rows,
    // so pass &row_ptr[1] and a_len = num_rows. To simplify indexing we create a small pointer:
    const unsigned *row_end_offsets = &row_ptr[1]; // length num_rows; entries in [0..nnz]

    unsigned num_merge_items = num_rows + nnz;
    unsigned items_per_thread = (num_merge_items + (unsigned)num_threads - 1u) / (unsigned)num_threads;

    // allocate carry arrays (one per thread)
    unsigned *row_carry_out = (unsigned*) malloc((size_t)num_threads * sizeof(unsigned));
    double   *value_carry_out = (double*) malloc((size_t)num_threads * sizeof(double));

    if (!row_carry_out || !value_carry_out) {
        fprintf(stderr, "Allocation failure in merge_path_csr_mv\n");
        if (row_carry_out) free(row_carry_out);
        if (value_carry_out) free(value_carry_out);
        return;
    }

    // Initialize y to 0.0
    for (unsigned r = 0; r < num_rows; ++r) y[r] = 0.0;

    #pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int tid = 0; tid < num_threads; ++tid)
    {
        unsigned diagonal = items_per_thread * (unsigned)tid;
        if (diagonal > num_merge_items) diagonal = num_merge_items;
        unsigned diagonal_end = diagonal + items_per_thread;
        if (diagonal_end > num_merge_items) diagonal_end = num_merge_items;

        // Find starting coordinate and ending coordinate
        MergeCoord coord = MergePathSearch(diagonal, row_end_offsets, num_rows, nnz);
        MergeCoord coord_end = MergePathSearch(diagonal_end, row_end_offsets, num_rows, nnz);

        unsigned local_row = coord.row;
        unsigned local_nz  = coord.nz;

        double running_total = 0.0;

        unsigned items_to_consume = diagonal_end - diagonal;
        // Use Algorithm 4 style: single loop for items_per_thread iterations
        for (unsigned it = 0; it < items_to_consume; ++it) {
            // Check bounds before accessing row_end_offsets
            // If local_row < num_rows and local_nz < row_end_offsets[local_row] consume nonzero
            if (local_row < num_rows && local_nz < row_end_offsets[local_row]) {
                double a = val[local_nz];
                unsigned col = col_ind[local_nz];
                running_total = fma_fallback(a, x[col], running_total);
                ++local_nz;
            } else {
                // Move right: flush row result and reset accumulator, advance row index
                if (local_row < num_rows && running_total != 0.0) {
                    // store the running_total as the row's value
                    // IMPORTANT: Only the thread that completes a full row writes y[row].
                    // If the row is shared across threads, the last thread to reach
                    // the row-end will write the final y[row]. This matches the paper.
                    y[local_row] = running_total;
                }
                running_total = 0.0;
                ++local_row;
            }
        }

        // Save carry-out: the row index that would be next and the partial running_total
        row_carry_out[tid] = local_row;
        value_carry_out[tid] = running_total;
    } // end parallel region

    // Carry-out fix-up: only threads 0 .. num_threads-2 may have partial sums for rows
    // that were completed by a later thread. The last thread's carry is the global end.
    for (int t = 0; t < num_threads - 1; ++t) {
        unsigned r = row_carry_out[t];
        if (r < num_rows) {
            // Serial addition (this is outside parallel region)
            y[r] += value_carry_out[t];
        }
    }

    free(row_carry_out);
    free(value_carry_out);
}

unsigned coo_count(Sparse_Coordinate *p){
    if (p == NULL || p->nnz == 0)
        return 0;

    unsigned i, n = p->nnz; 
    if (n == 0) return 0;
    unsigned count = 1;
    for (i=1; i<n; i++){
        if (p->row_indices[i-1] !=p->row_indices[i] || 
            p->col_indices[i-1] !=p->col_indices[i]){
            count++;
        }
    }
    return count;
}

Sparse_CSR *coo_to_csr_matrix(Sparse_Coordinate *p) {
    Sparse_CSR *q;
    unsigned count, i;
    unsigned r,c, ri, ci, cols, k, l, rows;
    unsigned *col_ind, *row_ptr, *prow_ind, *pcol_ind;
    double x, *val, *pval;
    unsigned n = p->nnz;
    coo_quicksort(p, 0, n);
    k = coo_count(p);
    rows = p->n_rows;
    prow_ind=p->row_indices;
    pcol_ind=p->col_indices;
    pval = p->values;
    q = surely_malloc(sizeof(Sparse_CSR));
    val = surely_malloc(k * sizeof(double));
    col_ind = surely_malloc(k * sizeof(unsigned));
    row_ptr = surely_malloc((rows+1) * sizeof(unsigned));
    r=-1;
    c=0; 
    l=0;
    /* partial_csr_0 */
    for (i=0; i<n; i++) {
        ri = prow_ind[i];
        ci = pcol_ind[i];
        x = pval[i];
        if (ri==r){
            if (ci==c)
                val[l-1] += x; /* partial_csr_duplicate */
            else {
                c=ci;
                col_ind[l] = ci;
                val[l] = x;
                l++;           /* partial_csr_newcol */
            }
        }
        else{
            while (r+1<=ri) row_ptr[++r]=l; /* partial_csr_skiprow */
            c= ci;
            col_ind[l] = ci;
            val[l] = x;
            l++;            /* partial_csr_newrow */
        }
    }
    cols = p->n_cols;
    while (r+1<=rows) row_ptr[++r]=l;  /* partial_csr_lastrows */
    q->values = val;
    q->col_ind = col_ind;
    q->row_ptr = row_ptr;
    q->n_rows = rows;
    q->n_cols = cols;
    return q;          /* partial_CSR_properties */
}

/*
For SpMV, focus on memory/cache optimizations first (reordering, 
blocking, prefetching, improve locality, reduce indirection) 
— they yield larger gains. Then focus on optimizing computation (vectorization,
parallelization)
*/
void csr_mv_multiply(Sparse_CSR *m, double *v, double *p) {
    unsigned i, rows = m->n_rows;
    double *val = m->values;
    unsigned *col_ind = m->col_ind;
    unsigned *row_ptr = m->row_ptr;
    unsigned next=row_ptr[0];

    // sequential implementation
    #pragma omp parallel for schedule(static)
    for (i = 0; i < rows; i++) {
        double s = 0.0; // private scope to each thread
        for (unsigned h = row_ptr[i]; h < row_ptr[i + 1]; h++) {
            double x = val[h];
            unsigned j = col_ind[h];
            s = fma_fallback(x, v[j], s);
        }
        p[i] = s;
    }

}



static void init_random_vector(double *v, unsigned n) {
    for (unsigned i = 0; i < n; ++i) {
        v[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    }
}

int validate_results(const double *a, const double *b, unsigned n,
                     double tol, int print_errors)
{
    int errors = 0;
    for (unsigned i = 0; i < n; ++i) {
        double diff = fabs(a[i] - b[i]);
        if (diff > tol) {
            errors++;
            if (print_errors && errors < 20) {
                printf("Mismatch at row %u: seq=%g, merge=%g (diff=%g)\n",
                       i, a[i], b[i], diff);
            }
        }
    }
    return errors;
}

void test_csr_merge_path(Sparse_CSR *csr)
{
    unsigned n = csr->n_cols;
    unsigned m = csr->n_rows;

    // allocate vectors
    double *x = malloc(n * sizeof(double));
    double *y_seq = malloc(m * sizeof(double));
    double *y_merge = malloc(m * sizeof(double));

    srand(0);
    init_random_vector(x, n);

    // --- sequential ---
    double t0 = omp_get_wtime();
    csr_mv_multiply(csr, x, y_seq);
    double t1 = omp_get_wtime();

    // --- merge path ---
    double t2 = omp_get_wtime();
    merge_path_csr_mv(csr, x, y_merge, 0);  // 0 → use OMP default number of threads
    double t3 = omp_get_wtime();

    printf("\nSequential CSR time      = %.6f s\n", t1 - t0);
    printf("Merge-Path CSR time      = %.6f s\n", t3 - t2);

    // validate
    int errors = validate_results(y_seq, y_merge, m, 1e-12, 1);
    if (errors == 0) {
        printf("\n Results match. Merge-Path CSR is correct.\n");
    } else {
        printf("\n ERROR: %d mismatches found.\n", errors);
    }

    free(x);
    free(y_seq);
    free(y_merge);
}


int main(int argc, char *argv[])
{
    
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   // M=rows, N=cols, nz=nonzeroes
    int i;
    unsigned *I, *J;
    double *val;

    // Initialize struct for sparse matrix 
    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else    
    { 
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    //  This is how one can screen matrix types if their application 
    //  only supports a subset of the Matrix Market data types.     
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    // find out size of sparse matrix .... 
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    // reserve memory for matrices 
    I = (unsigned *) surely_malloc(nz * sizeof(unsigned));
    J = (unsigned *) surely_malloc(nz * sizeof(unsigned));
    val = (double *) surely_malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    for (i=0; i<nz; i++)
    {
        int temp_i, temp_j;
        fscanf(f, "%d %d %lg\n", &temp_i, &temp_j, &val[i]);
        I[i] = (unsigned)(temp_i - 1);
        J[i] = (unsigned)(temp_j - 1);
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    // for (i=0; i<nz; i++){
    //     fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);
    // }

    // create struct with data read from .mtx file
    Sparse_Coordinate* struct_COO = initialize_COO((unsigned)M, (unsigned)N, (unsigned)nz, I, J, val);

    // INITIALIZE MATRIX VECTOR MULTIPLICATION
    double* res = surely_malloc(M * sizeof(double));
    double* vec = surely_malloc(N * sizeof(double));

    // INITIALIZE RANDOM VECTOR
    srand(0);
    for(int i = 0; i < N; i++){
        vec[i] = rand() % 10;
    }

    // INITIALIZE CSR MATRIX FROM COO
    Sparse_CSR* struct_CSR = coo_to_csr_matrix(struct_COO);
    double* res_csr = surely_malloc(M * sizeof(double));

   
    // test_csr_merge_path(struct_CSR); -> misleading for times, avoid it unless you want to check correctness

    // COMPUTE SpMV WITH CSR MREGED
    int N_ITERS = 20;
    double total = 0;

    for(int i = 0; i < N_ITERS; i++) {
        double t0 = omp_get_wtime();
        //csr_mv_multiply(struct_CSR, vec, res_csr);
        merge_path_csr_mv(struct_CSR, vec, res_csr, 0);
        double t1 = omp_get_wtime();
        total += (t1 - t0);
    }

    double avg_time = total/N_ITERS;
    printf("\nElapsed time: %.6f seconds\n", avg_time);

    // // COMPUTE SpMV WITH COO (for verification)
    // SpMV_COO(struct_COO, vec, res);

    // printf("\nResult (first 10 entries):\n");
    // for (int i = 0; i < M && i < 10; i++) {
    //     printf("res[%d] = %g\n", i, res[i]);
    // }

    // printf("\nCSR Result (first 10 entries):\n");
    // for (int i = 0; i < M && i < 10; i++) {
    //     printf("res_csr[%d] = %g\n", i, res_csr[i]);
    // }   

    free(I);
    free(J);
    free(val);
    free(vec);
    free(res);
    free(res_csr);
    free(struct_CSR);      
    free(struct_COO);
    
    return 0;
}
