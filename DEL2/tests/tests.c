/*
 * test_spmv_correctness.c
 *
 * Correctness verification for parallel SpMV across 1D, 2D block, and 2D cyclic.
 *
 * Builds a known 12x12 sparse matrix on every rank, computes the serial
 * reference y_ref = A * x, then runs each distribution's SpMV pipeline
 * and gathers the distributed result back to rank 0 for comparison.
 *
 * Compile (example):
 *   mpicc -O2 -o test_spmv test_spmv_correctness.c specifications.c \
 *         distributions.c communication.c communication_2d.c mmio.c -lm
 *
 * Run:
 *   mpirun -np 4 ./test_spmv
 *
 * The test matrix is 12x12 with nonzeros spread across all quadrants so
 * that every distribution must actually communicate to get correct results.
 * If ghost detection or communication is broken, mismatches will appear.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "mmio.h"
#include "specifications.h"
#include "distributions.h"
#include "communication.h"

/* ========================================================================
 * 1. Build a known 12x12 test matrix in COO format
 *
 * Structure: tridiagonal + some far-off-diagonal entries to force
 * communication in every distribution scheme.
 *
 *   Row i has: A[i][i] = 4.0   (diagonal)
 *              A[i][i-1] = -1.0 (sub-diagonal, if i>0)
 *              A[i][i+1] = -1.0 (super-diagonal, if i<11)
 *              A[i][(i+6)%12] = 0.5  (far off-diagonal, wraps around)
 *
 * This guarantees that in 2D block with a 2x2 grid (blocks of 6 rows/cols),
 * process (0,0) needs column 6..11 values and vice versa.
 * ======================================================================== */
#define N 12

static void build_test_matrix(Sparse_Coordinate *coo)
{
    /* Count nnz first */
    unsigned nnz = 0;
    for (int i = 0; i < N; i++) {
        nnz++;                          /* diagonal */
        if (i > 0)     nnz++;          /* sub-diagonal */
        if (i < N - 1) nnz++;          /* super-diagonal */
        nnz++;                          /* far off-diagonal */
    }

    coo->n_rows = N;
    coo->n_cols = N;
    coo->nnz    = nnz;
    coo->row_indices = (unsigned *)surely_malloc(nnz * sizeof(unsigned));
    coo->col_indices = (unsigned *)surely_malloc(nnz * sizeof(unsigned));
    coo->values      = (double *)surely_malloc(nnz * sizeof(double));

    unsigned k = 0;
    for (int i = 0; i < N; i++) {
        /* sub-diagonal */
        if (i > 0) {
            coo->row_indices[k] = i;
            coo->col_indices[k] = i - 1;
            coo->values[k]      = -1.0;
            k++;
        }
        /* diagonal */
        coo->row_indices[k] = i;
        coo->col_indices[k] = i;
        coo->values[k]      = 4.0;
        k++;
        /* super-diagonal */
        if (i < N - 1) {
            coo->row_indices[k] = i;
            coo->col_indices[k] = i + 1;
            coo->values[k]      = -1.0;
            k++;
        }
        /* far off-diagonal: wraps around to force cross-block communication */
        coo->row_indices[k] = i;
        coo->col_indices[k] = (i + 6) % N;
        coo->values[k]      = 0.5;
        k++;
    }
}

/* ========================================================================
 * 2. Build the known x vector: x[i] = i + 1
 * ======================================================================== */
static void build_test_vector(double *x)
{
    for (int i = 0; i < N; i++) {
        x[i] = (double)(i + 1);
    }
}

/* ========================================================================
 * 3. Compute serial reference: y_ref = A * x  (using COO SpMV)
 * ======================================================================== */
static void compute_reference(Sparse_Coordinate *coo, double *x, double *y_ref)
{
    SpMV_COO(coo, x, y_ref);
}

/* ========================================================================
 * 4. Compare two vectors, print mismatches
 * ======================================================================== */
static int compare_results(const char *label, double *y_ref, double *y_test, int n)
{
    double tol = 1e-10;
    int mismatches = 0;
    int max_report = 20;

    for (int i = 0; i < n; i++) {
        double diff = fabs(y_ref[i] - y_test[i]);
        double scale = fmax(fabs(y_ref[i]), fabs(y_test[i]));
        double allowed = tol + tol * scale;

        if (diff > allowed) {
            if (mismatches < max_report) {
                printf("  [%s] MISMATCH row %2d: ref=%.12f  got=%.12f  diff=%.2e\n",
                       label, i, y_ref[i], y_test[i], diff);
            }
            mismatches++;
        }
    }

    if (mismatches == 0) {
        printf("  [%s] PASS — all %d entries match reference (tol=%.0e)\n", label, n, tol);
    } else {
        printf("  [%s] FAIL — %d/%d mismatches\n", label, mismatches, n);
    }
    return (mismatches == 0);
}


/* ========================================================================
 * 5. Test 1D distribution SpMV
 *
 * 1D uses cyclic row distribution: row i -> process (i % P)
 * Vector element j -> process (j % P), local index (j / P)
 * After SpMV, local_res[k] corresponds to global row (k * P + rank).
 * ======================================================================== */
static void test_1d(Sparse_Coordinate *coo, double *x_global, double *y_ref,
                    int rank, int comm_size)
{
    Sparse_Coordinate local_matrix = {0};
    Sparse_CSR local_csr = {0};

    /* distribute */
    distribution_1D(coo, &local_matrix, &local_csr, rank, comm_size);

    /* build local x vector (cyclic: x_local[k] = x_global[rank + k*P]) */
    unsigned local_vec_size = N / comm_size;
    if (rank < (int)(N % comm_size)) local_vec_size++;

    double *local_x = (double *)surely_malloc(local_vec_size * sizeof(double));
    for (unsigned k = 0; k < local_vec_size; k++) {
        unsigned global_idx = rank + k * comm_size;
        local_x[k] = x_global[global_idx];
    }

    /* set up ghost communication */
    Comm_Pattern cp;
    memset(&cp, 0, sizeof(cp));
    find_ghost_vals(&local_csr, &cp, rank, comm_size, N);

    /* exchange + SpMV */
    double *local_res = (double *)surely_malloc(local_csr.n_rows * sizeof(double));
    exchange_ghost_vals(local_x, &cp, rank, comm_size);
    local_SpMV(&local_csr, local_x, local_res, &cp, rank, comm_size);

    /* gather results back to rank 0 */
    /* local_res[k] -> global row (rank + k * comm_size) */
    double *y_gathered = NULL;
    if (rank == 0) y_gathered = (double *)surely_malloc(N * sizeof(double));

    /* Each rank sends (local_row_index, value) pairs; simpler: use an
       intermediate full-size vector and Reduce. */
    double *y_full = (double *)calloc(N, sizeof(double));
    for (unsigned k = 0; k < local_csr.n_rows; k++) {
        unsigned global_row = rank + k * comm_size;
        y_full[global_row] = local_res[k];
    }
    /* MPI_Reduce with SUM: only one process contributes to each row, rest are 0 */
    MPI_Reduce(y_full, y_gathered, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        compare_results("1D", y_ref, y_gathered, N);
        free(y_gathered);
    }

    free(y_full);
    free(local_x);
    free(local_res);
    free_sparse(&local_matrix);
    free_csr(&local_csr);
    /* note: should also free cp fields, but for a test this is fine */
}


/* ========================================================================
 * 6. Test 2D block distribution SpMV
 *
 * 2D block: process grid pr x pc. Process (r,c) owns rows
 *   [r*rpp, (r+1)*rpp) and its column block is [c*cpp, (c+1)*cpp).
 * Vector x is distributed by column block.
 * After fold, y_result contains the final y for local rows.
 * Local row k -> global row (my_row_start + k).
 * ======================================================================== */
static void test_2d_block(Sparse_Coordinate *coo, double *x_global, double *y_ref,
                          int rank, int comm_size)
{
    Sparse_Coordinate local_matrix = {0};
    Sparse_CSR local_csr = {0};

    /* distribute (keeps global column indices) */
    distribution_2D(coo, &local_matrix, &local_csr, rank, comm_size);

    /* set up 2D comm pattern */
    Comm_Pattern_2D cp2d;
    memset(&cp2d, 0, sizeof(cp2d));
    setup_comm_pattern_2d(&local_csr, &cp2d, rank, comm_size, N, N);

    /* build local x from global x (column block) */
    unsigned local_vec_size = cp2d.my_col_end - cp2d.my_col_start;
    double *local_x = (double *)surely_malloc(local_vec_size * sizeof(double));
    for (unsigned k = 0; k < local_vec_size; k++) {
        local_x[k] = x_global[cp2d.my_col_start + k];
    }

    /* allocate result buffers */
    double *y_partial = (double *)surely_malloc(local_csr.n_rows * sizeof(double));
    double *y_result  = (double *)surely_malloc(local_csr.n_rows * sizeof(double));

    /* run 2D SpMV: expand + local compute + fold */
    spmv_2d(&local_csr, local_x, y_partial, y_result, &cp2d, rank, comm_size);

    /* Debug: print ghost info per rank */
    printf("  [2D block] Rank %d: grid(%d,%d) rows[%u,%u) cols[%u,%u) "
           "local_nnz=%u ghosts=%d local_nrows=%u\n",
           rank, cp2d.my_row_coord, cp2d.my_col_coord,
           cp2d.my_row_start, cp2d.my_row_end,
           cp2d.my_col_start, cp2d.my_col_end,
           local_matrix.nnz, cp2d.expand_ghost_count, local_csr.n_rows);

    /* gather y_result back to rank 0 */
    /* y_result[k] -> global row (my_row_start + k) */
    /* But ONLY the "fold root" (rank 0 in row_comm, i.e., proc-col 0) has
       the correct summed result. Other proc-cols in the same row also have it
       because we used Allreduce. So ALL processes can contribute, but we must
       avoid double-counting: only proc-col 0 contributes to the gather. */
    double *y_full = (double *)calloc(N, sizeof(double));
    if (cp2d.my_col_coord == 0) {
        /* Only the first process in each proc-row contributes */
        for (unsigned k = 0; k < local_csr.n_rows; k++) {
            unsigned global_row = cp2d.my_row_start + k;
            if (global_row < N) {
                y_full[global_row] = y_result[k];
            }
        }
    }

    double *y_gathered = NULL;
    if (rank == 0) y_gathered = (double *)surely_malloc(N * sizeof(double));
    MPI_Reduce(y_full, y_gathered, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        compare_results("2D block", y_ref, y_gathered, N);
        free(y_gathered);
    }

    free(y_full);
    free(local_x);
    free(y_partial);
    free(y_result);
    free_sparse(&local_matrix);
    free_csr(&local_csr);
    free_comm_pattern_2d(&cp2d);
}


/* ========================================================================
 * 7. Test 2D cyclic distribution SpMV
 *
 * 2D cyclic: process grid pr x pc.
 * Row i -> proc-row (i % pr), local row (i / pr)
 * Col j -> proc-col (j % pc)
 *
 * Vector x is distributed by column block in the 2D comm pattern,
 * but for cyclic the "column block" concept is different. We use the
 * same setup_comm_pattern_2d which computes block-based column ranges.
 *
 * IMPORTANT: For cyclic, the column ownership is different from block.
 * The current setup_comm_pattern_2d uses block column ranges, but the
 * cyclic distribution assigns columns cyclically. This mismatch needs
 * to be addressed. For this test, we check what happens.
 *
 * Actually, looking at the code more carefully:
 * - distribution_2D_cyclic keeps GLOBAL column indices in CSR
 * - setup_comm_pattern_2d uses block-based column ranges to decide
 *   what's local vs ghost
 * - The vector is distributed by block column range
 *
 * For cyclic 2D: row i is owned by proc-row (i % pr), but the VECTOR
 * x[j] should be owned by proc-col (j % pc) for consistency. However,
 * setup_comm_pattern_2d assigns x[j] to proc-col (j / cols_per_proc).
 *
 * This means 2D_cyclic with the current comm pattern uses BLOCK vector
 * distribution even though the MATRIX rows/cols are cyclic. This might
 * still work because the SpMV just needs x[j] from somewhere, and the
 * expand phase will fetch it regardless of which process owns it —
 * as long as ownership is consistent between vector init and expand.
 * ======================================================================== */
static void test_2d_cyclic(Sparse_Coordinate *coo, double *x_global, double *y_ref,
                           int rank, int comm_size)
{
    Sparse_Coordinate local_matrix = {0};
    Sparse_CSR local_csr = {0};

    /* distribute (keeps global column indices) */
    distribution_2D_cyclic(coo, &local_matrix, &local_csr, rank, comm_size);

    /* set up 2D comm pattern (uses block column ranges for vector ownership) */
    Comm_Pattern_2D cp2d;
    memset(&cp2d, 0, sizeof(cp2d));
    setup_comm_pattern_2d(&local_csr, &cp2d, rank, comm_size, N, N);

    /* build local x from global x (block column distribution) */
    unsigned local_vec_size = cp2d.my_col_end - cp2d.my_col_start;
    double *local_x = (double *)surely_malloc(local_vec_size * sizeof(double));
    for (unsigned k = 0; k < local_vec_size; k++) {
        local_x[k] = x_global[cp2d.my_col_start + k];
    }

    double *y_partial = (double *)surely_malloc(local_csr.n_rows * sizeof(double));
    double *y_result  = (double *)surely_malloc(local_csr.n_rows * sizeof(double));

    spmv_2d(&local_csr, local_x, y_partial, y_result, &cp2d, rank, comm_size);

    /* Debug */
    printf("  [2D cyclic] Rank %d: grid(%d,%d) local_rows=%u local_nnz=%u ghosts=%d\n",
           rank, cp2d.my_row_coord, cp2d.my_col_coord,
           local_csr.n_rows, local_matrix.nnz, cp2d.expand_ghost_count);

    /* Gather y_result back to rank 0.
     * For cyclic: local row k -> global row (my_row_coord + k * proc_rows)
     * But setup_comm_pattern_2d used BLOCK row ranges (my_row_start + k).
     *
     * This is a mismatch! The comm pattern thinks rows are block-distributed,
     * but distribution_2D_cyclic distributes them cyclically.
     * For now, we use the BLOCK mapping since that's what the fold operates on.
     * This will likely produce wrong results — which is useful diagnostic info.
     *
     * Actually: setup_comm_pattern_2d sets my_row_start/end based on block ranges.
     * The fold reduces across proc-cols in the same proc-row. If the matrix rows
     * are cyclic but the fold assumes block rows, the result may be scrambled.
     * Let's try both mappings and see.
     */

    /* Try block mapping (what setup_comm_pattern_2d thinks): */
    double *y_full_block = (double *)calloc(N, sizeof(double));
    if (cp2d.my_col_coord == 0) {
        for (unsigned k = 0; k < local_csr.n_rows; k++) {
            unsigned global_row = cp2d.my_row_start + k;
            if (global_row < N) {
                y_full_block[global_row] = y_result[k];
            }
        }
    }

    /* Try cyclic mapping (what distribution_2D_cyclic actually does): */
    double *y_full_cyclic = (double *)calloc(N, sizeof(double));
    if (cp2d.my_col_coord == 0) {
        for (unsigned k = 0; k < local_csr.n_rows; k++) {
            unsigned global_row = cp2d.my_row_coord + k * cp2d.proc_rows;
            if (global_row < N) {
                y_full_cyclic[global_row] = y_result[k];
            }
        }
    }

    double *y_gathered_block = NULL;
    double *y_gathered_cyclic = NULL;
    if (rank == 0) {
        y_gathered_block  = (double *)surely_malloc(N * sizeof(double));
        y_gathered_cyclic = (double *)surely_malloc(N * sizeof(double));
    }

    MPI_Reduce(y_full_block,  y_gathered_block,  N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(y_full_cyclic, y_gathered_cyclic, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("  [2D cyclic — block row mapping]:\n");
        compare_results("2D_cyc(block_map)", y_ref, y_gathered_block, N);
        printf("  [2D cyclic — cyclic row mapping]:\n");
        compare_results("2D_cyc(cyclic_map)", y_ref, y_gathered_cyclic, N);
        free(y_gathered_block);
        free(y_gathered_cyclic);
    }

    free(y_full_block);
    free(y_full_cyclic);
    free(local_x);
    free(y_partial);
    free(y_result);
    free_sparse(&local_matrix);
    free_csr(&local_csr);
    free_comm_pattern_2d(&cp2d);
}


/* ========================================================================
 * MAIN
 * ======================================================================== */
int main(int argc, char *argv[])
{
    int rank, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (rank == 0) {
        printf("==========================================================\n");
        printf("  SpMV Correctness Test — %d x %d matrix, %d processes\n", N, N, comm_size);
        printf("==========================================================\n\n");
    }

    /* Build matrix and vector on ALL ranks (small enough to replicate) */
    Sparse_Coordinate coo = {0};
    build_test_matrix(&coo);

    double x_global[N];
    build_test_vector(x_global);

    double y_ref[N];
    compute_reference(&coo, x_global, y_ref);

    if (rank == 0) {
        printf("Reference y = A*x (serial):\n");
        for (int i = 0; i < N; i++) {
            printf("  y[%2d] = %12.6f\n", i, y_ref[i]);
        }
        printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* ---- Test 1D ---- */
    if (rank == 0) printf("--- Testing 1D distribution ---\n");
    MPI_Barrier(MPI_COMM_WORLD);
    test_1d(&coo, x_global, y_ref, rank, comm_size);
    MPI_Barrier(MPI_COMM_WORLD);

    /* ---- Test 2D block ---- */
    if (rank == 0) printf("\n--- Testing 2D block distribution ---\n");
    MPI_Barrier(MPI_COMM_WORLD);
    test_2d_block(&coo, x_global, y_ref, rank, comm_size);
    MPI_Barrier(MPI_COMM_WORLD);

    /* ---- Test 2D cyclic ---- */
    if (rank == 0) printf("\n--- Testing 2D cyclic distribution ---\n");
    MPI_Barrier(MPI_COMM_WORLD);
    test_2d_cyclic(&coo, x_global, y_ref, rank, comm_size);
    MPI_Barrier(MPI_COMM_WORLD);

    /* ---- Print summary ---- */
    if (rank == 0) {
        printf("\n==========================================================\n");
        printf("  Test complete. Check PASS/FAIL above for each method.\n");
        printf("==========================================================\n");
    }

    free_sparse(&coo);
    MPI_Finalize();
    return 0;
}