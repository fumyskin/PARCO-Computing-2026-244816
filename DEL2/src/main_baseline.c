#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "mmio.h"
#include "specifications.h"
#include "distributions.h"
#include "communication.h"       
#include "tests.h"

/* format better later
 NOTE: THE FOLLOWING MAIN IS A PLAYGROUND THAT CONTAINS MANY FUNCTIONALITIES :
 - 1D CYCLIC DISTRIBUTION
 - 2D BLOCK DISTRIBUTION
 - 2D CYCLIC DISTRIBUTION
 - HYBRID MPI IMPLEMENTATION

 All of them can be used and tested, HOWEVER, I have decided for this deliverable to focus ONLY on 
 - 1D (baseline)
 - 2D BLOCK DISTRIBUTION (bonus 1)
 - PARALLEL I/O READING ATTEMPT (see parallel_bonus.c) (bonus 2)
 - EXTRA METRICS (bonus 3)

 Hence, the deliverable 2 content will explore only such options (Professor is ok with that)
*/

int main(int argc, char *argv[])
{
    int rank, comm_size;
    int num_threads = 1;

    #ifdef HYBRID
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
        if (provided < MPI_THREAD_FUNNELED) {
            fprintf(stderr, "MPI does not support MPI_THREAD_FUNNELED\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        num_threads = omp_get_max_threads();
    #else
        MPI_Init(&argc, &argv);
    #endif

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    srand(time(NULL) + rank * 17);

    double start_total = MPI_Wtime();

    if (argc < 4) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <matrix_file|dummy> <repeats> <distribution_type> [rows_per_proc] [nnz_per_row]\n", argv[0]);
            fprintf(stderr, "  distribution_type: 1D, 2D, or 2D_cyclic\n");
            fprintf(stderr, "  For dummy matrices: provide rows_per_proc and nnz_per_row\n");
            fprintf(stderr, "Error: Expected at least 3 arguments, but received %d.\n", argc - 1);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    char *matrix = argv[1];
    int repeats  = atoi(argv[2]);
    char *distr_type = argv[3];
    int dummy = (strcmp(matrix, "dummy") == 0);

    int use_2d = 0;
    int use_cyclic = 0;
    if (strcmp(distr_type, "1D") == 0) {
        use_2d = 0;
    } else if (strcmp(distr_type, "2D") == 0) {
        use_2d = 1;
        use_cyclic = 0;
    } else if (strcmp(distr_type, "2D_cyclic") == 0) {
        use_2d = 1;
        use_cyclic = 1;
    } else {
        if (rank == 0) {
            fprintf(stderr, "Error: Invalid distribution type '%s'. Use 1D, 2D, or 2D_cyclic\n", distr_type);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    Sparse_Coordinate coo_matrix = {0};
    Sparse_CSR csr_matrix = {0};
    unsigned matrix_rows, matrix_cols, matrix_nnz;
    unsigned rows_per_process = 0;
    unsigned nnz_per_row = 0;

    if (dummy) {
        if (argc < 6) {
            if (rank == 0) {
                fprintf(stderr, "Error: weak scaling simulation requires #rows per process and #nnz per row\n");
            }
            MPI_Finalize();
            return EXIT_FAILURE;
        }
        rows_per_process = atoi(argv[4]);
        nnz_per_row = atoi(argv[5]);
    }

    // ========================================================================
    // LOAD / GENERATE MATRIX
    // ========================================================================
    if (!dummy) {
        if (rank == 0) {
            if (load_matrix_market(matrix, &coo_matrix) != 0) {
                fprintf(stderr, "Failed to load matrix\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            matrix_rows = coo_matrix.n_rows;
            matrix_cols = coo_matrix.n_cols;
            matrix_nnz  = coo_matrix.nnz;

            Sparse_CSR *temp_csr = coo_to_csr_matrix(&coo_matrix);
            csr_matrix = *temp_csr;
            free(temp_csr);

            printf("Real Matrix: Rows=%d | Cols=%d | NNZ=%d | Procs=%d | Distribution=%s\n",
                   matrix_rows, matrix_cols, matrix_nnz, comm_size, distr_type);
        }

        MPI_Bcast(&matrix_rows, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(&matrix_cols, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(&matrix_nnz,  1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        if (rank != 0) {
            coo_matrix.n_rows = matrix_rows;
            coo_matrix.n_cols = matrix_cols;
            coo_matrix.nnz    = matrix_nnz;
            coo_matrix.row_indices = (unsigned *)surely_malloc(matrix_nnz * sizeof(unsigned));
            coo_matrix.col_indices = (unsigned *)surely_malloc(matrix_nnz * sizeof(unsigned));
            coo_matrix.values      = (double *)surely_malloc(matrix_nnz * sizeof(double));
        }

        MPI_Bcast(coo_matrix.row_indices, matrix_nnz, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(coo_matrix.col_indices, matrix_nnz, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(coo_matrix.values,      matrix_nnz, MPI_DOUBLE,   0, MPI_COMM_WORLD);

    } else {
        if (use_2d) {
            int proc_rows, proc_cols;
            create_2d_grid(comm_size, &proc_rows, &proc_cols);
            matrix_rows = rows_per_process * proc_rows;
            matrix_cols = matrix_rows;
        } else {
            matrix_rows = rows_per_process * comm_size;
            matrix_cols = matrix_rows;
        }
        matrix_nnz = 0;

        if (rank == 0) {
            printf("Synthetic Matrix (Weak Scaling): Rows/proc=%d | Total rows=%d | NNZ/row=~%d | Procs=%d\n",
                   rows_per_process, matrix_rows, nnz_per_row, comm_size);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // ========================================================================
    // DISTRIBUTE MATRIX
    // ========================================================================
    Sparse_Coordinate local_matrix = {0};
    Sparse_CSR local_csr = {0};

    if (!dummy) {
        if (use_2d) {
            if (use_cyclic) {
                distribution_2D_cyclic(&coo_matrix, &local_matrix, &local_csr, rank, comm_size);
            } else {
                distribution_2D(&coo_matrix, &local_matrix, &local_csr, rank, comm_size);
            }
            if (rank == 0) {
                printf("Applied 2D%s distribution (global column indices preserved)\n",
                       use_cyclic ? " cyclic" : " block");
            }
        } else {
            distribution_1D(&coo_matrix, &local_matrix, &local_csr, rank, comm_size);
            if (rank == 0) {
                printf("Applied 1D distribution\n");
            }
        }
    } else {
        if (use_2d) {
            if (use_cyclic) {
                generate_dummy_matrix_2d_cyclic(&local_matrix, rows_per_process, nnz_per_row, rank, comm_size);
            } else {
                generate_dummy_matrix_2d_block(&local_matrix, rows_per_process, nnz_per_row, rank, comm_size);
            }
        } else {
            generate_dummy_matrix(&local_matrix, rows_per_process, nnz_per_row, rank, comm_size);
        }

        Sparse_CSR *temp_csr = coo_to_csr_matrix(&local_matrix);
        local_csr = *temp_csr;
        free(temp_csr);

        unsigned local_nnz_count = local_matrix.nnz;
        MPI_Allreduce(&local_nnz_count, &matrix_nnz, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

        if (!use_2d) {
            matrix_rows = rows_per_process * comm_size;
            matrix_cols = matrix_rows;
        }

        if (rank == 0) {
            printf("Dummy Matrix Generated: %u Global Rows | Total NNZ: %u\n", matrix_rows, matrix_nnz);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ========================================================================
    // SET UP COMMUNICATION PATTERN + VECTOR
    // ========================================================================

    // Variables used by both paths
    double *local_vec = NULL;
    double *local_res = NULL;
    unsigned local_vec_size = 0;

    // 1D-specific
    Comm_Pattern comm_pattern;
    memset(&comm_pattern, 0, sizeof(comm_pattern));

    // 2D-specific
    Comm_Pattern_2D cp2d;
    memset(&cp2d, 0, sizeof(cp2d));
    double *y_partial_2d = NULL;  // partial result before fold
    double *y_result_2d  = NULL;  // final result after fold

    if (use_2d) {
        // ====== 2D PATH ======

        // Set up 2D communication pattern
        setup_comm_pattern_2d(&local_csr, &cp2d, rank, comm_size, matrix_rows, matrix_cols);

        // Vector distribution: each process owns x entries for its COLUMN BLOCK
        // local_vec_size = my_col_end - my_col_start
        local_vec_size = cp2d.my_col_end - cp2d.my_col_start;

        local_vec = (double *)surely_malloc(local_vec_size * sizeof(double));
        for (unsigned i = 0; i < local_vec_size; i++) {
            unsigned global_col = cp2d.my_col_start + i;
            local_vec[i] = ((double)rand() / RAND_MAX) * 8.0 - 4.0;
        }

        // Allocate result buffers
        y_partial_2d = (double *)surely_malloc(local_csr.n_rows * sizeof(double));
        y_result_2d  = (double *)surely_malloc(local_csr.n_rows * sizeof(double));
        local_res = y_result_2d;  // point to the final result for metrics

        if (rank == 0) {
            printf("2D: local vector size = %u (column block), local rows = %u\n",
                   local_vec_size, local_csr.n_rows);
        }

        // Cache warmup
        spmv_2d(&local_csr, local_vec, y_partial_2d, y_result_2d, &cp2d, rank, comm_size);
        MPI_Barrier(MPI_COMM_WORLD);

    } else {
        // ====== 1D PATH ======

        find_ghost_vals(&local_csr, &comm_pattern, rank, comm_size, matrix_cols);

        if (rank == 0) {
            printf("Ghost columns identified: %d total ghost values needed\n",
                   comm_pattern.num_ghost_cols);
        }

        // Cyclic vector distribution
        local_vec_size = matrix_cols / comm_size;
        if (rank < (int)(matrix_cols % comm_size)) {
            local_vec_size++;
        }

        local_vec = (double *)surely_malloc(local_vec_size * sizeof(double));
        local_res = (double *)surely_malloc(local_csr.n_rows * sizeof(double));

        for (unsigned i = 0; i < local_vec_size; i++) {
            int global_idx = rank + i * comm_size;
            local_vec[i] = ((double)rand() / RAND_MAX) * 8.0 - 4.0;
        }

        if (rank == 0) {
            printf("Local vector initialized (size=%d per process)\n", local_vec_size);
        }

        // Cache warmup
        exchange_ghost_vals(local_vec, &comm_pattern, rank, comm_size);
        local_SpMV(&local_csr, local_vec, local_res, &comm_pattern, rank, comm_size);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // ========================================================================
    // PERFORMANCE BENCHMARK
    // ========================================================================
    Performance_Metrics pmetrics;
    pmetrics.local_nnz     = local_matrix.nnz;
    pmetrics.ghost_entries  = use_2d ? cp2d.expand_ghost_count : comm_pattern.num_ghost_cols;
    pmetrics.local_flops    = 2LL * local_matrix.nnz;
    pmetrics.num_repeats    = repeats;
    pmetrics.elapsed_times  = (double *)surely_malloc(repeats * sizeof(double));
    pmetrics.comm_times     = (double *)surely_malloc(repeats * sizeof(double));

    for (int r = 0; r < repeats; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        double total_iter_start = MPI_Wtime();

        if (use_2d) {
            // --- 2D: expand + compute + fold ---
            double comm_start = MPI_Wtime();
            expand_x_2d(local_vec, &cp2d, rank, comm_size);
            double comm_mid = MPI_Wtime();

            local_spmv_2d(&local_csr, local_vec, y_partial_2d, &cp2d, rank, comm_size);

            double fold_start = MPI_Wtime();
            fold_y_2d(y_partial_2d, y_result_2d, &cp2d);
            double fold_end = MPI_Wtime();

            // comm_time = expand + fold
            pmetrics.comm_times[r] = (comm_mid - comm_start) + (fold_end - fold_start);
        } else {
            // --- 1D: exchange + compute ---
            double comm_start = MPI_Wtime();
            exchange_ghost_vals(local_vec, &comm_pattern, rank, comm_size);
            double comm_end = MPI_Wtime();
            pmetrics.comm_times[r] = comm_end - comm_start;

            local_SpMV(&local_csr, local_vec, local_res, &comm_pattern, rank, comm_size);
        }

        double total_iter_end = MPI_Wtime();
        pmetrics.elapsed_times[r] = total_iter_end - total_iter_start;
    }

    // ========================================================================
    // GATHER STATISTICS 
    // ========================================================================
    unsigned local_nnz = local_matrix.nnz;
    unsigned local_ghosts = use_2d ? (unsigned)cp2d.expand_ghost_count : comm_pattern.num_ghost_cols;

    unsigned nnz_min, nnz_max;
    double nnz_avg;
    MPI_Reduce(&local_nnz, &nnz_min, 1, MPI_UNSIGNED, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_nnz, &nnz_max, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
    unsigned long long nnz_sum = local_nnz, nnz_total;
    MPI_Reduce(&nnz_sum, &nnz_total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    nnz_avg = (double)nnz_total / comm_size;

    unsigned ghost_min, ghost_max;
    double ghost_avg;
    MPI_Reduce(&local_ghosts, &ghost_min, 1, MPI_UNSIGNED, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_ghosts, &ghost_max, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
    unsigned long long ghost_sum = local_ghosts, ghost_total;
    MPI_Reduce(&ghost_sum, &ghost_total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    ghost_avg = (double)ghost_total / comm_size;

    double recv_bytes = local_ghosts * sizeof(double);
    double recv_bytes_avg;
    MPI_Reduce(&recv_bytes, &recv_bytes_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    recv_bytes_avg /= comm_size;

    // Compute average times
    double total_time_local = 0.0, comm_time_local = 0.0;
    for (int r = 0; r < repeats; r++) {
        total_time_local += pmetrics.elapsed_times[r];
        comm_time_local  += pmetrics.comm_times[r];
    }
    total_time_local /= repeats;
    comm_time_local  /= repeats;
    double comp_time_local = total_time_local - comm_time_local;

    double total_time_max, comm_time_max, comp_time_max;
    MPI_Reduce(&total_time_local, &total_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time_local,  &comm_time_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comp_time_local,  &comp_time_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    long long total_flops = 2LL * matrix_nnz;
    double gflops_total = 0.0, gflops_comp = 0.0;
    if (total_time_max > 0) gflops_total = (total_flops / 1e9) / total_time_max;
    if (comp_time_max > 0)  gflops_comp  = (total_flops / 1e9) / comp_time_max;

    // ========================================================================
    // DISPLAY METRICS
    // ========================================================================
    if (rank == 0) {
        printf("\n========================================\n");
        printf("COMPREHENSIVE METRICS SUMMARY\n");
        printf("========================================\n");
        printf("PROCESSES: %d\n", comm_size);
        printf("DISTRIBUTION: %s\n", distr_type);
        if (use_2d) {
            printf("GRID: %d x %d (proc_rows x proc_cols)\n", cp2d.proc_rows, cp2d.proc_cols);
        }
        printf("ITERATIONS: %d\n", repeats);
        printf("Matrix NNZ: %d\n", matrix_nnz);
        printf("BALANCE: nnz_local   min=%8u  avg=%10.2f  max=%8u\n", nnz_min, nnz_avg, nnz_max);
        printf("COMM:    ghost       min=%8u  avg=%10.2f  max=%8u  | recv_double_avg=%10.2f B\n",
               ghost_min, ghost_avg, ghost_max, recv_bytes_avg);
        printf("METRICS: t_total     = %.3e  t_comp = %.3e  t_comm = %.3e\n",
               total_time_max, comp_time_max, comm_time_max);
        printf("         gflops_total= %.6e  gflops_comp= %.6e\n", gflops_total, gflops_comp);
        printf("========================================\n\n");
    }

    // RAW OUTPUT
    if (rank == 0) {
        for (int r = 0; r < repeats; r++) {
            printf("[RESULT] %d,%d,%s,%d,%.9f,%.9f,%d,%d,%lld\n",
                   rank, comm_size, distr_type, r,
                   pmetrics.elapsed_times[r], pmetrics.comm_times[r],
                   pmetrics.local_nnz, pmetrics.ghost_entries,
                   (long long)pmetrics.local_flops);
            fflush(stdout);
        }
    }

    double end_total = MPI_Wtime();
    if (rank == 0) {
        printf("Total execution time: %.3f seconds\n", end_total - start_total);
    }

    // ========================================================================
    // CLEANUP
    // ========================================================================
    free(pmetrics.elapsed_times);
    free(pmetrics.comm_times);
    free(local_vec);

    if (use_2d) {
        free(y_partial_2d);
        // y_result_2d == local_res, only free once
        free(y_result_2d);
        free_comm_pattern_2d(&cp2d);
    } else {
        free(local_res);
        // free 1D comm_pattern fields (if you have a free function for it)
    }

    free_sparse(&coo_matrix);
    free_sparse(&local_matrix);
    free_csr(&local_csr);

    MPI_Finalize();
    return 0;
}