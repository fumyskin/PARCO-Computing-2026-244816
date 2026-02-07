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

int main(int argc, char *argv[])
{
    int rank, comm_size; 
    int num_threads = 1;

    #ifdef HYBRID
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
        if(provided < MPI_THREAD_FUNNELED){
            fprintf(stderr, "MPI does not support MPI_THREAD_FUNNELED\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        num_threads = omp_get_max_threads();
    #else
        MPI_Init(&argc, &argv);
    #endif

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    srand(time(NULL) + rank*17); 

    double start_total = MPI_Wtime();

    // Usage: ./main_parallel <matrix_file|dummy> <repeats> <io_type> [rows_per_proc] [nnz_per_row]
    // io_type: sequential or parallel
    if (argc < 4) {
        if (rank == 0){
            fprintf(stderr, "Usage: %s <matrix_file|dummy> <repeats> <io_type> [rows_per_proc] [nnz_per_row]\n", argv[0]);
            fprintf(stderr, "  io_type: sequential or parallel\n");
            fprintf(stderr, "  For dummy matrices: provide rows_per_proc and nnz_per_row\n");
            fprintf(stderr, "Error: Expected at least 3 arguments, but received %d.\n", argc - 1);
        }
        MPI_Finalize();
        return EXIT_FAILURE; 
    }

    char* matrix = argv[1];
    int repeats = atoi(argv[2]);
    char* io_type = argv[3];
    int dummy = (strcmp(matrix, "dummy") == 0);

    // Check I/O type
    int use_parallel_io = 0;
    if (strcmp(io_type, "sequential") == 0){
        use_parallel_io = 0;
    }else if (strcmp(io_type, "parallel") == 0){
        use_parallel_io = 1;
    }else{
        if (rank == 0){
            fprintf(stderr, "Error: Invalid I/O type '%s'. Use 'sequential' or 'parallel'\n", io_type);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    Sparse_Coordinate coo_matrix = {0};
    Sparse_CSR csr_matrix = {0};
    unsigned matrix_rows, matrix_cols, matrix_nnz;
    unsigned rows_per_process = 0;
    unsigned nnz_per_row = 0;

    // Input check for dummy matrices
    if (dummy){
        if(argc < 6){
            if(rank == 0){
                fprintf(stderr, "Error: dummy matrix requires #rows per process and #nnz per row\n");
            }
            MPI_Finalize();
            return EXIT_FAILURE;
        }

        rows_per_process = atoi(argv[4]);
        nnz_per_row = atoi(argv[5]);
    }

    // ============================================================
    // MATRIX LOADING SECTION
    // ============================================================
    double io_start = MPI_Wtime();
    
    Sparse_Coordinate local_matrix = {0};
    Sparse_CSR local_csr = {0};

    if(!dummy){
        if (use_parallel_io) {
            // ====== PARALLEL I/O PATH ======
            if (rank == 0) {
                printf("Using PARALLEL I/O (MPI_File_read_at_all)\n");
            }
            
            if (load_mm_chunked(matrix, &local_matrix, rank, comm_size) != 0) {
                if (rank == 0) fprintf(stderr, "Failed to load matrix with parallel I/O\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            
            // Extract global dimensions from local_matrix (they're broadcast inside load_mm_chunked)
            matrix_rows = local_matrix.n_rows;
            matrix_cols = local_matrix.n_cols;
            
            // Get total nnz by summing local nnz counts
            unsigned local_nnz_count = local_matrix.nnz;
            MPI_Allreduce(&local_nnz_count, &matrix_nnz, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            
            // Convert local COO to CSR
            Sparse_CSR *temp_csr = coo_to_csr_matrix(&local_matrix);
            local_csr = *temp_csr;
            free(temp_csr);
            
            if (rank == 0) {
                printf("Parallel I/O Complete: Rows=%u | Cols=%u | Total NNZ=%u | Procs=%d\n", 
                       matrix_rows, matrix_cols, matrix_nnz, comm_size);
            }
            
        } else {
            // ====== SEQUENTIAL I/O PATH ======
            if (rank == 0) {
                printf("Using SEQUENTIAL I/O (rank 0 reads, then broadcasts)\n");
                
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
               
                printf("Sequential I/O Complete: Rows=%u | Cols=%u | NNZ=%u | Procs=%d\n", 
                       matrix_rows, matrix_cols, matrix_nnz, comm_size);
            }
            
            // Broadcast dimensions
            MPI_Bcast(&matrix_rows, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(&matrix_cols, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(&matrix_nnz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            // Allocate memory on non-root ranks
            if (rank != 0) {
                coo_matrix.n_rows = matrix_rows;
                coo_matrix.n_cols = matrix_cols;
                coo_matrix.nnz = matrix_nnz;
            
                coo_matrix.row_indices = (unsigned*)surely_malloc(matrix_nnz * sizeof(unsigned));
                coo_matrix.col_indices = (unsigned*)surely_malloc(matrix_nnz * sizeof(unsigned));
                coo_matrix.values      = (double*)surely_malloc(matrix_nnz * sizeof(double));
            }
            
            // Broadcast the actual data
            MPI_Bcast(coo_matrix.row_indices, matrix_nnz, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(coo_matrix.col_indices, matrix_nnz, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(coo_matrix.values, matrix_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            
            // Apply 1D distribution
            distribution_1D(&coo_matrix, &local_matrix, &local_csr, rank, comm_size);
            
            if (rank == 0){
                printf("Applied 1D distribution\n");
            }
        }
        
    } else {
        // ====== DUMMY MATRIX PATH ======
        matrix_rows = rows_per_process * comm_size;
        matrix_cols = matrix_rows;
        matrix_nnz = 0;

        if (rank == 0) {
            printf("Generating Synthetic Matrix: Rows/proc=%u | Total rows=%u | NNZ/row=~%u | Procs=%d\n", 
                   rows_per_process, matrix_rows, nnz_per_row, comm_size);
        }
        
        // Generate random local matrix for each process
        generate_dummy_matrix(&local_matrix, rows_per_process, nnz_per_row, rank, comm_size);
        
        // Convert to CSR
        Sparse_CSR *temp_csr = coo_to_csr_matrix(&local_matrix);
        local_csr = *temp_csr;
        free(temp_csr); 

        // Calculate total nnz
        unsigned local_nnz_count = local_matrix.nnz;
        MPI_Allreduce(&local_nnz_count, &matrix_nnz, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0) {
            printf("Dummy Matrix Generated: Total Rows=%u | Total NNZ=%u\n", matrix_rows, matrix_nnz);
        }
    }
    
    double io_end = MPI_Wtime();
    double io_time = io_end - io_start;
    double io_time_max;
    MPI_Reduce(&io_time, &io_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("I/O Time (max across ranks): %.6f seconds\n", io_time_max);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ============================================================
    // GHOST ELEMENTS IDENTIFICATION
    // ============================================================
    Comm_Pattern comm_pattern;
    find_ghost_vals(&local_csr, &comm_pattern, rank, comm_size, matrix_cols);

    if(rank == 0){
        printf("Ghost columns identified: %d total ghost values needed\n", comm_pattern.num_ghost_cols);
    }

    // ============================================================
    // VECTOR INITIALIZATION (cyclic distribution)
    // ============================================================
    unsigned local_vec_size = matrix_cols / comm_size;
    if(rank < (matrix_cols % comm_size)){
        local_vec_size++;
    }

    double *local_rand_vec = (double*)surely_malloc(local_vec_size * sizeof(double));
    double *local_res = (double*)surely_malloc(local_csr.n_rows * sizeof(double));

    for (unsigned i = 0; i < local_vec_size; i++) {
        local_rand_vec[i] = ((double)rand() / RAND_MAX) * 8.0 - 4.0;
    }
    
    if (rank == 0) {
        printf("Local vector initialized (size=%u per process)\n", local_vec_size);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // ============================================================
    // CACHE WARMUP
    // ============================================================
    exchange_ghost_vals(local_rand_vec, &comm_pattern, rank, comm_size);
    local_SpMV(&local_csr, local_rand_vec, local_res, &comm_pattern, rank, comm_size);
    MPI_Barrier(MPI_COMM_WORLD);

    // ============================================================
    // PERFORMANCE BENCHMARK
    // ============================================================
    Performance_Metrics pmetrics;
    pmetrics.local_nnz = local_matrix.nnz;
    pmetrics.ghost_entries = comm_pattern.num_ghost_cols;
    pmetrics.local_flops = 2LL * local_matrix.nnz;
    pmetrics.num_repeats = repeats;

    pmetrics.elapsed_times = (double*)surely_malloc(repeats * sizeof(double));
    pmetrics.comm_times = (double*)surely_malloc(repeats * sizeof(double));

    for (int r = 0; r < repeats; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        double total_iter_start = MPI_Wtime();   
        
        double comm_start = MPI_Wtime();       
        exchange_ghost_vals(local_rand_vec, &comm_pattern, rank, comm_size);
        double comm_end = MPI_Wtime();             
        pmetrics.comm_times[r] = comm_end - comm_start;
        
        local_SpMV(&local_csr, local_rand_vec, local_res, &comm_pattern, rank, comm_size);
        
        double total_iter_end = MPI_Wtime();  
        pmetrics.elapsed_times[r] = total_iter_end - total_iter_start;
    }

    // ============================================================
    // METRICS COLLECTION
    // ============================================================
    unsigned nnz = local_csr.row_ptr[local_csr.n_rows];
    
    // Calculate local memory usage
    size_t local_bytes = 0;
    local_bytes += (local_csr.n_rows + 1) * sizeof(unsigned);  // row_ptr
    local_bytes += nnz * sizeof(unsigned);                      // col_indices
    local_bytes += nnz * sizeof(double);                        // values
    local_bytes += local_vec_size * sizeof(double);             // local_rand_vec
    local_bytes += local_csr.n_rows * sizeof(double);           // local_res
    local_bytes += comm_pattern.num_ghost_cols * sizeof(double); // ghost values
    
    unsigned local_nnz = local_matrix.nnz;
    unsigned local_ghosts = comm_pattern.num_ghost_cols;
    
    // Load balance metrics
    unsigned nnz_min, nnz_max;
    double nnz_avg;
    MPI_Reduce(&local_nnz, &nnz_min, 1, MPI_UNSIGNED, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_nnz, &nnz_max, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
    unsigned long long nnz_sum = local_nnz;
    unsigned long long nnz_total;
    MPI_Reduce(&nnz_sum, &nnz_total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    nnz_avg = (double)nnz_total / comm_size;
    
    // Communication volume metrics
    unsigned ghost_min, ghost_max;
    double ghost_avg;
    MPI_Reduce(&local_ghosts, &ghost_min, 1, MPI_UNSIGNED, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_ghosts, &ghost_max, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
    unsigned long long ghost_sum = local_ghosts;
    unsigned long long ghost_total;
    MPI_Reduce(&ghost_sum, &ghost_total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    ghost_avg = (double)ghost_total / comm_size;
    
    double recv_bytes = local_ghosts * sizeof(double);
    double recv_bytes_avg;
    MPI_Reduce(&recv_bytes, &recv_bytes_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    recv_bytes_avg /= comm_size;
    
    // Memory metrics
    size_t mem_min, mem_max;
    double mem_avg;
    MPI_Reduce(&local_bytes, &mem_min, 1, MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_bytes, &mem_max, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    unsigned long long mem_sum = local_bytes;
    unsigned long long mem_total;
    MPI_Reduce(&mem_sum, &mem_total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mem_avg = (double)mem_total / comm_size;
    
    // Ghost ratio
    double ghost_ratio = (local_nnz > 0) ? (100.0 * local_ghosts / local_nnz) : 0.0;
    double ghost_ratio_min, ghost_ratio_max, ghost_ratio_avg;
    MPI_Reduce(&ghost_ratio, &ghost_ratio_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ghost_ratio, &ghost_ratio_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ghost_ratio, &ghost_ratio_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    ghost_ratio_avg /= comm_size;
    
    // Timing metrics
    double total_time_local = 0.0, comm_time_local = 0.0;
    for (int r = 0; r < repeats; r++) {
        total_time_local += pmetrics.elapsed_times[r];
        comm_time_local += pmetrics.comm_times[r];
    }
    total_time_local /= repeats;
    comm_time_local /= repeats;
    
    double comp_time_local = total_time_local - comm_time_local;
    
    double total_time_max, comm_time_max, comp_time_max;
    MPI_Reduce(&total_time_local, &total_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time_local, &comm_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comp_time_local, &comp_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    // GFLOPS calculations
    long long total_flops = 2LL * matrix_nnz;
    double gflops_total = 0.0, gflops_comp = 0.0;
    if (total_time_max > 0) {
        gflops_total = (total_flops / 1e9) / total_time_max;
    }
    if (comp_time_max > 0) {
        gflops_comp = (total_flops / 1e9) / comp_time_max;
    }
    
    // ============================================================
    // DISPLAY METRICS 
    // ============================================================
    if (rank == 0) {
        printf("\n");
        printf("========================================\n");
        printf("COMPREHENSIVE METRICS SUMMARY\n");
        printf("========================================\n");
        printf("I/O TYPE:     %s\n", use_parallel_io ? "PARALLEL (MPI_File_read_at_all)" : "SEQUENTIAL (Rank 0 + Bcast)");
        printf("PROCESSES:    %d\n", comm_size);
        printf("DISTRIBUTION: 1D (row-cyclic)\n");
        printf("ITERATIONS:   %d\n", repeats);
        printf("Matrix NNZ:   %u\n", matrix_nnz);
        printf("I/O TIME:     %.6f seconds\n", io_time_max);
        printf("----------------------------------------\n");
        printf("BALANCE: nnz_local   min=%8u  avg=%10.2f  max=%8u\n", 
               nnz_min, nnz_avg, nnz_max);
        printf("COMM:    ghost       min=%8u  avg=%10.2f  max=%8u  | recv_bytes_avg=%10.2f B\n",
               ghost_min, ghost_avg, ghost_max, recv_bytes_avg);
        printf("MEM:     bytes       min=%8zu  avg=%10.0f  max=%8zu\n",
               mem_min, mem_avg, mem_max);
        printf("STRUCT:  ghost_ratio min=%8.3f%%  avg=%10.3f%%  max=%8.3f%%\n",
               ghost_ratio_min, ghost_ratio_avg, ghost_ratio_max);
        printf("----------------------------------------\n");
        printf("TIMING:  t_total     = %.6e s\n", total_time_max);
        printf("         t_comp      = %.6e s\n", comp_time_max);
        printf("         t_comm      = %.6e s\n", comm_time_max);
        printf("PERF:    gflops_total= %.6f GFLOPS\n", gflops_total);
        printf("         gflops_comp = %.6f GFLOPS\n", gflops_comp);
        printf("========================================\n\n");
    }

    // RAW OUTPUT DATA
    if (rank == 0){
        for (int r = 0; r < repeats; r++) {
            printf("[RESULT] %d,%d,%s,%d,%.9f,%.9f,%u,%u,%d,%.9f\n",
                rank, 
                comm_size,
                io_type,
                r,
                pmetrics.elapsed_times[r], 
                pmetrics.comm_times[r],
                pmetrics.local_nnz,
                pmetrics.ghost_entries,
                pmetrics.local_flops,
                io_time_max);
            fflush(stdout); 
        }
    }

    double end_total = MPI_Wtime();
    if (rank == 0) {
        printf("Total execution time: %.3f seconds\n", end_total - start_total);
    }

    // Cleanup
    if (!use_parallel_io && !dummy) {
        free_sparse(&coo_matrix);
    }
    free_sparse(&local_matrix);
    free_csr(&local_csr);
    free(local_rand_vec);
    free(local_res);
    free(pmetrics.elapsed_times);
    free(pmetrics.comm_times);
    
    MPI_Finalize();
    return 0;
}