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
    // Implementation of the baseline
    // 1. Read a matrix in matrix market format
    // Rank 0 reads the entire file and distributes matrix entries to all other processes.
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

    // changed to consider also 2d 
    if (argc < 4) {
        if (rank == 0){
            // print a better message ?
            fprintf(stderr, "Usage: %s <matrix_file|dummy> <repeats> <distribution_type> [rows_per_proc] [nnz_per_row]\n", argv[0]);
            fprintf(stderr, "  distribution_type: 1D, 2D, or 2D_cyclic\n");
            fprintf(stderr, "  For dummy matrices: provide rows_per_proc and nnz_per_row\n");
            fprintf(stderr, "Error: Expected at least 3 arguments, but received %d.\n", argc - 1);;
        }
        MPI_Finalize();
        return EXIT_FAILURE; 
    }

    char* matrix = argv[1];
    int repeats = atoi(argv[2]);
    char* distr_type = argv[3];
    int dummy = (strcmp(matrix, "dummy") == 0);

    // check distribution type
    int use_2d = 0;
    int use_cyclic = 0;
    if (strcmp(distr_type, "1D") == 0){
        use_2d = 0;
    }else if (strcmp(distr_type, "2D") == 0){
        use_2d = 1;
        use_cyclic = 0;
    }else if (strcmp(distr_type, "2D_cyclic") == 0){
        use_2d = 1;
        use_cyclic = 1;
    }else{
        if (rank == 0){
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


    // input check
    if (dummy){
        if(argc < 6){
            if(rank == 0){
                fprintf(stderr, "Error: weak scaling simulation requires #rows per process and #nnz per row\n");
            }
            MPI_Finalize();
            return EXIT_FAILURE;
        }

        rows_per_process = atoi(argv[4]);
        nnz_per_row = atoi(argv[5]);
    }

    if(!dummy){
        if (rank == 0){
            if (load_matrix_market(matrix, &coo_matrix) != 0) {
                    fprintf(stderr, "Failed to load matrix\n");
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            matrix_rows = coo_matrix.n_rows;
            matrix_cols = coo_matrix.n_cols;
            matrix_nnz  = coo_matrix.nnz;

            // allocation for fields of csr happens in function
            Sparse_CSR *temp_csr = coo_to_csr_matrix(&coo_matrix);
            csr_matrix = *temp_csr; 
            free(temp_csr);
           
            printf("Real Matrix: Rows=%d | Cols=%d | NNZ=%d | Procs=%d | Distribution=%s\n", 
                   matrix_rows, matrix_cols, matrix_nnz, comm_size, distr_type);

        }
        MPI_Bcast(&matrix_rows, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(&matrix_cols, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(&matrix_nnz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        // allocate memory on non-root ranks
        if (rank != 0) {
            coo_matrix.n_rows = matrix_rows;
            coo_matrix.n_cols = matrix_cols;
            coo_matrix.nnz = matrix_nnz;
        
            coo_matrix.row_indices = (unsigned*)surely_malloc(matrix_nnz * sizeof(unsigned));
            coo_matrix.col_indices = (unsigned*)surely_malloc(matrix_nnz * sizeof(unsigned));
            coo_matrix.values      = (double*)surely_malloc(matrix_nnz * sizeof(double));

        }
        // broadcast the actual data
        MPI_Bcast(coo_matrix.row_indices, matrix_nnz, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(coo_matrix.col_indices, matrix_nnz, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(coo_matrix.values, matrix_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }else{
        // dummy matrix setup
        if (use_2d){
            int proc_rows, proc_cols;
            create_2d_grid(comm_size, &proc_rows, &proc_cols);
            matrix_rows = rows_per_process * proc_rows;
            matrix_cols = matrix_rows;
        }else{
            matrix_rows = rows_per_process * comm_size;
            matrix_cols = matrix_rows;
        }
        matrix_nnz = 0;

        if (rank == 0) {
            printf("Synthetic Matrix (Weak Scaling): Rows/proc=%d | Total rows=%d | NNZ/row=~%d | Procs=%d\n", rows_per_process, matrix_rows, nnz_per_row, comm_size);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // idea: sparse_coordinate matrix -> bcast to all processes-> call distr_1D for row selection->
    // -> get local sparse matrix and local csr matrix for each
    // 1D distribution
    Sparse_Coordinate local_matrix = {0};
    Sparse_CSR local_csr = {0};
    if(!dummy){
        if(use_2d){
            if(use_cyclic){
                distribution_2D_cyclic(&coo_matrix, &local_matrix, &local_csr, rank, comm_size);
            }else{
                distribution_2D(&coo_matrix, &local_matrix, &local_csr, rank, comm_size);
            }
            if(rank == 0){
                printf("Applied 2D%s distribution\n", use_cyclic ? "cyclic" : "");
            }
        }else{
            distribution_1D(&coo_matrix, &local_matrix, &local_csr, rank, comm_size);
            if (rank == 0){
                printf("Applied 1D distribution\n");
            }
        }    
    }else{
        if(use_2d){
            generate_dummy_matrix_2d(&local_matrix, rows_per_process, nnz_per_row, rank, comm_size);
        }else{
            //generate a random local matrix for each process
            generate_dummy_matrix(&local_matrix, rows_per_process, nnz_per_row, rank, comm_size);
        }
        
        // conversion in CSR
        Sparse_CSR *temp_csr = coo_to_csr_matrix(&local_matrix);
        local_csr = *temp_csr;
        free(temp_csr); 

        unsigned local_nnz_count = local_matrix.nnz;
        MPI_Allreduce(&local_nnz_count, &matrix_nnz, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        
        // Also update global rows/cols if you use those variables later
        matrix_rows = rows_per_process * comm_size;
        matrix_cols = matrix_rows;

        if (rank == 0) {
            printf("Dummy Matrix Generated: %llu Global Rows | Total NNZ: %u\n", (unsigned long long)matrix_rows, matrix_nnz);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ghost elements identification -> WE FIND GHOST COLUMNS
    Comm_Pattern comm_pattern;
    find_ghost_vals(&local_csr, &comm_pattern, rank, comm_size, matrix_cols);

    //verify_ghost_entries(&local_csr, &comm_pattern, rank, comm_size, matrix_cols);  //ok

    // verification at this point done
    if(rank == 0){
        printf("Ghost columns identified: %d total ghost values needed\n", comm_pattern.num_ghost_cols);
    }

    // RANDOM VECTOR DISTRIBUTION (cyclic)
    unsigned local_vec_size = matrix_cols/comm_size;
    if(rank < (matrix_cols % comm_size)){
        local_vec_size++;
    }

    double *local_rand_vec = (double*)surely_malloc(local_vec_size * sizeof(double));
    double *local_res = (double*)surely_malloc(local_csr.n_rows * sizeof(double)); //?

    for (int i = 0; i < local_vec_size; i++) {
        // global index with cyclic distribution
        int global_idx = rank + i * comm_size;
        local_rand_vec[i] = ((double)rand() / RAND_MAX) * 8.0 - 4.0; //!!!
    }
    
    if (rank == 0) {
        printf("Local vector initialized (size=%d per process)\n", local_vec_size);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // CACHE WARMUP 
    exchange_ghost_vals(local_rand_vec, &comm_pattern, rank, comm_size);
    local_SpMV(&local_csr, local_rand_vec, local_res, &comm_pattern, rank, comm_size);
    MPI_Barrier(MPI_COMM_WORLD);

    // PERFORMANCE BENCHMARK OPERATIONS 
    Performance_Metrics pmetrics;
    pmetrics.local_nnz = local_matrix.nnz;
    pmetrics.ghost_entries = comm_pattern.num_ghost_cols;
    pmetrics.local_flops = 2LL * local_matrix.nnz; // fix
    pmetrics.num_repeats = repeats;

    pmetrics.elapsed_times = (double*)surely_malloc(repeats * sizeof(double));
    pmetrics.comm_times = (double*)surely_malloc(repeats * sizeof(double));

    for (int r = 0; r < repeats; r++) {
        MPI_Barrier(MPI_COMM_WORLD);   // sync all processes
        double total_iter_start = MPI_Wtime();   
        
        double comm_start = MPI_Wtime();       
        exchange_ghost_vals(local_rand_vec, &comm_pattern, rank, comm_size);
        double comm_end = MPI_Wtime();             
        pmetrics.comm_times[r] = comm_end - comm_start;
        
        local_SpMV(&local_csr, local_rand_vec, local_res, &comm_pattern, rank, comm_size);
        
        double total_iter_end = MPI_Wtime();  
        pmetrics.elapsed_times[r] = total_iter_end - total_iter_start;
    }

    

    // RAW OUTPUT DATA
    // fix visualization bug in benchmark results
    for (int r = 0; r < repeats; r++) {
        printf("[RESULT] %d,%d,%s,%d,%.9f,%.9f,%d,%d,%d\n",
            rank, 
            comm_size,
            distr_type,
            r,
            pmetrics.elapsed_times[r], 
            pmetrics.comm_times[r],
            pmetrics.local_nnz,
            pmetrics.ghost_entries,
            pmetrics.local_flops);
        fflush(stdout); 
    }

    if (rank == 0) {
        printf("\nBenchmark Summary\n");
        printf("Processes: %d\n", comm_size);
        printf("Distribution: %s\n", distr_type);
        printf("Iterations: %d\n", repeats);
        printf("Matrix NNZ: %d\n", matrix_nnz);
    }

    double end_total = MPI_Wtime();
    if (rank == 0) {
        printf("Total execution time: %.3f seconds\n", end_total - start_total);
    }


    // removing from the heap
    free_sparse(&coo_matrix);
    free_sparse(&local_matrix);
    free_csr(&local_csr);
    
    MPI_Finalize();
    return 0;
}




