
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <omp.h> 
#include <stdbool.h>

#include "mmio.h"
#include "specifications.h"
#include "distributions.h"
#include "communication.h"


void create_test_matrix(const char* filename) {
    FILE *f = fopen(filename, "w");
    fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "4 4 8\n");
    // Row 0: columns 0, 1
    fprintf(f, "1 1 1.0\n");
    fprintf(f, "1 2 2.0\n");
    // Row 1: columns 2, 3
    fprintf(f, "2 3 3.0\n");
    fprintf(f, "2 4 4.0\n");
    // Row 2: columns 0, 2
    fprintf(f, "3 1 5.0\n");
    fprintf(f, "3 3 6.0\n");
    // Row 3: columns 1, 3
    fprintf(f, "4 2 7.0\n");
    fprintf(f, "4 4 8.0\n");
    fclose(f);
}

void verify_ghost_entries(Sparse_CSR *local_csr, Comm_Pattern *comm_pattern, 
                          int rank, int comm_size, int total_cols) {
    printf("\n[Rank %d] ===== VERIFICATION =====\n", rank);
    
    unsigned local_nnz = local_csr->row_ptr[local_csr->n_rows];
    bool *needs_ghost = (bool*)calloc(total_cols, sizeof(bool));
    int expected_ghosts = 0;
    
    // Find which columns we need but don't own
    for (unsigned i = 0; i < local_nnz; i++) {
        int col = local_csr->col_ind[i];
        int owner = col % comm_size;
        if (owner != rank && !needs_ghost[col]) {
            needs_ghost[col] = true;
            expected_ghosts++;
        }
    }
    
    // Verify count matches
    if (expected_ghosts == comm_pattern->num_ghost_cols) {
        printf("[Rank %d] ✓ Ghost count CORRECT: %d\n", rank, expected_ghosts);
    } else {
        printf("[Rank %d] ✗ Ghost count MISMATCH: expected %d, got %d\n", 
               rank, expected_ghosts, comm_pattern->num_ghost_cols);
    }
    
    // Verify all ghost columns are correct
    bool all_correct = true;
    for (int i = 0; i < comm_pattern->num_ghost_cols; i++) {
        int col = comm_pattern->ghost_col_indices[i];
        int owner = col % comm_size;
        
        if (owner == rank) {
            printf("[Rank %d] ✗ ERROR: Ghost col %d is owned by this rank!\n", rank, col);
            all_correct = false;
        }
        if (!needs_ghost[col]) {
            printf("[Rank %d] ✗ ERROR: Ghost col %d not in local matrix!\n", rank, col);
            all_correct = false;
        }
    }
    
    if (all_correct) {
        printf("[Rank %d] ✓ All ghost columns CORRECT\n", rank);
    }
    
    // Verify send/recv symmetry
    int *all_recv_totals = NULL;
    int *all_send_totals = NULL;
    if (rank == 0) {
        all_recv_totals = (int*)malloc(comm_size * sizeof(int));
        all_send_totals = (int*)malloc(comm_size * sizeof(int));
    }
    
    int my_recv_total = 0, my_send_total = 0;
    for (int p = 0; p < comm_size; p++) {
        my_recv_total += comm_pattern->recv_counts[p];
        my_send_total += comm_pattern->send_counts[p];
    }
    
    MPI_Gather(&my_recv_total, 1, MPI_INT, all_recv_totals, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&my_send_total, 1, MPI_INT, all_send_totals, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        int global_recv = 0, global_send = 0;
        for (int p = 0; p < comm_size; p++) {
            global_recv += all_recv_totals[p];
            global_send += all_send_totals[p];
        }
        
        if (global_recv == global_send) {
            printf("[Rank 0] ✓ Global send/recv symmetry CORRECT: %d == %d\n", 
                   global_send, global_recv);
        } else {
            printf("[Rank 0] ✗ Global send/recv MISMATCH: send=%d, recv=%d\n", 
                   global_send, global_recv);
        }
        
        free(all_recv_totals);
        free(all_send_totals);
    }
    
    // Verify ghost_to_local mapping
    int mapping_errors = 0;
    for (int i = 0; i < comm_pattern->num_ghost_cols; i++) {
        int col = comm_pattern->ghost_col_indices[i];
        int mapped_idx = comm_pattern->ghost_to_local[col];
        if (mapped_idx != i) {
            printf("[Rank %d] ✗ ghost_to_local[%d] = %d, expected %d\n", 
                   rank, col, mapped_idx, i);
            mapping_errors++;
        }
    }
    
    if (mapping_errors == 0) {
        printf("[Rank %d] ✓ ghost_to_local mapping CORRECT\n", rank);
    } else {
        printf("[Rank %d] ✗ ghost_to_local has %d errors\n", rank, mapping_errors);
    }
    
    free(needs_ghost);
    printf("[Rank %d] ===== END VERIFICATION =====\n\n", rank);
}


int main(int argc, char *argv[])
{
    // Implementation of the baseline
    // 1. Read a matrix in matrix market format
    // Rank 0 reads the entire file and distributes matrix entries to all other processes.
    int rank, comm_size; 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    srand(time(NULL) + rank*1000);

    double start_total = MPI_Wtime();

    // change
    if (argc < 3) {
        if (rank == 0){
            // print a better message
            fprintf(stderr, "Usage: %s <arg1> <arg2> <arg3>\n", argv[0]);
            fprintf(stderr, "Error: Expected 3 arguments, but received %d.\n", argc - 1);
        }
        MPI_Finalize();
        return EXIT_FAILURE; 
    }

    char* matrix = argv[1];
    int repeats = atoi(argv[2]);
    int dummy = (strcmp(matrix, "dummy") == 0);

    Sparse_Coordinate coo_matrix = {0};
    Sparse_CSR csr_matrix = {0};
    unsigned matrix_rows, matrix_cols, matrix_nnz;
    unsigned rows_per_process = 0;
    unsigned nnz_per_row = 0;


    // input check
    if (dummy){
        if(argc < 5){
            if(rank == 0){
                fprintf(stderr, "Error: weak scaling simulation requires #rows per process and #nnz per row\n");
            }
            MPI_Finalize();
            return EXIT_FAILURE;
        }

        rows_per_process = atoi(argv[3]);
        nnz_per_row = atoi(argv[4]);
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
           
            printf("Real Matrix (Strong Scaling): Rows=%d | Cols=%d | NNZ=%d | Procs=%d\n", matrix_rows, matrix_cols, matrix_nnz, comm_size);

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
        matrix_rows = rows_per_process*comm_size;
        matrix_cols = matrix_rows;
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
        distribution_1D(&coo_matrix, &local_matrix, &local_csr, rank, comm_size);
        if (rank == 0){
            printf("ok 1D distribtuion\n");
        }
    }else{
        //generate a random local matrix for each process
        generate_dummy_matrix(&local_matrix, rows_per_process, nnz_per_row, rank, comm_size);
        
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
            printf("Dummy Matrix Generated: %llu Global Rows | Total NNZ: %u\n", 
                   (unsigned long long)matrix_rows, matrix_nnz);
        }
    }

    // ghost elements identification
    Comm_Pattern comm_pattern;
    find_ghost_entries(&local_csr, &comm_pattern, rank, comm_size, matrix_cols);

   // Add verification
    verify_ghost_entries(&local_csr, &comm_pattern, rank, comm_size, matrix_cols);

    if(rank == 0){
        printf("Ghost columns ok\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // removing from the heap
    free_sparse(&coo_matrix);
    free_sparse(&local_matrix);
    free_csr(&local_csr);
    
    MPI_Finalize();
    return 0;
}




