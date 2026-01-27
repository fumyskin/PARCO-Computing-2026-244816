
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <omp.h> 

#include "mmio.h"
#include "specifications.h"
#include "distributions.h"


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

    if (argc != 3) {
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
    int dummy_matrix = (strcmp(matrix, "dummy") == 0);

    Sparse_Coordinate sparse_matrix;
    unsigned matrix_rows, matrix_cols, matrix_nnz;

    if (rank == 0){
        if (load_matrix_market(matrix, &sparse_matrix) != 0) {
                fprintf(stderr, "Failed to load matrix\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    
        matrix_rows = sparse_matrix.n_rows;
        matrix_cols = sparse_matrix.n_cols;
        matrix_nnz  = sparse_matrix.nnz;

        printf("Real Matrix (Strong Scaling): Rows=%d | Cols=%d | NNZ=%d | Procs=%d\n", matrix_rows, matrix_cols, matrix_nnz, comm_size);

    }

    MPI_Bcast(&matrix_rows, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&matrix_cols, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&matrix_nnz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    // allocate memory on non-root ranks
    if (rank != 0) {
        sparse_matrix.n_rows = matrix_rows;
        sparse_matrix.n_cols = matrix_cols;
        sparse_matrix.nnz = matrix_nnz;
    
        sparse_matrix.row_indices = (unsigned*)surely_malloc(matrix_nnz * sizeof(unsigned));
        sparse_matrix.col_indices = (unsigned*)surely_malloc(matrix_nnz * sizeof(unsigned));
        sparse_matrix.values      = (double*)surely_malloc(matrix_nnz * sizeof(double));
    }

    // broadcast the actual data
    MPI_Bcast(sparse_matrix.row_indices, matrix_nnz, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(sparse_matrix.col_indices, matrix_nnz, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(sparse_matrix.values, matrix_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // idea: sparse_coordinate matrix -> bcast to all processes-> call distr_1D for row selection->
    // -> get local sparse matrix and local csr matrix for each
    // 1D distribution
    Sparse_Coordinate local_matrix;
    Sparse_CSR local_csr;
    distribution_1D(&sparse_matrix, &local_matrix, &local_csr, rank, comm_size);

    if (rank == 0){
        printf("ok 1D distribtuion\n");
    }

    // removing from the heap
    free_sparse(&sparse_matrix);
    free_sparse(&local_matrix);
    free_csr(&local_csr);
    
    MPI_Finalize();
    return 0;
}


