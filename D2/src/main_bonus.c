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
    int rank, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    srand(time(NULL) + rank*1000);

    double start_total = MPI_Wtime();

    if (argc != 3) {
        if (rank == 0){
            fprintf(stderr, "Usage: %s <matrix_file> <repeats>\n", argv[0]);
            fprintf(stderr, "Error: Expected 2 arguments, but received %d.\n", argc - 1);
        }
        MPI_Finalize();
        return EXIT_FAILURE; 
    }

    char* matrix = argv[1];
    int repeats = atoi(argv[2]);

    // ALL processes call load_mm_chunked collectively
    Sparse_Coordinate local_matrix;
    if (load_mm_chunked(matrix, &local_matrix, rank, comm_size) != 0) {
        fprintf(stderr, "Rank %d: Failed to load matrix\n", rank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Each process now has its local portion in local_matrix
    // local_matrix.nnz is the LOCAL number of non-zeros for this rank
    
    if (rank == 0) {
        printf("Matrix loaded: Rows=%d | Cols=%d | Procs=%d\n", 
               local_matrix.n_rows, local_matrix.n_cols, comm_size);
    }

    // Now convert local coordinate format to CSR if needed
    Sparse_CSR* local_csr;
    local_csr = coo_to_csr_matrix(&local_matrix);

    if (rank == 0){
        printf("Matrix distributed and converted\n");
    }

    // COMPUTATIONS HERE FOR SPMV


    free_sparse(&local_matrix);
    free_csr(local_csr);
    
    MPI_Finalize();
    return 0;
}