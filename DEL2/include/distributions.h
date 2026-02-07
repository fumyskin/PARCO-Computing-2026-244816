#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "mmio.h"
#include "specifications.h"

// 1D DISTRIBUTION
void distribution_1D(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix, Sparse_CSR *local_csr, int rank, int comm_size);
void generate_dummy_matrix(Sparse_Coordinate* local_matrix, unsigned rows_per_process, unsigned nnz_per_row, int rank, int comm_size);

// 2D DISTRIBUTION
// create a 2D process grid (proc_rows x proc_cols)
void create_2d_grid(int comm_size, int *proc_rows, int *proc_cols);

// convert between rank and 2D coordinates
void rank_to_coords_2d(int rank, int proc_rows, int proc_cols, int *row_coord, int *col_coord);
int coords_to_rank_2d(int row_coord, int col_coord, int proc_cols);

// determine owner of a matrix element in 2D distribution
int owner_2d(unsigned row, unsigned col, unsigned matrix_rows, unsigned matrix_cols, int proc_rows, int proc_cols);
int owner_2d_cyclic(unsigned row, unsigned col, int proc_rows, int proc_cols);

// main 2D distribution functions
void distribution_2D(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix, Sparse_CSR *local_csr, int rank, int comm_size);
void distribution_2D_cyclic(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix, Sparse_CSR *local_csr, int rank, int comm_size);

// gnerate dummy matrix for 2D weak scaling tests
void generate_dummy_matrix_2d_block(Sparse_Coordinate *local_matrix, unsigned rows_per_proc, 
                                    unsigned nnz_per_row, int rank, int comm_size);
void generate_dummy_matrix_2d_cyclic(Sparse_Coordinate *local_matrix, unsigned rows_per_proc, 
                                     unsigned nnz_per_row, int rank, int comm_size);

#endif