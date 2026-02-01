#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "mmio.h"
#include "specifications.h"

void distribution_1D(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix, Sparse_CSR *local_csr, int rank, int comm_size);
void generate_dummy_matrix(Sparse_Coordinate* local_matrix, unsigned rows_per_process, unsigned nnz_per_row, int rank, int comm_size);

#endif