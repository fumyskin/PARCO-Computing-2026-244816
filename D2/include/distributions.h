#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "mmio.h"
#include "specifications.h"

void distribution_1D(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix, Sparse_CSR *local_csr, int rank, int comm_size);

#endif