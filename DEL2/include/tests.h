#ifndef TESTS_H
#define TESTS_H

#include "mmio.h"
#include "specifications.h"
#include "communication.h"

void create_test_matrix(const char* filename);
void verify_ghost_entries(Sparse_CSR *local_csr, Comm_Pattern *comm_pattern, 
                          int rank, int comm_size, int total_cols);

#endif
