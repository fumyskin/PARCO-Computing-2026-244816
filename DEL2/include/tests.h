#ifndef TEST_SPMV_CORRECTNESS_H
#define TEST_SPMV_CORRECTNESS_H

#include "mmio.h"
#include "specifications.h"
#include "communication.h"
#include "distributions.h"

/* Build a known 12x12 sparse test matrix (tridiagonal + far off-diagonal) */
void build_test_matrix(Sparse_Coordinate *coo);

/* Build test vector: x[i] = i + 1 */
void build_test_vector(double *x);

/* Serial reference: y_ref = A * x */
void compute_reference(Sparse_Coordinate *coo, double *x, double *y_ref);

/* Compare two vectors, return 1 if match, 0 if mismatch */
int compare_results(const char *label, double *y_ref, double *y_test, int n);

/* Test each distribution, compare against reference */
void test_1d(Sparse_Coordinate *coo, double *x_global, double *y_ref, int rank, int comm_size);
void test_2d_block(Sparse_Coordinate *coo, double *x_global, double *y_ref, int rank, int comm_size);
void test_2d_cyclic(Sparse_Coordinate *coo, double *x_global, double *y_ref, int rank, int comm_size);

#endif /* TEST_SPMV_CORRECTNESS_H */