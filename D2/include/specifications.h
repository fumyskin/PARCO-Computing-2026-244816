#ifndef SPECIFICATIONS_H
#define SPECIFICATIONS_H

#include<stddef.h>
#include "mmio.h"

typedef struct {
    unsigned n_rows;
    unsigned n_cols;
    unsigned nnz;
    unsigned *row_indices;
    unsigned *col_indices;
    double *values;
}Sparse_Coordinate;

typedef struct {
    unsigned n_rows;
    unsigned n_cols;
    unsigned* row_ptr;
    unsigned* col_ind;
    double* values;
}Sparse_CSR;

static inline double fma_fallback(double a, double b, double c);

void *surely_malloc(size_t size);
void coo_quicksort(Sparse_Coordinate *p, unsigned base, unsigned n);
int load_matrix_market(const char *filename, Sparse_Coordinate* matrix);
unsigned coo_count(Sparse_Coordinate *p);

Sparse_CSR *coo_to_csr_matrix(Sparse_Coordinate *p);
Sparse_Coordinate* initialize_COO(
    unsigned n_rows,
    unsigned n_cols,
    unsigned nnz,
    unsigned* row_indices,
    unsigned* col_indices,
    double* values
);

void SpMV_COO(Sparse_Coordinate* COO, double* vec, double* res);
void csr_mv_multiply(Sparse_CSR *m, double *v, double *p);
void free_sparse(Sparse_Coordinate * matrix);
void free_csr(Sparse_CSR *q);

#endif