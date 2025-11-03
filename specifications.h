#ifndef SPECIFICATIONS_H
#define SPECIFICATIONS_H

#include<stddef.h>

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

void *surely_malloc(size_t size);
void coo_quicksort(Sparse_Coordinate *p, unsigned base, unsigned n);
unsigned coo_count(Sparse_Coordinate *p);
Sparse_CSR *coo_to_csr_matrix(Sparse_Coordinate *p);


#endif