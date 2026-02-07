#ifndef SPECIFICATIONS_H
#define SPECIFICATIONS_H

#include<stddef.h>
#include <immintrin.h>
#include "mmio.h"

typedef struct Comm_Pattern Comm_Pattern;

// define fma block operation
#if defined(__x86_64__) && defined(__FMA__)
static inline double fma_fallback(double a, double b, double c) {
    __m128d A = _mm_set_sd(a);
    __m128d B = _mm_set_sd(b);
    __m128d C = _mm_set_sd(c);
    __m128d R = _mm_fmadd_sd(A, B, C);
    return _mm_cvtsd_f64(R);
}
#else
static inline double fma_fallback(double a, double b, double c) {
    return a * b + c; // may be fused by compiler with -O3 -ffp-contract=fast
}
#endif

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

typedef struct {
    unsigned local_nnz;
    unsigned ghost_entries;
    unsigned local_flops; //FLOPs: multiply + add per non-zero
    unsigned num_repeats;

    double* elapsed_times; // total time per iteration
    double* comm_times; // communication time per iteration 
}Performance_Metrics;

// OTHER UTILS
void *surely_malloc(size_t size);
void coo_quicksort(Sparse_Coordinate *p, unsigned base, unsigned n);
int load_matrix_market(const char *filename, Sparse_Coordinate* matrix);
int load_mm_chunked(const char *filename, Sparse_Coordinate* matrix, int rank, int comm_size);
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