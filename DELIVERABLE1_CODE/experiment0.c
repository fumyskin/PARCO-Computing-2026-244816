#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <pthread.h>
#include <omp.h>
#include "mmio.h"
#include "specifications.h"

//define fma block operation
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

// function to initialize a struct COO given the data extracted from .mtx file
Sparse_Coordinate* initialize_COO(
    int n_rows,
    int n_cols,
    int nnz,
    int* row_indices,
    int* col_indices,
    double* values
)
{
    Sparse_Coordinate* struct_COO = surely_malloc(sizeof(Sparse_Coordinate));
    struct_COO->n_rows = n_rows;
    struct_COO->n_cols = n_cols;
    struct_COO->nnz = nnz;
    struct_COO->row_indices = row_indices;
    struct_COO->col_indices = col_indices;
    struct_COO->values = values;

    return struct_COO;
}


// function to perform spmv on COO
void SpMV_COO(Sparse_Coordinate* COO, double* vec, double* res){
    for(int i = 0; i < COO->n_rows; i++){
        res[i] = 0;
    }

    for(ssize_t nnz_id = 0; nnz_id < COO->nnz; nnz_id++){
        ssize_t i = COO->row_indices[nnz_id];
        ssize_t j = COO->col_indices[nnz_id];
        double val = COO->values[nnz_id];

        res[i] += val * vec[j]; 
    } 

    return;
}


unsigned coo_count(Sparse_Coordinate *p){
    if (p == NULL || p->nnz == 0)
        return 0;

    unsigned i, n = p->nnz; 
    if (n == 0) return 0;
    unsigned count = 1;
    for (i=1; i<n; i++){
        if (p->row_indices[i-1] !=p->row_indices[i] || 
            p->col_indices[i-1] !=p->col_indices[i]){
            count++;
        }
    }
    return count;
}

Sparse_CSR *coo_to_csr_matrix(Sparse_Coordinate *p) {
    Sparse_CSR *q;
    unsigned count, i;
    unsigned r,c, ri, ci, cols, k, l, rows;
    unsigned *col_ind, *row_ptr, *prow_ind, *pcol_ind;
    double x, *val, *pval;
    unsigned n = p->nnz;
    coo_quicksort(p, 0, n);
    k = coo_count(p);
    rows = p->n_rows;
    prow_ind=p->row_indices;
    pcol_ind=p->col_indices;
    pval = p->values;
    q = surely_malloc(sizeof(Sparse_CSR));
    val = surely_malloc(k * sizeof(double));
    col_ind = surely_malloc(k * sizeof(unsigned));
    row_ptr = surely_malloc((rows+1) * sizeof(unsigned));
    r=-1;
    c=0; 
    l=0;
    /* partial_csr_0 */
    for (i=0; i<n; i++) {
        ri = prow_ind[i];
        ci = pcol_ind[i];
        x = pval[i];
        if (ri==r){
            if (ci==c)
                val[l-1] += x; /* partial_csr_duplicate */
            else {
                c=ci;
                col_ind[l] = ci;
                val[l] = x;
                l++;           /* partial_csr_newcol */
            }
        }
        else{
            while (r+1<=ri) row_ptr[++r]=l; /* partial_csr_skiprow */
            c= ci;
            col_ind[l] = ci;
            val[l] = x;
            l++;            /* partial_csr_newrow */
        }
    }
    cols = p->n_cols;
    while (r+1<=rows) row_ptr[++r]=l;  /* partial_csr_lastrows */
    q->values = val;
    q->col_ind = col_ind;
    q->row_ptr = row_ptr;
    q->n_rows = rows;
    q->n_cols = cols;
    return q;          /* partial_CSR_properties */
}

/*
For SpMV, focus on memory/cache optimizations first (reordering, 
blocking, prefetching, improve locality, reduce indirection) 
— they yield larger gains. Then focus on optimizing computation (vectorization,
parallelization)
*/
void csr_mv_multiply(Sparse_CSR *m, double *v, double *p) {
    unsigned i, rows = m->n_rows;
    double *val = m->values;
    unsigned *col_ind = m->col_ind;
    unsigned *row_ptr = m->row_ptr;
    unsigned next=row_ptr[0];

    // sequential implementation
    for (i = 0; i < rows; i++) {
        double s = 0.0; // private scope to each thread
        for (unsigned h = row_ptr[i]; h < row_ptr[i + 1]; h++) {
            double x = val[h];
            unsigned j = col_ind[h];
            s = fma_fallback(x, v[j], s);
        }
        p[i] = s;
    }

}


int main(int argc, char *argv[])
{
    //FOR NOW I'LL USE THE EXAMPLE GIVEN
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   
    int i, *I, *J;
    double *val;

    /*Initialize struct for sparse matrix */
    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else    
    { 
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    //  This is how one can screen matrix types if their application 
    //  only supports a subset of the Matrix Market data types.     
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    // find out size of sparse matrix .... 
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    // reseve memory for matrices 
    I = (int *) surely_malloc(nz * sizeof(int));
    J = (int *) surely_malloc(nz * sizeof(int));
    val = (double *) surely_malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    // for (i=0; i<nz; i++){
    //     fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);
    // }

    //create struct with data read from .mtx file
    Sparse_Coordinate* struct_COO = initialize_COO(M, N, nz, I, J, val);
    //convert COO to CSR
    //Sparse_CSR* struct_CSR = convert_COO_CSR(M, N, nz, &struct_COO);

    //INITIALIZE MATRIX VECTOR MULTIPLICATION
    double* res = surely_malloc(M * sizeof(double));
    double* vec = surely_malloc(N * sizeof(double));

    //INITIALIZE RANDOM VECTOR
    srand(0);
    for(int i = 0; i < N; i++){
        vec[i] = rand() % 10;
    }

    //compute SpMV with COO 
    SpMV_COO(struct_COO, vec, res);

    //INITIALIZE CSR MATRIX FROM COO
    Sparse_CSR* struct_CSR = coo_to_csr_matrix(struct_COO);
    double* res_csr = malloc(M * sizeof(double));

    //COMPUTE SpMV WITH CSR
    double start = omp_get_wtime();
    csr_mv_multiply(struct_CSR, vec, res_csr);
    double end = omp_get_wtime();
    printf("\nElapsed time: %g seconds\n", end - start);
    
    free(I);
    free(J);
    free(val);
    free(vec);
    free(res);
    free(res_csr);
    free(struct_CSR);      
    free(struct_COO);
    
    return 0;
}





/*
NOTES ON PARALLELIZATION
CACHE

- efficient parallel code we need to assure that SEQUENTIAL PARTS ARE PERFORMANT
- TRY TO COORDINATE SEQUENTIAL PARTS in such a way that exploit cache as much as possible (eg row major, ...)
- WARNING:
    - FALSE SHARING MAY BE A PROBLEM
    -> if you use cache, pay attention to cache lines: if data is invalidated, the entire content of cache line is
    -> FORCE VARIABLES WHICH ARE ACCESSED BY DIFFERENT THREADS TO BE ON DIFFERENT CACHE LINES
    
        struct alignTo64ByteCacheLine {
            int _onCacheLine1 __attribute__((aligned(64)))
            int _onCacheLine2 __attribute__((aligned(64)))
        }

    - BRANCH PREDICTION
    -> Random data leads to unpredictable branches, slowing execution
    -> SORTING data can improve branch prediction and speed up execution

*/