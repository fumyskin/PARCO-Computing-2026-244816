#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <pthread.h>
#include <omp.h> 
#include "mmio.h"
#include "specifications.h"


/*
*    PREMISE:
*    - in case one wants to try optimized version, run with:
*      gcc -O3 mmio.h specifications.c experiment2.c -fopenmp -ffp-contract=fast -o spmv_manual
*
*    - otherwise, withouth optimization of brach prediction, run with:
*      gcc mmio.h specifications.c experiment2.c -fopenmp -o spmv_manual
*
*/

// padding structure to avoid false sharing between threads
// not used in the end but a possible interesting solutions would have been
// to implement this struct to prevent also boundary false sharing 
// typedef struct alignTo64ByteCacheLine {
//     char pad[64];
// } alignTo64ByteCacheLine;


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

// function to initialize a struct COO given the data extracted from .mtx file
Sparse_Coordinate* initialize_COO(
    unsigned n_rows,
    unsigned n_cols,
    unsigned nnz,
    unsigned* row_indices,
    unsigned* col_indices,
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


// function to perform Spmv on COO
void SpMV_COO(Sparse_Coordinate* COO, double* vec, double* res){
    for(unsigned i = 0; i < COO->n_rows; i++){
        res[i] = 0;
    }

    for(unsigned nnz_id = 0; nnz_id < COO->nnz; nnz_id++){
        unsigned i = COO->row_indices[nnz_id];
        unsigned j = COO->col_indices[nnz_id];
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
    // partial_csr_0 
    for (i=0; i<n; i++) {
        ri = prow_ind[i];
        ci = pcol_ind[i];
        x = pval[i];
        if (ri==r){
            if (ci==c)
                val[l-1] += x; // partial_csr_duplicate
            else {
                c=ci;
                col_ind[l] = ci;
                val[l] = x;
                l++;           // partial_csr_newcol
            }
        }
        else{
            while (r+1<=ri) row_ptr[++r]=l; // partial_csr_skiprow 
            c= ci;
            col_ind[l] = ci;
            val[l] = x;
            l++;            // partial_csr_newrow 
        }
    }
    cols = p->n_cols;
    while (r+1<=rows) row_ptr[++r]=l;  // partial_csr_lastrows 
    q->values = val;
    q->col_ind = col_ind;
    q->row_ptr = row_ptr;
    q->n_rows = rows;
    q->n_cols = cols;
    return q;          // partial_CSR_properties 
}

/*
For SpMV, focus on memory/cache optimizations first (reordering, 
blocking, prefetching, improve locality, reduce indirection) 
— they yield larger gains. Then focus on optimizing computation (vectorization,
parallelization)
*/
void csr_mv_multiply(Sparse_CSR *m, double *v, double *p) {
    if (!m || !v || !p) return;  // null check

    unsigned rows = m->n_rows;
    double *val = m->values;
    unsigned *col_ind = m->col_ind;
    unsigned *row_ptr = m->row_ptr;

    //  manual static partitioning: each thread writes a contiguous block of rows
    //  (avoids different threads writing elements that share a cache line) 
    #pragma omp parallel            
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        unsigned start = (rows * (unsigned)tid) / (unsigned)nthreads;
        unsigned end   = (rows * (unsigned)(tid + 1)) / (unsigned)nthreads;

        for (unsigned i = start; i < end; ++i) {
            double s = 0.0;
            for (unsigned h = row_ptr[i]; h < row_ptr[i + 1]; ++h) {
                s = fma_fallback(val[h], v[col_ind[h]], s);
            }
            p[i] = s;
        }
    }
}



int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   // M=rows, N=cols, nz=nonzeroes
    int i;
    unsigned *I, *J;
    double *val;

    // Initialize struct for sparse matrix 
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


    // reserve memory for matrices 
    I = (unsigned *) surely_malloc(nz * sizeof(unsigned));
    J = (unsigned *) surely_malloc(nz * sizeof(unsigned));
    val = (double *) surely_malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    for (i=0; i<nz; i++)
    {
        int temp_i, temp_j;
        fscanf(f, "%d %d %lg\n", &temp_i, &temp_j, &val[i]);
        I[i] = (unsigned)(temp_i - 1);
        J[i] = (unsigned)(temp_j - 1);
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

    // create struct with data read from .mtx file
    Sparse_Coordinate* struct_COO = initialize_COO((unsigned)M, (unsigned)N, (unsigned)nz, I, J, val);

    // INITIALIZE MATRIX VECTOR MULTIPLICATION
    double* res;
    double* vec;

    //allocate aligned result/vector buffers to 64 bytes boundary  
    /*  posix_memaling from man:
    *
    *   int posix_memalign(void **memptr, size_t alignment, size_t size);
    *   
    *   posix_memalign() allocates size bytes and places the address of the allocated memory in *memptr.  The address of the allo‐
    *   cated memory will be a multiple of alignment, which must be a power of two and a multiple of sizeof(void *).  This address
    *   can  later  be  successfully passed to free(3).  If size is 0, then the value placed in *memptr is either NULL or a unique
    *   pointer value.
    *
    */
    if (posix_memalign((void**)&res, 64, M * sizeof(double)) != 0) { // allocate res pointer aligned to 64 bytes
        perror("posix_memalign res");
        exit(1);
    }
    if (posix_memalign((void**)&vec, 64, N * sizeof(double)) != 0) {
        perror("posix_memalign vec");
        exit(1);
    }

    // INITIALIZE RANDOM VECTOR
    srand(0);
    for(int i = 0; i < N; i++){
        vec[i] = rand() % 10;
    }
    

    //INITIALIZE CSR MATRIX FROM COO
    Sparse_CSR* struct_CSR = coo_to_csr_matrix(struct_COO);
    double *res_csr;
    if (posix_memalign((void**)&res_csr, 64, M * sizeof(double)) != 0) {
        perror("posix_memalign res_csr");
        exit(1);
    }

    int N_ITERS = 20;
    double total = 0;

    for(int i = 0; i < N_ITERS; i++) {
        double t0 = omp_get_wtime();
        csr_mv_multiply(struct_CSR, vec, res_csr);
        double t1 = omp_get_wtime();
        total += (t1 - t0);
    }

    double avg_time = total/N_ITERS;
    printf("\nElapsed time: %.6f seconds\n", avg_time);



    //COMPUTE SpMV WITH CSR
    // double start = omp_get_wtime();
    // csr_mv_multiply(struct_CSR, vec, res_csr);
    // double end = omp_get_wtime();
    // printf("\nElapsed time: %g seconds\n", end - start);

    // //COMPUTE SpMV WITH COO (for verification)
    // SpMV_COO(struct_COO, vec, res);

    // printf("\nResult (first 10 entries):\n");
    // for (int i = 0; i < M && i < 10; i++) {
    //     printf("res[%d] = %g\n", i, res[i]);
    // }

    // printf("\nCSR Result (first 10 entries):\n");
    // for (int i = 0; i < M && i < 10; i++) {
    //     printf("res_csr[%d] = %g\n", i, res_csr[i]);
    // } 
 
    
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

