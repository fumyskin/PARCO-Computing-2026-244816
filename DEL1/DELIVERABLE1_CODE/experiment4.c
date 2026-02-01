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
*      gcc -O3 mmio.h specifications.c experiment4.c -fopenmp -ffp-contract=fast -o spmv_bind
*
*    - otherwise, withouth optimization of brach prediction, run with:
*      gcc mmio.h specifications.c experiment4.c -fopenmp -o spmv_bind
*
*   - use compile1.sh script to run this code with proc binding
*
*/


/* ATTEMPT TO USE PROC_BINDING -> see compile1.sh */
/*
    NUMA INFORMATION (used -lscpu on cluster):
    Architecture:          x86_64
    CPU op-mode(s):        32-bit, 64-bit
    Byte Order:            Little Endian
    CPU(s):                96
    On-line CPU(s) list:   0-95
    Thread(s) per core:    1
    Core(s) per socket:    24
    Socket(s):             4
    NUMA node(s):          4
    Vendor ID:             GenuineIntel
    CPU family:            6
    Model:                 85
    Model name:            Intel(R) Xeon(R) Gold 6252N CPU @ 2.30GHz
    Stepping:              7
    CPU MHz:               2300.000
    BogoMIPS:              4600.00
    Virtualization:        VT-x
    L1d cache:             32K
    L1i cache:             32K
    L2 cache:              1024K
    L3 cache:              36608K
    NUMA node0 CPU(s):     0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92
    NUMA node1 CPU(s):     1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93
    NUMA node2 CPU(s):     2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90,94
    NUMA node3 CPU(s):     3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63,67,71,75,79,83,87,91,95
*/


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
// let's modify slightly the code
void csr_mv_multiply(Sparse_CSR *m, double *v, double *p) {
    unsigned i, rows = m->n_rows;
    double *val = m->values;
    unsigned *col_ind = m->col_ind;
    unsigned *row_ptr = m->row_ptr;
    unsigned next=row_ptr[0];

    // note1: private is not necessary and will throw a compilation error, since h, j and i are already private
    // note2: parallelisation happens only on outer loop -> each thread takes care of a cerain amount of rows 
    // (fixed row based partitioning)
    #pragma omp parallel for schedule(static)  
    for (i = 0; i < rows; i++) {
        double s = 0.0; // private scope to each thread
        for (unsigned h = row_ptr[i]; h < row_ptr[i + 1]; h++) {
            unsigned j = col_ind[h];
            s += val[h]*v[j];
        }
        p[i] = s;
    }
}


/*
Pthread implementation:
https://inria.hal.science/hal-03768726/document

*/

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

    //COMPUTE SpMV WITH CSR
    double start = omp_get_wtime();
    csr_mv_multiply(struct_CSR, vec, res_csr);
    double end = omp_get_wtime();
    printf("\nElapsed time: %g seconds\n", end - start);

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


