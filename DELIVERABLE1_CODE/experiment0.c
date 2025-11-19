#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <pthread.h>
#include <omp.h>
//#include <math.h>  
#include "mmio.h"
#include "specifications.h"


/*
*    PREMISE: 
*    - in case one wants to verify correctness of CSR COO SpMVM, uncomment the math.h library
*      and body in the main. Finally, run with the following flags: 
*      gcc mmio.h specifications.c experiment0.c -fopenmp -lm -o spmv_sequential
*      ./spmv_sequential MATRICES/<name_matrix_folder>/<name_matrix>.mtx
*
*    - in case one wants to try optimized version, run with:
*      gcc -O3 mmio.h specifications.c experiment0.c -fopenmp -ffp-contract=fast -o spmv_sequential
*
*    - otherwise, withouth optimization of brach prediction, run with:
*      gcc mmio.h specifications.c experiment0.c -fopenmp -o spmv_sequential
*
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
— they yield larger gains. Then focus on optimizing computation (eventual vectorization,
parallelization)
*/
void csr_mv_multiply(Sparse_CSR *m, double *v, double *p) {
    if (!m || !v || !p) return;  // null check

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

// /*
//  * Compare SpMV results computed with COO and CSR.
//  *
//  * Parameters:
//  *   coo   - pointer to Sparse_Coordinate (input COO matrix)
//  *   csr   - pointer to Sparse_CSR (CSR matrix converted from COO)
//  *   vec   - input vector (length = coo->n_cols)
//  *   tol   - tolerance for comparison (suggested 1e-12..1e-9)
//  *
//  * Returns:
//  *   1 if results match within tolerance, 0 otherwise.
//  *
//  * Notes:
//  *   - Uses a combined absolute+relative tolerance:
//  *       |a - b| <= tol_abs + tol_rel * max(|a|, |b|)
//  *   - Prints up to 10 mismatches for debugging.
//  */
// int compare_spmv_results(Sparse_Coordinate *coo, Sparse_CSR *csr, double *vec, double tol) {
//     if (coo == NULL || csr == NULL || vec == NULL) {
//         fprintf(stderr, "compare_spmv_results: NULL pointer input\n");
//         return 0;
//     }

//     int rows = coo->n_rows;
//     double *res_coo = (double *) surely_malloc(rows * sizeof(double));
//     double *res_csr = (double *) surely_malloc(rows * sizeof(double));
//     if (!res_coo || !res_csr) {
//         fprintf(stderr, "compare_spmv_results: allocation failed\n");
//         free(res_coo); free(res_csr);
//         return 0;
//     }

//     /* compute using COO */
//     SpMV_COO(coo, vec, res_coo);

//     /* compute using CSR */
//     csr_mv_multiply(csr, vec, res_csr);

//     /* compare */
//     const double tol_abs = tol * 1e-3; /* small absolute floor so very tiny values compare well */
//     int mismatches = 0;
//     int max_report = 10;
//     for (int i = 0; i < rows; ++i) {
//         double a = res_coo[i];
//         double b = res_csr[i];
//         double diff = a - b;
//         double adiff = fabs(diff);
//         double scale = fmax(fabs(a), fabs(b));
//         double allowed = tol_abs + tol * scale;

//         if (adiff > allowed) {
//             if (mismatches < max_report) {
//                 printf("Mismatch row %d: COO = %.17g, CSR = %.17g, diff = %.17g, allowed = %.17g\n",
//                        i, a, b, diff, allowed);
//             }
//             mismatches++;
//         }
//     }

//     if (mismatches == 0) {
//         printf("compare_spmv_results: OK — results match within tol = %g\n", tol);
//     } else {
//         printf("compare_spmv_results: FAILED — %d mismatches (showing up to %d)\n", mismatches, max_report);
//     }

//     free(res_coo);
//     free(res_csr);

//     return (mismatches == 0) ? 1 : 0;
// }




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
    double* res = surely_malloc(M * sizeof(double));
    double* vec = surely_malloc(N * sizeof(double));

    // INITIALIZE RANDOM VECTOR
    srand(0);
    for(int i = 0; i < N; i++){
        vec[i] = rand() % 10;
    }

    // INITIALIZE CSR MATRIX FROM COO
    Sparse_CSR* struct_CSR = coo_to_csr_matrix(struct_COO);
    double* res_csr = surely_malloc(M * sizeof(double));

    // COMPUTE SpMV WITH CSR
    double start = omp_get_wtime();
    csr_mv_multiply(struct_CSR, vec, res_csr);
    double end = omp_get_wtime();
    printf("\nElapsed time: %g seconds\n", end - start);
 

    // //verify CSR vs COO correctness -> VERIFIED AND CORRECT
    // double tol = 1e-12; // tolerance delta to compare similar results up to tol
    // int ok = compare_spmv_results(struct_COO, struct_CSR, vec, tol);
    // if (!ok) {
    //     fprintf(stderr, "Warning: CSR and COO SpMV differ!\n");
    //     // optional: you can exit non-zero here if you want immediate failure 
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