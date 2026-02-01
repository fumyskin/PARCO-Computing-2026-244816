#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include "mmio.h"
#include "specifications.h"


//implement surely_malloc function
void *surely_malloc(size_t size) {
    void *ptr = malloc(size);
    if (ptr == NULL) {
        fprintf(stderr, "Out of memory while allocating %zu bytes\n", size);
        //MPI aborting
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    return ptr;
}

// quicksort taken from Appel paper and modified for COO struct : qsort3.c; 
// https://github.com/cverified/cbench-vst/blob/master/qsort/qsort3.c

int coo_less (Sparse_Coordinate *p, unsigned a, unsigned b) {
    unsigned ra = p->row_indices[a], rb = p->row_indices[b];
    if (ra<rb) return 1;
    if (ra>rb) return 0;
    return p->col_indices[a] < p->col_indices[b];
}

void swap(Sparse_Coordinate *p, unsigned a, unsigned b) {
    unsigned i,j; double x;
    i=p->row_indices[a];
    j=p->col_indices[a];
    x=p->values[a];
    p->row_indices[a]=p->row_indices[b];
    p->col_indices[a]=p->col_indices[b];
    p->values[a]=p->values[b];
    p->row_indices[b]=i;
    p->col_indices[b]=j;
    p->values[b]=x;
}

/* Lexicographic quicksort by (row, col) */
void coo_quicksort(Sparse_Coordinate *p, unsigned base, unsigned n)
{
    unsigned lo, hi, left, right, mid;

    if (n == 0)
    return;
    lo = base;
    hi = lo + n - 1;
    while (lo < hi) {
    mid = lo + ((hi - lo) >> 1);

    if (coo_less(p,mid,lo))
        swap(p, mid, lo);
    if (coo_less(p,hi,mid)) {
        swap(p, mid, hi);
        if (coo_less(p,mid,lo))
        swap(p, mid, lo);
    }
    left = lo + 1;
    right = hi - 1;
    do {
        while (coo_less(p,left,mid))
        left++;
        while (coo_less(p,mid,right))
        right--;
        if (left < right) {
    swap(p, left, right);
        if (mid == left)
            mid = right;
        else if (mid == right)
            mid = left;
        left++;
        right--;
        } else if (left == right) {
        left++;
        right--;
        break;
        }
    } while (left <= right);
    if (right - lo > hi - left) {
        coo_quicksort(p, left, hi - left + 1);
        hi = right;
    } else {
        coo_quicksort(p, lo, right - lo + 1);
        lo = left;
    }
    }
}

int load_matrix_market(const char *filename, Sparse_Coordinate* matrix)
{
    MM_typecode matcode;
    FILE *f;
    int ret_code;
    int M, N, nz;   // M=rows, N=cols, nz=nonzeroes
    int i;

    int temp_i, temp_j;

    printf("Retrieving the matrix from '%s'...\n", filename);

    if ((f = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "Error: could not open file %s\n", filename);
        return -1;
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
   
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
    matrix->n_rows = (unsigned)M;
    matrix->n_cols = (unsigned)N;
    matrix->nnz = (unsigned)nz;

    // allocate memory for struct members
    matrix->row_indices = (unsigned *) malloc(nz * sizeof(unsigned));
    matrix->col_indices = (unsigned *) malloc(nz * sizeof(unsigned));
    matrix->values = (double *) malloc(nz * sizeof(double));


    if (!matrix->row_indices || !matrix->col_indices || !matrix->values) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        if (f) fclose(f);
        return -1;
    }

    for (i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &temp_i, &temp_j, &matrix->values[i]);
        matrix->row_indices[i] = (unsigned)(temp_i - 1);
        matrix->col_indices[i] = (unsigned)(temp_j - 1);
    }

    if (f != stdin) fclose(f);

    // matrix write out -> ok
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, matrix->n_rows, matrix->n_cols, matrix->nnz);
    // for (i=0; i < matrix->nnz; i++){
    //     fprintf(stdout, "%d %d %20.19g\n", 
    //         matrix->row_indices[i]+1, matrix->col_indices[i]+1, matrix->values[i]);
    // }

    return 0;
}

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




void free_sparse(Sparse_Coordinate * matrix){

    free(matrix->values);
    free(matrix->col_indices);
    free(matrix->row_indices);
}

void free_csr(Sparse_CSR *q)
{
    free(q->values);
    free(q->col_ind);
    free(q->row_ptr);
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


