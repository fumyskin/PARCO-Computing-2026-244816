#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "mmio.h"
#include "specifications.h"

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

int load_mm_chunked(const char *filename, Sparse_Coordinate* matrix, int rank, int comm_size){

    /*
    int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype, MPI_Status * status)

    MPI_Offset = file offset
    MPI_File_read_at_all is a collective routine that attempts to read from 
    the file associated with fh (at the offset position) a total number of 
    count data items having datatype type into the user’s buffer buf. 
    The offset is in etype units relative to the current view. 
    That is, holes are not counted when locating an offset. The data 
    is taken out of those parts of the file specified by the current view.
     MPI_File_read_at_all stores the number of datatype elements 
     actually read in status. All other fields of status 
     are undefined. It is erroneous to call this function if 
     MPI_MODE_SEQUENTIAL mode was specified when the file was opened.
    */

    MPI_File fh;
    MPI_Offset file_size;
    MPI_Status status;
    int ret;
    // every process in MPI_COMM_WORLD open the filename in RDMODE
    ret = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (ret != MPI_SUCCESS) {
        if (rank == 0) fprintf(stderr, "Error opening file: %s\n", filename);
        return -1;
    }

    // how many bytes each process reads
    MPI_File_get_size(fh, &file_size);

    MPI_Offset chunk_size = file_size / comm_size;
    MPI_Offset my_offset = rank * chunk_size;

    //  check for remainder: last one takes the rest -> mhhh is it the best approach?
    if (rank == comm_size - 1) {
        chunk_size = file_size - my_offset;
    }


    // allocate buffer (add padding for alignment/overlap!)
    int padding = 1024;
    MPI_Offset bytes_to_read = chunk_size + padding;
    if (my_offset + bytes_to_read > file_size) {
        bytes_to_read = file_size - my_offset;
    }

    char *buffer = (char*)surely_malloc(bytes_to_read+1); // +1 for null terminator

    // everyone reads their chunk simultaneously + checking we're not going out
    ret = MPI_File_read_at_all(fh, my_offset, buffer, bytes_to_read, MPI_CHAR, &status);
    if (ret != MPI_SUCCESS) {
        if (rank == 0) fprintf(stderr, "Error reading file\n");
        free(buffer);
        MPI_File_close(&fh);
        return -1;
    }
    
    buffer[bytes_to_read] = '\0'; // Null terminate for safety

    char* cursor = buffer;
    char* end_of_buffer = buffer + bytes_to_read;

    int n_rows = 0, n_cols = 0, nnz = 0;

    if (rank == 0){
        // deal with header ONLY ON RANK 0
        while (cursor < end_of_buffer && *cursor == '%') {
            char* next = memchr(cursor, '\n', end_of_buffer - cursor);
            if (!next) break;
            cursor = next + 1;
        }
        // extract n_rows, n_cols, nnz
        int args_found = sscanf(cursor, "%d %d %d", &n_rows, &n_cols, &nnz);
        if (args_found != 3){
            fprintf(stderr, "Error parsing matrix header\n");
            free(buffer);
            MPI_File_close(&fh);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }

        char* next = memchr(cursor, '\n', end_of_buffer - cursor);
        if (next) cursor = next + 1;
    }else{
        if (my_offset > 0) {
            char* eol = memchr(cursor, '\n', bytes_to_read);
            if (eol != NULL) {
                cursor = eol + 1;
            } else {
                // No complete line in this chunk
                cursor = end_of_buffer;
            }
        }
    }

    MPI_Bcast(&n_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int initial_capacity = (nnz / comm_size) + (nnz / (10 * comm_size)) + 100;
    unsigned local_count = 0;
    unsigned local_capacity = initial_capacity;

    unsigned* local_rows = (unsigned*)surely_malloc(local_capacity * sizeof(unsigned));
    unsigned* local_cols = (unsigned*)surely_malloc(local_capacity * sizeof(unsigned));
    double* local_vals = (double*)surely_malloc(local_capacity * sizeof(double));
    
    int* send_counts = (int*)calloc(comm_size, sizeof(int));
    int* send_capacities = (int*)surely_malloc(comm_size*sizeof(int));

    // 2d array
    unsigned** send_rows = (unsigned**)surely_malloc(comm_size*sizeof(unsigned*));
    unsigned** send_cols = (unsigned**)surely_malloc(comm_size*sizeof(unsigned*));
    double** send_vals   = (double**)surely_malloc(comm_size*sizeof(double*));

    int default_send_cap = 1000;
    for (int i = 0; i < comm_size; i++) {
        send_capacities[i] = default_send_cap;
        send_rows[i] = (unsigned*)surely_malloc(default_send_cap * sizeof(unsigned));
        send_cols[i] = (unsigned*)surely_malloc(default_send_cap * sizeof(unsigned));
        send_vals[i] = (double*)surely_malloc(default_send_cap * sizeof(double));
    }

    char* chunk_end = buffer + chunk_size;
    if (chunk_end > end_of_buffer) chunk_end = end_of_buffer;
    
    // parsing -> FOR ALL PROCESSES
    while(cursor < chunk_end){
        // next current line check
        char* next_newline = memchr(cursor, '\n', end_of_buffer - cursor);
        if (!next_newline) break;

        int curr_row, curr_col;
        double curr_val;
        int args_found = sscanf(cursor, "%d %d %lf", &curr_row, &curr_col, &curr_val);

        if (args_found == 3) {
            // adjust to 0 based ! mm is 1 based
            // curr_row--, curr_col-- 
            curr_row--;
            curr_col--;

            int destination_rank = curr_row % comm_size;

            if(destination_rank == rank){
                 if (local_count >= local_capacity) {
                    // grow local storage
                    local_capacity *= 2;
                    unsigned* new_rows = (unsigned*)realloc(local_rows, local_capacity * sizeof(unsigned));
                    unsigned* new_cols = (unsigned*)realloc(local_cols, local_capacity * sizeof(unsigned));
                    double* new_vals = (double*)realloc(local_vals, local_capacity * sizeof(double));
                    
                    if (!new_rows || !new_cols || !new_vals) {
                        fprintf(stderr, "Rank %d: Memory allocation failed\n", rank);
            
                        free(buffer);
                        MPI_File_close(&fh);
                        MPI_Abort(MPI_COMM_WORLD, 1);
                        return -1;
                    }
                    
                    local_rows = new_rows;
                    local_cols = new_cols;
                    local_vals = new_vals;
                }
                local_rows[local_count] = curr_row;
                local_cols[local_count] = curr_col;
                local_vals[local_count] = curr_val;
                local_count++;
            }else{
                int cnt = send_counts[destination_rank];
                int cap = send_capacities[destination_rank];

                if(cnt >= cap){
                    cap *= 2; // again, not sure if this makes sense
                    send_capacities[destination_rank] = cap;
                   
                    unsigned* new_rows = (unsigned*)realloc(send_rows[destination_rank], cap * sizeof(unsigned));
                    unsigned* new_cols = (unsigned*)realloc(send_cols[destination_rank], cap * sizeof(unsigned));
                    double* new_vals = (double*)realloc(send_vals[destination_rank], cap * sizeof(double));

                    if (!new_rows || !new_cols || !new_vals) {
                        fprintf(stderr, "Rank %d: Memory allocation failed for send buffer\n", rank);
                        free(buffer);
                        MPI_File_close(&fh);
                        MPI_Abort(MPI_COMM_WORLD, 1);
                        return -1;
                    }

                    send_rows[destination_rank] = new_rows;
                    send_cols[destination_rank] = new_cols;
                    send_vals[destination_rank] = new_vals;
                }

                send_rows[destination_rank][cnt] = curr_row;
                send_cols[destination_rank][cnt] = curr_col;
                send_vals[destination_rank][cnt] = curr_val;
                send_counts[destination_rank]++;
            }
        }
        cursor = next_newline + 1;
    }

    free(buffer);

    /*
    int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
        void *recvbuf, int recvcount, MPI_Datatype recvtype,
        MPI_Comm comm)
    */
    
    // EXCHANGE 
    int* recv_counts = (int*)surely_malloc(comm_size * sizeof(int));
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

    // displacements aaa
    int* sdispls = (int*)surely_malloc(comm_size * sizeof(int));
    int* rdispls = (int*)surely_malloc(comm_size * sizeof(int));
    int total_send = 0;
    int total_recv = 0;

    for (int i = 0; i < comm_size; i++) {
        sdispls[i] = total_send;
        rdispls[i] = total_recv;
        total_send += send_counts[i];
        total_recv += recv_counts[i];
    }

    unsigned* packed_s_rows = (unsigned*)surely_malloc(total_send * sizeof(unsigned));
    unsigned* packed_s_cols = (unsigned*)surely_malloc(total_send * sizeof(unsigned));
    double* packed_s_vals = (double*)surely_malloc(total_send * sizeof(double));

    for (int r = 0; r < comm_size; r++) {
        // copy each bucket into the big contiguous array at the right offset
        memcpy(&packed_s_rows[sdispls[r]], send_rows[r], send_counts[r] * sizeof(unsigned));
        memcpy(&packed_s_cols[sdispls[r]], send_cols[r], send_counts[r] * sizeof(unsigned));
        memcpy(&packed_s_vals[sdispls[r]], send_vals[r], send_counts[r] * sizeof(double));
        
        free(send_rows[r]);
        free(send_cols[r]);
        free(send_vals[r]);
    }

    unsigned* recv_rows_buf = (unsigned*)surely_malloc(total_recv * sizeof(unsigned));
    unsigned* recv_cols_buf = (unsigned*)surely_malloc(total_recv * sizeof(unsigned));
    double* recv_vals_buf = (double*)surely_malloc(total_recv * sizeof(double));


    // execute transfer 
    MPI_Alltoallv(packed_s_rows, send_counts, sdispls, MPI_UNSIGNED,
                  recv_rows_buf, recv_counts, rdispls, MPI_UNSIGNED, MPI_COMM_WORLD);
    
    MPI_Alltoallv(packed_s_cols, send_counts, sdispls, MPI_UNSIGNED,
                  recv_cols_buf, recv_counts, rdispls, MPI_UNSIGNED, MPI_COMM_WORLD);
    
    MPI_Alltoallv(packed_s_vals, send_counts, sdispls, MPI_DOUBLE,
                  recv_vals_buf, recv_counts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
    

    // MERGE 
    int new_total_nnz = local_count + total_recv;

    unsigned* final_rows = (unsigned*)realloc(local_rows, new_total_nnz * sizeof(unsigned));
    unsigned* final_cols = (unsigned*)realloc(local_cols, new_total_nnz * sizeof(unsigned));
    double* final_vals = (double*)realloc(local_vals, new_total_nnz * sizeof(double));

    if (!final_rows || !final_cols || !final_vals) {
        fprintf(stderr, "Rank %d: Final memory allocation failed\n", rank);
        free(buffer);
        MPI_File_close(&fh);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }

    local_rows = final_rows;
    local_cols = final_cols;
    local_vals = final_vals;

    memcpy(&local_rows[local_count], recv_rows_buf, total_recv * sizeof(unsigned));
    memcpy(&local_cols[local_count], recv_cols_buf, total_recv * sizeof(unsigned));
    memcpy(&local_vals[local_count], recv_vals_buf, total_recv * sizeof(double));


    matrix->n_rows = n_rows;
    matrix->n_cols = n_cols;
    matrix->nnz = new_total_nnz;  // local nnz for this rank
    matrix->row_indices = local_rows;
    matrix->col_indices = local_cols;
    matrix->values = final_vals;

    free(send_counts);
    free(recv_counts);
    free(sdispls);
    free(rdispls);
    free(send_capacities);
    free(send_rows);
    free(send_cols);
    free(send_vals);
    free(packed_s_rows);
    free(packed_s_cols);
    free(packed_s_vals);
    free(recv_rows_buf);
    free(recv_cols_buf);
    free(recv_vals_buf);
    
    MPI_File_close(&fh);
    
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


