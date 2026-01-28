#include <stdio.h>
#include <stdlib.h>

#include "mmio.h"
#include "specifications.h"
#include "distributions.h"

// to implement yet:
// - local to global mapping
// - 2D (save as last)

void distribution_1D(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix, Sparse_CSR* local_csr, int rank, int comm_size){

    // 1D : start with all sparse coo matrices, then convert to csr
    int rows_per_process = matrix->n_rows/comm_size;
    // in case of nonzero remainder !! CHECK AGAIN
    if (rank < (matrix->n_rows % comm_size)) {
        rows_per_process++; //?
    }

    local_matrix->n_rows = rows_per_process;
    local_matrix->n_cols = matrix->n_cols;

    // find nnz: owner (i) = imod(P)
    local_matrix->nnz = 0;
    for(unsigned i = 0; i < matrix->nnz; i++){
        if(rank == matrix->row_indices[i] % comm_size){
            local_matrix->nnz++;
        }
    }
    // find row_indices and col_indices
    if (local_matrix->nnz != 0){
        local_matrix->row_indices = (unsigned*)surely_malloc((local_matrix->nnz)*sizeof(unsigned));
        local_matrix->col_indices = (unsigned*)surely_malloc((local_matrix->nnz)*sizeof(unsigned));
        local_matrix->values      = (double*)surely_malloc((local_matrix->nnz)*sizeof(double));
    }else{
        local_matrix->row_indices = NULL;
        local_matrix->col_indices = NULL;
        local_matrix->values      = NULL;
    }
    
    int temp_i = 0;
    for(int i = 0; i < matrix->nnz; i++){
        int assigned_row = matrix->row_indices[i];
        int owner_row = assigned_row % comm_size;

        if(owner_row == rank){
            local_matrix->row_indices[temp_i] = assigned_row/comm_size; // reframe index to local
            local_matrix->col_indices[temp_i] = matrix->col_indices[i];
            local_matrix->values[temp_i] = matrix->values[i];
            temp_i++;
        }
    }
    
    // local CSR conversion
    if (local_matrix->nnz != 0){
        Sparse_CSR *temp_csr = coo_to_csr_matrix(local_matrix);
       
        *local_csr = *temp_csr;
        free(temp_csr);
    } else {
        local_csr->n_rows = local_matrix->n_rows;
        local_csr->n_cols = local_matrix->n_cols;
        
        local_csr->row_ptr = (unsigned*)calloc(local_matrix->n_rows + 1, sizeof(unsigned));
        local_csr->col_ind = NULL;
        local_csr->values  = NULL;
    }

}

void generate_dummy_matrix(Sparse_Coordinate* local_matrix, unsigned rows_per_process, unsigned nnz_per_row, int rank, int comm_size) {
    
    unsigned long long global_matrix_rows = (unsigned long long)rows_per_process * comm_size;
    unsigned long long global_matrix_cols = global_matrix_rows; // ASSUMING SQUARE MATRIX !
    unsigned long long local_nnz = (unsigned long long)rows_per_process * nnz_per_row;
    
    local_matrix->row_indices = (unsigned*)surely_malloc(local_nnz * sizeof(unsigned));
    local_matrix->col_indices = (unsigned*)surely_malloc(local_nnz * sizeof(unsigned));
    local_matrix->values      = (double*)surely_malloc(local_nnz * sizeof(double));

    local_matrix->n_rows = global_matrix_rows;
    local_matrix->n_cols = global_matrix_cols; 
    local_matrix->nnz    = local_nnz;

    // seed random based on rank so every process generates unique, reproducible numbers.
    srand(rank * 12345 + 1); 

    unsigned current_idx = 0;

    for (unsigned k = 0; k < rows_per_process; k++) {
        // 1D CYCLIC RULE: 
        // The k-th row in my local memory corresponds to Global Row: (k * P) + rank
        // Example (P=4, Rank=1):
        // k=0 -> Global Row 1
        // k=1 -> Global Row 5
        // k=2 -> Global Row 9
        unsigned global_row_idx = (k * comm_size) + rank;
        
        for (unsigned n = 0; n < nnz_per_row; n++) {
    
            unsigned rand_col = rand() % global_matrix_cols;
            
            local_matrix->row_indices[current_idx] = global_row_idx;
            local_matrix->col_indices[current_idx] = rand_col;
            local_matrix->values[current_idx]      = 1.0; // use 1.0 for easy verification later
            
            current_idx++;
        }
    }
}