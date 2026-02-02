#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

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



// 2D DISTRIBUTION
// helper function to create a 2D process grid
// returns the process grid dimensions (proc_rows x proc_cols)
void create_2d_grid(int comm_size, int *proc_rows, int *proc_cols) {
    // try to create a square-ish grid for better balance
    // find the closest factors to sqrt(comm_size)
    int sqrt_p = (int)sqrt((double)comm_size);
    
    for (int i = sqrt_p; i >= 1; i--) {
        if (comm_size % i == 0) {
            *proc_rows = i;
            *proc_cols = comm_size / i;
            return;
        }
    }
    
    // fallback (should never happen)
    *proc_rows = 1;
    *proc_cols = comm_size;
}

// get 2D coordinates from rank
void rank_to_coords_2d(int rank, int proc_rows, int proc_cols, int *row_coord, int *col_coord) {
    *row_coord = rank / proc_cols;  
    *col_coord = rank % proc_cols;
}

// get rank from 2D coordinates
int coords_to_rank_2d(int row_coord, int col_coord, int proc_cols) {
    return row_coord * proc_cols + col_coord;
}

// determine which process owns a matrix element (i, j) in 2D block distribution
int owner_2d(unsigned row, unsigned col, unsigned matrix_rows, unsigned matrix_cols, 
             int proc_rows, int proc_cols) {
    // block distribution: divide matrix into blocks
    unsigned rows_per_proc = (matrix_rows + proc_rows - 1) / proc_rows;
    unsigned cols_per_proc = (matrix_cols + proc_cols - 1) / proc_cols;
    
    int owner_row = row / rows_per_proc;
    int owner_col = col / cols_per_proc;
    
    if (owner_row >= proc_rows) owner_row = proc_rows - 1;
    if (owner_col >= proc_cols) owner_col = proc_cols - 1;
    
    return coords_to_rank_2d(owner_row, owner_col, proc_cols);
}

// Alternative: Cyclic 2D distribution (better load balance for irregular matrices)
int owner_2d_cyclic(unsigned row, unsigned col, int proc_rows, int proc_cols) {
    int owner_row = row % proc_rows;
    int owner_col = col % proc_cols;
    return coords_to_rank_2d(owner_row, owner_col, proc_cols);
}

// main
void distribution_2D(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix, 
                     Sparse_CSR *local_csr, int rank, int comm_size) {
    
    // 1) create 2D process grid
    int proc_rows, proc_cols;
    create_2d_grid(comm_size, &proc_rows, &proc_cols);
   
    int my_row_coord, my_col_coord;
    rank_to_coords_2d(rank, proc_rows, proc_cols, &my_row_coord, &my_col_coord);
    
    if (rank == 0) {
        printf("2D Process Grid: %d x %d (rows x cols)\n", proc_rows, proc_cols);
    }
    
    // 2) determine local matrix dimensions -> use block distribution
    unsigned rows_per_proc = (matrix->n_rows + proc_rows - 1) / proc_rows;
    unsigned cols_per_proc = (matrix->n_cols + proc_cols - 1) / proc_cols;
    
    unsigned my_start_row = my_row_coord * rows_per_proc;
    unsigned my_end_row = my_start_row + rows_per_proc;
    if (my_end_row > matrix->n_rows) {
        my_end_row = matrix->n_rows;
    }
    
    unsigned my_start_col = my_col_coord * cols_per_proc;
    unsigned my_end_col = my_start_col + cols_per_proc;
    if (my_end_col > matrix->n_cols) {
        my_end_col = matrix->n_cols;
    }
    
    local_matrix->n_rows = my_end_row - my_start_row;
    local_matrix->n_cols = matrix->n_cols;  // keep full column space for SpMV ?
    
    // 3) count local non-zeros (first pass)
    local_matrix->nnz = 0;
    for (unsigned i = 0; i < matrix->nnz; i++) {
        unsigned global_row = matrix->row_indices[i];
        unsigned global_col = matrix->col_indices[i];
        
        // check if this element belongs to this process using block distribution
        int element_owner = owner_2d(global_row, global_col, matrix->n_rows, matrix->n_cols, proc_rows, proc_cols);
        
        if (element_owner == rank) {
            local_matrix->nnz++;
        }
    }
   
    // 4) allocation
    if (local_matrix->nnz > 0) {
        local_matrix->row_indices = (unsigned*)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->col_indices = (unsigned*)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->values = (double*)surely_malloc(local_matrix->nnz * sizeof(double));
    } else {
        local_matrix->row_indices = NULL;
        local_matrix->col_indices = NULL;
        local_matrix->values = NULL;
    }
    
    // 5) fill local arrays (second pass)
    unsigned local_idx = 0;
    for (unsigned i = 0; i < matrix->nnz; i++) {
        unsigned global_row = matrix->row_indices[i];
        unsigned global_col = matrix->col_indices[i];
        
        int element_owner = owner_2d(global_row, global_col, matrix->n_rows, 
                                     matrix->n_cols, proc_rows, proc_cols);
        
        if (element_owner == rank) {
            // Convert to local row index (subtract the starting row)
            local_matrix->row_indices[local_idx] = global_row - my_start_row;
            local_matrix->col_indices[local_idx] = global_col;  // Keep global column index
            local_matrix->values[local_idx] = matrix->values[i];
            local_idx++;
        }
    }
    
    // 6) convert to CSR format
    if (local_matrix->nnz > 0) {
        Sparse_CSR *temp_csr = coo_to_csr_matrix(local_matrix);
        *local_csr = *temp_csr;
        free(temp_csr);
    } else {
        local_csr->n_rows = local_matrix->n_rows;
        local_csr->n_cols = local_matrix->n_cols;
        local_csr->row_ptr = (unsigned*)calloc(local_matrix->n_rows + 1, sizeof(unsigned));
        local_csr->col_ind = NULL;
        local_csr->values = NULL;
    }
    
    // debug output
    printf("Rank %d: Grid pos (%d,%d) | Rows [%u-%u) | Cols [%u-%u) | Local NNZ: %u\n",
           rank, my_row_coord, my_col_coord, my_start_row, my_end_row, 
           my_start_col, my_end_col, local_matrix->nnz);
}




// generate dummy matrix for 2D distribution (weak scaling)
void generate_dummy_matrix_2d(Sparse_Coordinate *local_matrix, unsigned rows_per_proc,
                              unsigned nnz_per_row, int rank, int comm_size) {
    
    int proc_rows, proc_cols;
    create_2d_grid(comm_size, &proc_rows, &proc_cols);
    
    int my_row_coord, my_col_coord;
    rank_to_coords_2d(rank, proc_rows, proc_cols, &my_row_coord, &my_col_coord);
    
    unsigned total_rows = rows_per_proc * proc_rows;
    unsigned total_cols = total_rows;  // Square matrix
    
    local_matrix->n_rows = rows_per_proc;
    local_matrix->n_cols = total_cols;
    
    // Generate approximately nnz_per_row non-zeros per row
    unsigned estimated_nnz = rows_per_proc * nnz_per_row;
    
    // Allocate temporary arrays (we'll trim later)
    unsigned *temp_rows = (unsigned*)surely_malloc(estimated_nnz * 2 * sizeof(unsigned));
    unsigned *temp_cols = (unsigned*)surely_malloc(estimated_nnz * 2 * sizeof(unsigned));
    double *temp_vals = (double*)surely_malloc(estimated_nnz * 2 * sizeof(double));
    
    unsigned actual_nnz = 0;
    
    // Generate non-zeros for each local row
    for (unsigned local_row = 0; local_row < rows_per_proc; local_row++) {
        unsigned global_row = my_row_coord + local_row * proc_rows;  // Cyclic mapping
        
        // Add diagonal element
        if (global_row < total_cols) {
            temp_rows[actual_nnz] = local_row;
            temp_cols[actual_nnz] = global_row;
            temp_vals[actual_nnz] = 4.0;
            actual_nnz++;
        }
        
        // Add random off-diagonal elements
        for (unsigned k = 1; k < nnz_per_row && actual_nnz < estimated_nnz * 2; k++) {
            unsigned global_col = rand() % total_cols;
            
            if (global_col != global_row) {
                temp_rows[actual_nnz] = local_row;
                temp_cols[actual_nnz] = global_col;
                temp_vals[actual_nnz] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
                actual_nnz++;
            }
        }
    }
    
    // Copy to final arrays
    local_matrix->nnz = actual_nnz;
    local_matrix->row_indices = (unsigned*)surely_malloc(actual_nnz * sizeof(unsigned));
    local_matrix->col_indices = (unsigned*)surely_malloc(actual_nnz * sizeof(unsigned));
    local_matrix->values = (double*)surely_malloc(actual_nnz * sizeof(double));
    
    for (unsigned i = 0; i < actual_nnz; i++) {
        local_matrix->row_indices[i] = temp_rows[i];
        local_matrix->col_indices[i] = temp_cols[i];
        local_matrix->values[i] = temp_vals[i];
    }
    
    free(temp_rows);
    free(temp_cols);
    free(temp_vals);
}

/////////////////////////////////////////////////////////////////////////////////////////////

// AAAA
// alternative: 2D Cyclic distribution (better for load balancing)
void distribution_2D_cyclic(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix, 
                            Sparse_CSR *local_csr, int rank, int comm_size) {
    
    // Step 1: Create 2D process grid
    int proc_rows, proc_cols;
    create_2d_grid(comm_size, &proc_rows, &proc_cols);
    
    int my_row_coord, my_col_coord;
    rank_to_coords_2d(rank, proc_rows, proc_cols, &my_row_coord, &my_col_coord);
    
    if (rank == 0) {
        printf("2D Process Grid (Cyclic): %d x %d\n", proc_rows, proc_cols);
    }
    
    // Step 2: Count rows owned by this process
    // In cyclic distribution: row i belongs to process (i % proc_rows)
    unsigned local_rows = 0;
    for (unsigned i = 0; i < matrix->n_rows; i++) {
        if ((int)(i % proc_rows) == my_row_coord) {
            local_rows++;
        }
    }
    
    local_matrix->n_rows = local_rows;
    local_matrix->n_cols = matrix->n_cols;
    
    // Step 3: Count local non-zeros
    local_matrix->nnz = 0;
    for (unsigned i = 0; i < matrix->nnz; i++) {
        unsigned global_row = matrix->row_indices[i];
        unsigned global_col = matrix->col_indices[i];
        
        int element_owner = owner_2d_cyclic(global_row, global_col, proc_rows, proc_cols);
        
        if (element_owner == rank) {
            local_matrix->nnz++;
        }
    }
    
    // Step 4: Allocate
    if (local_matrix->nnz > 0) {
        local_matrix->row_indices = (unsigned*)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->col_indices = (unsigned*)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->values = (double*)surely_malloc(local_matrix->nnz * sizeof(double));
    } else {
        local_matrix->row_indices = NULL;
        local_matrix->col_indices = NULL;
        local_matrix->values = NULL;
    }
    
    // Step 5: Fill arrays with local-to-global mapping
    unsigned local_idx = 0;
    for (unsigned i = 0; i < matrix->nnz; i++) {
        unsigned global_row = matrix->row_indices[i];
        unsigned global_col = matrix->col_indices[i];
        
        int element_owner = owner_2d_cyclic(global_row, global_col, proc_rows, proc_cols);
        
        if (element_owner == rank) {
            // Map global row to local row index
            // Local row index = global_row / proc_rows
            local_matrix->row_indices[local_idx] = global_row / proc_rows;
            local_matrix->col_indices[local_idx] = global_col;
            local_matrix->values[local_idx] = matrix->values[i];
            local_idx++;
        }
    }
    
    // Step 6: Convert to CSR
    if (local_matrix->nnz > 0) {
        Sparse_CSR *temp_csr = coo_to_csr_matrix(local_matrix);
        *local_csr = *temp_csr;
        free(temp_csr);
    } else {
        local_csr->n_rows = local_matrix->n_rows;
        local_csr->n_cols = local_matrix->n_cols;
        local_csr->row_ptr = (unsigned*)calloc(local_matrix->n_rows + 1, sizeof(unsigned));
        local_csr->col_ind = NULL;
        local_csr->values = NULL;
    }
    
    printf("Rank %d: Grid pos (%d,%d) | Local rows: %u | Local NNZ: %u\n",
           rank, my_row_coord, my_col_coord, local_rows, local_matrix->nnz);
}

