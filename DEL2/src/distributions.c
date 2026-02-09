#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "mmio.h"
#include "specifications.h"
#include "distributions.h"

// ============================================================================
// 1D CYCLIC DISTRIBUTION 
// ============================================================================
void distribution_1D(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix,
                     Sparse_CSR *local_csr, int rank, int comm_size)
{
    int rows_per_process = matrix->n_rows / comm_size;
    if (rank < (matrix->n_rows % comm_size)) {
        rows_per_process++;
    }

    local_matrix->n_rows = rows_per_process;
    local_matrix->n_cols = matrix->n_cols;

    // count local nnz: owner(i) = i % P
    local_matrix->nnz = 0;
    for (unsigned i = 0; i < matrix->nnz; i++) {
        if (rank == (int)(matrix->row_indices[i] % comm_size)) {
            local_matrix->nnz++;
        }
    }

    if (local_matrix->nnz > 0) {
        local_matrix->row_indices = (unsigned *)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->col_indices = (unsigned *)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->values      = (double *)surely_malloc(local_matrix->nnz * sizeof(double));
    } else {
        local_matrix->row_indices = NULL;
        local_matrix->col_indices = NULL;
        local_matrix->values      = NULL;
    }

    int temp_i = 0;
    for (int i = 0; i < (int)matrix->nnz; i++) {
        int assigned_row = matrix->row_indices[i];
        int owner_row = assigned_row % comm_size;
        if (owner_row == rank) {
            local_matrix->row_indices[temp_i] = assigned_row / comm_size;  // local row
            local_matrix->col_indices[temp_i] = matrix->col_indices[i];    // global col (1D uses global cols)
            local_matrix->values[temp_i]      = matrix->values[i];
            temp_i++;
        }
    }

    if (local_matrix->nnz != 0) {
        Sparse_CSR *temp_csr = coo_to_csr_matrix(local_matrix);
        *local_csr = *temp_csr;
        free(temp_csr);
    } else {
        local_csr->n_rows  = local_matrix->n_rows;
        local_csr->n_cols  = local_matrix->n_cols;
        local_csr->row_ptr = (unsigned *)calloc(local_matrix->n_rows + 1, sizeof(unsigned));
        local_csr->col_ind = NULL;
        local_csr->values  = NULL;
    }
}

// ============================================================================
// DUMMY MATRIX GENERATION (1D cyclic, unchanged)
// ============================================================================
void generate_dummy_matrix(Sparse_Coordinate *local_matrix, unsigned rows_per_process,
                           unsigned nnz_per_row, int rank, int comm_size)
{
    unsigned long long global_matrix_rows = (unsigned long long)rows_per_process * comm_size;
    unsigned long long global_matrix_cols = global_matrix_rows;
    unsigned long long local_nnz = (unsigned long long)rows_per_process * nnz_per_row;

    local_matrix->row_indices = (unsigned *)surely_malloc(local_nnz * sizeof(unsigned));
    local_matrix->col_indices = (unsigned *)surely_malloc(local_nnz * sizeof(unsigned));
    local_matrix->values      = (double *)surely_malloc(local_nnz * sizeof(double));

    local_matrix->n_rows = global_matrix_rows;
    local_matrix->n_cols = global_matrix_cols;
    local_matrix->nnz    = local_nnz;

    srand(rank * 12345 + 1);
    unsigned current_idx = 0;

    for (unsigned k = 0; k < rows_per_process; k++) {
        unsigned global_row_idx = (k * comm_size) + rank;
        for (unsigned n = 0; n < nnz_per_row; n++) {
            unsigned rand_col = rand() % global_matrix_cols;
            local_matrix->row_indices[current_idx] = global_row_idx;
            local_matrix->col_indices[current_idx] = rand_col;
            local_matrix->values[current_idx]      = 1.0;
            current_idx++;
        }
    }
}


// ============================================================================
// 2D GRID HELPERS
// ============================================================================

void create_2d_grid(int comm_size, int *proc_rows, int *proc_cols)
{
    int sqrt_p = (int)sqrt((double)comm_size);
    for (int i = sqrt_p; i >= 1; i--) {
        if (comm_size % i == 0) {
            *proc_rows = i;
            *proc_cols = comm_size / i;
            return;
        }
    }
    *proc_rows = 1;
    *proc_cols = comm_size;
}

void rank_to_coords_2d(int rank, int proc_rows, int proc_cols, int *row_coord, int *col_coord)
{
    *row_coord = rank / proc_cols;
    *col_coord = rank % proc_cols;
}

int coords_to_rank_2d(int row_coord, int col_coord, int proc_cols)
{
    return row_coord * proc_cols + col_coord;
}


// ============================================================================
// 2D BLOCK DISTRIBUTION
//
// KEY CHANGE: Column indices remain GLOBAL in the local CSR.
// Only row indices are remapped to local (0-based).
// This is essential for the 2D communication pattern to work.
// ============================================================================

static int owner_2d_block(unsigned row, unsigned col,
                          unsigned rows_per_proc, unsigned cols_per_proc,
                          int proc_rows, int proc_cols)
{
    int owner_row = (int)(row / rows_per_proc);
    int owner_col = (int)(col / cols_per_proc);
    if (owner_row >= proc_rows) owner_row = proc_rows - 1;
    if (owner_col >= proc_cols) owner_col = proc_cols - 1;
    return coords_to_rank_2d(owner_row, owner_col, proc_cols);
}

void distribution_2D(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix,
                     Sparse_CSR *local_csr, int rank, int comm_size)
{
    // 1) Create 2D process grid
    int proc_rows, proc_cols;
    create_2d_grid(comm_size, &proc_rows, &proc_cols);

    int my_row_coord, my_col_coord;
    rank_to_coords_2d(rank, proc_rows, proc_cols, &my_row_coord, &my_col_coord);

    if (rank == 0) {
        printf("2D Process Grid (Block): %d x %d (rows x cols)\n", proc_rows, proc_cols);
    }

    // 2) Determine block ranges
    unsigned rows_per_proc = (matrix->n_rows + proc_rows - 1) / proc_rows;
    unsigned cols_per_proc = (matrix->n_cols + proc_cols - 1) / proc_cols;

    unsigned my_start_row = my_row_coord * rows_per_proc;
    unsigned my_end_row   = my_start_row + rows_per_proc;
    if (my_end_row > matrix->n_rows) my_end_row = matrix->n_rows;

    unsigned my_start_col = my_col_coord * cols_per_proc;
    unsigned my_end_col   = my_start_col + cols_per_proc;
    if (my_end_col > matrix->n_cols) my_end_col = matrix->n_cols;

    // Local dimensions
    // n_rows = number of rows in my block (for CSR structure)
    // n_cols = GLOBAL number of columns (since col indices stay global)
    local_matrix->n_rows = my_end_row - my_start_row;
    local_matrix->n_cols = matrix->n_cols;  // *** KEEP GLOBAL ***

    // 3) Count local non-zeros
    local_matrix->nnz = 0;
    for (unsigned i = 0; i < matrix->nnz; i++) {
        int element_owner = owner_2d_block(matrix->row_indices[i], matrix->col_indices[i],
                                           rows_per_proc, cols_per_proc,
                                           proc_rows, proc_cols);
        if (element_owner == rank) {
            local_matrix->nnz++;
        }
    }

    // 4) Allocate
    if (local_matrix->nnz > 0) {
        local_matrix->row_indices = (unsigned *)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->col_indices = (unsigned *)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->values      = (double *)surely_malloc(local_matrix->nnz * sizeof(double));
    } else {
        local_matrix->row_indices = NULL;
        local_matrix->col_indices = NULL;
        local_matrix->values      = NULL;
    }

    // 5) Fill local arrays
    unsigned local_idx = 0;
    for (unsigned i = 0; i < matrix->nnz; i++) {
        unsigned global_row = matrix->row_indices[i];
        unsigned global_col = matrix->col_indices[i];

        int element_owner = owner_2d_block(global_row, global_col,
                                           rows_per_proc, cols_per_proc,
                                           proc_rows, proc_cols);
        if (element_owner == rank) {
            // Row: remap to local (0-based within this block)
            local_matrix->row_indices[local_idx] = global_row - my_start_row;
            // Column: KEEP GLOBAL — this is the critical change!
            local_matrix->col_indices[local_idx] = global_col;
            local_matrix->values[local_idx]      = matrix->values[i];
            local_idx++;
        }
    }

    // 6) Convert to CSR
    if (local_matrix->nnz > 0) {
        Sparse_CSR *temp_csr = coo_to_csr_matrix(local_matrix);
        *local_csr = *temp_csr;
        free(temp_csr);
    } else {
        local_csr->n_rows  = local_matrix->n_rows;
        local_csr->n_cols  = local_matrix->n_cols;
        local_csr->row_ptr = (unsigned *)calloc(local_matrix->n_rows + 1, sizeof(unsigned));
        local_csr->col_ind = NULL;
        local_csr->values  = NULL;
    }

    printf("Rank %d: Grid (%d,%d) | Rows [%u,%u) | Cols [%u,%u) | LocalRows=%u | NNZ=%u | ColsGlobal\n",
           rank, my_row_coord, my_col_coord, my_start_row, my_end_row,
           my_start_col, my_end_col, local_matrix->n_rows, local_matrix->nnz);
}


// ============================================================================
// 2D CYCLIC DISTRIBUTION
//
// Same fix: column indices remain GLOBAL.
// ============================================================================

static int owner_2d_cyclic_func(unsigned row, unsigned col, int proc_rows, int proc_cols)
{
    int owner_row = (int)(row % proc_rows);
    int owner_col = (int)(col % proc_cols);
    return coords_to_rank_2d(owner_row, owner_col, proc_cols);
}

void distribution_2D_cyclic(Sparse_Coordinate *matrix, Sparse_Coordinate *local_matrix,
                            Sparse_CSR *local_csr, int rank, int comm_size)
{
    int proc_rows, proc_cols;
    create_2d_grid(comm_size, &proc_rows, &proc_cols);

    int my_row_coord, my_col_coord;
    rank_to_coords_2d(rank, proc_rows, proc_cols, &my_row_coord, &my_col_coord);

    if (rank == 0) {
        printf("2D Process Grid (Cyclic): %d x %d\n", proc_rows, proc_cols);
    }

    // Count local rows (cyclic assignment)
    unsigned local_rows = 0;
    for (unsigned i = 0; i < matrix->n_rows; i++) {
        if ((int)(i % proc_rows) == my_row_coord) {
            local_rows++;
        }
    }

    local_matrix->n_rows = local_rows;
    local_matrix->n_cols = matrix->n_cols;  // *** KEEP GLOBAL ***

    // Count local nnz
    local_matrix->nnz = 0;
    for (unsigned i = 0; i < matrix->nnz; i++) {
        int element_owner = owner_2d_cyclic_func(matrix->row_indices[i], matrix->col_indices[i],
                                                  proc_rows, proc_cols);
        if (element_owner == rank) {
            local_matrix->nnz++;
        }
    }

    // Allocate
    if (local_matrix->nnz > 0) {
        local_matrix->row_indices = (unsigned *)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->col_indices = (unsigned *)surely_malloc(local_matrix->nnz * sizeof(unsigned));
        local_matrix->values      = (double *)surely_malloc(local_matrix->nnz * sizeof(double));
    } else {
        local_matrix->row_indices = NULL;
        local_matrix->col_indices = NULL;
        local_matrix->values      = NULL;
    }

    // Fill arrays: local row index, GLOBAL column index
    unsigned local_idx = 0;
    for (unsigned i = 0; i < matrix->nnz; i++) {
        unsigned global_row = matrix->row_indices[i];
        unsigned global_col = matrix->col_indices[i];

        int element_owner = owner_2d_cyclic_func(global_row, global_col, proc_rows, proc_cols);
        if (element_owner == rank) {
            // Row: cyclic local index
            local_matrix->row_indices[local_idx] = global_row / proc_rows;
            // Column: KEEP GLOBAL
            local_matrix->col_indices[local_idx] = global_col;
            local_matrix->values[local_idx]      = matrix->values[i];
            local_idx++;
        }
    }

    // Convert to CSR
    if (local_matrix->nnz > 0) {
        Sparse_CSR *temp_csr = coo_to_csr_matrix(local_matrix);
        *local_csr = *temp_csr;
        free(temp_csr);
    } else {
        local_csr->n_rows  = local_matrix->n_rows;
        local_csr->n_cols  = local_matrix->n_cols;
        local_csr->row_ptr = (unsigned *)calloc(local_matrix->n_rows + 1, sizeof(unsigned));
        local_csr->col_ind = NULL;
        local_csr->values  = NULL;
    }

    printf("Rank %d: Grid (%d,%d) | Local rows=%u | NNZ=%u | ColsGlobal\n",
           rank, my_row_coord, my_col_coord, local_rows, local_matrix->nnz);
}


// ============================================================================
// DUMMY MATRIX GENERATORS FOR 2D (updated to keep global columns)
// ============================================================================

void generate_dummy_matrix_2d_block(Sparse_Coordinate *local_matrix, unsigned rows_per_proc,
                                    unsigned nnz_per_row, int rank, int comm_size)
{
    int proc_rows, proc_cols;
    create_2d_grid(comm_size, &proc_rows, &proc_cols);

    int my_row_coord, my_col_coord;
    rank_to_coords_2d(rank, proc_rows, proc_cols, &my_row_coord, &my_col_coord);

    unsigned total_rows = rows_per_proc * proc_rows;
    unsigned total_cols = total_rows;

    local_matrix->n_rows = rows_per_proc;
    local_matrix->n_cols = total_cols;  // global columns

    unsigned estimated_nnz = rows_per_proc * nnz_per_row;

    unsigned *temp_rows = (unsigned *)surely_malloc(estimated_nnz * 2 * sizeof(unsigned));
    unsigned *temp_cols = (unsigned *)surely_malloc(estimated_nnz * 2 * sizeof(unsigned));
    double   *temp_vals = (double *)surely_malloc(estimated_nnz * 2 * sizeof(double));

    unsigned actual_nnz = 0;

    for (unsigned local_row = 0; local_row < rows_per_proc; local_row++) {
        unsigned global_row = my_row_coord * rows_per_proc + local_row;
        if (global_row >= total_rows) continue;

        if (global_row < total_cols) {
            temp_rows[actual_nnz] = local_row;
            temp_cols[actual_nnz] = global_row;  // global column (diagonal)
            temp_vals[actual_nnz] = 4.0;
            actual_nnz++;
        }

        for (unsigned k = 1; k < nnz_per_row && actual_nnz < estimated_nnz * 2; k++) {
            unsigned global_col = rand() % total_cols;
            if (global_col != global_row) {
                temp_rows[actual_nnz] = local_row;
                temp_cols[actual_nnz] = global_col;  // global column
                temp_vals[actual_nnz] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
                actual_nnz++;
            }
        }
    }

    local_matrix->nnz = actual_nnz;
    local_matrix->row_indices = (unsigned *)surely_malloc(actual_nnz * sizeof(unsigned));
    local_matrix->col_indices = (unsigned *)surely_malloc(actual_nnz * sizeof(unsigned));
    local_matrix->values      = (double *)surely_malloc(actual_nnz * sizeof(double));

    for (unsigned i = 0; i < actual_nnz; i++) {
        local_matrix->row_indices[i] = temp_rows[i];
        local_matrix->col_indices[i] = temp_cols[i];
        local_matrix->values[i]      = temp_vals[i];
    }

    free(temp_rows);
    free(temp_cols);
    free(temp_vals);
}

void generate_dummy_matrix_2d_cyclic(Sparse_Coordinate *local_matrix, unsigned rows_per_proc,
                                     unsigned nnz_per_row, int rank, int comm_size)
{
    int proc_rows, proc_cols;
    create_2d_grid(comm_size, &proc_rows, &proc_cols);

    int my_row_coord, my_col_coord;
    rank_to_coords_2d(rank, proc_rows, proc_cols, &my_row_coord, &my_col_coord);

    unsigned total_rows = rows_per_proc * proc_rows;
    unsigned total_cols = total_rows;

    local_matrix->n_rows = rows_per_proc;
    local_matrix->n_cols = total_cols;  // global columns

    unsigned estimated_nnz = rows_per_proc * nnz_per_row;

    unsigned *temp_rows = (unsigned *)surely_malloc(estimated_nnz * 2 * sizeof(unsigned));
    unsigned *temp_cols = (unsigned *)surely_malloc(estimated_nnz * 2 * sizeof(unsigned));
    double   *temp_vals = (double *)surely_malloc(estimated_nnz * 2 * sizeof(double));

    unsigned actual_nnz = 0;

    for (unsigned local_row = 0; local_row < rows_per_proc; local_row++) {
        unsigned global_row = my_row_coord + local_row * proc_rows;
        if (global_row >= total_rows) continue;

        if (global_row < total_cols) {
            temp_rows[actual_nnz] = local_row;
            temp_cols[actual_nnz] = global_row;  // global
            temp_vals[actual_nnz] = 4.0;
            actual_nnz++;
        }

        for (unsigned k = 1; k < nnz_per_row && actual_nnz < estimated_nnz * 2; k++) {
            unsigned global_col = rand() % total_cols;
            if (global_col != global_row) {
                temp_rows[actual_nnz] = local_row;
                temp_cols[actual_nnz] = global_col;  // global
                temp_vals[actual_nnz] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
                actual_nnz++;
            }
        }
    }

    local_matrix->nnz = actual_nnz;
    local_matrix->row_indices = (unsigned *)surely_malloc(actual_nnz * sizeof(unsigned));
    local_matrix->col_indices = (unsigned *)surely_malloc(actual_nnz * sizeof(unsigned));
    local_matrix->values      = (double *)surely_malloc(actual_nnz * sizeof(double));

    for (unsigned i = 0; i < actual_nnz; i++) {
        local_matrix->row_indices[i] = temp_rows[i];
        local_matrix->col_indices[i] = temp_cols[i];
        local_matrix->values[i]      = temp_vals[i];
    }

    free(temp_rows);
    free(temp_cols);
    free(temp_vals);
}