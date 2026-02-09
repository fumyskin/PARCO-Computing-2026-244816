#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include "mmio.h"
#include "specifications.h"

typedef struct Comm_Pattern{
    // ghost data
    int num_ghost_cols;           // ghost columns count
    int *ghost_col_indices;       // ghost global indices
    double *ghost_values;         
    int *ghost_to_local;          // global->local mapping (size = total_cols)
    
    int *send_counts;             
    int *recv_counts;             
    int *send_displs;             
    int *recv_displs;             

    int total_cols;
    int *send_indices;      
    int total_to_send;
} Comm_Pattern;


void find_ghost_vals(Sparse_CSR *local_csr,
                     Comm_Pattern *cp,
                     int rank,
                     int comm_size,
                     unsigned matrix_cols);

void exchange_ghost_vals(double *local_vec,
                         Comm_Pattern *cp,
                         int rank,
                         int comm_size);
// void find_ghost_vals(Sparse_CSR *local_matrix, Comm_Pattern *comm_pattern, int rank, int comm_size, int total_cols);
// void exchange_ghost_vals(double *local_x, Comm_Pattern *comm_pattern, int rank, int comm_size);
void local_SpMV(Sparse_CSR *local_csr, double *local_vec, double *local_result, Comm_Pattern *comm_pattern, int rank, int comm_size);


// ============================================================================
// 2D Communication Pattern
//
// In 2D SpMV (y = A*x), the matrix is distributed on a pr x pc process grid.
// Each process (r,c) owns a block of rows and a block of columns.
//
// The SpMV has two communication phases:
//   EXPAND:  Gather x-vector entries needed for local multiply.
//            Process (r,c) needs x[j] for every column j in its local nonzeros.
//            x[j] is owned by the process responsible for that column block.
//   FOLD:    After local multiply, partial y-results must be sent to the
//            process that owns the corresponding row of y.
//            Process (r,c) produces partial y[i] for its local rows, but
//            y[i] is owned by the process responsible for that row block.
//
// For BLOCK 2D distribution:
//   - Rows [r*rows_per_proc .. (r+1)*rows_per_proc) belong to proc-row r
//   - Cols [c*cols_per_proc .. (c+1)*cols_per_proc) belong to proc-col c
//   - x[j] is owned by process (j / cols_per_proc, j / cols_per_proc) -- NO!
//     We need a 1D vector owner. Convention: x[j] is owned by the process
//     on the DIAGONAL or by a designated "vector owner" row.
//
// SIMPLER APPROACH (used here):
//   - x[j] is owned by process-column (j / cols_per_proc), specifically
//     by the process whose proc-row == proc-col (diagonal), or row 0.
//     We use: x_owner(j) = coords_to_rank(0, j's_proc_col, pc)
//     i.e., row 0 of each processor column owns the x entries.
//   - y[i] is owned by process-row (i / rows_per_proc), specifically
//     by the process whose proc-col == 0 (or the diagonal).
//     We use: y_owner(i) = coords_to_rank(i's_proc_row, 0, pc)
//     i.e., column 0 of each processor row owns the y entries.
//
// EVEN SIMPLER (what we actually implement):
//   Each process in a processor column shares the SAME x-segment after expand.
//   Each process in a processor row contributes partial y that gets summed.
//
// We implement expand as communication WITHIN processor columns,
// and fold as communication WITHIN processor rows, using MPI subcommunicators.
// ============================================================================

typedef struct {
    // --- Grid info ---
    int proc_rows;          // number of process rows (pr)
    int proc_cols;          // number of process columns (pc)
    int my_row_coord;       // this process's row in the grid
    int my_col_coord;       // this process's column in the grid

    // --- Matrix dimensions (global) ---
    unsigned global_rows;
    unsigned global_cols;

    // --- Block ranges for this process ---
    unsigned my_row_start;  // first global row owned
    unsigned my_row_end;    // one past last global row owned
    unsigned my_col_start;  // first global col in this proc-column's block
    unsigned my_col_end;    // one past last global col

    // --- Block sizes (for mapping global -> proc coordinate) ---
    unsigned rows_per_proc; // ceil(global_rows / pr)
    unsigned cols_per_proc; // ceil(global_cols / pc)

    // --- EXPAND phase (gather x values) ---
    // We need x[j] for all columns j referenced by local nonzeros.
    // Columns in [my_col_start, my_col_end) are "local" to our proc-column.
    // Columns outside that range are "remote" and must be fetched.
    //
    // expand_ghost_count: number of unique remote columns needed
    // expand_ghost_global_cols: the global column indices of those ghosts
    // expand_ghost_values: buffer to receive ghost x-values after expand
    // expand_ghost_to_local: maps global_col -> index into expand_ghost_values
    //                        -1 if the column is local to our proc-column
    //
    // Communication is with ALL processes (Alltoallv), since a remote column
    // could be owned by any process. We identify the vector owner of column j
    // as the process at grid position (my_row_coord, j / cols_per_proc).
    // NOTE: this means each proc-row has its OWN copy of x for its column block.
    // The "vector owner" for column j from the perspective of proc (r,c) is
    // process (r, col_block(j)) -- i.e., same proc-row, different proc-col.
    // This way, expand is communication WITHIN processor rows.
    
    int expand_ghost_count;
    int *expand_ghost_global_cols;  // [expand_ghost_count] global col indices
    double *expand_ghost_values;    // [expand_ghost_count] received x values
    int *expand_col_to_ghost;       // [global_cols] -> ghost index or -1

    // Alltoallv arrays for expand (comm within proc row, or global)
    int *expand_recv_counts;   // [comm_size] how many x values we receive from each proc
    int *expand_recv_displs;   // [comm_size]
    int *expand_send_counts;   // [comm_size] how many x values we send to each proc
    int *expand_send_displs;   // [comm_size]
    int expand_total_to_send;
    int *expand_send_indices;  // [expand_total_to_send] local x indices to pack

    // --- FOLD phase (reduce partial y values) ---
    // After local SpMV, process (r,c) has partial y[i] for rows in its local block.
    // The owner of y[i] is process (row_block(i), my_col_coord) -- same proc-col.
    // But wait: in our 2D block distribution, process (r,c) only has nonzeros for
    // rows in proc-row r's range. So the partial y[i] for row i (which is in our
    // row range) needs to go to the process that "accumulates" y[i].
    //
    // Convention: y[i] is owned by process at (row_block(i), 0), i.e., first
    // proc-col in that proc-row. All processes in the same proc-row that have
    // partial contributions to the same y[i] must send them there.
    //
    // Actually, in BLOCK 2D: process (r,c) owns rows [r*rpp, (r+1)*rpp).
    // ALL processes in proc-row r have the SAME row range. Each produces a
    // partial y for those rows. They need to be summed (reduced).
    // => fold = MPI_Reduce (or Allreduce) within the processor row.
    //
    // We use an MPI sub-communicator for each processor row.
    
    MPI_Comm row_comm;     // communicator for processes in the same proc-row
    MPI_Comm col_comm;     // communicator for processes in the same proc-col
    int row_comm_rank;     // my rank within row_comm
    int row_comm_size;     // size of row_comm (== proc_cols)
    int col_comm_rank;     // my rank within col_comm
    int col_comm_size;     // size of col_comm (== proc_rows)

    // local result size (number of local rows)
    unsigned local_nrows;
    
} Comm_Pattern_2D;


// Set up the 2D communication pattern.
// Call AFTER distribution_2D has been applied and local_csr has global column indices.
void setup_comm_pattern_2d(
    Sparse_CSR *local_csr,
    Comm_Pattern_2D *cp,
    int rank, int comm_size,
    unsigned global_rows, unsigned global_cols
);

// EXPAND: gather remote x-values needed for local SpMV.
// local_x: this process's portion of the x vector (entries for my column block).
// After this call, cp->expand_ghost_values is filled.
void expand_x_2d(
    double *local_x,
    Comm_Pattern_2D *cp,
    int rank, int comm_size
);

// Local SpMV: y_partial = A_local * x  (using local + ghost x values)
void local_spmv_2d(
    Sparse_CSR *local_csr,
    double *local_x,
    double *y_partial,
    Comm_Pattern_2D *cp,
    int rank, int comm_size
);

// FOLD: reduce partial y across processes in the same processor row.
// y_partial (input):  this process's partial result for local rows
// y_result  (output): the summed result (valid on the "owner" process, i.e., rank 0 in row_comm)
// If you want all processes in the row to have the result, use Allreduce variant.
void fold_y_2d(
    double *y_partial,
    double *y_result,
    Comm_Pattern_2D *cp
);

// Combined: expand + spmv + fold
void spmv_2d(
    Sparse_CSR *local_csr,
    double *local_x,
    double *y_partial,
    double *y_result,
    Comm_Pattern_2D *cp,
    int rank, int comm_size
);

// Cleanup
void free_comm_pattern_2d(Comm_Pattern_2D *cp);




#endif