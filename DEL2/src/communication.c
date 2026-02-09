#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <stdbool.h>
#include "specifications.h"
#include "communication.h"
#include "distributions.h"

#ifdef HYBRID
#include <omp.h>
#endif
// ---------------------------------------------------------------------------
// find_ghost_vals
// EXPLANATION:
// Cyclic vector ownership: global column j is owned by process (j % comm_size)
// at local vector index (j / comm_size).
//
// Steps:
//   1. Scan local_csr, collect unique ghost columns (owner != rank).
//   2. Sort ghosts by owner — this gives the contiguous layout that
//      recv_displs expects (all ghosts from rank 0, then rank 1, ...).
//   3. Build recv_counts / recv_displs from the sorted order.
//   4. Fill ghost_to_local: dense map, global_col -> index in ghost_values.
//      -1 means "not a ghost" (i.e. local).
//   5. Alltoall recv_counts -> we learn send_counts (how many values each
//      other process wants FROM us). Build send_displs, total_to_send.
//   6. Alltoallv to exchange the actual remote indices each process needs.
//      We send our ghost_remote_idx (what WE need, grouped by owner).
//      We receive into send_indices (what OTHERS need from us).
// ---------------------------------------------------------------------------
void find_ghost_vals(Sparse_CSR *local_csr, Comm_Pattern *cp, int rank, int comm_size, unsigned matrix_cols){
    cp->total_cols = (int)matrix_cols;
    unsigned nnz = local_csr->row_ptr[local_csr->n_rows];

    // collect unique ghost columns, use boolean array
    // note: for large matrices an hash-set could be even more handy for speed
    // assumption: all matrices in the set fit all in memory
    bool *seen = (bool*)calloc(matrix_cols, sizeof(bool));

    int ghost_count = 0;
    for (unsigned i = 0; i < nnz; i++) {
        unsigned col = local_csr->col_ind[i];
        int owner = (int)(col % comm_size);
        if (owner != rank && !seen[col]) {
            seen[col] = true;
            ghost_count++;
        }
    }

    cp->num_ghost_cols = ghost_count;

    // temporary flat list of ghost columns (unsorted)
    int *raw_ghosts = (int *)surely_malloc(ghost_count * sizeof(int));
    int idx = 0;
    for (unsigned col = 0; col < matrix_cols; col++) {
        if (seen[col]) {
            raw_ghosts[idx++] = (int)col;
        }
    }
    free(seen);

    // CRITICAL ASSUMPTION: COLUMNS MUST BE SORTED
    // sort ghost columns by owner rank 
    // final result: [all ghosts owned by 0 | owned by 1 | ... | owned by comm_size-1]
    cp->recv_counts = (int *)calloc(comm_size, sizeof(int));
    for (int i = 0; i < ghost_count; i++)
        cp->recv_counts[raw_ghosts[i] % comm_size]++;

    // prefix sum -> recv_displs
    cp->recv_displs = (int *)surely_malloc(comm_size * sizeof(int));
    cp->recv_displs[0] = 0;
    for (int p = 1; p < comm_size; p++)
        cp->recv_displs[p] = cp->recv_displs[p-1] + cp->recv_counts[p-1];

    // scatter into sorted order using a running cursor per owner
    cp->ghost_col_indices = (int *)surely_malloc(ghost_count * sizeof(int));
    int *cursor = (int *)calloc(comm_size, sizeof(int));

    for (int i = 0; i < ghost_count; i++) {
        int col   = raw_ghosts[i];
        int owner = col % comm_size;
        int pos   = cp->recv_displs[owner] + cursor[owner];
        cp->ghost_col_indices[pos] = col;
        cursor[owner]++;
    }
    free(raw_ghosts);
    free(cursor);

    // allocate the receive buffer for ghost vector values
    cp->ghost_values = (double *)surely_malloc(ghost_count * sizeof(double));

    // ghost_to_local: dense reverse map, global_col -> ghost index
    cp->ghost_to_local = (int *)surely_malloc(matrix_cols * sizeof(int));
    memset(cp->ghost_to_local, 0xFF, matrix_cols * sizeof(int));   // fill with -1 to indicate column is owned locally!

    for (int i = 0; i < ghost_count; i++)
        cp->ghost_to_local[cp->ghost_col_indices[i]] = i;

    // Alltoall: learn send_counts
    // Every process shares its recv_counts row.  Process p reads
    // column `rank` to find out how many values each p wants from us.
    // Equivalently: MPI_Alltoall transposes the matrix of counts.
    // recv_counts[p]  = "I want this many from p"
    // After Alltoall the result in send_counts is:
    // send_counts[p]  = "p wants this many from me"
    cp->send_counts = (int *)surely_malloc(comm_size * sizeof(int));
    MPI_Alltoall(cp->recv_counts, 1, MPI_INT,
                 cp->send_counts, 1, MPI_INT,
                 MPI_COMM_WORLD);

    // prefix sum -> send_displs, total_to_send
    cp->send_displs = (int *)surely_malloc(comm_size * sizeof(int));
    cp->send_displs[0] = 0;
    for (int p = 1; p < comm_size; p++)
        cp->send_displs[p] = cp->send_displs[p-1] + cp->send_counts[p-1];

    cp->total_to_send = cp->send_displs[comm_size - 1] + cp->send_counts[comm_size - 1];

    // Alltoallv: exchange remote indices:
    //     Each process SENDS its ghost_remote_idx — i.e., for each ghost
    //     column it needs, the local vector index on the owning process.
    //     ghost_remote_idx for ghost_col_indices[i] = ghost_col_indices[i] / comm_size
    //     These are already grouped by owner (= destination) thanks to the
    //     sort in step 2, so send displacement = recv_displs, send count = recv_counts.
    //     Each process RECEIVES into send_indices — the indices that OTHER
    //     processes need from our local vector.  recv displacement = send_displs,
    //     recv count = send_counts.
    // build the flat send buffer: remote index for each ghost, in sorted order
    int *ghost_remote_idx = (int *)surely_malloc(ghost_count * sizeof(int));
    for (int i = 0; i < ghost_count; i++)
        ghost_remote_idx[i] = cp->ghost_col_indices[i] / comm_size;

    cp->send_indices = (int *)surely_malloc(cp->total_to_send * sizeof(int));

    MPI_Alltoallv(ghost_remote_idx,      
                    cp->recv_counts, 
                    cp->recv_displs, 
                    MPI_INT,
                    cp->send_indices,       
                    cp->send_counts, 
                    cp->send_displs, 
                    MPI_INT,
                    MPI_COMM_WORLD);
    // send_indices now contains, for each value we must send, the local
    // index into local_rand_vec to read from.
    free(ghost_remote_idx);
}


// ---------------------------------------------------------------------------
// exchange_ghost_vals
//
// Pack -> Alltoallv -> done.
//
// send_buf[send_displs[p] .. +send_counts[p]-1] are the values process p
// requested from us (indices given by send_indices).
// They arrive in ghost_values[recv_displs[p] .. +recv_counts[p]-1].
// ---------------------------------------------------------------------------
void exchange_ghost_vals(double *local_vec, Comm_Pattern *cp, int rank, int comm_size){
    // Pack outgoing values
    double *send_buf = (double *)surely_malloc(cp->total_to_send * sizeof(double));
    #ifdef HYBRID
        #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < cp->total_to_send; i++)
        send_buf[i] = local_vec[cp->send_indices[i]];

    // single collective: sends and receives everything in one call
    MPI_Alltoallv(send_buf, 
                    cp->send_counts, 
                    cp->send_displs, 
                    MPI_DOUBLE,
                    cp->ghost_values,    
                    cp->recv_counts, 
                    cp->recv_displs, 
                    MPI_DOUBLE,
                    MPI_COMM_WORLD);

    free(send_buf);
}

void local_SpMV(Sparse_CSR *local_csr, double *local_vec, double *local_result, Comm_Pattern *comm_pattern, int rank, int comm_size) {
    
    
    #ifdef HYBRID
        #pragma omp parallel for schedule(dynamic, 64)
    #endif
    for (unsigned i = 0; i < local_csr->n_rows; i++) {
        local_result[i] = 0.0;
    }
    
    // perform SpMV row by row
    unsigned i, j;
    for (i = 0; i < local_csr->n_rows; i++) {
        double sum = 0.0;
        unsigned row_start = local_csr->row_ptr[i];
        unsigned row_end = local_csr->row_ptr[i + 1];
        
        for (j = row_start; j < row_end; j++) {
            unsigned col = local_csr->col_ind[j];
            double val = local_csr->values[j];
            
            // determine where to get the vector value
            int col_owner = col % comm_size;
            
            if (col_owner == rank) {
                // local column - use local_vec
                unsigned local_col_idx = col / comm_size;
                sum += val * local_vec[local_col_idx];
            } else {
                // ghost column - use ghost_values
                int ghost_idx = comm_pattern->ghost_to_local[col];
                if (ghost_idx >= 0) {
                    sum += val * comm_pattern->ghost_values[ghost_idx];
                } else {
                    fprintf(stderr, "[Rank %d] ERROR: Column %d not found in ghost mapping!\n", rank, col);
                }
            }
        }
        
        local_result[i] = sum;
    }
}


void setup_comm_pattern_2d(
    Sparse_CSR *local_csr,
    Comm_Pattern_2D *cp,
    int rank, int comm_size,
    unsigned global_rows, unsigned global_cols)
{
    // --- 1. Grid setup ---
    create_2d_grid(comm_size, &cp->proc_rows, &cp->proc_cols);
    rank_to_coords_2d(rank, cp->proc_rows, cp->proc_cols,
                      &cp->my_row_coord, &cp->my_col_coord);

    cp->global_rows = global_rows;
    cp->global_cols = global_cols;

    cp->rows_per_proc = (global_rows + cp->proc_rows - 1) / cp->proc_rows;
    cp->cols_per_proc = (global_cols + cp->proc_cols - 1) / cp->proc_cols;

    // My block ranges
    cp->my_row_start = cp->my_row_coord * cp->rows_per_proc;
    cp->my_row_end   = cp->my_row_start + cp->rows_per_proc;
    if (cp->my_row_end > global_rows) cp->my_row_end = global_rows;

    cp->my_col_start = cp->my_col_coord * cp->cols_per_proc;
    cp->my_col_end   = cp->my_col_start + cp->cols_per_proc;
    if (cp->my_col_end > global_cols) cp->my_col_end = global_cols;

    cp->local_nrows = local_csr->n_rows;  // == my_row_end - my_row_start

    // --- 2. Create sub-communicators ---
    // row_comm: all processes in the same processor ROW (same my_row_coord)
    //           These share the same row range, communicate for expand + fold.
    // col_comm: all processes in the same processor COLUMN (same my_col_coord)
    //           Could be used for broadcasting x within a column if needed.
    MPI_Comm_split(MPI_COMM_WORLD, cp->my_row_coord, cp->my_col_coord, &cp->row_comm);
    MPI_Comm_rank(cp->row_comm, &cp->row_comm_rank);
    MPI_Comm_size(cp->row_comm, &cp->row_comm_size);

    MPI_Comm_split(MPI_COMM_WORLD, cp->my_col_coord, cp->my_row_coord, &cp->col_comm);
    MPI_Comm_rank(cp->col_comm, &cp->col_comm_rank);
    MPI_Comm_size(cp->col_comm, &cp->col_comm_size);

    // --- 3. Identify ghost columns for EXPAND phase ---
    // Scan all column indices in local_csr. If a column j is NOT in
    // [my_col_start, my_col_end), it's a "ghost" that we need from another process.
    // The owner of x[j] (within our proc-row) is the process at
    // grid position (my_row_coord, j / cols_per_proc).
    // In terms of row_comm rank, that's (j / cols_per_proc).

    unsigned nnz = local_csr->row_ptr[local_csr->n_rows];

    // Use a boolean array to find unique ghost columns
    bool *is_ghost = (bool *)calloc(global_cols, sizeof(bool));
    int ghost_count = 0;

    for (unsigned k = 0; k < nnz; k++) {
        unsigned col = local_csr->col_ind[k];
        if (col < cp->my_col_start || col >= cp->my_col_end) {
            if (!is_ghost[col]) {
                is_ghost[col] = true;
                ghost_count++;
            }
        }
    }

    cp->expand_ghost_count = ghost_count;

    // Build sorted ghost list grouped by owner (proc-col coordinate)
    // Owner in row_comm for column j is: j / cols_per_proc
    // (clamped to proc_cols - 1)
    cp->expand_recv_counts = (int *)calloc(cp->row_comm_size, sizeof(int));

    // First pass: count ghosts per owner
    for (unsigned col = 0; col < global_cols; col++) {
        if (is_ghost[col]) {
            int owner_col_coord = col / cp->cols_per_proc;
            if (owner_col_coord >= cp->proc_cols) owner_col_coord = cp->proc_cols - 1;
            cp->expand_recv_counts[owner_col_coord]++;
        }
    }

    // Prefix sum for recv_displs
    cp->expand_recv_displs = (int *)surely_malloc(cp->row_comm_size * sizeof(int));
    cp->expand_recv_displs[0] = 0;
    for (int p = 1; p < cp->row_comm_size; p++) {
        cp->expand_recv_displs[p] = cp->expand_recv_displs[p - 1] + cp->expand_recv_counts[p - 1];
    }

    // Build sorted ghost column list
    cp->expand_ghost_global_cols = (int *)surely_malloc((ghost_count > 0 ? ghost_count : 1) * sizeof(int));
    int *cursor = (int *)calloc(cp->row_comm_size, sizeof(int));

    for (unsigned col = 0; col < global_cols; col++) {
        if (is_ghost[col]) {
            int owner = col / cp->cols_per_proc;
            if (owner >= cp->proc_cols) owner = cp->proc_cols - 1;
            int pos = cp->expand_recv_displs[owner] + cursor[owner];
            cp->expand_ghost_global_cols[pos] = (int)col;
            cursor[owner]++;
        }
    }
    free(cursor);
    free(is_ghost);

    // Allocate ghost values buffer
    cp->expand_ghost_values = (double *)surely_malloc((ghost_count > 0 ? ghost_count : 1) * sizeof(double));

    // Build reverse map: global_col -> ghost index (-1 if local)
    cp->expand_col_to_ghost = (int *)surely_malloc(global_cols * sizeof(int));
    memset(cp->expand_col_to_ghost, 0xFF, global_cols * sizeof(int)); // fill with -1

    for (int i = 0; i < ghost_count; i++) {
        cp->expand_col_to_ghost[cp->expand_ghost_global_cols[i]] = i;
    }

    // --- 4. Exchange: learn what other processes need from us (within row_comm) ---
    cp->expand_send_counts = (int *)surely_malloc(cp->row_comm_size * sizeof(int));
    MPI_Alltoall(cp->expand_recv_counts, 1, MPI_INT,
                 cp->expand_send_counts, 1, MPI_INT,
                 cp->row_comm);

    cp->expand_send_displs = (int *)surely_malloc(cp->row_comm_size * sizeof(int));
    cp->expand_send_displs[0] = 0;
    for (int p = 1; p < cp->row_comm_size; p++) {
        cp->expand_send_displs[p] = cp->expand_send_displs[p - 1] + cp->expand_send_counts[p - 1];
    }
    cp->expand_total_to_send = cp->expand_send_displs[cp->row_comm_size - 1]
                             + cp->expand_send_counts[cp->row_comm_size - 1];

    // --- 5. Exchange the actual column indices so senders know what to pack ---
    // We send our ghost_global_cols (what we need), converted to LOCAL indices
    // on the owning process. Local index of global col j on owner process:
    //   j - owner_col_start, where owner_col_start = owner_col_coord * cols_per_proc
    // But since the ghost list is grouped by owner, and within each owner group
    // the columns all belong to that owner's block, we can convert:
    int *ghost_remote_local_idx = (int *)surely_malloc((ghost_count > 0 ? ghost_count : 1) * sizeof(int));
    for (int i = 0; i < ghost_count; i++) {
        int gcol = cp->expand_ghost_global_cols[i];
        int owner = gcol / cp->cols_per_proc;
        if (owner >= cp->proc_cols) owner = cp->proc_cols - 1;
        unsigned owner_col_start = owner * cp->cols_per_proc;
        ghost_remote_local_idx[i] = (int)(gcol - owner_col_start);
    }

    cp->expand_send_indices = (int *)surely_malloc((cp->expand_total_to_send > 0 ? cp->expand_total_to_send : 1) * sizeof(int));

    MPI_Alltoallv(ghost_remote_local_idx,
                  cp->expand_recv_counts, cp->expand_recv_displs, MPI_INT,
                  cp->expand_send_indices,
                  cp->expand_send_counts, cp->expand_send_displs, MPI_INT,
                  cp->row_comm);

    free(ghost_remote_local_idx);

    if (rank == 0) {
        printf("2D Comm Pattern: grid %d x %d | expand ghosts on rank 0: %d\n",
               cp->proc_rows, cp->proc_cols, ghost_count);
    }
}


// ============================================================================
// EXPAND: Gather remote x values needed for local SpMV
//
// local_x: the portion of x owned by this process, indexed [0 .. local_x_size)
//          where local_x_size = my_col_end - my_col_start
//
// Communication happens within row_comm (processes in the same processor row).
// ============================================================================
void expand_x_2d(double *local_x, Comm_Pattern_2D *cp, int rank, int comm_size)
{
    (void)rank; (void)comm_size;  // unused, we use row_comm

    // Pack outgoing values: other processes in my row_comm want certain
    // entries from my local_x.
    double *send_buf = (double *)surely_malloc(
        (cp->expand_total_to_send > 0 ? cp->expand_total_to_send : 1) * sizeof(double));

    #ifdef HYBRID
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < cp->expand_total_to_send; i++) {
        send_buf[i] = local_x[cp->expand_send_indices[i]];
    }

    // Alltoallv within row_comm
    MPI_Alltoallv(send_buf,
                  cp->expand_send_counts, cp->expand_send_displs, MPI_DOUBLE,
                  cp->expand_ghost_values,
                  cp->expand_recv_counts, cp->expand_recv_displs, MPI_DOUBLE,
                  cp->row_comm);

    free(send_buf);
}


// ============================================================================
// LOCAL SpMV for 2D distribution
//
// local_csr has GLOBAL column indices!!
// For each column j:
//   - If j is in [my_col_start, my_col_end), use local_x[j - my_col_start]
//   - Otherwise, use cp->expand_ghost_values[cp->expand_col_to_ghost[j]]
// ============================================================================
void local_spmv_2d(
    Sparse_CSR *local_csr,
    double *local_x,
    double *y_partial,
    Comm_Pattern_2D *cp,
    int rank, int comm_size)
{
    (void)rank; (void)comm_size;

    unsigned nrows = local_csr->n_rows;

    #ifdef HYBRID
    #pragma omp parallel for schedule(dynamic, 64)
    #endif
    for (unsigned i = 0; i < nrows; i++) {
        double sum = 0.0;
        unsigned row_start = local_csr->row_ptr[i];
        unsigned row_end   = local_csr->row_ptr[i + 1];

        for (unsigned k = row_start; k < row_end; k++) {
            unsigned col = local_csr->col_ind[k];
            double   val = local_csr->values[k];

            if (col >= cp->my_col_start && col < cp->my_col_end) {
                // Local column: direct access
                sum += val * local_x[col - cp->my_col_start];
            } else {
                // Ghost column: look up in expand buffer
                int ghost_idx = cp->expand_col_to_ghost[col];
                if (ghost_idx >= 0) {
                    sum += val * cp->expand_ghost_values[ghost_idx];
                } else {
                    fprintf(stderr, "[Rank %d] ERROR in local_spmv_2d: col %u not found "
                            "(not local [%u,%u) and not ghost)\n",
                            rank, col, cp->my_col_start, cp->my_col_end);
                }
            }
        }
        y_partial[i] = sum;
    }
}


// ============================================================================
// FOLD: Reduce partial y values across processes in the same processor row
//
// All processes in the same proc-row (r, *) have the same set of local rows
// (global rows [r*rpp, (r+1)*rpp)). Each has computed a partial y for those rows.
// We sum them up. The result is valid on ALL processes (Allreduce) or on
// one designated process (Reduce to root=0 in row_comm).
//
// We use MPI_Allreduce so every process in the row has the final y.
// (You could use MPI_Reduce to rank 0 of row_comm if only one needs it.)
// ============================================================================
void fold_y_2d(double *y_partial, double *y_result, Comm_Pattern_2D *cp){
    MPI_Allreduce(y_partial, y_result, (int)cp->local_nrows, MPI_DOUBLE, MPI_SUM, cp->row_comm);
}


// ============================================================================
// COMBINED: expand + local SpMV + fold -> i don't think I'll really use it
// ============================================================================
void spmv_2d(
    Sparse_CSR *local_csr,
    double *local_x,
    double *y_partial,
    double *y_result,
    Comm_Pattern_2D *cp,
    int rank, int comm_size)
{
    // Phase 1: Expand - gather needed x values
    expand_x_2d(local_x, cp, rank, comm_size);

    // Phase 2: Local compute
    local_spmv_2d(local_csr, local_x, y_partial, cp, rank, comm_size);

    // Phase 3+4: Fold - reduce partial y across processor row
    fold_y_2d(y_partial, y_result, cp);
}


// ============================================================================
// CLEANUP
// ============================================================================
void free_comm_pattern_2d(Comm_Pattern_2D *cp)
{
    if (cp->expand_ghost_global_cols) free(cp->expand_ghost_global_cols);
    if (cp->expand_ghost_values)      free(cp->expand_ghost_values);
    if (cp->expand_col_to_ghost)      free(cp->expand_col_to_ghost);
    if (cp->expand_recv_counts)       free(cp->expand_recv_counts);
    if (cp->expand_recv_displs)       free(cp->expand_recv_displs);
    if (cp->expand_send_counts)       free(cp->expand_send_counts);
    if (cp->expand_send_displs)       free(cp->expand_send_displs);
    if (cp->expand_send_indices)      free(cp->expand_send_indices);

    if (cp->row_comm != MPI_COMM_NULL) MPI_Comm_free(&cp->row_comm);
    if (cp->col_comm != MPI_COMM_NULL) MPI_Comm_free(&cp->col_comm);
}
