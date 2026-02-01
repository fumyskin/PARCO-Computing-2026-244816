#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <stdbool.h>
#include "specifications.h"
#include "communication.h"

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


// void find_ghost_vals(Sparse_CSR *local_matrix, Comm_Pattern *comm_pattern, int rank, int comm_size, int total_cols) {
    
//     unsigned i, j;
//     unsigned nnz = local_matrix->row_ptr[local_matrix->n_rows];
    
//     // count unique columns: form vector of vector elements long
//     char *seen = (char*)calloc(total_cols, sizeof(char));
//     if (!seen) {
//         fprintf(stderr, "[Rank %d] ERROR: Cannot allocate seen array\n", rank);
//         MPI_Abort(MPI_COMM_WORLD, 1);
//     }

//     // identify all columns referenced by LOCAL rows
//     for (j = 0; j < local_matrix->n_rows; j++) {
//         for (i = local_matrix->row_ptr[j]; i < local_matrix->row_ptr[j + 1]; i++) {
//             seen[local_matrix->col_ind[i]] = 1;
//         }
//     }
    
//     // communication arrays: 
//     comm_pattern->send_counts = (int*)calloc(comm_size, sizeof(int));
//     comm_pattern->recv_counts = (int*)calloc(comm_size, sizeof(int));
//     comm_pattern->send_displs = (int*)calloc(comm_size, sizeof(int));
//     comm_pattern->recv_displs = (int*)calloc(comm_size, sizeof(int));
    
//     if (!comm_pattern->send_counts || !comm_pattern->recv_counts || 
//         !comm_pattern->send_displs || !comm_pattern->recv_displs) {
//         fprintf(stderr, "[Rank %d] ERROR: Cannot allocate comm arrays\n", rank);
//         MPI_Abort(MPI_COMM_WORLD, 1);
//     }
    
//     // count ghost columns
//     int ghost_count = 0;
//     for (i = 0; i < total_cols; i++) {
//         if (!seen[i]) continue;

//         int owner = i % comm_size;
//         // identify ghosts (ie, NON LOCAL DEPENDENCIES)
//         if (owner != rank) {
//             comm_pattern->recv_counts[owner]++;
//             ghost_count++;
//         }
//     }
    
//     comm_pattern->num_ghost_cols = ghost_count;
//     comm_pattern->total_cols = total_cols;
    
//     // receive displacements -> for each i finds how many cells each recv occupy
//     comm_pattern->recv_displs[0] = 0;
//     for (i = 1; i < comm_size; i++) {
//         comm_pattern->recv_displs[i] = comm_pattern->recv_displs[i-1] + comm_pattern->recv_counts[i-1];
//     }
    
//     comm_pattern->ghost_col_indices = (int*)surely_malloc((ghost_count > 0 ? ghost_count : 1) * sizeof(int));
//     comm_pattern->ghost_values = (double*)surely_malloc((ghost_count > 0 ? ghost_count : 1) * sizeof(double));
//     comm_pattern->ghost_to_local = (int*)surely_malloc(total_cols * sizeof(int));
    
//     for (i = 0; i < total_cols; i++) {
//         comm_pattern->ghost_to_local[i] = -1;
//     }
    
//     int *current_pos = (int*)surely_malloc(comm_size * sizeof(int));
//     for (i = 0; i < comm_size; i++) {
//         current_pos[i] = comm_pattern->recv_displs[i];
//     }
    
//     for (i = 0; i < total_cols; i++) {
//         if (!seen[i]) continue;

//         int owner = i % comm_size;
//         if (owner != rank) {
//             int pos = current_pos[owner]++;
//             // store ghost metadata 
//             comm_pattern->ghost_col_indices[pos] = i;
//             comm_pattern->ghost_to_local[i] = pos;
//         }
//     }
    
//     free(current_pos);
//     free(seen);
    
//     // exchange counts using Alltoall -> each rank now knows how many and to whom
//     MPI_Alltoall(comm_pattern->recv_counts, 1, MPI_INT, comm_pattern->send_counts, 1, MPI_INT, MPI_COMM_WORLD);
    
//     // calc send displacements
//     comm_pattern->send_displs[0] = 0;
//     for (i = 1; i < comm_size; i++) {
//         comm_pattern->send_displs[i] = comm_pattern->send_displs[i-1] + comm_pattern->send_counts[i-1];
//     }
    
//     comm_pattern->total_to_send = comm_pattern->send_displs[comm_size-1] + comm_pattern->send_counts[comm_size-1];
//     comm_pattern->send_indices = (int*)malloc((comm_pattern->total_to_send > 0 ? comm_pattern->total_to_send : 1) * sizeof(int));
    
    
//     if (comm_pattern->send_indices == NULL) {
//         fprintf(stderr, "[Rank %d] ERROR: Cannot allocate send_indices\n", rank);
//         MPI_Abort(MPI_COMM_WORLD, 1);
//     }
    
//     // exchange column indices
//     MPI_Alltoallv(comm_pattern->ghost_col_indices,  // which indices must receive
//                   comm_pattern->recv_counts,        // how many indices must be received
//                   comm_pattern->recv_displs,        // from whom must be received
//                   MPI_INT,
//                   comm_pattern->send_indices,       // which local x entries to send
//                   comm_pattern->send_counts,        // how many entries must be sent
//                   comm_pattern->send_displs,        // to whom must be sent
//                   MPI_INT,
//                   MPI_COMM_WORLD);
    
//     // convert GLOBAL indices to LOCAL indices
//     for (i = 0; i < comm_pattern->total_to_send; i++) {
//         comm_pattern->send_indices[i] /= comm_size;
//     }
// }


// void exchange_ghost_vals(double *local_x, Comm_Pattern *comm_pattern, int rank, int comm_size) {
//     // allocate send buffer
//     double *send_buf = (double*)surely_malloc((comm_pattern->total_to_send > 0 ? comm_pattern->total_to_send : 1) * sizeof(double));
    
//     // pack outgoing vector entries
//     unsigned i;
//     for (i = 0; i < comm_pattern->total_to_send; i++) {
//         int local_idx = comm_pattern->send_indices[i];
//         send_buf[i] = local_x[local_idx];
//     }
    
//     // exchange values
//     MPI_Alltoallv(send_buf,
//                   comm_pattern->send_counts, 
//                   comm_pattern->send_displs, 
//                   MPI_DOUBLE,
//                   comm_pattern->ghost_values,
//                   comm_pattern->recv_counts, 
//                   comm_pattern->recv_displs, 
//                   MPI_DOUBLE,
//                   MPI_COMM_WORLD);
    
//     free(send_buf);
// }

