#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "communication.h"

void find_ghost_entries(Sparse_CSR *local_matrix, Comm_Pattern *comm_pattern, int rank, int comm_size, int total_cols) {
    
    unsigned i, r;
    unsigned nnz = local_matrix->row_ptr[local_matrix->n_rows];
    
    // count unique columns: form vector of vector elements long
    char *seen = (char*)calloc(total_cols, sizeof(char));
    if (!seen) {
        fprintf(stderr, "[Rank %d] ERROR: Cannot allocate seen array\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (r = 0; r < local_matrix->n_rows; r++) {
        for (i = local_matrix->row_ptr[r];
             i < local_matrix->row_ptr[r + 1]; i++) {
            seen[ local_matrix->col_ind[i] ] = 1;
        }
    }
    
    // communication arrays: 
    comm_pattern->send_counts = (int*)calloc(comm_size, sizeof(int));
    comm_pattern->recv_counts = (int*)calloc(comm_size, sizeof(int));
    comm_pattern->send_displs = (int*)calloc(comm_size, sizeof(int));
    comm_pattern->recv_displs = (int*)calloc(comm_size, sizeof(int));
    
    if (!comm_pattern->send_counts || !comm_pattern->recv_counts || 
        !comm_pattern->send_displs || !comm_pattern->recv_displs) {
        fprintf(stderr, "[Rank %d] ERROR: Cannot allocate comm arrays\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // count ghost columns
    int ghost_count = 0;
    for (i = 0; i < total_cols; i++) {
        if (!seen[i]) continue;

        int owner = i % comm_size;
        if (owner != rank) {
            comm_pattern->recv_counts[owner]++;
            ghost_count++;
        }
    }
    
    comm_pattern->num_ghost_cols = ghost_count;
    comm_pattern->total_cols = total_cols;
    
    // receive displacements -> for each i finds how many cells each recv occupy
    comm_pattern->recv_displs[0] = 0;
    for (i = 1; i < comm_size; i++) {
        comm_pattern->recv_displs[i] =
            comm_pattern->recv_displs[i - 1] +
            comm_pattern->recv_counts[i - 1];
    }
    
    comm_pattern->ghost_col_indices = (int*)surely_malloc((ghost_count > 0 ? ghost_count : 1) * sizeof(int));
    comm_pattern->ghost_values = (double*)surely_malloc((ghost_count > 0 ? ghost_count : 1) * sizeof(double));
    comm_pattern->ghost_to_local = (int*)surely_malloc(total_cols * sizeof(int));
    
    for (i = 0; i < total_cols; i++) {
        comm_pattern->ghost_to_local[i] = -1;
    }
    
    int *current_pos = (int*)surely_malloc(comm_size * sizeof(int));
    for (i = 0; i < comm_size; i++) {
        current_pos[i] = comm_pattern->recv_displs[i];
    }
    
    for (i = 0; i < total_cols; i++) {
        if (!seen[i]) continue;

        int owner = i % comm_size;
        if (owner != rank) {
            int pos = current_pos[owner]++;
            comm_pattern->ghost_col_indices[pos] = i;
            comm_pattern->ghost_to_local[i] = pos;
        }
    }
    
    free(current_pos);
    free(seen);
    
    // exchange counts
    MPI_Alltoall(comm_pattern->recv_counts, 1, MPI_INT, comm_pattern->send_counts, 1, MPI_INT, MPI_COMM_WORLD);
    
    // calc send displacements
    comm_pattern->send_displs[0] = 0;
    for (i = 1; i < comm_size; i++) {
        comm_pattern->send_displs[i] = comm_pattern->send_displs[i-1] + comm_pattern->send_counts[i-1];
    }
    
    comm_pattern->total_to_send = comm_pattern->send_displs[comm_size-1] + comm_pattern->send_counts[comm_size-1];
    comm_pattern->send_indices = (int*)malloc((comm_pattern->total_to_send > 0 ? comm_pattern->total_to_send : 1) * sizeof(int));
    
    
    if (comm_pattern->send_indices == NULL) {
        fprintf(stderr, "[Rank %d] ERROR: Cannot allocate send_indices\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // exchange column indices
    MPI_Alltoallv(comm_pattern->ghost_col_indices,
                  comm_pattern->recv_counts, comm_pattern->recv_displs, MPI_INT,
                  comm_pattern->send_indices,
                  comm_pattern->send_counts, comm_pattern->send_displs, MPI_INT,
                  MPI_COMM_WORLD);
    
    // convert global indices to local indices
    for (i = 0; i < comm_pattern->total_to_send; i++) {
        comm_pattern->send_indices[i] /= comm_size;
    }
}


void exchange_ghost_vals(double *local_x, Comm_Pattern *comm_pattern, int rank, int comm_size) {
    double *send_buf = (double*)surely_malloc((comm_pattern->total_to_send > 0 ? comm_pattern->total_to_send : 1) * sizeof(double));
    
    unsigned i;
    for (i = 0; i < comm_pattern->total_to_send; i++) {
        int local_idx = comm_pattern->send_indices[i];
        send_buf[i] = local_x[local_idx];
    }
    
    MPI_Alltoallv(send_buf,
                  comm_pattern->send_counts, comm_pattern->send_displs, MPI_DOUBLE,
                  comm_pattern->ghost_values,
                  comm_pattern->recv_counts, comm_pattern->recv_displs, MPI_DOUBLE,
                  MPI_COMM_WORLD);
    
    free(send_buf);
}