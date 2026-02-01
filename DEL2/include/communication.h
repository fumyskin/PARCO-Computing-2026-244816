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

#endif