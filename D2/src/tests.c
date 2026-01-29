#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>

#include "mmio.h"
#include "specifications.h"
#include "communication.h"
// bunch of test functions to check bit by bit if what built was correct

// TESTING OF GHOST ENTRIES
void create_test_matrix(const char* filename) {
    FILE *f = fopen(filename, "w");
    fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "4 4 8\n");
    // Row 0: columns 0, 1
    fprintf(f, "1 1 1.0\n");
    fprintf(f, "1 2 2.0\n");
    // Row 1: columns 2, 3
    fprintf(f, "2 3 3.0\n");
    fprintf(f, "2 4 4.0\n");
    // Row 2: columns 0, 2
    fprintf(f, "3 1 5.0\n");
    fprintf(f, "3 3 6.0\n");
    // Row 3: columns 1, 3
    fprintf(f, "4 2 7.0\n");
    fprintf(f, "4 4 8.0\n");
    fclose(f);
}

void verify_ghost_entries(Sparse_CSR *local_csr, Comm_Pattern *comm_pattern, 
                          int rank, int comm_size, int total_cols) {
    printf("\n[Rank %d] ===== VERIFICATION =====\n", rank);
    
    unsigned local_nnz = local_csr->row_ptr[local_csr->n_rows];
    bool *needs_ghost = (bool*)calloc(total_cols, sizeof(bool));
    int expected_ghosts = 0;
    
    // Find which columns we need but don't own
    for (unsigned i = 0; i < local_nnz; i++) {
        int col = local_csr->col_ind[i];
        int owner = col % comm_size;
        if (owner != rank && !needs_ghost[col]) {
            needs_ghost[col] = true;
            expected_ghosts++;
        }
    }
    
    // Verify count matches
    if (expected_ghosts == comm_pattern->num_ghost_cols) {
        printf("[Rank %d] ✓ Ghost count CORRECT: %d\n", rank, expected_ghosts);
    } else {
        printf("[Rank %d] ✗ Ghost count MISMATCH: expected %d, got %d\n", 
               rank, expected_ghosts, comm_pattern->num_ghost_cols);
    }
    
    // Verify all ghost columns are correct
    bool all_correct = true;
    for (int i = 0; i < comm_pattern->num_ghost_cols; i++) {
        int col = comm_pattern->ghost_col_indices[i];
        int owner = col % comm_size;
        
        if (owner == rank) {
            printf("[Rank %d] ✗ ERROR: Ghost col %d is owned by this rank!\n", rank, col);
            all_correct = false;
        }
        if (!needs_ghost[col]) {
            printf("[Rank %d] ✗ ERROR: Ghost col %d not in local matrix!\n", rank, col);
            all_correct = false;
        }
    }
    
    if (all_correct) {
        printf("[Rank %d] ✓ All ghost columns CORRECT\n", rank);
    }
    
    // Verify send/recv symmetry
    int *all_recv_totals = NULL;
    int *all_send_totals = NULL;
    if (rank == 0) {
        all_recv_totals = (int*)malloc(comm_size * sizeof(int));
        all_send_totals = (int*)malloc(comm_size * sizeof(int));
    }
    
    int my_recv_total = 0, my_send_total = 0;
    for (int p = 0; p < comm_size; p++) {
        my_recv_total += comm_pattern->recv_counts[p];
        my_send_total += comm_pattern->send_counts[p];
    }
    
    MPI_Gather(&my_recv_total, 1, MPI_INT, all_recv_totals, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&my_send_total, 1, MPI_INT, all_send_totals, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        int global_recv = 0, global_send = 0;
        for (int p = 0; p < comm_size; p++) {
            global_recv += all_recv_totals[p];
            global_send += all_send_totals[p];
        }
        
        if (global_recv == global_send) {
            printf("[Rank 0] ✓ Global send/recv symmetry CORRECT: %d == %d\n", 
                   global_send, global_recv);
        } else {
            printf("[Rank 0] ✗ Global send/recv MISMATCH: send=%d, recv=%d\n", 
                   global_send, global_recv);
        }
        
        free(all_recv_totals);
        free(all_send_totals);
    }
    
    // Verify ghost_to_local mapping
    int mapping_errors = 0;
    for (int i = 0; i < comm_pattern->num_ghost_cols; i++) {
        int col = comm_pattern->ghost_col_indices[i];
        int mapped_idx = comm_pattern->ghost_to_local[col];
        if (mapped_idx != i) {
            printf("[Rank %d] ✗ ghost_to_local[%d] = %d, expected %d\n", 
                   rank, col, mapped_idx, i);
            mapping_errors++;
        }
    }
    
    if (mapping_errors == 0) {
        printf("[Rank %d] ✓ ghost_to_local mapping CORRECT\n", rank);
    } else {
        printf("[Rank %d] ✗ ghost_to_local has %d errors\n", rank, mapping_errors);
    }
    
    free(needs_ghost);
    printf("[Rank %d] ===== END VERIFICATION =====\n\n", rank);
}
