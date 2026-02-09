#ifndef PARALLEL_READING_H
#define PARALLEL_READING_H

#include "mmio.h"
#include "specifications.h"

typedef struct {
    unsigned** rows;
    unsigned** cols;
    double** vals;
    int* counts;
    int* capacities;
} SendBuffers;

// BONUS: PARALLEL READING
void skip_partial_line(char* buffer, MPI_Offset bytes_read, char** cursor_out);
void free_send_buffers(SendBuffers* sb, int comm_size);

char* read_file_chunk(MPI_File fh, int rank, int comm_size, MPI_Offset file_size, MPI_Offset* bytes_read_out);
SendBuffers* create_send_buffers(int comm_size, int default_capacity);

int parse_mm_header(char* buffer, char* buffer_end, char** cursor_out, int* n_rows, int* n_cols, int* nnz);
int grow_array(void** array, unsigned* capacity, size_t elem_size);
int add_to_send_buffer(SendBuffers* sb, int dest_rank, unsigned row, unsigned col, double val);
int parse_and_partition(char* buffer, MPI_Offset chunk_size, 
                        MPI_Offset bytes_read, int rank, int comm_size,
                        unsigned** local_rows, unsigned** local_cols, 
                        double** local_vals, unsigned* local_count,
                        unsigned* local_capacity, SendBuffers* send_bufs,
                        unsigned n_rows);
int pack_send_data(SendBuffers* sb, int comm_size,
                        unsigned** packed_rows, unsigned** packed_cols, 
                        double** packed_vals, int** sdispls, int* total_send);
int exchange_data(int* send_counts, int comm_size,
                        unsigned* packed_s_rows, unsigned* packed_s_cols, 
                        double* packed_s_vals, int* sdispls,
                        unsigned** recv_rows, unsigned** recv_cols, 
                        double** recv_vals, int* total_recv);
int merge_data(unsigned** local_rows, unsigned** local_cols, 
                        double** local_vals, unsigned local_count,
                        unsigned* recv_rows, unsigned* recv_cols, 
                        double* recv_vals, int total_recv,
                        unsigned* final_nnz);


#endif