#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include "mmio.h"
#include "specifications.h"
#include "bonus_utils.h"


// NEEDS TO BE FIXED PLEASE FIX IT
// skip comment lines starting with '%' in matrix market file
// parse matrix dimensions
// only called by rank 0 !
int parse_mm_header(char* buffer, char* buffer_end, char** cursor_out, int* n_rows, int* n_cols, int* nnz){
    char* cursor = buffer;

    while(cursor < buffer_end && *cursor == '%'){
        char* next = memchr(cursor, '\n', buffer_end - cursor);
        if (!next) return -1;
        cursor = next + 1;
    }

    // parse dimensions
    int args_found = sscanf(cursor, "%d %d %d", n_rows, n_cols, nnz);
    if (args_found != 3) {
        return -1;
    }

    // move past header line
    char* next = memchr(cursor, '\n', buffer_end - cursor);
    if (next) {
        cursor = next + 1;
    } else {
        return -1; // no newline after dimensions is an error
    }
    
    *cursor_out = cursor;
    return 0;
}

// adjust cursor for non-zero ranks to skip partial first line
void skip_partial_line(char* buffer, MPI_Offset bytes_read, char** cursor_out){
    if (my_offset > 0) {
        char* eol = memchr(buffer, '\n', bytes_read);
        if (eol != NULL) {
            *cursor_out = eol + 1;
        } else {
            *cursor_out = buffer + bytes_read;
        }
    }
}

// grow dynamic array (safe realloc pattern ?)
int grow_array(void** array, size_t* capacity, size_t elem_size){
    *capacity *= 2;
    void* new_ptr = realloc(*array, (*capacity) * elem_size);
    if (!new_ptr) {
        return -1;
    }
    *array = new_ptr;
    return 0;
}

// READ FILE CHUNKS
char* read_file_chunk(MPI_File fh, int rank, int comm_size, MPI_Offset file_size, MPI_Offset* bytes_read_out) {
    // calculate chunk boundaries
    MPI_Offset chunk_size = file_size / comm_size;
    MPI_Offset my_offset = rank * chunk_size;
    
    // last rank takes remainder
    if (rank == comm_size - 1) {
        chunk_size = file_size - my_offset;
    }
    
    // add padding to ensure we can complete partial lines
    int padding = 1024;
    MPI_Offset bytes_to_read = chunk_size + padding;
    if (my_offset + bytes_to_read > file_size) {
        bytes_to_read = file_size - my_offset;
    }
    
    char* buffer = (char*)surely_malloc(bytes_to_read + 1);
    if (!buffer) return NULL;
    
    MPI_Status status;
    int ret = MPI_File_read_at_all(fh, my_offset, buffer, bytes_to_read, MPI_CHAR, &status);
    if (ret != MPI_SUCCESS) {
        free(buffer);
        return NULL;
    }
    
    buffer[bytes_to_read] = '\0';
    *bytes_read_out = bytes_to_read;
    
    return buffer;
}

// PARSE AND PARITION
SendBuffers* create_send_buffers(int comm_size, int default_capacity) {
    SendBuffers* sb = (SendBuffers*)surely_malloc(sizeof(SendBuffers));
    if (!sb) return NULL;
    
    sb->rows = (unsigned**)surely_malloc(comm_size * sizeof(unsigned*));
    sb->cols = (unsigned**)surely_malloc(comm_size * sizeof(unsigned*));
    sb->vals = (double**)surely_malloc(comm_size * sizeof(double*));
    sb->counts = (int*)calloc(comm_size, sizeof(int));
    sb->capacities = (int*)surely_malloc(comm_size * sizeof(int));
    
    if (!sb->rows || !sb->cols || !sb->vals || !sb->counts || !sb->capacities) {
        free(sb);
        return NULL;
    }
    
    for (int i = 0; i < comm_size; i++) {
        sb->capacities[i] = default_capacity;
        sb->rows[i] = (unsigned*)surely_malloc(default_capacity * sizeof(unsigned));
        sb->cols[i] = (unsigned*)surely_malloc(default_capacity * sizeof(unsigned));
        sb->vals[i] = (double*)surely_malloc(default_capacity * sizeof(double));
        
        if (!sb->rows[i] || !sb->cols[i] || !sb->vals[i]) {
            // Cleanup and return NULL
            for (int j = 0; j <= i; j++) {
                free(sb->rows[j]);
                free(sb->cols[j]);
                free(sb->vals[j]);
            }
            free(sb);
            return NULL;
        }
    }
    
    return sb;
}

int add_to_send_buffer(SendBuffers* sb, int dest_rank, unsigned row, unsigned col, double val) {
    int cnt = sb->counts[dest_rank];
    int cap = sb->capacities[dest_rank];
    
    // Grow if needed
    if (cnt >= cap) {
        sb->capacities[dest_rank] *= 2;
        cap = sb->capacities[dest_rank];
        
        unsigned* new_rows = (unsigned*)realloc(sb->rows[dest_rank], cap * sizeof(unsigned));
        unsigned* new_cols = (unsigned*)realloc(sb->cols[dest_rank], cap * sizeof(unsigned));
        double* new_vals = (double*)realloc(sb->vals[dest_rank], cap * sizeof(double));
        
        if (!new_rows || !new_cols || !new_vals) {
            return -1;
        }
        
        sb->rows[dest_rank] = new_rows;
        sb->cols[dest_rank] = new_cols;
        sb->vals[dest_rank] = new_vals;
    }
    
    sb->rows[dest_rank][cnt] = row;
    sb->cols[dest_rank][cnt] = col;
    sb->vals[dest_rank][cnt] = val;
    sb->counts[dest_rank]++;
    
    return 0;
}


void free_send_buffers(SendBuffers* sb, int comm_size) {
    if (!sb) return;
    
    for (int i = 0; i < comm_size; i++) {
        free(sb->rows[i]);
        free(sb->cols[i]);
        free(sb->vals[i]);
    }
    free(sb->rows);
    free(sb->cols);
    free(sb->vals);
    free(sb->counts);
    free(sb->capacities);
    free(sb);
}

int parse_and_partition(char* buffer, MPI_Offset chunk_size, 
                               MPI_Offset bytes_read, int rank, int comm_size,
                               unsigned** local_rows, unsigned** local_cols, 
                               double** local_vals, unsigned* local_count,
                               unsigned* local_capacity, SendBuffers* send_bufs) {
    char* cursor = buffer;
    char* end_of_buffer = buffer + bytes_read;
    char* chunk_end = buffer + chunk_size;
    if (chunk_end > end_of_buffer) chunk_end = end_of_buffer;
    
    while (cursor < chunk_end) {
        char* next_newline = memchr(cursor, '\n', end_of_buffer - cursor);
        if (!next_newline) break;
        
        int curr_row, curr_col;
        double curr_val;
        int args_found = sscanf(cursor, "%d %d %lf", &curr_row, &curr_col, &curr_val);
        
        if (args_found == 3) {
            // Convert to 0-based indexing
            curr_row--;
            curr_col--;
            
            // Determine owner
            int destination_rank = curr_row % comm_size;
            
            if (destination_rank == rank) {
                // Store locally
                if (*local_count >= *local_capacity) {
                    if (grow_array((void**)local_rows, local_capacity, sizeof(unsigned)) != 0 ||
                        grow_array((void**)local_cols, local_capacity, sizeof(unsigned)) != 0 ||
                        grow_array((void**)local_vals, local_capacity, sizeof(double)) != 0) {
                        return -1;
                    }
                }
                
                (*local_rows)[*local_count] = curr_row;
                (*local_cols)[*local_count] = curr_col;
                (*local_vals)[*local_count] = curr_val;
                (*local_count)++;
            } else {
                // Add to send buffer
                if (add_to_send_buffer(send_bufs, destination_rank, 
                                      curr_row, curr_col, curr_val) != 0) {
                    return -1;
                }
            }
        }
        
        cursor = next_newline + 1;
    }
    
    return 0;
}

int pack_send_data(SendBuffers* sb, int comm_size,
                         unsigned** packed_rows, unsigned** packed_cols, 
                         double** packed_vals, int** sdispls, int* total_send) {
    *sdispls = (int*)surely_malloc(comm_size * sizeof(int));
    if (!*sdispls) return -1;
    
    *total_send = 0;
    for (int i = 0; i < comm_size; i++) {
        (*sdispls)[i] = *total_send;
        *total_send += sb->counts[i];
    }
    
    *packed_rows = (unsigned*)surely_malloc(*total_send * sizeof(unsigned));
    *packed_cols = (unsigned*)surely_malloc(*total_send * sizeof(unsigned));
    *packed_vals = (double*)surely_malloc(*total_send * sizeof(double));
    
    if (!*packed_rows || !*packed_cols || !*packed_vals) {
        free(*sdispls);
        return -1;
    }
    
    for (int r = 0; r < comm_size; r++) {
        memcpy(&(*packed_rows)[(*sdispls)[r]], sb->rows[r], 
               sb->counts[r] * sizeof(unsigned));
        memcpy(&(*packed_cols)[(*sdispls)[r]], sb->cols[r], 
               sb->counts[r] * sizeof(unsigned));
        memcpy(&(*packed_vals)[(*sdispls)[r]], sb->vals[r], 
               sb->counts[r] * sizeof(double));
    }
    
    return 0;
}

int exchange_data(int* send_counts, int comm_size,
                        unsigned* packed_s_rows, unsigned* packed_s_cols, 
                        double* packed_s_vals, int* sdispls,
                        unsigned** recv_rows, unsigned** recv_cols, 
                        double** recv_vals, int* total_recv) {
    // Exchange counts
    int* recv_counts = (int*)surely_malloc(comm_size * sizeof(int));
    if (!recv_counts) return -1;
    
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);
    
    // Calculate receive displacements
    int* rdispls = (int*)surely_malloc(comm_size * sizeof(int));
    if (!rdispls) {
        free(recv_counts);
        return -1;
    }
    
    *total_recv = 0;
    for (int i = 0; i < comm_size; i++) {
        rdispls[i] = *total_recv;
        *total_recv += recv_counts[i];
    }
    
    // Allocate receive buffers
    *recv_rows = (unsigned*)surely_malloc(*total_recv * sizeof(unsigned));
    *recv_cols = (unsigned*)surely_malloc(*total_recv * sizeof(unsigned));
    *recv_vals = (double*)surely_malloc(*total_recv * sizeof(double));
    
    if (!*recv_rows || !*recv_cols || !*recv_vals) {
        free(recv_counts);
        free(rdispls);
        return -1;
    }
    
    // Execute exchanges
    MPI_Alltoallv(packed_s_rows, send_counts, sdispls, MPI_UNSIGNED,
                  *recv_rows, recv_counts, rdispls, MPI_UNSIGNED, MPI_COMM_WORLD);
    
    MPI_Alltoallv(packed_s_cols, send_counts, sdispls, MPI_UNSIGNED,
                  *recv_cols, recv_counts, rdispls, MPI_UNSIGNED, MPI_COMM_WORLD);
    
    MPI_Alltoallv(packed_s_vals, send_counts, sdispls, MPI_DOUBLE,
                  *recv_vals, recv_counts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
    
    free(recv_counts);
    free(rdispls);
    return 0;
}

// MERGE DATA
int merge_data(unsigned** local_rows, unsigned** local_cols, 
                     double** local_vals, unsigned local_count,
                     unsigned* recv_rows, unsigned* recv_cols, 
                     double* recv_vals, int total_recv,
                     unsigned* final_nnz) {
    *final_nnz = local_count + total_recv;
    
    unsigned* final_rows = (unsigned*)realloc(*local_rows, *final_nnz * sizeof(unsigned));
    unsigned* final_cols = (unsigned*)realloc(*local_cols, *final_nnz * sizeof(unsigned));
    double* final_vals = (double*)realloc(*local_vals, *final_nnz * sizeof(double));
    
    if (!final_rows || !final_cols || !final_vals) {
        return -1;
    }
    
    *local_rows = final_rows;
    *local_cols = final_cols;
    *local_vals = final_vals;
    
    // Append received data
    memcpy(&(*local_rows)[local_count], recv_rows, total_recv * sizeof(unsigned));
    memcpy(&(*local_cols)[local_count], recv_cols, total_recv * sizeof(unsigned));
    memcpy(&(*local_vals)[local_count], recv_vals, total_recv * sizeof(double));
    
    return 0;
}

// MAIN FUNCTION
int load_mm_chunked(const char* filename, Sparse_Coordinate* matrix, 
                    int rank, int comm_size) {
    MPI_File fh;
    int ret;
    
    // open file collectively
    ret = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (ret != MPI_SUCCESS) {
        if (rank == 0) fprintf(stderr, "Error opening file: %s\n", filename);
        return -1;
    }
    
    MPI_Offset file_size;
    MPI_File_get_size(fh, &file_size);
    
    // PHASE 1: read file chunks
    MPI_Offset bytes_read;
    char* buffer = read_file_chunk(fh, rank, comm_size, file_size, &bytes_read);
    if (!buffer) {
        MPI_File_close(&fh);
        return -1;
    }
    
    MPI_Offset chunk_size = file_size / comm_size;
    MPI_Offset my_offset = rank * chunk_size;
    if (rank == comm_size - 1) {
        chunk_size = file_size - my_offset;
    }
    
    // parse header, adjust cursor
    char* cursor = buffer;
    char* end_of_buffer = buffer + bytes_read;
    int n_rows = 0, n_cols = 0, nnz = 0;
    
    if (rank == 0) {
        if (parse_mm_header(buffer, end_of_buffer, &cursor, 
                           &n_rows, &n_cols, &nnz) != 0) {
            fprintf(stderr, "Error parsing matrix header\n");
            free(buffer);
            MPI_File_close(&fh);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }
    } else {
        skip_partial_line(buffer, bytes_read, &cursor);
    }
    
    MPI_Bcast(&n_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    unsigned local_capacity = (nnz / comm_size) + (nnz / (10 * comm_size)) + 100;
    unsigned local_count = 0;
    unsigned* local_rows = (unsigned*)surely_malloc(local_capacity * sizeof(unsigned));
    unsigned* local_cols = (unsigned*)surely_malloc(local_capacity * sizeof(unsigned));
    double* local_vals = (double*)surely_malloc(local_capacity * sizeof(double));

    // PHASE 2: parse and partition
    SendBuffers* send_bufs = create_send_buffers(comm_size, 1000);
    if (!send_bufs) {
        free(local_rows);
        free(local_cols);
        free(local_vals);
        free(buffer);
        MPI_File_close(&fh);
        return -1;
    }
    
    // Adjust cursor position in buffer
    MPI_Offset cursor_offset = cursor - buffer;
    MPI_Offset bytes_remaining = bytes_read - cursor_offset;

    // FIXED: For rank 0, adjust chunk_size to account for header
    MPI_Offset effective_chunk_size = chunk_size;
    if (rank == 0) {
        effective_chunk_size = chunk_size - cursor_offset;
    }

    ret = parse_and_partition(cursor, effective_chunk_size, 
                             bytes_remaining, rank, comm_size,
                             &local_rows, &local_cols, &local_vals, 
                             &local_count, &local_capacity, send_bufs);
    
    free(buffer);
    
    if (ret != 0) {
        free_send_buffers(send_bufs, comm_size);
        free(local_rows);
        free(local_cols);
        free(local_vals);
        MPI_File_close(&fh);
        return -1;
    }
    
    // PHASE 3: Exchange data
    unsigned* packed_s_rows;
    unsigned* packed_s_cols;
    double* packed_s_vals;
    int* sdispls;
    int total_send;
    
    if (pack_send_data(send_bufs, comm_size, &packed_s_rows, &packed_s_cols, 
                      &packed_s_vals, &sdispls, &total_send) != 0) {
        free_send_buffers(send_bufs, comm_size);
        free(local_rows);
        free(local_cols);
        free(local_vals);
        MPI_File_close(&fh);
        return -1;
    }
    
    unsigned* recv_rows;
    unsigned* recv_cols;
    double* recv_vals;
    int total_recv;
    
    ret = exchange_data(send_bufs->counts, comm_size, 
                       packed_s_rows, packed_s_cols, packed_s_vals, sdispls,
                       &recv_rows, &recv_cols, &recv_vals, &total_recv);
    
    free_send_buffers(send_bufs, comm_size);
    free(packed_s_rows);
    free(packed_s_cols);
    free(packed_s_vals);
    free(sdispls);
    
    if (ret != 0) {
        free(local_rows);
        free(local_cols);
        free(local_vals);
        MPI_File_close(&fh);
        return -1;
    }
    
    // PHASE 4: merge 
    unsigned final_nnz;
    ret = merge_data(&local_rows, &local_cols, &local_vals, local_count,
                    recv_rows, recv_cols, recv_vals, total_recv, &final_nnz);
    
    free(recv_rows);
    free(recv_cols);
    free(recv_vals);
    
    if (ret != 0) {
        free(local_rows);
        free(local_cols);
        free(local_vals);
        MPI_File_close(&fh);
        return -1;
    }
    
    // populate output
    matrix->n_rows = n_rows;
    matrix->n_cols = n_cols;
    matrix->nnz = final_nnz;
    matrix->row_indices = local_rows;
    matrix->col_indices = local_cols;
    matrix->values = local_vals;
    
    MPI_File_close(&fh);
    return 0;
}


// int load_mm_chunked(const char *filename, Sparse_Coordinate* matrix, int rank, int comm_size){

//     /*
//     int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
//                          int count, MPI_Datatype datatype, MPI_Status * status)

//     MPI_Offset = file offset
//     MPI_File_read_at_all is a collective routine that attempts to read from 
//     the file associated with fh (at the offset position) a total number of 
//     count data items having datatype type into the user’s buffer buf. 
//     The offset is in etype units relative to the current view. 
//     That is, holes are not counted when locating an offset. The data 
//     is taken out of those parts of the file specified by the current view.
//      MPI_File_read_at_all stores the number of datatype elements 
//      actually read in status. All other fields of status 
//      are undefined. It is erroneous to call this function if 
//      MPI_MODE_SEQUENTIAL mode was specified when the file was opened.
//     */

//     MPI_File fh;
//     MPI_Offset file_size;
//     MPI_Status status;
//     int ret;
//     // every process in MPI_COMM_WORLD open the filename in RDMODE
//     ret = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
//     if (ret != MPI_SUCCESS) {
//         if (rank == 0) fprintf(stderr, "Error opening file: %s\n", filename);
//         return -1;
//     }

//     // how many bytes each process reads
//     MPI_File_get_size(fh, &file_size);

//     MPI_Offset chunk_size = file_size / comm_size;
//     MPI_Offset my_offset = rank * chunk_size;

//     //  check for remainder: last one takes the rest -> mhhh is it the best approach?
//     if (rank == comm_size - 1) {
//         chunk_size = file_size - my_offset;
//     }


//     // allocate buffer (add padding for alignment/overlap!)
//     int padding = 1024;
//     MPI_Offset bytes_to_read = chunk_size + padding;
//     if (my_offset + bytes_to_read > file_size) {
//         bytes_to_read = file_size - my_offset;
//     }

//     char *buffer = (char*)surely_malloc(bytes_to_read+1); // +1 for null terminator

//     // everyone reads their chunk simultaneously + checking we're not going out
//     ret = MPI_File_read_at_all(fh, my_offset, buffer, bytes_to_read, MPI_CHAR, &status);
//     if (ret != MPI_SUCCESS) {
//         if (rank == 0) fprintf(stderr, "Error reading file\n");
//         free(buffer);
//         MPI_File_close(&fh);
//         return -1;
//     }
    
//     buffer[bytes_to_read] = '\0'; // Null terminate for safety

//     char* cursor = buffer;
//     char* end_of_buffer = buffer + bytes_to_read;

//     int n_rows = 0, n_cols = 0, nnz = 0;

//     if (rank == 0){
//         // deal with header ONLY ON RANK 0
//         while (cursor < end_of_buffer && *cursor == '%') {
//             char* next = memchr(cursor, '\n', end_of_buffer - cursor);
//             if (!next) break;
//             cursor = next + 1;
//         }
//         // extract n_rows, n_cols, nnz
//         int args_found = sscanf(cursor, "%d %d %d", &n_rows, &n_cols, &nnz);
//         if (args_found != 3){
//             fprintf(stderr, "Error parsing matrix header\n");
//             free(buffer);
//             MPI_File_close(&fh);
//             MPI_Abort(MPI_COMM_WORLD, 1);
//             return -1;
//         }

//         char* next = memchr(cursor, '\n', end_of_buffer - cursor);
//         if (next) cursor = next + 1;
//     }else{
//         if (my_offset > 0) {
//             char* eol = memchr(cursor, '\n', bytes_to_read);
//             if (eol != NULL) {
//                 cursor = eol + 1;
//             } else {
//                 // No complete line in this chunk
//                 cursor = end_of_buffer;
//             }
//         }
//     }

//     MPI_Bcast(&n_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     MPI_Bcast(&n_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

//     int initial_capacity = (nnz / comm_size) + (nnz / (10 * comm_size)) + 100;
//     unsigned local_count = 0;
//     unsigned local_capacity = initial_capacity;

//     unsigned* local_rows = (unsigned*)surely_malloc(local_capacity * sizeof(unsigned));
//     unsigned* local_cols = (unsigned*)surely_malloc(local_capacity * sizeof(unsigned));
//     double* local_vals = (double*)surely_malloc(local_capacity * sizeof(double));
    
//     int* send_counts = (int*)calloc(comm_size, sizeof(int));
//     int* send_capacities = (int*)surely_malloc(comm_size*sizeof(int));

//     // 2d array
//     unsigned** send_rows = (unsigned**)surely_malloc(comm_size*sizeof(unsigned*));
//     unsigned** send_cols = (unsigned**)surely_malloc(comm_size*sizeof(unsigned*));
//     double** send_vals   = (double**)surely_malloc(comm_size*sizeof(double*));

//     int default_send_cap = 1000;
//     for (int i = 0; i < comm_size; i++) {
//         send_capacities[i] = default_send_cap;
//         send_rows[i] = (unsigned*)surely_malloc(default_send_cap * sizeof(unsigned));
//         send_cols[i] = (unsigned*)surely_malloc(default_send_cap * sizeof(unsigned));
//         send_vals[i] = (double*)surely_malloc(default_send_cap * sizeof(double));
//     }

//     char* chunk_end = buffer + chunk_size;
//     if (chunk_end > end_of_buffer) chunk_end = end_of_buffer;
    
//     // parsing -> FOR ALL PROCESSES
//     while(cursor < chunk_end){
//         // next current line check
//         char* next_newline = memchr(cursor, '\n', end_of_buffer - cursor);
//         if (!next_newline) break;

//         int curr_row, curr_col;
//         double curr_val;
//         int args_found = sscanf(cursor, "%d %d %lf", &curr_row, &curr_col, &curr_val);

//         if (args_found == 3) {
//             // adjust to 0 based ! mm is 1 based
//             // curr_row--, curr_col-- 
//             curr_row--;
//             curr_col--;

//             int destination_rank = curr_row % comm_size;

//             if(destination_rank == rank){
//                  if (local_count >= local_capacity) {
//                     // grow local storage
//                     local_capacity *= 2;
//                     unsigned* new_rows = (unsigned*)realloc(local_rows, local_capacity * sizeof(unsigned));
//                     unsigned* new_cols = (unsigned*)realloc(local_cols, local_capacity * sizeof(unsigned));
//                     double* new_vals = (double*)realloc(local_vals, local_capacity * sizeof(double));
                    
//                     if (!new_rows || !new_cols || !new_vals) {
//                         fprintf(stderr, "Rank %d: Memory allocation failed\n", rank);
            
//                         free(buffer);
//                         MPI_File_close(&fh);
//                         MPI_Abort(MPI_COMM_WORLD, 1);
//                         return -1;
//                     }
                    
//                     local_rows = new_rows;
//                     local_cols = new_cols;
//                     local_vals = new_vals;
//                 }
//                 local_rows[local_count] = curr_row;
//                 local_cols[local_count] = curr_col;
//                 local_vals[local_count] = curr_val;
//                 local_count++;
//             }else{
//                 int cnt = send_counts[destination_rank];
//                 int cap = send_capacities[destination_rank];

//                 if(cnt >= cap){
//                     cap *= 2; // again, not sure if this makes sense
//                     send_capacities[destination_rank] = cap;
                   
//                     unsigned* new_rows = (unsigned*)realloc(send_rows[destination_rank], cap * sizeof(unsigned));
//                     unsigned* new_cols = (unsigned*)realloc(send_cols[destination_rank], cap * sizeof(unsigned));
//                     double* new_vals = (double*)realloc(send_vals[destination_rank], cap * sizeof(double));

//                     if (!new_rows || !new_cols || !new_vals) {
//                         fprintf(stderr, "Rank %d: Memory allocation failed for send buffer\n", rank);
//                         free(buffer);
//                         MPI_File_close(&fh);
//                         MPI_Abort(MPI_COMM_WORLD, 1);
//                         return -1;
//                     }

//                     send_rows[destination_rank] = new_rows;
//                     send_cols[destination_rank] = new_cols;
//                     send_vals[destination_rank] = new_vals;
//                 }

//                 send_rows[destination_rank][cnt] = curr_row;
//                 send_cols[destination_rank][cnt] = curr_col;
//                 send_vals[destination_rank][cnt] = curr_val;
//                 send_counts[destination_rank]++;
//             }
//         }
//         cursor = next_newline + 1;
//     }

//     free(buffer);

//     /*
//     int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
//         void *recvbuf, int recvcount, MPI_Datatype recvtype,
//         MPI_Comm comm)
//     */
    
//     // EXCHANGE 
//     int* recv_counts = (int*)surely_malloc(comm_size * sizeof(int));
//     MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

//     // displacements aaa
//     int* sdispls = (int*)surely_malloc(comm_size * sizeof(int));
//     int* rdispls = (int*)surely_malloc(comm_size * sizeof(int));
//     int total_send = 0;
//     int total_recv = 0;

//     for (int i = 0; i < comm_size; i++) {
//         sdispls[i] = total_send;
//         rdispls[i] = total_recv;
//         total_send += send_counts[i];
//         total_recv += recv_counts[i];
//     }

//     unsigned* packed_s_rows = (unsigned*)surely_malloc(total_send * sizeof(unsigned));
//     unsigned* packed_s_cols = (unsigned*)surely_malloc(total_send * sizeof(unsigned));
//     double* packed_s_vals = (double*)surely_malloc(total_send * sizeof(double));

//     for (int r = 0; r < comm_size; r++) {
//         // copy each bucket into the big contiguous array at the right offset
//         memcpy(&packed_s_rows[sdispls[r]], send_rows[r], send_counts[r] * sizeof(unsigned));
//         memcpy(&packed_s_cols[sdispls[r]], send_cols[r], send_counts[r] * sizeof(unsigned));
//         memcpy(&packed_s_vals[sdispls[r]], send_vals[r], send_counts[r] * sizeof(double));
        
//         free(send_rows[r]);
//         free(send_cols[r]);
//         free(send_vals[r]);
//     }

//     unsigned* recv_rows_buf = (unsigned*)surely_malloc(total_recv * sizeof(unsigned));
//     unsigned* recv_cols_buf = (unsigned*)surely_malloc(total_recv * sizeof(unsigned));
//     double* recv_vals_buf = (double*)surely_malloc(total_recv * sizeof(double));


//     // execute transfer 
//     MPI_Alltoallv(packed_s_rows, send_counts, sdispls, MPI_UNSIGNED,
//                   recv_rows_buf, recv_counts, rdispls, MPI_UNSIGNED, MPI_COMM_WORLD);
    
//     MPI_Alltoallv(packed_s_cols, send_counts, sdispls, MPI_UNSIGNED,
//                   recv_cols_buf, recv_counts, rdispls, MPI_UNSIGNED, MPI_COMM_WORLD);
    
//     MPI_Alltoallv(packed_s_vals, send_counts, sdispls, MPI_DOUBLE,
//                   recv_vals_buf, recv_counts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
    

//     // MERGE 
//     int new_total_nnz = local_count + total_recv;

//     unsigned* final_rows = (unsigned*)realloc(local_rows, new_total_nnz * sizeof(unsigned));
//     unsigned* final_cols = (unsigned*)realloc(local_cols, new_total_nnz * sizeof(unsigned));
//     double* final_vals = (double*)realloc(local_vals, new_total_nnz * sizeof(double));

//     if (!final_rows || !final_cols || !final_vals) {
//         fprintf(stderr, "Rank %d: Final memory allocation failed\n", rank);
//         free(buffer);
//         MPI_File_close(&fh);
//         MPI_Abort(MPI_COMM_WORLD, 1);
//         return -1;
//     }

//     local_rows = final_rows;
//     local_cols = final_cols;
//     local_vals = final_vals;

//     memcpy(&local_rows[local_count], recv_rows_buf, total_recv * sizeof(unsigned));
//     memcpy(&local_cols[local_count], recv_cols_buf, total_recv * sizeof(unsigned));
//     memcpy(&local_vals[local_count], recv_vals_buf, total_recv * sizeof(double));


//     matrix->n_rows = n_rows;
//     matrix->n_cols = n_cols;
//     matrix->nnz = new_total_nnz;  // local nnz for this rank
//     matrix->row_indices = local_rows;
//     matrix->col_indices = local_cols;
//     matrix->values = final_vals;

//     free(send_counts);
//     free(recv_counts);
//     free(sdispls);
//     free(rdispls);
//     free(send_capacities);
//     free(send_rows);
//     free(send_cols);
//     free(send_vals);
//     free(packed_s_rows);
//     free(packed_s_cols);
//     free(packed_s_vals);
//     free(recv_rows_buf);
//     free(recv_cols_buf);
//     free(recv_vals_buf);
    
//     MPI_File_close(&fh);
    
//     return 0;
// }