
#ifndef SPMAT_SPMVM_UTIL_MERGE_H
#define SPMAT_SPMVM_UTIL_MERGE_H

#include <stddef.h>
#include <stdint.h>

/**
 * @brief Calculate which coordinate in both lists the merge path will have on the given diagonal
 * 
 * Merge list B is implicitly given by the natural numbers.
 * 
 * @param diagonal diagonal on which coordinate must lie
 * @param list_A First merge list
 * @param nor Number of rows of corresponding matrix
 * @param nnz Number of nonzero elements of corresponding matrix
 * @param list_a_coord Outputs the coordinate in list A 
 * @param list_b_coord Outputs the coordinate in list B
 */
void searchPathOnDiag(int32_t diagonal, const int32_t* list_A, int32_t nor, int32_t nnz,
                      int32_t* list_a_coord, int32_t* list_b_coord) {
    // Search range for A list
    int32_t a_coord_min = 0;
    if (diagonal > nnz) {
        a_coord_min = diagonal - nnz;
    }
    int32_t a_coord_max = (diagonal < nor) ? diagonal : nor;

    // Binary-search along the diagonal
    int32_t pivot;
    while (a_coord_min < a_coord_max) {
        pivot = (a_coord_min + a_coord_max) >> 1; // Division by shift for speed

        if (list_A[pivot] <= diagonal - pivot - 1) {
            a_coord_min = pivot + 1; // Keep top-right of diag
        } else {
            a_coord_max = pivot; // Keep bottom-left of diag
        }
    }

    *list_a_coord = (a_coord_min < nor) ? a_coord_min : nor;
    *list_b_coord = diagonal - a_coord_min;
}

#endif // SPMAT_SPMVM_UTIL_MERGE_H
