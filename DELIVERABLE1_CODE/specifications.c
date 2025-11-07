#include <stdlib.h>
#include <stdio.h>
#include "specifications.h"


//implement surely_malloc function
void *surely_malloc(size_t size) {
    void *ptr = malloc(size);
    if (ptr == NULL) {
        fprintf(stderr, "Out of memory while allocating %zu bytes\n", size);
        exit(EXIT_FAILURE);
    }
    return ptr;
}


// quicksort taken from Appel paper and modified for COO struct : qsort3.c; 
// https://github.com/cverified/cbench-vst/blob/master/qsort/qsort3.c

int coo_less (Sparse_Coordinate *p, unsigned a, unsigned b) {
    unsigned ra = p->row_indices[a], rb = p->row_indices[b];
    if (ra<rb) return 1;
    if (ra>rb) return 0;
    return p->col_indices[a] < p->col_indices[b];
}

void swap(Sparse_Coordinate *p, unsigned a, unsigned b) {
    unsigned i,j; double x;
    i=p->row_indices[a];
    j=p->col_indices[a];
    x=p->values[a];
    p->row_indices[a]=p->row_indices[b];
    p->col_indices[a]=p->col_indices[b];
    p->values[a]=p->values[b];
    p->row_indices[b]=i;
    p->col_indices[b]=j;
    p->values[b]=x;
}

/* Lexicographic quicksort by (row, col) */
void coo_quicksort(Sparse_Coordinate *p, unsigned base, unsigned n)
{
    unsigned lo, hi, left, right, mid;

    if (n == 0)
    return;
    lo = base;
    hi = lo + n - 1;
    while (lo < hi) {
    mid = lo + ((hi - lo) >> 1);

    if (coo_less(p,mid,lo))
        swap(p, mid, lo);
    if (coo_less(p,hi,mid)) {
        swap(p, mid, hi);
        if (coo_less(p,mid,lo))
        swap(p, mid, lo);
    }
    left = lo + 1;
    right = hi - 1;
    do {
        while (coo_less(p,left,mid))
        left++;
        while (coo_less(p,mid,right))
        right--;
        if (left < right) {
    swap(p, left, right);
        if (mid == left)
            mid = right;
        else if (mid == right)
            mid = left;
        left++;
        right--;
        } else if (left == right) {
        left++;
        right--;
        break;
        }
    } while (left <= right);
    if (right - lo > hi - left) {
        coo_quicksort(p, left, hi - left + 1);
        hi = right;
    } else {
        coo_quicksort(p, lo, right - lo + 1);
        lo = left;
    }
    }
}
