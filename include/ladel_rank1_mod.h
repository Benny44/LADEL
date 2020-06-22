#ifndef LADEL_RANK1_MOD_H
#define LADEL_RANK1_MOD_H

#include "ladel_types.h"

ladel_set *ladel_init_set(ladel_int *set_vals, ladel_int size_set, ladel_int max_size_set);

ladel_col *ladel_init_col(ladel_int *i, ladel_double *x, ladel_int nz , ladel_int nzmax);

void ladel_set_col(ladel_col* col,ladel_int *i, ladel_double *x, ladel_int nz , ladel_int nzmax);

ladel_int ladel_add_nonzero_pattern_to_col_of_L(ladel_sparse_matrix* L, ladel_int col, ladel_set *col_set, ladel_set *set, ladel_set *difference, ladel_int* offset, ladel_int* insertions);

ladel_int ladel_set_union(ladel_set *first_set, ladel_set *second_set, ladel_set *difference, ladel_int *offset, ladel_int *insertions, ladel_int minimum_index);

#endif