#ifndef LADEL_RANK1_MOD_H
#define LADEL_RANK1_MOD_H

#include "types.h"

typedef struct ladel_set_struct {
    ladel_int *set;
    ladel_int size_set;
    ladel_int max_size_set;
} ladel_set;

ladel_set *ladel_init_set(ladel_int *set_vals, ladel_int size_set, ladel_int max_size_set);

ladel_int ladel_add_nonzero_pattern_to_L(ladel_int *Li, ladel_double *Lx, ladel_int *nz_L, ladel_int nzmax_L, ladel_int *set, ladel_int size_set);

ladel_int ladel_set_union(ladel_set *first_set, ladel_set *second_set);

ladel_int ladel_rank1_update(ladel_factor **LD, ladel_symbolics **sym, ladel_sparse_matrix *W, ladel_int col_in_W);

#endif