#ifndef LADEL_RANK1_MOD_H
#define LADEL_RANK1_MOD_H

#include "types.h"

ladel_int ladel_set_addition(ladel_int *set1, ladel_int *set2, ladel_int *size_set1, ladel_int max_size_set1, ladel_int size_set2);

ladel_int ladel_rank1_update(ladel_factor **LD, ladel_symbolics **sym, ladel_sparse_matrix *W, ladel_int col_in_W);

#endif