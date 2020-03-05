#ifndef LADEL_ROW_MOD_H
#define LADEL_ROW_MOD_H

#include "types.h"

ladel_int ladel_row_add(ladel_factor *LD, ladel_symbolics *sym, ladel_int row_in_L, ladel_sparse_matrix *W, ladel_int col_in_W, ladel_double diag, ladel_work *work);

ladel_int ladel_row_del(ladel_factor *LD, ladel_symbolics *sym, ladel_int row_in_L, ladel_work *work);

#endif /*LADEL_ROW_MOD_H*/