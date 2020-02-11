#ifndef LADEL_PATTERN_H
#define LADEL_PATTERN_H

#include "types.h"

ladel_int ladel_nonzero_pattern_of_row_in_L(ladel_sparse_matrix *M, 
                                                ladel_symbolics *sym,
                                                ladel_int row);

#endif /*LADEL_PATTERN_H*/