#ifndef LADEL_PATTERN_H
#define LADEL_PATTERN_H

#include "types.h"

ladel_int ladel_nonzero_pattern_of_row_in_L(ladel_sparse_matrix *M, 
                                                ladel_symbolics *sym,
                                                ladel_int row);

ladel_int ladel_etree_dfs(ladel_sparse_matrix *W, 
                            ladel_symbolics *sym,
                            ladel_int row,
                            ladel_int maximum_row);

#endif /*LADEL_PATTERN_H*/