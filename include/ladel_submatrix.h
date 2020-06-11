#ifndef LADEL_SUBMATRIX_H
#define LADEL_SUBMATRIX_H

#include "ladel_types.h"

ladel_sparse_matrix *ladel_column_submatrix(ladel_sparse_matrix *M, ladel_int* cols, ladel_int nb_cols);
#endif /* LADEL_SUBMATRIX_H */