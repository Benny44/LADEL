#ifndef LADEL_COPY_H
#define LADEL_COPY_H

#include "ladel_types.h"

ladel_sparse_matrix *ladel_sparse_allocate_and_copy(ladel_sparse_matrix *M);

void ladel_sparse_copy(ladel_sparse_matrix *M, ladel_sparse_matrix *M_copy);

void ladel_int_vector_copy(ladel_int *x, ladel_int size, ladel_int *y);

void ladel_double_vector_copy(ladel_double *x, ladel_int size, ladel_double *y);

#endif /*LADEL_COPY_H*/