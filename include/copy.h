#ifndef LADEL_COPY_H
#define LADEL_COPY_H

#include "types.h"
void ladel_sparse_copy(ladel_sparse_matrix *M, ladel_sparse_matrix *M_copy);

void ladel_int_vector_copy(ladel_int *x, ladel_int size, ladel_int *y);

#endif /*LADEL_COPY_H*/