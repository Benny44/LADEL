#ifndef LADEL_PERMUTATION_H
#define LADEL_PERMUTATION_H

#include "types.h"

void ladel_permute_vector(ladel_double *x, ladel_int *p, ladel_int size, ladel_double *y);

void ladel_inverse_permute_vector(ladel_double *x, ladel_int *pinv, ladel_int size, ladel_double *y);

void ladel_permute_symmetric_matrix(ladel_sparse_matrix *M, ladel_int *p, ladel_sparse_matrix *Mpp);

#endif /*LADEL_PERMUTATION_H*/