#ifndef LADEL_H
#define LADEL_H

#include "types.h"

ladel_int ladel_factorize(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor *LD);

void ladel_dense_solve(const ladel_factor *LD, const ladel_double *rhs, ladel_double *y);

#endif /*LADEL_H*/