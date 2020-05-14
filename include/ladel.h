#ifndef LADEL_H
#define LADEL_H

#include "ladel_types.h"

ladel_int ladel_factorize(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_work* work);

ladel_int ladel_factorize_advanced(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_sparse_matrix *Mbasis, ladel_work* work);

ladel_int ladel_factorize_with_prior_basis(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work);

void ladel_dense_solve(const ladel_factor *LD, const ladel_double *rhs, ladel_double *y, ladel_work* work);

#endif /*LADEL_H*/