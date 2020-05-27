#ifndef LADEL_MATVEC_H
#define LADEL_MATVEC_H

#include "ladel_global.h"
#include "ladel_types.h"

void ladel_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset);

void ladel_tpose_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset);

void ladel_symmetric_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset);

#endif /* LADEL_MATVEC_H */