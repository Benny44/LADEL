#ifndef LADEL_SCALE_H
#define LADEL_SCALE_H

#include "ladel_types.h"

void ladel_scale_columns(ladel_sparse_matrix *M, ladel_double* S);

void ladel_scale_rows(ladel_sparse_matrix *M, ladel_double* S);

void ladel_scale_scalar(ladel_sparse_matrix *M, ladel_double s);

void ladel_infinity_norm_columns(ladel_sparse_matrix *M, ladel_double *norms);

void ladel_infinity_norm_rows(ladel_sparse_matrix *M, ladel_double *norms);


#endif /*LADEL_SCALE_H*/