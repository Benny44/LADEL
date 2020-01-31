#ifndef LADEL_SCALE_H
#define LADEL_SCALE_H

#include "types.h"

void ladel_scale_columns(ladel_sparse_matrix* M, ladel_double* S);

void ladel_scale_rows(ladel_sparse_matrix* M, ladel_double* S);

void ladel_scale_scalar(ladel_sparse_matrix* M, ladel_double s);


#endif /*LADEL_SCALE_H*/