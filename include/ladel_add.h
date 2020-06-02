#ifndef LADEL_ADD_H
#define LADEL_ADD_H

#include "ladel_types.h"

ladel_sparse_matrix *ladel_add_matrices(ladel_double alpha, ladel_sparse_matrix* A, ladel_double beta, ladel_sparse_matrix *B, ladel_work *work);
#endif /* LADEL_ADD_H */