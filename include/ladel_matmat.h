#ifndef LADEL_MATMAT_H
#define LADEL_MATMAT_H

#include "ladel_types.h"

ladel_sparse_matrix *ladel_mat_mat_transpose(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_work *work);

ladel_sparse_matrix *ladel_mat_diag_mat_transpose(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_double *diag, ladel_work *work);

#endif /* LADEL_MATMAT_H */