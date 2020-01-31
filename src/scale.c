#include "types.h"
#include "scale.h"


void ladel_scale_columns(ladel_sparse_matrix* M, ladel_double* S)
{

}

void ladel_scale_rows(ladel_sparse_matrix* M, ladel_double* S)
{
    ladel_int index;
    for (index = 0; index < M->nzmax; index++) {
        M->x[index] *= S[M->i[index]];
    }
}

void ladel_scale_scalar(ladel_sparse_matrix* M, ladel_double s)
{

}