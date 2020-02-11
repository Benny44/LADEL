#include "types.h"
#include "global.h"

void ladel_sparse_copy(ladel_sparse_matrix *M, ladel_sparse_matrix *M_copy)
{
    if (!M || !M_copy)
    {
        M_copy = NULL;
    } else
    {
        M_copy->ncol = M->ncol;
        M_copy->nrow = M->nrow;
        M_copy->nzmax = M->nzmax;
        M_copy->symmetry = M->symmetry;
        M_copy->values = M->values;
        ladel_int index;
        for (index = 0; index < M->ncol+1; index++) M_copy->p[index] = M->p[index];
        for (index = 0; index < M->nzmax; index++)
        {
            M_copy->i[index] = M->i[index];
            if (M->values) M_copy->x[index] = M->x[index];
        }
    }
}

void ladel_int_vector_copy(ladel_int *x, ladel_int size, ladel_int *y)
{
    ladel_int index;
    for (index = 0; index < size; index++)
    {
        y[index] = x[index];
    }
}