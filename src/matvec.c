#include "global.h"
#include "types.h"
#include "matvec.h"

void ladel_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset)
{
    ladel_int index, col;

    if (reset) for (index = 0; index < M->nrow; index++) y[index] = 0;

    for (col = 0; col < M->ncol; col++)
    {
        for (index = M->p[col]; index < M->p[col+1]; index++)
        {
            y[M->i[index]] += M->x[index] * x[col];
        }
    }
}

void ladel_tpose_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset)
{
    ladel_int index, col;

    if (reset) for (index = 0; index < M->ncol; index++) y[index] = 0; 
    
    for (col = 0; col < M->ncol; col++)
    {
        for (index = M->p[col]; index < M->p[col+1]; index++)
        {
            y[col] += M->x[index] * x[M->i[index]];
        }
    }
}