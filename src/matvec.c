#include "global.h"
#include "types.h"
#include "matvec.h"

void ladel_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset)
{
    ladel_int index, col, row;
    if (reset)
    {
        for (index = 0; index < M->nrow; index++) y[index] = 0;
    }
    for (col = 0; col < M->ncol; col++)
    {
        for (index = M->p[col]; index < M->p[col+1]; index++)
        {
            row = M->i[index];
            y[row] += M->x[index] * x[col];
        }
    }
}

void ladel_tpose_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset)
{
    
}