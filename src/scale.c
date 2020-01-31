#include "global.h"
#include "types.h"
#include "scale.h"


void ladel_scale_columns(ladel_sparse_matrix* M, ladel_double* S)
{
    ladel_int col, index;
    for (col = 0; col < M->ncol; col++)
    {
        for (index = M->p[col]; index < M->p[col+1]; index++)
        {
            M->x[index] *= S[col];
        }    
    }
}

void ladel_scale_rows(ladel_sparse_matrix* M, ladel_double* S)
{
    ladel_int index;
    for (index = 0; index < M->nzmax; index++) 
    {
        M->x[index] *= S[M->i[index]];
    }
}

void ladel_scale_scalar(ladel_sparse_matrix* M, ladel_double s)
{
    ladel_int index;
    for (index = 0; index < M->nzmax; index++) 
    {
        M->x[index] *= s;
    }
}

void ladel_infinity_norm_columns(ladel_sparse_matrix *M, ladel_double *norms)
{
    ladel_int index, col;
    
    for (col = 0; col < M->ncol; col++)
    {
        norms[col] = 0;
        for (index = M->p[col]; index < M->p[col+1]; index++)
        {
            norms[col] = LADEL_MAX(norms[col], LADEL_ABS(M->x[index]));
        }
    }
}

void ladel_infinity_norm_rows(ladel_sparse_matrix *M, ladel_double *norms)
{
    ladel_int index, row;
    for (row = 0; row < M->nrow; row++) norms[row] = 0;
    for (index = 0; index < M->nzmax; index++)
    {
        row = M->i[index];
        norms[row] = LADEL_MAX(norms[row], LADEL_ABS(M->x[index]));
    }
}