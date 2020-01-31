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

}