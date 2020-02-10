#include "types.h"
#include "constants.h"

void ladel_permute_vector(ladel_double *x, ladel_int *p, ladel_int size, ladel_double *y)
{
    ladel_int index;
    for (index = 0; index < size; index++) y[index] = x[p[index]];
}

void ladel_inverse_permute_vector(ladel_double *x, ladel_int *pinv, ladel_int size, ladel_double *y)
{
    ladel_int index;
    for (index = 0; index < size; index++) y[pinv[index]] = x[index];
}

void ladel_permute_symmetric_matrix(ladel_sparse_matrix *M, ladel_int *p, ladel_sparse_matrix *Mpp)
{
    if (!M || !Mpp) return;
    if (!p) 
    {
        ladel_sparse_copy(M, Mpp);
    } else
    {

    }
    
    
    
}