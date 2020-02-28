#include "types.h"
#include "constants.h"
#include "global.h"

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

void ladel_permute_symmetric_matrix(ladel_sparse_matrix *M, ladel_int *p, ladel_sparse_matrix *Mpp, ladel_work* work)
{
    if (!M || !Mpp) return;
    if (!p) 
    {
        ladel_sparse_copy(M, Mpp);
    } else
    {
        ladel_int col, pcol, row, prow, index, pindex, prev_col_count, ncol = M->ncol;
        ladel_int *col_counts = work->array_int_ncol1, *pinv = work->array_int_ncol2;
        for (index = 0; index < ncol; index++) col_counts[index] = 0;
        for (col = 0; col < ncol; col++) pinv[p[col]] = col;
        for (col = 0; col < ncol; col++)
        {
            pcol = pinv[col];
            for (index = M->p[col]; index < M->p[col+1]; index++)
            {
                prow = pinv[M->i[index]];
                col_counts[LADEL_MAX(pcol, prow)]++;
            }
        }
        Mpp->p[0] = 0;
        for (col = 1; col < ncol; col++)
        {
            prev_col_count = col_counts[col-1];
            Mpp->p[col] = prev_col_count;
            col_counts[col] += prev_col_count; 
            col_counts[col-1] = Mpp->p[col-1]; 
        }
        Mpp->p[ncol] = col_counts[ncol-1];
        col_counts[ncol-1] = Mpp->p[ncol-1];

        for (col = 0; col < ncol; col++)
        {
            pcol = pinv[col];
            for (index = M->p[col]; index < M->p[col+1]; index++)
            {
                prow = pinv[M->i[index]];
                if (pcol < prow)
                {
                    pindex = col_counts[prow]++;
                    Mpp->i[pindex] = pcol;
                } else 
                {
                    pindex = col_counts[pcol]++;
                    Mpp->i[pindex] = prow;
                }
                if (M->values) Mpp->x[pindex] = M->x[index]; 
            }
        }
    }
    
    
    
}