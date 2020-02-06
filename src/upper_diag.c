#include "types.h"
#include "global.h"
#include "upper_diag.h"


void ladel_to_upper_diag(ladel_sparse_matrix *M)
{
    ladel_int index, row, col, Mptemp, nzM = 0;

    for (col = 0; col < M->ncol; col++)
    {
        Mptemp = M->p[col];
        M->p[col] = nzM;
        for (index = Mptemp; index < M->p[col+1]; index++)
        {
            row = M->i[index];
            if (row <= col)
            {
                M->i[nzM] = row;
                if (M->values) M->x[nzM] = M->x[index];
                nzM++;
            }
        }
    }
    M->p[M->ncol] = nzM;
    ladel_sparse_realloc(M, nzM);
}