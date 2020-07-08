#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_upper_diag.h"
#include "ladel_transpose.h"
#include "ladel_debug_print.h"

void ladel_to_upper_diag(ladel_sparse_matrix *M)
{
    ladel_int index, row, col, Mptemp, nzM = 0;
    ladel_sparse_matrix *Mt;

    switch (M->symmetry)
    {
    case UPPER:
        /* do nothing */
        __attribute__ ((fallthrough));
    case LOWER:
        /* first transpose the matrix */
        Mt = ladel_transpose(M, TRUE, NULL);
        ladel_sparse_copy(Mt, M);
        ladel_sparse_free(Mt);
        __attribute__ ((fallthrough));
    default:
        for (col = 0; col < M->ncol; col++)
        {
            Mptemp = M->p[col];
            M->p[col] = nzM;
            for (index = Mptemp; index < ((M->nz) ? Mptemp + M->nz[col] : M->p[col+1]); index++)
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
        M->symmetry = UPPER;
        break;
    }   
}