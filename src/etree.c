#include "types.h"
#include "global.h"
#include "constants.h"

ladel_int ladel_etree(ladel_sparse_matrix *M, ladel_int *etree)
{
    ladel_int index, row, col, next;
    ladel_int *ancestor = (ladel_int *) ladel_malloc(M->ncol, sizeof(ladel_int));
    if (!ancestor) return FAIL;

    for (col = 0; col < M->ncol; col++)
    {
        etree[col] = NONE;
        ancestor[col] = NONE;
        for (index = M->p[col]; index < M->p[col+1]; index++)
        {
            row = M->i[index];
            for (; row < col; row = next)
            {
                next = ancestor[row];
                ancestor[row] = col;
                if (next == NONE) 
                {
                    etree[row] = col;
                    break;
                }  
            }
        }
    }

    ladel_free(ancestor);
    return SUCCESS;
}