#include "types.h"
#include "global.h"
#include "constants.h"

ladel_sparse_matrix *ladel_transpose(ladel_sparse_matrix *M, int values)
{
    ladel_sparse_matrix *M_transpose = ladel_sparse_alloc(M->ncol, M->nrow, M->nzmax, -M->symmetry, values && M->values);
    ladel_int *col_pointers = ladel_calloc(M->nrow, sizeof(ladel_int));
    if (!M_transpose || !col_pointers) 
    {
        ladel_sparse_free(M_transpose);
        ladel_free(col_pointers);
        return NULL;
    }

    ladel_int index, col, new_index, prev_col_count;
    for (index = 0; index < M->nzmax; index++) col_pointers[M->i[index]]++;
    
    M_transpose->p[0] = 0;
    for (col = 1; col < M_transpose->ncol; col++)
    {
        prev_col_count = col_pointers[col-1];
        col_pointers[col] += prev_col_count;
        M_transpose->p[col] = prev_col_count;
        col_pointers[col-1] = M_transpose->p[col-1];
    } 
    M_transpose->p[M_transpose->ncol] = col_pointers[M_transpose->ncol-1];
    col_pointers[M_transpose->ncol-1] = M_transpose->p[M_transpose->ncol-1];

    for (col = 0; col < M->ncol; col++)
    {
        for (index = M->p[col]; index < M->p[col+1]; index++)
        {
            new_index = col_pointers[M->i[index]]++;
            M_transpose->i[new_index] = col;
            if (M_transpose->values) M_transpose->x[new_index] = M->x[index];
        }
    }
    ladel_free(col_pointers);
    return M_transpose;
}
