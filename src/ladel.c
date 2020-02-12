#include "types.h"
#include "constants.h"
#include "global.h"
#include "ldl_symbolic.h"
#include "ldl_numeric.h"

ladel_int ladel_factorize(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor *LD)
{
    ladel_int ok_symbolic, ok_numeric;
    ladel_sparse_matrix *Mpp = ladel_sparse_alloc(M->nrow, M->ncol, M->nzmax, M->symmetry, M->values);
    if (!Mpp) return FAIL;
    ok_symbolic = ladel_ldl_symbolic(M, sym, ordering_method, Mpp);
    ok_numeric = ladel_ldl_numeric(Mpp, sym, LD);
    ladel_sparse_free(Mpp);
    if (ok_symbolic && ok_numeric) return SUCCESS;
    else return FAIL;
}

void ladel_dense_solve(const ladel_factor *LD, const ladel_double *rhs, ladel_double *y)
{
    ladel_sparse_matrix *L = LD->L;
    ladel_double *Dinv = LD->Dinv;
    ladel_int index, row, ncol = LD->L->ncol;
    for (row = 0; row < ncol; row++) y[row] = rhs[row];

    for (row = 0; row < ncol; row++)
    {
        for (index = L->p[row]+1; index < L->p[row+1]; index++)
        {
            y[L->i[index]] -= L->x[index]*y[row];
        }
    }
    for (row = 0; row < ncol; row++) y[row] *= Dinv[row];
    for (row = ncol-1; row >= 0; row++)
    {
        for (index = L->p[row]+1; index < L->p[row+1]; index++)
        {
            y[row] -= L->x[index]*y[L->i[index]];
        }
    } 
}