#include "types.h"
#include "global.h"
#include "constants.h"
#include "pattern.h"
#include "rank1_mod.h"
#include "row_mod.h"

ladel_int ladel_row_add(ladel_factor *LD, ladel_symbolics *sym, ladel_int row_in_L, ladel_sparse_matrix *W, ladel_int col_in_W, ladel_double diag, ladel_work *work)
{
    if (!LD || !sym || !W || !work) return FAIL;
    
    ladel_int start, index_in_pattern, ncol = sym->ncol, row, index, index2, status;
    ladel_sparse_matrix* L = LD->L;
    ladel_double *Dinv = LD->Dinv;
    ladel_double temp, d22 = diag;
    ladel_double *l12 = work->array_double_all_zeros_ncol1;

    ladel_set *l31_pattern = work->set_preallocated1;
    l31_pattern->size_set = 0;
    ladel_set *set_L31 = work->set_unallocated_values2;
    ladel_set *difference = work->set_preallocated2;
    ladel_int *offset = work->array_int_ncol1;
    ladel_int *insertions = work->array_int_ncol2;
    
    /* 1. Solve lower triangular system L11*D11*l12 = W12 */
    for (index = W->p[col_in_W]; index < W->p[col_in_W+1]; index++) 
        l12[W->i[index]] = W->x[index];

    start = ladel_etree_dfs(W, sym, col_in_W, row_in_L);
    for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++)
    {
        row = sym->pattern[index_in_pattern];
        temp = l12[row];
        /* 2. d22 = c22 - l12^T D11*l12 */
        d22 -= temp*temp*Dinv[row];
        
        l12[row] *= Dinv[row];
        /* Gaussian elimination */  
        for (index = L->p[row]; index < (L->p[row] + L->nz[row]) && L->i[index] < row_in_L; index++)
        {
            l12[L->i[index]] -= L->x[index]*temp;
        }
        /* 3. l32 = (c32 - L31*D11*l12)/d22 */
        ladel_set_set(set_L31, L->i + index, (L->p[row] + L->nz[row])-index, ncol);
        ladel_set_union(l31_pattern, set_L31, difference, offset, insertions, row_in_L);
        for (index2 = (L->p[row] + L->nz[row] - 1); index2 >= index; index2--)
        {
            l12[L->i[index2]] -= L->x[index2]*temp;
            /* Shift the columns down by one to make room for l12*/
            L->i[index2+1] = L->i[index2];
            L->x[index2+1] = L->x[index2]; 
        }
        /* Insert l12[row] */
        L->i[index] = row_in_L;
        L->x[index] = l12[row];
        l12[row] = 0;
        L->nz[row]++;
    }

    /* Insert l31 */
    d22 = Dinv[row_in_L] = 1/d22;
    L->nz[row_in_L] = l31_pattern->size_set;
    for (index = L->p[row_in_L]; index < (L->p[row_in_L] + L->nz[row_in_L]); index++)
    {
        L->i[index] = row = l31_pattern->set[index];
        L->x[index] = l12[row];
        l12[row] = 0;
    } 
    
    /* 4. w = l32*sqrt(abs(d22)) */
    /* 5. Update or downdate L33*D33*L33^T = L33*D33*L33^T - sign(d22)*w*w^T */
    status = ladel_rank1_update(LD, sym, L, row_in_L, 1/sqrt(d22), d22 < 0, work);
    return status;
}

ladel_int ladel_row_del(ladel_factor *LD, ladel_symbolics *sym, ladel_int row_in_L, ladel_work *work)
{

}