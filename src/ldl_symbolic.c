#include "types.h"
#include "constants.h"
#include "global.h"
#include "permutation.h"
#include "etree.h"
#include "postorder.h"
#include "col_counts.h"

#ifdef DAMD
#include "amd.h"
#endif /*DAMD*/

ladel_int ladel_ldl_symbolic(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method)
{
    ladel_sparse_matrix *Mpp = M;
    if (ordering_method == AMD)
    {
        #ifdef DAMD
        ladel_int status;
        double Info [AMD_INFO];
        #ifdef DLONG
        status = amd_l_order(M->ncol, M->p, M->i, sym->p, NULL, Info);
        #else /*DLONG*/
        status = amd_order(M->ncol, M->p, M->i, sym->p, NULL, Info);
        #endif

        #else /*ifdef DAMD*/
        sym->p = NULL;
        #endif
    } else
    {
        sym->p = NULL;
    }
    if (sym->p)
    {
        Mpp = ladel_sparse_alloc(M->nrow, M->ncol, M->nzmax, M->symmetry, M->values);
        ladel_permute_symmetric_matrix(M, sym->p, Mpp);
    }

    #ifdef SIMPLE_COL_COUNTS
    ladel_etree_and_col_counts(M, sym);
    #else
    ladel_etree(M, sym);
    ladel_postorder(M, sym);
    ladel_col_counts(M, sym);
    #endif /*SIMPLE_COL_COUNTS*/

    if (sym->p)
    {
        ladel_sparse_free(Mpp);
    }
}