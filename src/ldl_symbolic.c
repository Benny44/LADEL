#include "types.h"
#include "constants.h"
#include "global.h"
#include "permutation.h"
#include "etree.h"
#include "postorder.h"
#include "col_counts.h"

#define DAMD
#ifdef DAMD
#include "amd.h"
#endif /*DAMD*/

ladel_int ladel_ldl_symbolic(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_sparse_matrix *Mpp)
{
    ladel_sparse_matrix *Mwork = M;
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
        if (status != AMD_OK) return FAIL;

        #else /*DAMD*/
        sym->p = ladel_free(sym->p);
        #endif
    } else
    {
        sym->p = ladel_free(sym->p);
    }
    if (sym->p)
    {
        ladel_permute_symmetric_matrix(M, sym->p, Mpp);
        Mwork = Mpp;
    }

    #ifdef SIMPLE_COL_COUNTS
    ladel_etree_and_col_counts(Mwork, sym);
    #else
    ladel_etree(Mwork, sym);
    ladel_postorder(Mwork, sym);
    ladel_col_counts(Mwork, sym);
    #endif /*SIMPLE_COL_COUNTS*/

    return SUCCESS;
}