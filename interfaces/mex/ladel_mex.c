#include "mex.h"
#include <string.h>
#include "ladel.h"
#include "global.h"
#include "types.h"
#include "row_mod.h"
#include "upper_diag.h"
#include "debug_print.h"

/* Modes of operation */
#define MODE_INIT "init"
#define MODE_FACTORIZE "factorize"
#define MODE_FACTORIZE_ADVANCED "factorize_advanced"
#define MODE_ROW_MOD "rowmod"
#define MODE_DENSE_SOLVE "solve"
#define MODE_DELETE "delete"

/* LADEL work identifier */
static ladel_work* work = NULL;
static ladel_symbolics *sym = NULL;
static ladel_factor *LD = NULL;

/* Mex calls this when it closes unexpectedly, freeing the workspace */
void exitFcn() {
  if (work != NULL) {
      ladel_workspace_free(work);
      work = NULL;
  }  
}

ladel_sparse_matrix *ladel_get_sparse_from_matlab(const mxArray *M_mex, ladel_sparse_matrix *M, ladel_int symmetry)
{
    // ladel_sparse_matrix *M;
    M->nrow = mxGetM(M_mex);
    M->ncol = mxGetN(M_mex);
    M->p = (ladel_int *) mxGetJc(M_mex);
    M->i = (ladel_int *) mxGetIr(M_mex);
    M->x = (ladel_double *) mxGetPr(M_mex);
    M->values = TRUE;
    M->symmetry = symmetry;
    M->nzmax = M->p[M->ncol];
    M->nz = NULL;
    return M;
}

/**
 * The gateway function to LADEL
 *
 * Usages:
 * ladel_mex('init', ncol);
 * ladel_mex('factorize', M);
 * ladel_mex('factorize', M, ordering);
 * ladel_mex('factorize_advanced', M, Mbasis);
 * ladel_mex('factorize_advanced', M, Mbasis, ordering);
 * ladel_mex('rowmod', row_index);
 * ladel_mex('rowmod', row_index, row);
 * y = ladel_mex('solve', x);
 * ladel_mex('delete');
 *
 * @param nlhs Number of output arguments
 * @param plhs Array of output argument pointers
 * @param nrhs Number of input arguments
 * @param prhs Array of input argument pointers
 */
void mexFunction(int nlhs, mxArray * plhs [], int nrhs, const mxArray * prhs []) {
    
    /* Set function to call when mex closes unexpectedly */
    mexAtExit(exitFcn);

    /* Get the command string */
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    if (strcmp(cmd, MODE_INIT) == 0) 
    {
        /* Warn if other commands were ignored */
        if (nrhs != 2)
            mexErrMsgTxt("Wrong number of input arguments for mode init.");
        
        if (work != NULL)
            mexErrMsgTxt("Work is already initialized.");

        ladel_int ncol = (ladel_int) *mxGetPr(prhs[1]);
        work = ladel_workspace_allocate(ncol);
        sym = ladel_symbolics_alloc(ncol);
        return;

    } 
    else if (strcmp(cmd, MODE_DELETE) == 0) 
    {    
        /* clean up the problem workspace */
        if(work != NULL){
            work = ladel_workspace_free(work);
            sym = ladel_symbolics_free(sym);
            LD = ladel_factor_free(LD);
        }
        /* Warn if other commands were ignored */
        if (nlhs != 0 || nrhs != 1)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;

    } 
    else if (strcmp(cmd, MODE_FACTORIZE) == 0)
    {
        if (nlhs != 0 || (nrhs != 2 && nrhs != 3))
            mexErrMsgTxt("Wrong number of input or output arguments for mode factorize.");

        if (LD != NULL) LD = ladel_factor_free(LD);

        ladel_sparse_matrix Mmatlab;
        ladel_sparse_matrix *M = ladel_get_sparse_from_matlab(prhs[1], &Mmatlab, UPPER);
        ladel_int ordering;
        if (nrhs == 3)
            ordering = (ladel_int) *mxGetPr(prhs[2]);
        else
            ordering = NO_ORDERING;
        
        ladel_int status = ladel_factorize(M, sym, ordering, &LD, work);
        if (status != SUCCESS)
            mexErrMsgTxt("Factorize: Something went wrong in the factorization.");
    }
    else if (strcmp(cmd, MODE_DENSE_SOLVE) == 0)
    {
        if (nlhs != 1 || nrhs != 2 )
            mexErrMsgTxt("Wrong number of input or output arguments for mode solve.");

        plhs[0] = mxCreateDoubleMatrix(LD->ncol,1,mxREAL);
        ladel_double *y = mxGetPr(plhs[0]); 
        ladel_double *x = mxGetPr(prhs[1]); 
        ladel_dense_solve(LD, x, y, work);
    }
    else if (strcmp(cmd, MODE_FACTORIZE_ADVANCED) == 0)
    {
        if (nlhs != 0 || (nrhs != 3 && nrhs != 4))
            mexErrMsgTxt("Wrong number of input or output arguments for mode factorize_advanced.");

        if (LD != NULL) LD = ladel_factor_free(LD);

        ladel_sparse_matrix Mmatlab;
        ladel_sparse_matrix *M = ladel_get_sparse_from_matlab(prhs[1], &Mmatlab, UPPER);

        ladel_sparse_matrix Mbasismatlab;
        ladel_sparse_matrix *Mbasis = ladel_get_sparse_from_matlab(prhs[2], &Mbasismatlab, UPPER);

        ladel_int ordering;
        if (nrhs == 4)
            ordering = (ladel_int) *mxGetPr(prhs[3]);
        else
            ordering = NO_ORDERING;
        

        ladel_int status = ladel_factorize_advanced(M, sym, ordering, &LD, Mbasis, work);
        if (status != SUCCESS)
            mexErrMsgTxt("Factorize_advanced: Something went wrong in the factorization.");
    }
    else 
    {
        mexErrMsgTxt("Invalid LADEL mode");
    }
}