#include "mex.h"
#include <string.h>
#include "ladel.h"
#include "global.h"
#include "types.h"
#include "row_mod.h"

/* Modes of operation */
#define MODE_INIT "init"
#define MODE_FACTORIZE "factorize"
#define MODE_FACTORIZE_ADVANCED "factorize_advanced"
#define MODE_ROW_MOD "rowmod"
#define MODE_DENSE_SOLVE "solve"
#define MODE_DELETE "delete"

/* LADEL work identifier */
static ladel_work* work = NULL;

/* Mex calls this when it closes unexpectedly, freeing the workspace */
void exitFcn() {
  if (work != NULL) {
      ladel_workspace_free(work);
      work = NULL;
  }  
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

    /* report the default settings */
    if (strcmp(cmd, MODE_INIT) == 0) {
        /* Warn if other commands were ignored */
        if (nrhs != 2)
            mexErrMsgTxt("Wrong number of input arguments for mode init.");
        
        if (work != NULL)
            mexErrMsgTxt("Work is already initialized.");

        ladel_int ncol = (int)*mxGetPr(prhs[1]);
        work = ladel_workspace_allocate(ncol);
        return;

    } else if (strcmp(cmd, MODE_DELETE) == 0) {    
        /* clean up the problem workspace */
        if(work != NULL){
            ladel_workspace_free(work);
            work = NULL;
        }
        /* Warn if other commands were ignored */
        if (nlhs != 0 || nrhs != 1)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;

    } else {
        mexErrMsgTxt("Invalid LADEL mode");
    }
}