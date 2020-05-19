#ifndef LADEL_GLOBAL_H
#define LADEL_GLOBAL_H

#include "ladel_types.h"
#include "ladel_constants.h"
#include "ladel_copy.h"
#include <stdlib.h>

#define LADEL_MAX(a, b) ((a) > (b) ? (a) : (b))
#define LADEL_MIN(a, b) ((a) > (b) ? (b) : (a))
#define LADEL_ABS(a) ((a) < 0 ? -(a) : (a))

#define LADEL_FOR(index, M, col) for(index = (M)->p[(col)]; index < (((M)->nz) ? (M)->p[(col)] + (M)->nz[(col)] : (M)->p[(col)+1]); index++)

void *ladel_malloc(ladel_int n, size_t size);

void *ladel_calloc(ladel_int n, size_t size);

void *ladel_free(void* p);

void *ladel_realloc(void *p, ladel_int n, size_t size, ladel_int *status);

#ifdef MATLAB
#include "mex.h"
#define ladel_print mexPrintf
#else
#define ladel_print printf
#endif

ladel_sparse_matrix *ladel_sparse_free(ladel_sparse_matrix *M);

ladel_sparse_matrix *ladel_sparse_alloc(ladel_int nrow, ladel_int ncol, 
                                                ladel_int nzmax, ladel_int symmetry,
                                                ladel_int values, ladel_int nz);

ladel_int ladel_sparse_realloc(ladel_sparse_matrix* M, ladel_int nzmax);

ladel_symbolics *ladel_symbolics_free(ladel_symbolics *sym);

ladel_symbolics *ladel_symbolics_alloc(ladel_int ncol);

ladel_factor *ladel_factor_free(ladel_factor *LD);

ladel_factor *ladel_factor_allocate(ladel_symbolics *sym);

ladel_set *ladel_set_free(ladel_set *set);

ladel_set *ladel_set_allocate(ladel_int ncol);

void ladel_set_set(ladel_set *set, ladel_int *set_vals, ladel_int size_set, ladel_int max_size_set);

ladel_work *ladel_workspace_free(ladel_work* work);

ladel_work *ladel_workspace_allocate(ladel_int ncol);


#endif /*LADEL_GLOBAL_H*/