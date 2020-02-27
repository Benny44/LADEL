#ifndef LADEL_GLOBAL_H
#define LADEL_GLOBAL_H

#include "types.h"
#include "constants.h"
#include "copy.h"
#include "stdlib.h"

#define LADEL_MAX(a, b) ((a) > (b) ? (a) : (b))
#define LADEL_MIN(a, b) ((a) > (b) ? (b) : (a))
#define LADEL_ABS(a) ((a) < 0 ? -(a) : (a))

static void *ladel_malloc(ladel_int n, size_t size) 
{
    return (malloc(LADEL_MAX(n, 1) * size));
}

static void *ladel_calloc(ladel_int n, size_t size)
{
    return (calloc(LADEL_MAX(n, 1), size));
}

static void *ladel_free(void* p) 
{
    if (p) free(p);
    return NULL;
}

static void *ladel_realloc(void *p, ladel_int n, size_t size, ladel_int *status)
{
    void *p_new;
    p_new = realloc(p, LADEL_MAX(n, 1) * size);
    *status = (p_new != NULL);
    return ((*status) ? p_new : p);
}

static ladel_sparse_matrix *ladel_sparse_free(ladel_sparse_matrix *M)
{
    if (!M) return NULL;
    ladel_free(M->p);
    ladel_free(M->i);
    ladel_free(M->x);
    ladel_free(M->nz);
    return ((ladel_sparse_matrix *) ladel_free(M));
}

static ladel_sparse_matrix *ladel_sparse_alloc(ladel_int nrow, ladel_int ncol, 
                                                ladel_int nzmax, ladel_int symmetry,
                                                ladel_int values)
{
    ladel_sparse_matrix *M = (ladel_sparse_matrix *) ladel_calloc(1, sizeof(ladel_sparse_matrix));
    if (!M) return NULL;
    M->nrow = nrow;
    M->ncol = ncol;
    M->nzmax = nzmax;
    M->values = values;
    M->symmetry = symmetry;
    M->p = (ladel_int *) ladel_malloc(ncol+1, sizeof(ladel_int));
    M->i = (ladel_int *) ladel_malloc(nzmax, sizeof(ladel_int));
    M->x = values ? (ladel_double *) ladel_malloc(nzmax, sizeof(ladel_double)) : NULL;
    M->nz = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    if (!M->p || !M->i || (values && !M->x) || !M->nz) M = ladel_sparse_free(M);
    return M;
}

static ladel_int ladel_sparse_realloc(ladel_sparse_matrix* M, ladel_int nzmax)
{
    ladel_int status_i, status, status_x = SUCCESS;
    if (!M) return NULL;
    if (nzmax <= 0) nzmax = M->p[M->ncol];
    M->i = (ladel_int *) ladel_realloc(M->i, nzmax, sizeof(ladel_int), &status_i);
    if (M->values) M->x = (ladel_double *) ladel_realloc(M->x, nzmax, sizeof(ladel_double), &status_x);
    status = status_i && status_x;
    if (status == SUCCESS) M->nzmax = nzmax;
    return status;
}

static ladel_symbolics *ladel_symbolics_free(ladel_symbolics *sym)
{
    if (!sym) return NULL;
    ladel_free(sym->etree);
    ladel_free(sym->postorder);
    ladel_free(sym->col_counts);
    ladel_free(sym->p);
    ladel_free(sym->pattern);
    ladel_free(sym->nodes);
    return (ladel_symbolics *) ladel_free(sym);
}

static ladel_symbolics *ladel_symbolics_alloc(ladel_int ncol)
{
    ladel_symbolics *sym = (ladel_symbolics *) ladel_calloc(1, sizeof(ladel_symbolics));
    if (!sym) return NULL;
    sym->ncol = ncol;
    sym->etree = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->postorder = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->col_counts = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->p = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->pattern = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->nodes = (ladel_int *) ladel_calloc(ncol, sizeof(ladel_int));
    if (!sym->etree || !sym->postorder || !sym->col_counts || !sym->pattern || !sym->nodes) 
        sym = ladel_symbolics_free(sym);
    return sym;
}

static ladel_factor *ladel_factor_free(ladel_factor *LD)
{
    if (!LD) return NULL;
    ladel_sparse_free(LD->L);
    ladel_free(LD->D);
    ladel_free(LD->Dinv);
    ladel_free(LD->p);
    return (ladel_factor *) ladel_free(LD);
}

static ladel_factor *ladel_factor_allocate(ladel_symbolics *sym)
{
    ladel_factor *LD = (ladel_factor *) ladel_calloc(1, sizeof(ladel_factor));
    if (!LD || !sym) return NULL;
    ladel_int ncol = LD->ncol = sym->ncol;
    LD->L = ladel_sparse_alloc(ncol, ncol, sym->col_counts[ncol-1], UNSYMMETRIC, TRUE);
    LD->D = ladel_malloc(ncol, sizeof(ladel_double));
    LD->Dinv = ladel_malloc(ncol, sizeof(ladel_double));
    if (!LD->L || !LD->D || !LD->Dinv)
    {
        ladel_factor_free(LD);
        return NULL;
    }
    if (sym->p)
    {
        LD->p = ladel_malloc(ncol, sizeof(ladel_int));
        if (!LD->p)
        {
            ladel_factor_free(LD);
            return NULL;
        }
        ladel_int_vector_copy(sym->p, ncol, LD->p);
    } else 
    {
        LD->p = NULL;
    }
    return LD; 
}

#endif /*LADEL_GLOBAL_H*/