#include "global.h"
#include "types.h"
#include "debug_print.h"
#include <stdio.h>

/* Print a sparse matrix so the output can be entered into matlab */
void ladel_print_sparse_matrix_matlab(ladel_sparse_matrix *M) {
    ladel_print("M = sparse(%ld, %ld);", M->nrow, M->ncol);
    ladel_int col, index = 0, row;
    ladel_double *Mx = M->x;
    ladel_int *Mi = M->i;
    ladel_int *Mp = M->p;

    for (col = 1; col <= M->ncol; col++) {
        for (row = Mp[col-1]; row < Mp[col]; row++) {
            ladel_print("M(%ld, %ld) = %.16le;", Mi[index]+1, col, Mx[index]);
            index++;
        }
    }
    ladel_print("\n");
}

void ladel_print_factor_matlab(ladel_factor *LD) {
    ladel_print("L = sparse(%ld, %ld);", LD->L->nrow, LD->L->ncol);
    ladel_int col, index = 0;
    ladel_double *Lx = LD->L->x;
    ladel_int *Li = LD->L->i;
    ladel_int *Lp = LD->L->p;
    ladel_int *nz = LD->L->nz;

    for (col = 1; col <= LD->L->ncol; col++) {
        for (index = Lp[col-1]; index < Lp[col-1] + nz[col-1]; index++) {
            ladel_print("L(%ld, %ld) = %.16le;", Li[index]+1, col, Lx[index]);
        }
    }
    ladel_print_dense_vector_matlab(LD->Dinv, LD->ncol);
}

void ladel_print_dense_vector_matlab(ladel_double* x, size_t len) {
    size_t k;
    ladel_print("x = zeros(%lu, 1);", len);
    for (k = 0; k < len; k++) {
        ladel_print("x(%lu) = %.16le;", k+1, x[k]);
    }
    ladel_print("\n");
}

void ladel_print_set(ladel_set *set)
{
    ladel_print("Size set %d (max %d)\n", set->size_set, set->max_size_set);
    ladel_print("Elements: ");
    ladel_int index;
    for (index = 0; index < set->size_set; index++) ladel_print("%d ", set->set[index]);
    ladel_print("\n");
}