#include "global.h"
#include "types.h"
#include "debug_print.h"
#include <stdio.h>

/* Print a sparse matrix so the output can be entered into matlab */
void ladel_print_sparse_matrix_matlab(ladel_sparse_matrix *M) {
    printf("M = sparse(%ld, %ld);", M->nrow, M->ncol);
    ladel_int col, index = 0, row;
    ladel_double *Mx = M->x;
    ladel_int *Mi = M->i;
    ladel_int *Mp = M->p;

    for (col = 1; col <= M->ncol; col++) {
        for (row = Mp[col-1]; row < Mp[col]; row++) {
            printf("M(%ld, %ld) = %.16le;", Mi[index]+1, col, Mx[index]);
            index++;
        }
    }
    printf("\n");
}

void ladel_print_factor_matlab(ladel_factor *LD) {
    printf("L = sparse(%ld, %ld);", LD->L->nrow, LD->L->ncol);
    ladel_int col, index = 0;
    ladel_double *Lx = LD->L->x;
    ladel_int *Li = LD->L->i;
    ladel_int *Lp = LD->L->p;
    ladel_int *nz = LD->L->nz;

    for (col = 1; col <= LD->L->ncol; col++) {
        for (index = Lp[col-1]; index < Lp[col-1] + nz[col-1]; index++) {
            printf("L(%ld, %ld) = %.16le;", Li[index]+1, col, Lx[index]);
        }
    }
    ladel_print_dense_vector_matlab(LD->Dinv, LD->ncol);
}

void ladel_print_dense_vector_matlab(ladel_double* x, size_t len) {
    size_t k;
    printf("x = zeros(%lu, 1);", len);
    for (k = 0; k < len; k++) {
        printf("x(%lu) = %.16le;", k+1, x[k]);
    }
    printf("\n");
}