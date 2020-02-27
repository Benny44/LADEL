#include "global.h"
#include "types.h"

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

void ladel_print_dense_vector_matlab(ladel_double* x, size_t len) {
    size_t k;
    printf("x = zeros(%lu, 1);", len);
    for (k = 0; k < len; k++) {
        printf("x(%lu) = %.16le;", k+1, x[k]);
    }
    printf("\n");
}