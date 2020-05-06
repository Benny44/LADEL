#ifndef LADEL_DEBUG_PRINT_H
#define LADEL_DEBUG_PRINT_H

#include "global.h"
#include "types.h"

void ladel_print_sparse_matrix_matlab(ladel_sparse_matrix *M);

void ladel_print_dense_vector_matlab(ladel_double* x, size_t len);

void ladel_print_dense_int_vector_matlab(ladel_int* x, size_t len);

void ladel_print_factor_matlab(ladel_factor *LD);

void ladel_print_set(ladel_set *set);

#endif