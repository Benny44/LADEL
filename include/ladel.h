#ifndef LADEL_H
#define LADEL_H

#include "types.h"

void ladel_factorize(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor *LD);

#endif /*LADEL_H*/