#ifndef LADEL_LDL_SYMBOLIC_H
#define LADEL_LDL_SYMBOLIC_H

#include "types.h"
ladel_int ladel_ldl_symbolic(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_sparse_matrix *Mpp, ladel_work* work);

#endif /*LADEL_LDL_SYMBOLIC_H*/
