#ifndef LDL_NUMERIC_H
#define LDL_NUMERIC_H

#include "ladel_types.h"

ladel_int ladel_ldl_numeric_with_diag(ladel_sparse_matrix *Mpp, ladel_diag d, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work);

ladel_int ladel_ldl_numeric(ladel_sparse_matrix *Mpp, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work);

#endif /*LDL_NUMERIC_H*/