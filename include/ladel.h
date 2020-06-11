#ifndef LADEL_H
#define LADEL_H

#include "ladel_col_counts.h"
#include "ladel_constants.h"
#include "ladel_copy.h"
#include "ladel_debug_print.h"
#include "ladel_etree.h"
#include "ladel_global.h"
#include "ladel_ldl_numeric.h"
#include "ladel_ldl_symbolic.h"
#include "ladel_matvec.h"
#include "ladel_matmat.h"
#include "ladel_pattern.h"
#include "ladel_permutation.h"
#include "ladel_postorder.h"
#include "ladel_rank1_mod.h"
#include "ladel_row_mod.h"
#include "ladel_scale.h"
#include "ladel_transpose.h"
#include "ladel_types.h"
#include "ladel_upper_diag.h"
#include "ladel_add.h"
#include "ladel_submatrix.h"

ladel_int ladel_factorize(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_work* work);

ladel_int ladel_factorize_with_diag(ladel_sparse_matrix *M, ladel_diag d, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_work* work);

ladel_int ladel_factorize_advanced(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_sparse_matrix *Mbasis, ladel_work* work);

ladel_int ladel_factorize_advanced_with_diag(ladel_sparse_matrix *M, ladel_diag d, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_sparse_matrix *Mbasis, ladel_work* work);

ladel_int ladel_factorize_with_prior_basis_with_diag(ladel_sparse_matrix *M, ladel_diag d, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work);

ladel_int ladel_dense_solve(const ladel_factor *LD, const ladel_double *rhs, ladel_double *y, ladel_work* work);

#endif /*LADEL_H*/