#ifndef LADEL_ETREE_H
#define LADEL_ETREE_H

#include "types.h"

ladel_int ladel_etree(ladel_sparse_matrix *M, ladel_int *etree);

#ifdef SIMPLE_COL_COUNTS
ladel_int ladel_etree_and_col_counts(ladel_sparse_matrix *M, ladel_int *etree, ladel_int *col_counts);
#endif

#endif /*LADEL_ETREE_H*/