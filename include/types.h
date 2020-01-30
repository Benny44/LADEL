#ifndef LADEL_TYPES_H
#define LADEL_TYPES_H

#ifdef DFLOAT
    typedef float ladel_double;
#else
    typedef double ladel_double;
#endif

#ifdef DINT
    typedef int ladel_int;
#else
    typedef long ladel_int;
#endif

typedef struct compressed_column_sparse_matrix 
{
    ladel_int nzmax; /**< number of nonzeros */
    ladel_int nrow; /**< number of rows */
    ladel_int ncol; /**< number of columns */

    ladel_int* p; /**< column pointers (size ncol+1) */
    ladel_int* i; /**< row pointers (size nzmax) */
    ladel_double* x; /**< numerical values (size nzmax) */

    ladel_int values; /**< has numerical values */
    ladel_int symmetry; /**< type of symmetry */    

} ladel_sparse_matrix;

#endif /*LADEL_TYPES_H*/