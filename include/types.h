#ifndef LADEL_TYPES_H
#define LADEL_TYPES_H

#ifdef DFLOAT
    typedef float ladel_double;
#else
    typedef double ladel_double;
#endif

#ifdef DLONG
    typedef long ladel_int;
#else
    typedef int ladel_int;
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

typedef struct symbolic_cholesky_information
{
    ladel_int ncol; /**<  number of columns in the analyzed matrix */
    ladel_int *etree; /**< eliminations tree*/
    ladel_int *postorder; /**< postordiring of the elimination tree */
    ladel_int *col_counts; /** < column counts, stored as column pointers */
    ladel_int *p; /** < fill-reducing ordering (AMD) */
    ladel_int *pattern; /** < stores the nonzero pattern of a row of L */ 
    ladel_int *nodes; /** < keeps track of which nodes have been marked */
} ladel_symbolics;

typedef struct ldl_factors
{
    ladel_int ncol;         /**< number of columns in the analyzed matrix */
    ladel_sparse_matrix *L; /**< L in LDL' factorization */
    ladel_double *D;        /**< D in LDL' factorization (stored as vector) */
    ladel_double *Dinv;     /**< D^-1 in LDL' factorization (stored as vector) */
    ladel_int *p;           /**< permutation vector */
} ladel_factor;

#endif /*LADEL_TYPES_H*/