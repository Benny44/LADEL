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

    ladel_int *p; /**< column pointers (size ncol+1) */
    ladel_int *i; /**< row pointers (size nzmax) */
    ladel_double *x; /**< numerical values (size nzmax) */

    ladel_int *nz; /** < number of elements in each column (size ncol) */

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
    ladel_int *pinv;        /**< inverse permutation vector */
} ladel_factor;

typedef struct ladel_set_struct {
    ladel_int *set;
    ladel_int size_set;
    ladel_int max_size_set;
} ladel_set;

typedef struct ladel_col_struct {
    ladel_int *i;
    ladel_double *x;
    ladel_int nz;
    ladel_int nzmax;
} ladel_col;

/* Workspace needed for the factorization and updates */
typedef struct workspace
{
    ladel_set *set_preallocated1;
    ladel_set *set_preallocated2;
    ladel_set *set_preallocated3;
    ladel_set *set_unallocated_values1;
    ladel_set *set_unallocated_values2;
    ladel_set *set_unallocated_values3;
    ladel_int *array_int_ncol1;
    ladel_int *array_int_ncol2;
    ladel_int *array_int_ncol3;
    ladel_int *array_int_ncol4;
    ladel_double *array_double_all_zeros_ncol1;
} ladel_work;

#endif /*LADEL_TYPES_H*/