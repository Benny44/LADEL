/**
 * @file ladel_constants.h
 * @author Ben Hermans
 * @brief Constants (and a few simple macros) used in LADEL.
 */

#ifndef LADEL_CONSTANTS_H
#define LADEL_CONSTANTS_H

#define TRUE 1 /**< For booleans */
#define FALSE 0 /**< For booleans */

#define SUCCESS 1 /**< For status returns */
#define FAIL -1 /**< For status returns */

#define NONE -1 /**< Indicates a root (a node without parent), or an untouched node */

#define UNSYMMETRIC 0 /**< No symmetry is assumed in the matrix */
#define UPPER 1 /**< Use only upper part of matrix */
#define LOWER -1 /**< Use only lower part of matrix */

#define AMD 1 /**< Ordering method during the symbolic part of the factorization */
#define NO_ORDERING 0 /**< No ordering is performed during the symbolic part of the factorization */
#define GIVEN_ORDERING 2 /**< The ordering was computed previously and is already stored in sym->p */

#define IS_ROOT(col, etree) ((etree)[(col)] == NONE) /**< Check whether a node is a root (i.e. has no parent) */

#define MARKED 1 
#define UNMARKED 0
#define MARK(nodes, k) (nodes[(k)] = MARKED)
#define UNMARK(nodes, k) (nodes[(k)] = UNMARKED)
#define IS_MARKED(nodes, k) (nodes[(k)] == MARKED)

#define SET_HAS_CHANGED TRUE /**< Possible return value of ladel_set_union indicating the set has changed */
#define SET_HAS_NOT_CHANGED FALSE /**< Possible return value of ladel_set_union indicating the set has not changed */
#define MAX_SET_SIZE_EXCEEDED NONE /**< Possible return value of ladel_set_union indicating the set has grown beyond the maximum size */

#define UPDATE TRUE /**< Flag in rank1_update to perform an update */
#define DOWNDATE FALSE /**< Flag in rank1_update to perform a downdate */

#endif /*LADEL_CONSTANTS_H*/