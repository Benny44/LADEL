#ifndef LADEL_CONSTANTS_H
#define LADEL_CONSTANTS_H

#define TRUE 1
#define FALSE 0

#define SUCCESS 1
#define FAIL -1

#define NONE -1

#define UNSYMMETRIC 0
#define UPPER 1 /**< Use only upper part of matrix*/
#define LOWER -1 /**< Use only lower part of matrix*/

#define AMD 1
#define NO_ORDERING 0

#define IS_ROOT(col, etree) ((etree)[(col)] == NONE)

#define MARKED 1
#define UNMARKED 0
#define MARK(nodes, k) (nodes[(k)] = MARKED)
#define UNMARK(nodes, k) (nodes[(k)] = UNMARKED)
#define IS_MARKED(nodes, k) (nodes[(k)] == MARKED)

#endif /*LADEL_CONSTANTS_H*/