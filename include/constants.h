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

#endif /*LADEL_CONSTANTS_H*/