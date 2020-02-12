#include "minunit.h"
#include "types.h"
#include "global.h"
#include "constants.h"
#include "ladel.h"

#define NCOL 5
#define NROW 5
#define NZMAX 10
#define TOL 1e-10

ladel_sparse_matrix *M;
ladel_symbolics *sym;
ladel_factor *LD;

void ldl_suite_setup(void)
{
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE);
    M->p[0] = 0; M->p[1] = 1; M->p[2] = 3; M->p[3] = 4; M->p[4] = 7; M->p[5] = 10;
    M->i[0] = 0; M->i[1] = 0; M->i[2] = 1; M->i[3] = 2; M->i[4] = 1; M->i[5] = 2; M->i[6] = 3; M->i[7] = 0; M->i[8] = 3; M->i[9] = 4;  
    M->x[0] = 1; M->x[1] = 10; M->x[2] = 2; M->x[3] = -3; M->x[4] = 11; M->x[5] = 12; M->x[6] = 4; M->x[7] = -3; M->x[8] = 2; M->x[9] = -5;

    
}

void ldl_suite_teardown(void)
{
    ladel_sparse_free(M);
}

void ldl_test_setup(void)
{
    sym = ladel_symbolics_alloc(NCOL);  
}

void ldl_test_teardown(void)
{
    ladel_symbolics_free(sym);
    ladel_factor_free(LD);
}

MU_TEST(test_simple_ldl)
{
    ladel_double x[NCOL] = {1, 2, 3, 4, 5};
    ladel_double y[NCOL], y_ref[NCOL] = {-1.738103756708408e-01, -1.081932021466905e-01, 4.379964221824687e-01, 3.594991055456172e-01, -7.519141323792485e-01};

    ladel_int status = ladel_factorize(M, sym, NO_ORDERING, &LD);
    
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, x, y);
    ladel_int index;

    for (index = 0; index < NCOL; index++)
        mu_assert_double_eq(y[index], y_ref[index], TOL);
    
}

#ifdef DAMD
MU_TEST(test_simple_ldl_with_amd)
{
    ladel_double x[NCOL] = {1, 2, 3, 4, 5};
    ladel_double y[NCOL], y_ref[NCOL] = {-1.738103756708408e-01, -1.081932021466905e-01, 4.379964221824687e-01, 3.594991055456172e-01, -7.519141323792485e-01};

    ladel_int status = ladel_factorize(M, sym, AMD, &LD);
    
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, x, y);
    ladel_int index;

    for (index = 0; index < NCOL; index++)
        mu_assert_double_eq(y[index], y_ref[index], TOL);
    
}
#endif
MU_TEST_SUITE(suite_ldl)
{
    MU_SUITE_CONFIGURE(ldl_suite_setup, ldl_suite_teardown, ldl_test_setup, ldl_test_teardown);
    MU_RUN_TEST(test_simple_ldl);
    #ifdef DAMD
    MU_RUN_TEST(test_simple_ldl_with_amd);
    #endif
}