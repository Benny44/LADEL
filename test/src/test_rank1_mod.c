#include "minunit.h"
#include "global.h"
#include "constants.h"
#include "types.h"
#include "ladel.h"
#include "rank1_mod.h"

#define NCOL 8
#define NROW 8
#define NZMAX 24
#define TOL 1e-8

ladel_sparse_matrix *M;
ladel_sparse_matrix *W;
ladel_factor *LD;
ladel_symbolics *sym;

void rank1_mod_suite_setup(void)
{
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE);
    M->p[0] = 0; M->p[1] = 1; M->p[2] = 2; M->p[3] = 4; M->p[4] = 7; M->p[5] = 8; M->p[6] = 10; M->p[7] = 14; M->p[8] = 16;
    M->i[0] = 0; M->i[1] = 1; M->i[2] = 0; M->i[3] = 2; M->i[4] = 1; M->i[5] = 2; M->i[6] = 3; M->i[7] = 4;  
    M->i[8] = 4; M->i[9] = 5; M->i[10] = 2; M->i[11] = 3; M->i[12] = 5; M->i[13] = 6; M->i[14] = 5; M->i[15] = 7; 
    M->x[0] = 1.484233718441293e+00; M->x[1] = 8.503811060900474e-01; M->x[2] = 7.946847833340259e-02; 
    M->x[3] = 4.093923468942169e-01; M->x[4] = 6.944674301359678e-02; M->x[5] = 1.761109237578858e-01; 
    M->x[6] = 3.766406539600695e-01; M->x[7] = 9.168193399034030e-01; M->x[8] = 9.238873678854942e-01; 
    M->x[9] = 1.897902116880709e+00; M->x[10] = 6.052733699027846e-01; M->x[11] = 2.665692902440686e-01; 
    M->x[12] = 7.650155176644616e-02; M->x[13] = 1.151758087181832e+00; M->x[14] = 7.767464464874704e-01; 
    M->x[15] = 6.404488288848778e-01; 

    W = ladel_sparse_alloc(NROW, 1, 3, UNSYMMETRIC, TRUE);
    W->p[0] = 0; W->p[1] = 3;
    W->i[0] = 3; W->i[1] = 5; W->i[2] = 7;
    W->x[0] = 1.418863386272153e-01; W->x[1] = 4.217612826262750e-01; W->x[2] = 9.157355251890671e-01;  

    sym = ladel_symbolics_alloc(NCOL);   
}

void rank1_mod_suite_teardown(void)
{
    ladel_sparse_free(M);
    ladel_sparse_free(W);
    ladel_symbolics_free(sym);
}

void rank1_mod_test_setup(void)
{

}

void rank1_mod_test_teardown(void)
{
    ladel_factor_free(LD);
}

MU_TEST(test_rank1_update)
{
    ladel_int status;
    status = ladel_factorize(M, sym, AMD, &LD);
    mu_assert_long_eq(status, SUCCESS);

    ladel_double rhs[8] = {7.922073295595544e-01, 9.594924263929030e-01, 6.557406991565868e-01,
                            3.571167857418955e-02, 8.491293058687771e-01, 9.339932477575505e-01,
                            6.787351548577735e-01, 7.577401305783334e-01};
    ladel_double sol[8] = {1.626739458964396e+01, 1.404806681196643e+00, -2.938574982347923e+02,
                            -3.385740249326764e+00, 6.516216869838241e+02, -6.457174936667400e+02,
                            1.986908367068955e+02, 7.843195055030803e+02};

    ladel_double x[8];
    ladel_dense_solve(LD, rhs, x);
    ladel_int index;
    for (index = 0; index < NCOL; index++) mu_assert_double_eq(x[index], sol[index], TOL);

    ladel_double sol_mod[8] = {4.021591066933866e-01, 1.257689300408343e+00, 2.457694262215502e+00,
                                -1.584275766315217e+00, 3.426556805737886e+00, -2.481259429010251e+00, 
                                -1.707843620820627e-01, 2.602541252030759e+00};
    
    status = ladel_rank1_update(&LD, &sym, W, 0);
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, rhs, x);
    for (index = 0; index < NCOL; index++) mu_assert_double_eq(x[index], sol[index], TOL);
}


MU_TEST_SUITE(suite_rank1_mod)
{
    MU_SUITE_CONFIGURE(rank1_mod_suite_setup, rank1_mod_suite_teardown, rank1_mod_test_setup, rank1_mod_test_teardown);
    MU_RUN_TEST(test_rank1_update);
}