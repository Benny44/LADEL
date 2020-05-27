#include "minunit.h"
#include "ladel.h"

#define NROW 4
#define NCOL 5
#define NZMAX 9
#define TOL 1e-8

static ladel_work *work;
static ladel_sparse_matrix *M;
/*
M = [1.2  0  0    0    0.5;
     0   -2  1.1  0    0;
     3.6  0  1.5  0  0;
     0   -3  0    1.7 -0.5;]
*/

void matmat_suite_setup(void)
{
    work = ladel_workspace_allocate(NCOL);
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UNSYMMETRIC, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 2; M->p[2] = 4; M->p[3] = 6; M->p[4] = 7; M->p[5] = 9;
    M->i[0] = 0; M->i[1] = 2; M->i[2] = 1; M->i[3] = 3; M->i[4] = 1; M->i[5] = 2; M->i[6] = 3; M->i[7] = 0; M->i[8] = 3;
    M->x[0] = 1.2; M->x[1] = 3.6; M->x[2] = -2; M->x[3] = -3; M->x[4] = 1.1; M->x[5] = 1.5; M->x[6] = 1.7; M->x[7] = 0.5; M->x[8] = -0.5;
}

void matmat_suite_teardown(void)
{
    ladel_sparse_free(M);
    ladel_workspace_free(work);
}

MU_TEST(test_mat_diag_mat_transpose)
{
    ladel_sparse_matrix *M_transpose = ladel_transpose(M, TRUE, work);
    ladel_double diag[NCOL] = {1, 2, 3, 4, 5};
    ladel_sparse_matrix *MMt = ladel_mat_diag_mat_transpose(M, M_transpose, diag, work);
    mu_assert_true(MMt != NULL);

    ladel_int p_sol[NROW+1] = {0, 1, 2, 5, 8};
    ladel_int i_sol[12] = {0, 1, 0, 2, 1, 1, 3, 0}; /* Unsorted! */
    ladel_double x_sol[12] = {2.69, 11.63, 4.32, 19.71, 4.95, 12.0, 30.81, -1.25};
    ladel_int index;
    for (index = 0; index < NROW+1; index++)
        mu_assert_long_eq(MMt->p[index], p_sol[index]);

    for (index = 0; index < 8; index++)
    {
        mu_assert_long_eq(MMt->i[index], i_sol[index]);
        mu_assert_double_eq(MMt->x[index], x_sol[index], TOL);
    }

    ladel_sparse_free(M_transpose);
    ladel_sparse_free(MMt);
}

MU_TEST_SUITE(suite_matmat)
{
    MU_SUITE_CONFIGURE(matmat_suite_setup, matmat_suite_teardown, NULL, NULL);
    MU_RUN_TEST(test_mat_diag_mat_transpose);
}