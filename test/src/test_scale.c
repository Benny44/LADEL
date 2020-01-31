#include "minunit.h"
#include "types.h"
#include "scale.h"

#define NROW 4
#define NCOL 5

ladel_sparse_matrix *M;
ladel_double col_scale[NCOL];
ladel_double row_scale[NROW];

void scale_suite_setup(void)
{

}

void scale_suite_teardown(void)
{

}

MU_TEST(test_scale_rows)
{

}

MU_TEST(test_scale_rows)
{

}

MU_TEST(test_scale_rows)
{

}

MU_TEST(test_scale_columns)
{

}

MU_TEST(test_scale_scalar)
{

}

MU_TEST_SUITE(suite_scale) 
{
    MU_SUITE_CONFIGURE(scale_suite_setup, scale_suite_teardown, NULL, NULL);

    MU_RUN_TEST(test_scale_rows);
    MU_RUN_TEST(test_scale_columns);
    MU_RUN_TEST(test_scale_scalar);
}