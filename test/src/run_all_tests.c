#include "minunit.h"
#include "test_scale.h"
#include "test_matvec.h"
#include "test_upper_diag.h"
#include "test_etree.h"
#include "test_postorder.h"
#include "test_col_counts.h"

int main(){
    MU_INITIALIZE();
    MU_RUN_SUITE(suite_scale);
    MU_RUN_SUITE(suite_matvec);
    MU_RUN_SUITE(suite_upper_diag);
    MU_RUN_SUITE(suite_etree);
    MU_RUN_SUITE(suite_postorder);
    MU_RUN_SUITE(suite_col_counts);
    MU_REPORT();
    return minunit_fail;
}