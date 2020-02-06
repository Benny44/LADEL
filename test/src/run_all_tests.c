#include "minunit.h"
#include "test_scale.h"
#include "test_matvec.h"
#include "test_upper_diag.h"

int main(){
    MU_INITIALIZE();
    MU_RUN_SUITE(suite_scale);
    MU_RUN_SUITE(suite_matvec);
    MU_RUN_SUITE(suite_upper_diag);
    MU_REPORT();
    return minunit_fail;
}