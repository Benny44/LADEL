#include "minunit.h"
#include "test_scale.h"
#include "test_matvec.h"

int main(){
    MU_INITIALIZE();
    MU_RUN_SUITE(suite_scale);
    MU_RUN_SUITE(suite_matvec);
    MU_REPORT();
    return minunit_fail;
}