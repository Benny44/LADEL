#include "minunit.h"
#include "test_scale.h"

int main(){
    MU_INITIALIZE();
    MU_RUN_SUITE(suite_scale);
    MU_REPORT();
    return minunit_fail;
}