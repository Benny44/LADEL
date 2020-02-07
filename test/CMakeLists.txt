cmake_minimum_required(VERSION 3.0.0)
project(LADELtest VERSION 0.1.0)

SET(TEST_SRC
    src/run_all_tests.c
    src/test_scale.c
    src/test_matvec.c
    src/test_upper_diag.c
    src/test_etree.c
    src/test_postorder.c
    src/test_transpose.c
    src/test_col_counts.c)

add_executable(run_all_tests
    ${TEST_SRC}
)

target_include_directories(run_all_tests
    PRIVATE include
)

target_link_libraries(run_all_tests ladel)

if(UNIX)
target_link_libraries(run_all_tests m)
endif(UNIX)