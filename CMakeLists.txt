cmake_minimum_required(VERSION 3.0.2)
project(LADEL VERSION 0.1.0)

############################################################
################## Configure directories ###################
############################################################
SET(LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/../lib)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

############################################################
################### Compilation options ####################
############################################################
option (UNITTESTS "Perform unit testing" ON)
option (SIMPLE_COL_COUNTS "Simplify symbolic analysis (worse asymptotic time complexity)" OFF)

if (SIMPLE_COL_COUNTS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSIMPLE_COL_COUNTS")
endif (SIMPLE_COL_COUNTS) 

if (UNIX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wconversion")
endif (UNIX)
############################################################
#################### Build the library #####################
############################################################

SET(SRC
    src/scale.c
    src/matvec.c
    src/upper_diag.c
    src/etree.c
    src/postorder.c
    src/transpose.c
    src/col_counts.c)

add_library(ladel ${SRC})
target_include_directories(ladel
    PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)

############################################################
#################### Build executables #####################
############################################################
if (UNITTESTS)
    add_subdirectory(test)
    include(CTest)
    enable_testing()
    add_test(NAME ladel_tester COMMAND $<TARGET_FILE:run_all_tests>)
endif()