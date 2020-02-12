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

option (LONG "Use long integers for indexing" ON)
if (LONG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDLONG")
endif (LONG)

if (FLOAT)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDFLOAT")
endif (FLOAT)

if (UNIX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wconversion")
endif (UNIX)

option (AMD "Use approximate minimum degree ordering" ON)
if (AMD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDAMD")
endif (AMD)
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
        src/col_counts.c
        src/ldl_symbolic.c
        src/permutation.c
        src/copy.c
        src/ldl_numeric.c
        src/pattern.c
        src/ladel.c)


if (AMD)
    SET(AMD_SRC
            amd/Source/amd_1.c    
            amd/Source/amd_control.c   
            amd/Source/amd_global.c  
            amd/Source/amd_postorder.c   
            amd/Source/amd_valid.c
            amd/Source/amd_2.c    
            amd/Source/amd_defaults.c  
            amd/Source/amd_info.c    
            amd/Source/amd_post_tree.c
            amd/Source/amd_aat.c  
            amd/Source/amd_dump.c      
            amd/Source/amd_order.c   
            amd/Source/amd_preprocess.c
            amd/Source/SuiteSparse_config.c)
    SET(SRC ${SRC} ${AMD_SRC})
endif (AMD)

add_library(ladel ${SRC})

if (AMD)
    target_include_directories(ladel
        PUBLIC  ${CMAKE_CURRENT_LIST_DIR}/include 
                ${CMAKE_CURRENT_LIST_DIR}/amd/Include)
else ()
    target_include_directories(ladel
        PUBLIC  ${CMAKE_CURRENT_LIST_DIR}/include)

endif (AMD)
############################################################
#################### Build executables #####################
############################################################
if (UNITTESTS)
    add_subdirectory(test)
    include(CTest)
    enable_testing()
    add_test(NAME ladel_tester COMMAND $<TARGET_FILE:run_all_tests>)
endif()