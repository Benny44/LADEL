cmake_minimum_required(VERSION 3.0.2)
project(LADEL VERSION 0.1.0)

############################################################
################## Configure directories ###################
############################################################
SET(LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/../lib)

############################################################
################### Compilation options ####################
############################################################
option (UNITTESTS "Perform unit testing" ON)


############################################################
#################### Build the library #####################
############################################################

SET(SRC
    src/scale.c)

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