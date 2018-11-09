cmake_minimum_required(VERSION 3.12)
project(BiCGSTAB)

set(CMAKE_CXX_STANDARD 11)

SET(GCC_COVERAGE_COMPILE_FLAGS "-fopenmp")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

add_executable(BiCGSTAB
        main_polus.cpp)
