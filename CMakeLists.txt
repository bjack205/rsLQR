
cmake_minimum_required(VERSION 3.17)
list(APPEND CMAKE_MESSAGE_CONTEXT rsLQR)
project(rsLQR VERSION 0.1 LANGUAGES C CXX)

include(CMakePrintHelpers)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

##############################
# Options 
##############################

# Handle default build type
set(RSLQR_DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified. Setting CMAKE_BUILD_TYPE to ${RSLQR_DEFAULT_BUILD_TYPE}.")
  set(CMAKE_BUILD_TYPE ${RSLQR_DEFAULT_BUILD_TYPE} CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Pick the linear algebra library
set(RSLQR_LINALG_LIBRARY "InternalRoutines" CACHE STRING "Linear algebra library.")
set_property(CACHE RSLQR_LINALG_LIBRARY PROPERTY STRINGS BLAS MKL Eigen InternalRoutines)

# Enable testing
option(RSLQR_BUILD_TESTS "Build tests for rsLQR" "ON") 

# Code Coverage
option(RSLQR_CODE_COVERAGE "Compile rsLQR with Code Coverage." OFF)
if(RSLQR_CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  message(STATUS "Compiling rsLQR with coverage info.")
  add_compile_options(-O0 -g -fprofile-arcs -ftest-coverage)
  add_link_options(-fprofile-arcs -ftest-coverage)
endif()

# Test threads
set(RSLQR_NTHREADS_TEST 2 CACHE STRING "Number of threads to use in tests.")
set_property(CACHE RSLQR_NTHREADS_TEST PROPERTY STRINGS 1 2 4 6 16)

# Run full test suite
option(RSLQR_RUN_FULL_TEST "Run the full, computationally intensive test suite." ON)

##############################
# Dependencies
##############################
add_subdirectory(deps)

# Find required packages
find_package(OpenMP REQUIRED)
include(FindLinearAlgebra)


##############################
# Build
##############################

include_directories(${rsLQR_SOURCE_DIR}/src)
add_subdirectory(src)

##############################
# Testing
##############################
if (RSLQR_BUILD_TESTS)
  enable_testing()
  include_directories(${rsLQR_SOURCE_DIR})
  add_subdirectory(test)
endif()