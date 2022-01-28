
cmake_minimum_required(VERSION 3.17)
list(APPEND CMAKE_MESSAGE_CONTEXT rsLQR)
project(rsLQR VERSION 0.1 LANGUAGES C CXX)

include(CMakePrintHelpers)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

##############################
# Options 
##############################
set(RSLQR_DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified. Setting CMAKE_BUILD_TYPE to ${RSLQR_DEFAULT_BUILD_TYPE}.")
  set(CMAKE_BUILD_TYPE ${RSLQR_DEFAULT_BUILD_TYPE} CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(RSLQR_LINALG_LIBRARY "InternalRoutines" CACHE STRING "Linear algebra library.")
set_property(CACHE RSLQR_LINALG_LIBRARY PROPERTY STRINGS BLAS MKL Eigen InternalRoutines)
option(RSLQR_USE_INTERNAL_LINALG "Use the internal linear algebra routines. Takes precedence over BLAS (default)." "OFF")
option(RSLQR_BUILD_TESTS "Build tests for ndLQR" "ON") 

if (RSLQR_LINALG_LIBRARY STREQUAL "MKL")
  set(RSLQR_USE_MKL ON) 
  message(STATUS "Using the Intel MKL library")
else()
  set(RSLQR_USE_MKL OFF) 
endif()

##############################
# Dependencies
##############################
add_subdirectory(deps)

# Find required packages
find_package(OpenMP REQUIRED)
include(FindLinearAlgebra)

enable_testing()

##############################
# Build
##############################
include_directories(${rsLQR_SOURCE_DIR}/src)
add_subdirectory(src)