
if (RSLQR_LINALG_LIBRARY STREQUAL "MKL")
  find_package(MKL CONFIG REQUIRED)
endif()

if (RSLQR_LINALG_LIBRARY STREQUAL "BLAS")
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)  # supresses warning from LAPACK library
  set(CBLAS "ON")
  set(LAPACKE "ON")
  list(APPEND CMAKE_MESSAGE_CONTEXT LAPACK)
  FetchContent_MakeAvailable(lapackelib)
  list(POP_BACK CMAKE_MESSAGE_CONTEXT LAPACK)
  add_library(LAPACK::CBLAS ALIAS cblas)
  add_library(LAPACK::LAPACKE ALIAS lapacke)
endif()

find_package(Eigen3 3.3 NO_MODULE)
if (Eigen3_FOUND)
  message(STATUS "Found Eigen3 library.")
elseif(RSLQR_LINALG_LIBRARY STREQUAL "Eigen")
  message(FATAL_ERROR "Couldn't find Eigen library, but the user specified it.")
endif()