# Add Compile options
add_compile_options(-Wall -Wextra -pedantic -Werror -Wno-error=unknown-pragmas)
add_compile_options(-fPIC -mavx2 -mfma)

if (Eigen3_FOUND)
  add_subdirectory(eigen_c)
endif()

# Matrix / Linear Algebra library
add_library(matrix STATIC
  matrix.h
  matrix.c

  linalg.h
  linalg.c

  linalg_utils.h
  linalg_utils.c

  linalg_custom.h
  linalg_custom.c
)
target_link_libraries(matrix
  PUBLIC
  m  # math library
    
  PRIVATE
  OpenMP::OpenMP_C
)
if (RSLQR_LINALG_LIBRARY STREQUAL "Eigen")
  message(STATUS "Using C++ Eigen Library for linear algebra.")
  target_compile_definitions(matrix PUBLIC USE_EIGEN=1)
  target_link_libraries(matrix
    PRIVATE
    eigen_c
  )
elseif(RSLQR_LINALG_LIBRARY STREQUAL "InternalRoutines")
  message(STATUS "Using internal routines for linear algebra.")
  target_compile_definitions(matrix PUBLIC USE_CLAP=1)

elseif(RSLQR_LINALG_LIBRARY STREQUAL "MKL")
  message(STATUS "Using Intel MKL library for linear algebra.")
  target_compile_options(matrix 
    PUBLIC 
    $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>
  )
  target_include_directories(matrix 
    PUBLIC 
    $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>
  )
  target_link_libraries(matrix
    PUBLIC
    $<LINK_ONLY:MKL::MKL>
  )
  target_compile_definitions(matrix
    PUBLIC
    USE_MKL=1
  )
elseif(RSLQR_LINALG_LIBRARY STREQUAL "BLAS")
  message(STATUS "Using open-source BLAS for linear algebra.")
  target_link_libraries(matrix
    PRIVATE
    LAPACK::LAPACKE
    LAPACK::CBLAS
  )
  target_compile_definitions(matrix
    PUBLIC
    USE_BLAS=1
  )
endif()

# Prep for install 
add_library(rsLQR::matrix ALIAS matrix)
set_target_properties(matrix PROPERTIES
  OUTPUT_NAME matrix  # sets the name of the installed library
  EXPORT_NAME matrix  # sets the name of the imported target
)
install(TARGETS matrix
  EXPORT rsLQR_export   # add to the export set
)


add_library(ndlqr 
  ndlqr.h

  binary_tree.h
  binary_tree.c

  # matmul.h
  # matmul.c

  nddata.h
  nddata.c

  solver.h
  solver.c

  cholesky_factors.h
  cholesky_factors.c

  nested_dissection.h
  nested_dissection.c

  solve.h
  solve.c

  lqr_data.h
  lqr_data.c

  lqr_problem.h
  lqr_problem.c

  json_utils.h
  json_utils.c

  riccati_solver.h
  riccati_solver.c

  riccati_solve.h
  riccati_solve.c

  utils.h
  utils.c
)
target_link_libraries(ndlqr
  PUBLIC
  matrix
  CJSON::CJSON

  PRIVATE
  OpenMP::OpenMP_C
)
set_property(TARGET ndlqr PROPERTY C_STANDARD 11)
target_include_directories(ndlqr
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include/rslqr>
)

get_target_property(NDLQR_INCLUDEDIR ndlqr INCLUDE_DIRECTORIES)
cmake_print_variables(NDLQR_INCLUDEDIR)

# Prep for install
add_library(rsLQR::rslqr ALIAS ndlqr)
set_target_properties(ndlqr PROPERTIES
  OUTPUT_NAME rslqr  # sets the name of the installed library
  EXPORT_NAME rslqr  # sets the name of the imported target
)
install(TARGETS ndlqr
  EXPORT rsLQR_export   # add to the export set
)