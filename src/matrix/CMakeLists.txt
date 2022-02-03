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