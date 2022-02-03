
add_library(rslqr 
  binary_tree.h
  binary_tree.c

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

  utils.h
  utils.c
)
target_link_libraries(rslqr
  PUBLIC
  matrix
  CJSON::CJSON
  lqr
  riccati

  PRIVATE
  OpenMP::OpenMP_C
)
set_property(TARGET rslqr PROPERTY C_STANDARD 11)
target_include_directories(rslqr
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include/rslqr>
)

get_target_property(NDLQR_INCLUDEDIR rslqr INCLUDE_DIRECTORIES)
cmake_print_variables(NDLQR_INCLUDEDIR)

# Prep for install
add_library(rsLQR::rslqr ALIAS rslqr)
set_target_properties(rslqr PROPERTIES
  OUTPUT_NAME rslqr  # sets the name of the installed library
  EXPORT_NAME rslqr  # sets the name of the imported target
)
install(TARGETS rslqr
  EXPORT rsLQR_export   # add to the export set
)