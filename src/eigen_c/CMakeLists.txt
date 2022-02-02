add_library(eigen_c SHARED
  eigen_c.h
  eigen_c.cpp
)
target_compile_definitions(eigen_c
  PRIVATE
  EIGEN_DONT_PARALLELIZE
)
target_link_libraries(eigen_c
  PRIVATE
  Eigen3::Eigen
)
add_library(rsLQR::eigen_c ALIAS eigen_c)
set_target_properties(eigen_c PROPERTIES
  OUTPUT_NAME eigen_c
  EXPORT_NAME eigen_c
)
install(TARGETS eigen_c
  EXPORT rsLQR_export
)