add_library(eigen_c SHARED
  eigen_c.h
  eigen_c.cpp
)
target_compile_definitions(eigen_c
  PRIVATE
  EIGEN_DONT_PARALLELIZE
)