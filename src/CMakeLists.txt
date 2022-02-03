if (Eigen3_FOUND)
  add_subdirectory(eigen_c)
endif()


add_subdirectory(matrix)
add_subdirectory(riccati)
add_subdirectory(lqr)
add_subdirectory(rslqr)