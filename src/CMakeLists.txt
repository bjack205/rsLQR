if (Eigen3_FOUND)
  add_subdirectory(eigen_c)
endif()


add_subdirectory(matrix)
add_subdirectory(riccati)
add_subdirectory(lqr)
add_subdirectory(rslqr)

add_library(rslqr_interface INTERFACE)
target_link_libraries(rslqr_interface
  INTERFACE
  rslqr
  riccati
  lqr
  matrix
)
add_library(rsLQR::rsLQR ALIAS rslqr_interface)
set_target_properties(rslqr_interface PROPERTIES
  EXPORT_NAME rsLQR
)
install(TARGETS rslqr_interface
  EXPORT rsLQR_export
)