add_library(lqr
  lqr_data.h
  lqr_data.c
  
  lqr_problem.h
  lqr_problem.c

  json_utils.h
  json_utils.c
)
target_link_libraries(lqr
  PUBLIC
  matrix

  PRIVATE
  CJSON::CJSON
)

# Prep for install
add_library(rsLQR::lqr ALIAS lqr)
set_target_properties(lqr PROPERTIES
  OUTPUT_NAME lqr
  EXPORT_NAME lqr
)
install(TARGETS lqr
  EXPORT rsLQR_export
)