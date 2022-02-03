add_library(riccati
  riccati_solver.h
  riccati_solver.c

  riccati_solve.h
  riccati_solve.c
)
target_link_libraries(riccati
  PUBLIC
  matrix
)

# Prep for install
add_library(rsLQR::riccati ALIAS riccati)
set_target_properties(riccati PROPERTIES
  OUTPUT_NAME riccati
  EXPORT_NAME riccati
)
install(TARGETS riccati
  EXPORT rsLQR_export
)