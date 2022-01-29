add_library(minunit
    minunit.h
    minunit.c
)

if (RSLQR_RUN_FULL_TEST)
  set(FULLTEST 1)
else()
  set(FULLTEST 0)
endif()

function (add_ndlqr_test name)
  set(TEST_NAME ${name}_test)
  add_executable(${TEST_NAME}
    ${TEST_NAME}.c

    test_problem.h
    test_problem.c
  )
  target_link_libraries(${TEST_NAME}
    PRIVATE
    ndlqr
    minunit
  )
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
  if (Eigen3_FOUND)
    target_link_libraries(${TEST_NAME}
      PRIVATE
      eigen_c
    )
    target_compile_definitions(${TEST_NAME}
      PRIVATE
      USE_EIGEN=1
    )
  endif()

  target_compile_definitions(${TEST_NAME}
    PRIVATE
    LQRDATAFILE="${rsLQR_SOURCE_DIR}/lqr_data.json"
    LQRPROBFILE="${rsLQR_SOURCE_DIR}/lqr_prob.json"
    LQRPROB256FILE="${rsLQR_SOURCE_DIR}/lqr_prob_256.json"
    SAMPLEPROBFILE="${rsLQR_SOURCE_DIR}/sample_problem.json"
    NTHREADS=${RSLQR_NTHREADS_TEST}
    FULLTEST=${FULLTEST}
  )
endfunction()

add_ndlqr_test(matrix)
add_ndlqr_test(linalg)
add_ndlqr_test(linalg_custom)
add_ndlqr_test(utils)
add_ndlqr_test(binarytree)
add_ndlqr_test(lqrdata)
add_ndlqr_test(nddata)
add_ndlqr_test(solver)
add_ndlqr_test(nested_dissection)
add_ndlqr_test(sample_problem)
add_ndlqr_test(riccati_solver)

# add_ndlqr_test(matmul)
# target_link_libraries(matmul_test
#   PRIVATE
#   OpenMP::OpenMP_C
# )

if (Eigen3_FOUND)
  add_ndlqr_test(cholesky)
  target_link_libraries(cholesky_test
    PRIVATE
    OpenMP::OpenMP_C
  )
endif()

add_ndlqr_test(parallel)
target_link_libraries(parallel_test
  PRIVATE
  OpenMP::OpenMP_C
)