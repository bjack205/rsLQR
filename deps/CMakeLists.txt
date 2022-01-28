include(FetchContent)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

if (RSLQR_LINALG_LIBRARY STREQUAL "BLAS")
  # Download LAPACK library
  FetchContent_Declare(lapackelib
    GIT_REPOSITORY https://github.com/Reference-LAPACK/lapack.git
    GIT_TAG c8cb980bb58d77738121d541acf063f89e5c5480 # version 3.10.0, commit with CMake fixes (by bjack205)
  )
endif()

# Add the CJSON library
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
FetchContent_Declare(cjson_lib
  GIT_REPOSITORY https://github.com/DaveGamble/cJSON.git
  GIT_TAG d348621ca93571343a56862df7de4ff3bc9b5667 # version 1.7.15
)
set(ENABLE_CJSON_TEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(cjson_lib)

# Copy the header file into build/include and add as a dependency to the target
set(CJSON_HEADER_FILE ${cjson_lib_SOURCE_DIR}/cJSON.h)
configure_file(${CJSON_HEADER_FILE} ${CMAKE_BINARY_DIR}/include/cjson/cJSON.h)
target_include_directories(cjson
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
)
add_library(CJSON::CJSON ALIAS cjson)