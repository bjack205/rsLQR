cmake_minimum_required(VERSION 3.17)
project(InstallExample VERSION 0.1 LANGUAGES C)

include(CMakePrintHelpers)

##############################
# Find the package 
##############################
# Note you can use CMAKE_PREFIX_PATH to provide a search path for the library
#   if it's not installed in a standard location
# Alternatively, you can provide guesses via the HINTS or PATHS arguments to find_package
find_package(rsLQR)
if (rsLQR_FOUND)
    message(STATUS "rsLQR package found!")
endif()


##############################
# Build the executable 
##############################
add_executable(installexample
    main.c
)
# Link against the imported library target
target_link_libraries(installexample
    rsLQR::rslqr
)
# Add a compile definition to provide the location of the JSON file
target_compile_definitions(installexample
    PRIVATE
    LQRPROBFILE="${PROJECT_SOURCE_DIR}/../../lqr_prob.json"
)