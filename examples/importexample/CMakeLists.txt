cmake_minimum_required(VERSION 3.17)
project(ImportExample VERSION 0.1 LANGUAGES C)

include(CMakePrintHelpers)

##############################
# Import the package 
##############################
include(FetchContent)
# Clones the repo in your build folder
FetchContent_Declare(librslqr
    GIT_REPOSITORY https://github.com/bjack205/rsLQR.git
    # NOTE: Replace this with whichever commit you want to use 
    GIT_TAG cda47a585079f8fed5311738aac9db766d21e135 
)
# Includes the root CMakeLists.txt of the library, importing all of it's targets
# to your build system
FetchContent_MakeAvailable(librslqr)

##############################
# Build the executable 
##############################
add_executable(importexample
    main.c
)
# Link against the imported library target
target_link_libraries(importexample
    rsLQR::rslqr
)

# Add a compile definition to provide the location of the JSON file
target_compile_definitions(importexample
    PRIVATE
    LQRPROBFILE="${PROJECT_SOURCE_DIR}/../../lqr_prob.json"
)