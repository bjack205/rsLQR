include(CMakePackageConfigHelpers)

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/src
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rslqr
    FILES_MATCHING
    PATTERN "*.h"
)