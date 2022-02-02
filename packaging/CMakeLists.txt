include(CMakePackageConfigHelpers)

# Add dependencies to export set 
install(
    TARGETS 
    cjson 
    EXPORT rsLQR_export
)

# Install targets
install(EXPORT rsLQR_export
    DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake
    NAMESPACE rsLQR::
    FILE rsLQR.cmake
    COMPONENT rsLQR
)

# Install header files
install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/src/   # final backslash copies contents of src/, not src/ itself
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rslqr
    FILES_MATCHING
    PATTERN "*.h"
)

# Write config and version file
install(FILES rsLQRConfig.cmake
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/rsLQRConfigVersion.cmake
  COMPATIBILITY AnyNewerVersion
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/rsLQRConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/rsLQRConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake
)