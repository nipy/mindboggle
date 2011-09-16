# The full VTK version number.
SET(PACKAGE_VERSION "5.6.1")

# This version is compatible only with matching major.minor versions.
IF("5.6" VERSION_EQUAL "${PACKAGE_FIND_VERSION_MAJOR}.${PACKAGE_FIND_VERSION_MINOR}")
  # This version is compatible with equal or lesser patch versions.
  IF(NOT "1" VERSION_LESS "${PACKAGE_FIND_VERSION_PATCH}")
    SET(PACKAGE_VERSION_COMPATIBLE 1)
    IF("1" VERSION_EQUAL "${PACKAGE_FIND_VERSION_PATCH}")
      SET(PACKAGE_VERSION_EXACT 1)
    ENDIF()
  ENDIF()
ENDIF()
