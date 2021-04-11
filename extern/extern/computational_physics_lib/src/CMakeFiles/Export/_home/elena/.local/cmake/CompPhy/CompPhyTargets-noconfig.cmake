#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "compphy-lib" for configuration ""
set_property(TARGET compphy-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(compphy-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "C"
  IMPORTED_LOCATION_NOCONFIG "/home/elena/.local/lib/libcompphy-lib.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS compphy-lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_compphy-lib "/home/elena/.local/lib/libcompphy-lib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
