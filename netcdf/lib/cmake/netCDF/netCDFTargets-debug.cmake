#----------------------------------------------------------------
# Generated CMake target import file for configuration "DEBUG".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "netcdf" for configuration "DEBUG"
set_property(TARGET netcdf APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(netcdf PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libnetcdf.15.dylib"
  IMPORTED_SONAME_DEBUG "@rpath/libnetcdf.15.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS netcdf )
list(APPEND _IMPORT_CHECK_FILES_FOR_netcdf "${_IMPORT_PREFIX}/lib/libnetcdf.15.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
