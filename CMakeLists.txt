cmake_minimum_required(VERSION 3.3)

project(ccpp_physics
        VERSION 5.0.0
        LANGUAGES Fortran)

#------------------------------------------------------------------------------
set(PACKAGE "ccpp-physics")
set(AUTHORS "Grant Firl" "Dom Heinzeller" "Man Zhang" "Mike Kavulich" "Chunxi Zhang")

#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  find_package(OpenMP REQUIRED)
endif()

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# Pass debug/release flag to Fortran files for preprocessor
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-DDEBUG)
endif()

#------------------------------------------------------------------------------
# Request a static build
option(BUILD_SHARED_LIBS "Build a shared library" OFF)

#------------------------------------------------------------------------------
# Set the sources: physics type definitions
set(TYPEDEFS $ENV{CCPP_TYPEDEFS})
if(TYPEDEFS)
  message(STATUS "Got CCPP TYPEDEFS from environment variable")
else(TYPEDEFS)
  include(${CMAKE_CURRENT_BINARY_DIR}/CCPP_TYPEDEFS.cmake)
  message(STATUS "Got CCPP TYPEDEFS from cmakefile include file")
endif(TYPEDEFS)
list(REMOVE_DUPLICATES TYPEDEFS)

# Generate list of Fortran modules from the CCPP type
# definitions that need need to be installed
foreach(typedef_module ${TYPEDEFS})
    list(APPEND MODULES_F90 ${CMAKE_CURRENT_BINARY_DIR}/${typedef_module})
endforeach()

#------------------------------------------------------------------------------
# Set the sources: physics schemes
set(SCHEMES $ENV{CCPP_SCHEMES})
if(SCHEMES)
  message(STATUS "Got CCPP SCHEMES from environment variable")
else(SCHEMES)
  include(${CMAKE_CURRENT_BINARY_DIR}/CCPP_SCHEMES.cmake)
  message(STATUS "Got CCPP SCHEMES from cmakefile include file")
endif(SCHEMES)
list(REMOVE_DUPLICATES SCHEMES)

# Set the sources: physics scheme caps
set(CAPS $ENV{CCPP_CAPS})
if(CAPS)
  message(STATUS "Got CCPP CAPS from environment variable")
else(CAPS)
  include(${CMAKE_CURRENT_BINARY_DIR}/CCPP_CAPS.cmake)
  message(STATUS "Got CCPP CAPS from cmakefile include file")
endif(CAPS)
list(REMOVE_DUPLICATES CAPS)

# Schemes and caps from the CCPP code generator use full paths with symlinks
# resolved, we need to do the same here for the below logic to work
get_filename_component(FULL_PATH_TO_CMAKELISTS CMakeLists.txt REALPATH BASE_DIR ${LOCAL_CURRENT_SOURCE_DIR})
get_filename_component(LOCAL_CURRENT_SOURCE_DIR ${FULL_PATH_TO_CMAKELISTS} DIRECTORY)

#------------------------------------------------------------------------------

# List of files that need to be compiled without OpenMP
set(SCHEMES_OPENMP_OFF ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_gas_optics.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_rrtmgp_constants.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_gas_concentrations.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_rrtmgp_util_string.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/kernels/mo_gas_optics_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_gas_optics_rrtmgp.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_rrtmgp_clr_all_sky.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_fluxes_byband.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/solar_variability/mo_solar_variability.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_heating_rates.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_fluxes_bygpoint.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_compute_bc.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/cloud_optics/mo_cloud_sampling.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/cloud_optics/mo_cloud_optics.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_config.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_source_functions.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_sw.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_fluxes.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_lw.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_util_array.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/kernels/mo_rte_solver_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/kernels/mo_optical_props_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/kernels/mo_fluxes_broadband_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_kind.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_optical_props.F90)

# List of files that need to be compiled with different precision
set(SCHEMES_DYNAMICS)

if(${LOCAL_CURRENT_SOURCE_DIR}/physics/fv_sat_adj.F90 IN_LIST SCHEMES)
  list(APPEND SCHEMES_DYNAMICS ${LOCAL_CURRENT_SOURCE_DIR}/physics/fv_sat_adj.F90)
endif()

# Remove files that need to be compiled with different precision
# of files with standard compiler flags, and assign OpenMP flags
if(SCHEMES_DYNAMICS)
  SET_PROPERTY(SOURCE ${SCHEMES_DYNAMICS}
               APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_DYNAMICS} ${OpenMP_Fortran_FLAGS}")
  list(REMOVE_ITEM SCHEMES ${SCHEMES_DYNAMICS})
endif()

# Remove files that need to be compiled without OpenMP from list
# of files with standard compiler flags, and assign no-OpenMP flags
if(SCHEMES_OPENMP_OFF)
  SET_PROPERTY(SOURCE ${SCHEMES_OPENMP_OFF}
               APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_PHYSICS}")
  list(REMOVE_ITEM SCHEMES ${SCHEMES_OPENMP_OFF})
endif()

# Assign standard floating point precision flags to all remaining schemes and caps
SET_PROPERTY(SOURCE ${SCHEMES} ${CAPS}
             APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_PHYSICS} ${OpenMP_Fortran_FLAGS}")

# Lower optimization for certain schemes when compiling with Intel in Release mode
if(CMAKE_BUILD_TYPE STREQUAL "Release" AND ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  # Define a list of schemes that need lower optimization with Intel in Release mode
  set(SCHEME_NAMES_LOWER_OPTIMIZATION module_sf_mynn.F90)
  foreach(SCHEME_NAME IN LISTS SCHEME_NAMES_LOWER_OPTIMIZATION)
    set(SCHEMES_TMP ${SCHEMES})
    # Need to determine the name of the scheme with its path
    list(FILTER SCHEMES_TMP INCLUDE REGEX ".*${SCHEME_NAME}$")
    SET_SOURCE_FILES_PROPERTIES(${SCHEMES_TMP}
                                APPEND_STRING PROPERTY COMPILE_FLAGS
                                " ${CMAKE_Fortran_FLAGS_PHYSICS} ${OpenMP_Fortran_FLAGS} -O1")
  endforeach()
endif()

# No optimization for certain schemes when compiling with Intel in Release mode
if(CMAKE_BUILD_TYPE STREQUAL "Release" AND ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  # Define a list of schemes that can't be optimized with Intel in Release mode
  set(SCHEME_NAMES_NO_OPTIMIZATION GFS_typedefs.F90)
  foreach(SCHEME_NAME IN LISTS SCHEME_NAMES_NO_OPTIMIZATION)
    set(SCHEMES_TMP ${SCHEMES})
    # Need to determine the name of the scheme with its path
    list(FILTER SCHEMES_TMP INCLUDE REGEX ".*${SCHEME_NAME}$")
    SET_SOURCE_FILES_PROPERTIES(${SCHEMES_TMP}
                                APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_PHYSICS} ${OpenMP_Fortran_FLAGS} -O0")
  endforeach()
endif()

# Reduce optimization for mo_gas_optics_kernels.F90 (to avoid an apparent compiler bug with Intel 19+)
if(${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/kernels/mo_gas_optics_kernels.F90 IN_LIST SCHEMES_OPENMP_OFF AND
   CMAKE_BUILD_TYPE STREQUAL "Release" AND ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/kernels/mo_gas_optics_kernels.F90
                              APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_PHYSICS} -O1")
endif()

#------------------------------------------------------------------------------

add_library(ccpp_physics STATIC ${SCHEMES} ${SCHEMES_OPENMP_OFF} ${SCHEMES_DYNAMICS} ${CAPS})
# Generate list of Fortran modules from defined sources
foreach(source_f90 ${CAPS})
    get_filename_component(tmp_source_f90 ${source_f90} NAME)
    string(REGEX REPLACE ".F90" ".mod" tmp_module_f90 ${tmp_source_f90})
    string(TOLOWER ${tmp_module_f90} module_f90)
    list(APPEND MODULES_F90 ${CMAKE_CURRENT_BINARY_DIR}/${module_f90})
endforeach()

set_target_properties(ccpp_physics PROPERTIES VERSION ${PROJECT_VERSION}
                                     SOVERSION ${PROJECT_VERSION_MAJOR})

target_include_directories(ccpp_physics PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

target_link_libraries(ccpp_physics PUBLIC w3emc::w3emc_d
                                          sp::sp_d
                                          NetCDF::NetCDF_Fortran)

# Define where to install the library
install(TARGETS ccpp_physics
        EXPORT ccpp_physics-targets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION lib
)
# Export our configuration
install(EXPORT ccpp_physics-targets
        FILE ccpp_physics-config.cmake
        DESTINATION lib/cmake
)
# Define where to install the C headers and Fortran modules
#install(FILES ${HEADERS_C} DESTINATION include)
install(FILES ${MODULES_F90} DESTINATION include)
