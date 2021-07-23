# Set default project to unknown
if(NOT PROJECT)
  message(STATUS "Setting CCPP project to 'unknown' as none was specified.")
  set(PROJECT "Unknown")
endif (NOT PROJECT)

#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.0)

project(ccpp_physics
        VERSION 4.0.0
        LANGUAGES C CXX Fortran)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif(POLICY CMP0042)

# CMP0057: Support new IN_LIST if() operator
if(POLICY CMP0057)
    cmake_policy(SET CMP0057 NEW)
endif(POLICY CMP0057)

#------------------------------------------------------------------------------
set(PACKAGE "ccpp-physics")
set(AUTHORS "Grant Firl" "Dom Heinzeller" "Man Zhang" "Laurie Carson")

#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  #set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  #set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Enable OpenMP support for C/C++/Fortran compiler")
else (OPENMP)
  message (STATUS "Disable OpenMP support for C/C++/Fortran compiler")
endif()

#------------------------------------------------------------------------------
# The Fortran compiler/linker flag inserted by cmake to create shared libraries
# with the Intel compiler is deprecated (-i_dynamic), correct here.
# CMAKE_Fortran_COMPILER_ID = {"Intel", "PGI", "GNU", "Clang", "MSVC", ...}
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")
endif()

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Bitforbit" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# Request a static build
option(BUILD_SHARED_LIBS "Build a shared library" OFF)

#------------------------------------------------------------------------------
# Set the sources: physics type definitions
set(TYPEDEFS $ENV{CCPP_TYPEDEFS})
if(TYPEDEFS)
  message(STATUS "Got CCPP TYPEDEFS from environment variable: ${TYPEDEFS}")
else(TYPEDEFS)
  include(${CMAKE_CURRENT_BINARY_DIR}/CCPP_TYPEDEFS.cmake)
  message(STATUS "Got CCPP TYPEDEFS from cmakefile include file: ${TYPEDEFS}")
endif(TYPEDEFS)

# Generate list of Fortran modules from the CCPP type
# definitions that need need to be installed
foreach(typedef_module ${TYPEDEFS})
    list(APPEND MODULES_F90 ${CMAKE_CURRENT_BINARY_DIR}/${typedef_module})
endforeach()

#------------------------------------------------------------------------------
# Set the sources: physics schemes
set(SCHEMES $ENV{CCPP_SCHEMES})
if(SCHEMES)
  message(STATUS "Got CCPP SCHEMES from environment variable: ${SCHEMES}")
else(SCHEMES)
  include(${CMAKE_CURRENT_BINARY_DIR}/CCPP_SCHEMES.cmake)
  message(STATUS "Got CCPP SCHEMES from cmakefile include file: ${SCHEMES}")
endif(SCHEMES)

# Set the sources: physics scheme caps
set(CAPS $ENV{CCPP_CAPS})
if(CAPS)
  message(STATUS "Got CCPP CAPS from environment variable: ${CAPS}")
else(CAPS)
  include(${CMAKE_CURRENT_BINARY_DIR}/CCPP_CAPS.cmake)
  message(STATUS "Got CCPP CAPS from cmakefile include file: ${CAPS}")
endif(CAPS)

# Create empty lists for schemes with special compiler optimization flags
set(SCHEMES_SFX_OPT "")
# Create empty lists for schemes with special floating point precision flags
set(SCHEMES_SFX_PREC "")
# Create a duplicate of the SCHEMES list for handling floating point precision flags
set(SCHEMES2 ${SCHEMES})

# Schemes and caps from the CCPP code generator use full paths with symlinks
# resolved, we need to do the same here for the below logic to work
get_filename_component(FULL_PATH_TO_CMAKELISTS CMakeLists.txt REALPATH BASE_DIR ${LOCAL_CURRENT_SOURCE_DIR})
get_filename_component(LOCAL_CURRENT_SOURCE_DIR ${FULL_PATH_TO_CMAKELISTS} DIRECTORY)

# List of files that need to be compiled without OpenMP
set(SCHEMES_OPENMP_OFF ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_gas_optics.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_rrtmgp_constants.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_rrtmgp_util_reorder.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_gas_concentrations.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/kernels-openacc/mo_gas_optics_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_rrtmgp_util_string.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/kernels/mo_gas_optics_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rrtmgp/mo_gas_optics_rrtmgp.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/tests/mo_testing_io.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/tests/clear_sky_regression.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_rrtmgp_clr_all_sky.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_fluxes_byband_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_fluxes_byband.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/solar_variability/mo_solar_variability.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_heating_rates.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_fluxes_bygpoint.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/mo_compute_bc.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/cloud_optics/mo_cloud_sampling.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/extensions/cloud_optics/mo_cloud_optics.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/examples/mo_load_coefficients.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/examples/rfmip-clear-sky/rrtmgp_rfmip_sw.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/examples/rfmip-clear-sky/mo_rfmip_io.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/examples/rfmip-clear-sky/rrtmgp_rfmip_lw.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/examples/mo_simple_netcdf.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/examples/all-sky/rrtmgp_allsky.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/examples/all-sky/mo_load_cloud_coefficients.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/examples/all-sky/mo_garand_atmos_io.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_config.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_source_functions.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_sw.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_fluxes.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_lw.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/kernels-openacc/mo_rte_solver_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/kernels-openacc/mo_optical_props_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_util_array.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/kernels/mo_rte_solver_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/kernels/mo_optical_props_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/kernels/mo_fluxes_broadband_kernels.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_rte_kind.F90
                       ${LOCAL_CURRENT_SOURCE_DIR}/physics/rte-rrtmgp/rte/mo_optical_props.F90)

#------------------------------------------------------------------------------
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

  SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/module_bfmicrophysics.f
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/sflx.f
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/sfc_diff.f
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/sfc_diag.f
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/module_nst_model.f90
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/calpreciptype.f90
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/mersenne_twister.f
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/module_nst_water_prop.f90
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/aer_cloud.F
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/wv_saturation.F
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/cldwat2m_micro.F
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/surface_perturbation.F90
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/module_mp_thompson_make_number_concentrations.F90
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/module_SF_JSFC.F90
                              ${LOCAL_CURRENT_SOURCE_DIR}/physics/module_BL_MYJPBL.F90
                              PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fdefault-double-8 ${OpenMP_Fortran_FLAGS}")

  if (PROJECT STREQUAL "CCPP-FV3")
    # Set 32-bit floating point precision flags for certain files
    # that are executed in the dynamics (fast physics part)
    if (DYN32)
      if (${LOCAL_CURRENT_SOURCE_DIR}/physics/gfdl_fv_sat_adj.F90 IN_LIST SCHEMES)
      # Reduce floating point precision from 64-bit to 32-bit, if necessary
      set(CMAKE_Fortran_FLAGS_PREC32 ${CMAKE_Fortran_FLAGS_DEFAULT_PREC})
      string(REPLACE "-fdefault-real-8" ""
             CMAKE_Fortran_FLAGS_PREC32 "${CMAKE_Fortran_FLAGS_PREC32}")
      string(REPLACE "-fdefault-double-8" ""
             CMAKE_Fortran_FLAGS_PREC32 "${CMAKE_Fortran_FLAGS_PREC32}")
      SET_PROPERTY(SOURCE ${LOCAL_CURRENT_SOURCE_DIR}/physics/gfdl_fv_sat_adj.F90
                   APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_PREC32} ${OpenMP_Fortran_FLAGS} ")
      # Add all of the above files to the list of schemes with special floating point precision flags
      list(APPEND SCHEMES_SFX_PREC ${LOCAL_CURRENT_SOURCE_DIR}/physics/gfdl_fv_sat_adj.F90)
      endif()
    endif()

    # Remove files with special floating point precision flags from list
    # of files with standard floating point precision flags
    if (SCHEMES_SFX_PREC)
      list(REMOVE_ITEM SCHEMES2 ${SCHEMES_SFX_PREC})
    endif ()

    if (PROJECT STREQUAL "CCPP-FV3")
      # Remove files that need to be compiled without OpenMP from list
      # of files with standard compiler flags, and assign no-OpenMP flags
      if (SCHEMES_OPENMP_OFF)
        list(REMOVE_ITEM SCHEMES2 ${SCHEMES_OPENMP_OFF})
      endif ()
      # Assign standard floating point precision flags to all remaining schemes and caps
      SET_PROPERTY(SOURCE ${SCHEMES_OPENMP_OFF}
                   APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_DEFAULT_PREC} ")
    endif()

    # Assign standard floating point precision flags to all remaining schemes and caps
    SET_PROPERTY(SOURCE ${SCHEMES2} ${CAPS}
                 APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_DEFAULT_PREC} ${OpenMP_Fortran_FLAGS} ")

  endif (PROJECT STREQUAL "CCPP-FV3")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  # Adjust settings for bit-for-bit reproducibility of NEMSfv3gfs
  if (PROJECT STREQUAL "CCPP-FV3")

    if (${LOCAL_CURRENT_SOURCE_DIR}/physics/module_sf_mynn.F90 IN_LIST SCHEMES)
      # Reduce optimization for module_sf_mynn.F90 (to avoid an apparent compiler bug with Intel 18 on Hera)
      SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/module_sf_mynn.F90
                                  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS_OPT} -O1")
      list(APPEND SCHEMES_SFX_OPT ${LOCAL_CURRENT_SOURCE_DIR}/physics/module_sf_mynn.F90)
    endif()

    if (${LOCAL_CURRENT_SOURCE_DIR}/physics/radiation_aerosols.f IN_LIST SCHEMES)
      # Replace -xHost or -xCORE-AVX2 with -xCORE-AVX-I for certain files
      set(CMAKE_Fortran_FLAGS_LOPT1 ${CMAKE_Fortran_FLAGS_OPT})
      string(REPLACE "-xHOST" "-xCORE-AVX-I"
             CMAKE_Fortran_FLAGS_LOPT1
             "${CMAKE_Fortran_FLAGS_LOPT1}")
      string(REPLACE "-xCORE-AVX2" "-xCORE-AVX-I"
             CMAKE_Fortran_FLAGS_LOPT1
             "${CMAKE_Fortran_FLAGS_LOPT1}")
      string(REPLACE "-axSSE4.2,CORE-AVX2" "-axSSE4.2,CORE-AVX-I"
             CMAKE_Fortran_FLAGS_LOPT1
             "${CMAKE_Fortran_FLAGS_LOPT1}")
      SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/radiation_aerosols.f
                                  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS_LOPT1}")
      # Add all of the above files to the list of schemes with special compiler flags
      list(APPEND SCHEMES_SFX_OPT ${LOCAL_CURRENT_SOURCE_DIR}/physics/radiation_aerosols.f)
    endif()

    # Remove files with special compiler flags from list of files with standard compiler flags
    if (SCHEMES_SFX_OPT)
      list(REMOVE_ITEM SCHEMES ${SCHEMES_SFX_OPT})
    endif(SCHEMES_SFX_OPT)
    # Assign standard compiler flags to all remaining schemes and caps
    SET_SOURCE_FILES_PROPERTIES(${SCHEMES} ${CAPS}
                                PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS_OPT}")

    # Set 32-bit floating point precision flags for certain files
    # that are executed in the dynamics (fast physics part)
    if (DYN32)
      if (${LOCAL_CURRENT_SOURCE_DIR}/physics/gfdl_fv_sat_adj.F90 IN_LIST SCHEMES)
        # Reduce floating point precision from 64-bit to 32-bit, if necessary
        set(CMAKE_Fortran_FLAGS_PREC32 ${CMAKE_Fortran_FLAGS_DEFAULT_PREC})
        string(REPLACE "-real-size 64" "-real-size 32"
               CMAKE_Fortran_FLAGS_PREC32 "${CMAKE_Fortran_FLAGS_PREC32}")
        SET_PROPERTY(SOURCE ${LOCAL_CURRENT_SOURCE_DIR}/physics/gfdl_fv_sat_adj.F90
                     APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_PREC32} ${OpenMP_Fortran_FLAGS} ")
        # Add all of the above files to the list of schemes with special floating point precision flags
        list(APPEND SCHEMES_SFX_PREC ${LOCAL_CURRENT_SOURCE_DIR}/physics/gfdl_fv_sat_adj.F90)
      endif()
    endif()

    # Remove files with special floating point precision flags from list
    # of files with standard floating point precision flags flags
    if (SCHEMES_SFX_PREC)
      list(REMOVE_ITEM SCHEMES2 ${SCHEMES_SFX_PREC})
    endif (SCHEMES_SFX_PREC)

    # Remove files that need to be compiled without OpenMP from list
    # of files with standard compiler flags, and assign no-OpenMP flags
    if (SCHEMES_OPENMP_OFF)
      list(REMOVE_ITEM SCHEMES2 ${SCHEMES_OPENMP_OFF})
      # Assign standard floating point precision flags to all remaining schemes and caps
      SET_PROPERTY(SOURCE ${SCHEMES_OPENMP_OFF}
                   APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_DEFAULT_PREC} ")
    endif ()

    # Assign standard floating point precision flags to all remaining schemes and caps
    SET_PROPERTY(SOURCE ${SCHEMES2} ${CAPS}
                 APPEND_STRING PROPERTY COMPILE_FLAGS " ${CMAKE_Fortran_FLAGS_DEFAULT_PREC} ${OpenMP_Fortran_FLAGS} ")

  else (PROJECT STREQUAL "CCPP-FV3")
    SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/module_bfmicrophysics.f
                                ${LOCAL_CURRENT_SOURCE_DIR}/physics/sflx.f
                                ${LOCAL_CURRENT_SOURCE_DIR}/physics/sfc_diff.f
                                ${LOCAL_CURRENT_SOURCE_DIR}/physics/sfc_diag.f
                                PROPERTIES COMPILE_FLAGS "-r8 ${OpenMP_Fortran_FLAGS} ")
    SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/module_nst_model.f90
                                ${LOCAL_CURRENT_SOURCE_DIR}/physics/calpreciptype.f90
                                PROPERTIES COMPILE_FLAGS "-r8 -free ${OpenMP_Fortran_FLAGS} ")
    SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/mersenne_twister.f
                                PROPERTIES COMPILE_FLAGS "-r8 -ftz ${OpenMP_Fortran_FLAGS} ")
    SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/module_nst_water_prop.f90
                                PROPERTIES COMPILE_FLAGS "-extend-source 132 -r8 -free ${OpenMP_Fortran_FLAGS} ")
    SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/aer_cloud.F
                                ${LOCAL_CURRENT_SOURCE_DIR}/physics/wv_saturation.F
                                ${LOCAL_CURRENT_SOURCE_DIR}/physics/cldwat2m_micro.F
                                ${LOCAL_CURRENT_SOURCE_DIR}/physics/surface_perturbation.F90
                                PROPERTIES COMPILE_FLAGS "-r8 ${OpenMP_Fortran_FLAGS} ")
    SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/module_mp_thompson_make_number_concentrations.F90
                                PROPERTIES COMPILE_FLAGS "-r8 ${OpenMP_Fortran_FLAGS} ")
    SET_SOURCE_FILES_PROPERTIES(${LOCAL_CURRENT_SOURCE_DIR}/physics/module_SF_JSFC.F90
                                ${LOCAL_CURRENT_SOURCE_DIR}/physics/module_BL_MYJPBL.F90
                                PROPERTIES COMPILE_FLAGS "-r8 ${OpenMP_Fortran_FLAGS} ")
  endif (PROJECT STREQUAL "CCPP-FV3")

else()
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${GFS_PHYS_SRC}/CMakeLists.txt")
endif()

# The auto-generated caps can contain calls to physics schemes in
# which some of the arguments (pointers, arrays) are not associated/allocated.
# This is on purpose to avoid allocating fields that are not used inside the
# scheme if, for example, certain conditions are not met. To avoid
# Fortran runtime errors, it is necessary to remove checks for pointers
# that are not associated and for array bounds from the caps ONLY. For the
# physics schemes, these checks can and should remain enabled. Overwriting
# the pointer check flags explicitly works for Intel and GNU, but not for PGI.
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set_property(SOURCE ${CAPS} APPEND_STRING PROPERTY COMPILE_FLAGS " -fcheck=no-pointer,no-bounds ")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set_property(SOURCE ${CAPS} APPEND_STRING PROPERTY COMPILE_FLAGS " -check nopointers,nobounds ")
endif ()

#------------------------------------------------------------------------------
add_library(ccpp_physics STATIC ${SCHEMES} ${SCHEMES_SFX_OPT} ${CAPS})
# Generate list of Fortran modules from defined sources
foreach(source_f90 ${CAPS})
    get_filename_component(tmp_source_f90 ${source_f90} NAME)
    string(REGEX REPLACE ".F90" ".mod" tmp_module_f90 ${tmp_source_f90})
    string(TOLOWER ${tmp_module_f90} module_f90)
    list(APPEND MODULES_F90 ${CMAKE_CURRENT_BINARY_DIR}/${module_f90})
endforeach()

target_link_libraries(ccpp_physics PUBLIC NetCDF::NetCDF_Fortran)
    
set_target_properties(ccpp_physics PROPERTIES VERSION ${PROJECT_VERSION}
                                     SOVERSION ${PROJECT_VERSION_MAJOR})

if (PROJECT STREQUAL "CCPP-FV3")
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
endif (PROJECT STREQUAL "CCPP-FV3")
