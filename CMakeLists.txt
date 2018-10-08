# Set default project to unknown
if(NOT PROJECT)
  message(STATUS "Setting CCPP project to 'unknown' as none was specified.")
  set(PROJECT "Unknown")
endif (NOT PROJECT)

#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.11)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
    project(ccppphys VERSION 2.0.0)
else(POLICY CMP0048)
    project(ccppphys)
    set(PROJECT_VERSION 2.0.0)
    set(PROJECT_VERSION_MAJOR 2)
    set(PROJECT_VERSION_MINOR 0)
    set(PROJECT_VERSION_PATCH 0)
endif(POLICY CMP0048)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif(POLICY CMP0042)

#------------------------------------------------------------------------------
set(PACKAGE "ccpp-physics")
set(AUTHORS  "Grant J. Firl" "Dom Heinzeller")

#------------------------------------------------------------------------------
# Enable Fortran
enable_language(Fortran)

if (PROJECT STREQUAL "CCPP-SCM")
  #------------------------------------------------------------------------------
  # CMake Modules
  # Set the CMake module path
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../ccpp-framework/cmake")
endif (PROJECT STREQUAL "CCPP-SCM")

#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
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
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# By default we want a shared library (unless a static build is requested)
if(STATIC)
  option(BUILD_SHARED_LIBS "Build a static library" OFF)
else(STATIC)
  option(BUILD_SHARED_LIBS "Build a shared library" ON)
endif(STATIC)

#------------------------------------------------------------------------------
# Add the CCPP include/module directory
set(CCPP_INCLUDE_DIRS "" CACHE FILEPATH "Path to ccpp includes")
set_property(DIRECTORY PROPERTY INCLUDE_DIRECTORIES ${CCPP_INCLUDE_DIRS})

#------------------------------------------------------------------------------
# Add the CCPP library
set(CCPP_LIB_DIRS "" CACHE FILEPATH "Path to ccpp library")
link_directories(${CCPP_LIB_DIRS})
list(APPEND LIBS "ccpp")

#------------------------------------------------------------------------------
# Set the sources: physics schemes
include(./CCPP_SCHEMES.cmake)
# Set the sources: physics scheme caps
include(./CCPP_CAPS.cmake)

#------------------------------------------------------------------------------
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  SET_SOURCE_FILES_PROPERTIES(./physics/module_bfmicrophysics.f ./physics/rascnvv2.f ./physics/sflx.f ./physics/sfc_diff.f ./physics/sfc_diag.f PROPERTIES COMPILE_FLAGS -fdefault-real-8)
  SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_model.f90 ./physics/calpreciptype.f90 PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -ffree-form")
  SET_SOURCE_FILES_PROPERTIES(./physics/mersenne_twister.f PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fno-range-check")
  SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_water_prop.f90 PROPERTIES COMPILE_FLAGS "-ffree-line-length-none -fdefault-real-8 -ffree-form")
  SET_SOURCE_FILES_PROPERTIES(./physics/aer_cloud.F ./physics/wv_saturation.F ./physics/cldwat2m_micro.F ./physics/surface_perturbation.F90 PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fdefault-double-8")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  # Adjust settings for bit-for-bit reproducibility of NEMSfv3gfs
  if (PROJECT STREQUAL "CCPP-FV3")
    SET_SOURCE_FILES_PROPERTIES(./physics/module_bfmicrophysics.f
                                ./physics/rascnvv2.f
                                ./physics/sflx.f
                                ./physics/sfc_diff.f
                                ./physics/sfc_diag.f
                                ./physics/module_nst_model.f90
                                ./physics/calpreciptype.f90
                                ./physics/mersenne_twister.f
                                ./physics/module_nst_water_prop.f90
                                ./physics/aer_cloud.F
                                ./physics/wv_saturation.F
                                ./physics/cldwat2m_micro.F
                                ./physics/surface_perturbation.F90
                                PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -ftz")
    # Replace -xHost or -xCORE-AVX2 with -xCORE-AVX-I for certain files (following FV3/gfsphysics/makefile)
    set(CMAKE_Fortran_FLAGS_LOPT ${CMAKE_Fortran_FLAGS})
    string(REPLACE "xHost" "xCORE-AVX-I"
           CMAKE_Fortran_FLAGS_LOPT
           "${CMAKE_Fortran_FLAGS_LOPT}")
    string(REPLACE "xCORE-AVX2" "xCORE-AVX-I"
           CMAKE_Fortran_FLAGS_LOPT
           "${CMAKE_Fortran_FLAGS_LOPT}")
    SET_SOURCE_FILES_PROPERTIES(./physics/radiation_aerosols.f
                                PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS_LOPT} -r8 -ftz")
  else (PROJECT STREQUAL "CCPP-FV3")
    SET_SOURCE_FILES_PROPERTIES(./physics/module_bfmicrophysics.f ./physics/rascnvv2.f ./physics/sflx.f ./physics/sfc_diff.f ./physics/sfc_diag.f PROPERTIES COMPILE_FLAGS -r8)
    SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_model.f90 ./physics/calpreciptype.f90 PROPERTIES COMPILE_FLAGS "-r8 -free")
    SET_SOURCE_FILES_PROPERTIES(./physics/mersenne_twister.f PROPERTIES COMPILE_FLAGS "-r8 -ftz")
    SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_water_prop.f90 PROPERTIES COMPILE_FLAGS "-extend-source 132 -r8 -free")
    SET_SOURCE_FILES_PROPERTIES(./physics/aer_cloud.F ./physics/wv_saturation.F ./physics/cldwat2m_micro.F ./physics/surface_perturbation.F90 PROPERTIES COMPILE_FLAGS "-r8")
  endif (PROJECT STREQUAL "CCPP-FV3")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
  SET_SOURCE_FILES_PROPERTIES(./physics/module_bfmicrophysics.f ./physics/rascnvv2.f ./physics/sflx.f ./physics/sfc_diff.f ./physics/sfc_diag.f PROPERTIES COMPILE_FLAGS -r8)
  SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_model.f90 ./physics/calpreciptype.f90 PROPERTIES COMPILE_FLAGS "-r8 -Mfree")
  SET_SOURCE_FILES_PROPERTIES(./physics/mersenne_twister.f PROPERTIES COMPILE_FLAGS "-r8 -Mnofptrap")
  SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_water_prop.f90 PROPERTIES COMPILE_FLAGS "-r8 -Mfree")
  SET_SOURCE_FILES_PROPERTIES(./physics/aer_cloud.F ./physics/wv_saturation.F ./physics/cldwat2m_micro.F ./physics/surface_perturbation.F90 PROPERTIES COMPILE_FLAGS "-r8")
else (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${GFS_PHYS_SRC}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

# The auto-generated caps can contain calls to physics schemes in
# which some of the arguments (pointers) are not associated. This is
# on purpose to avoid allocating fields that are not used inside the
# scheme if, for example, certain conditions are not met. To avoid
# Fortran runtime errors, it is necessary to remove checks for pointers
# that are not associated from the caps ONLY. For the physics schemes,
# these checks can and should remain enabled. Overwriting the check flags
# explicitly works for Intel and GNU, but not for PGI.
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set_property(SOURCE ${CAPS} PROPERTY COMPILE_FLAGS "-fcheck=no-pointer")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set_property(SOURCE ${CAPS} PROPERTY COMPILE_FLAGS "-check nopointers")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
  if (CMAKE_Fortran_FLAGS MATCHES ".*chkptr.*")
    message (FATAL_ERROR "PGI compiler option chkptr cannot be used for CCPP physics")
  endif (CMAKE_Fortran_FLAGS MATCHES ".*chkptr.*")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

if (PROJECT STREQUAL "CCPP-FV3")
  link_directories(${NCEPLIBS_DIR}/lib)
elseif (PROJECT STREQUAL "CCPP-SCM")
  SET(W3LIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../external/w3nco/v2.0.6/src)
  SET(BACIOLIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../external/bacio/v2.0.1/src)
  SET(SPLIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/../external/sp/v2.0.2/src)

  #add "sibling" directories (must specify the build directory too)
  ADD_SUBDIRECTORY(${W3LIB_SRC} ${CMAKE_BINARY_DIR}/w3nco)
  ADD_SUBDIRECTORY(${BACIOLIB_SRC} ${CMAKE_BINARY_DIR}/bacio)
  ADD_SUBDIRECTORY(${SPLIB_SRC} ${CMAKE_BINARY_DIR}/sp)

  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/w3nco)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/sp)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/bacio)

  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp-framework/src)
endif (PROJECT STREQUAL "CCPP-FV3")

#------------------------------------------------------------------------------

if(STATIC)
  add_library(ccppphys STATIC ${SCHEMES} ${CAPS})
  # Generate list of Fortran modules from defined sources
  foreach(source_f90 ${CAPS})
      string(REGEX REPLACE ".F90" ".mod" tmp_module_f90 ${source_f90})
      string(TOLOWER ${tmp_module_f90} module_f90)
      list(APPEND MODULES_F90 ${CMAKE_CURRENT_BINARY_DIR}/../${module_f90})
  endforeach()
else(STATIC)
  add_library(ccppphys ${SCHEMES} ${CAPS})
endif(STATIC)

if (PROJECT STREQUAL "CCPP-FV3")
  target_link_libraries(ccppphys LINK_PUBLIC ${LIBS} sp_v2.0.2_d bacio_4 w3nco_d)
elseif (PROJECT STREQUAL "CCPP-SCM")
  target_link_libraries(ccppphys LINK_PUBLIC ${LIBS} w3 sp bacio)
endif (PROJECT STREQUAL "CCPP-FV3")
set_target_properties(ccppphys PROPERTIES VERSION ${PROJECT_VERSION}
                                     SOVERSION ${PROJECT_VERSION_MAJOR}
                                     COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                                     LINK_FLAGS "${CMAKE_Fortran_FLAGS}")

# DH* hack for PGI compiler: rename objects in scheme cap object files for ISO_C compliancy
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
set(CAPOBJS)
foreach(cap ${CAPS})
    string(REPLACE "_cap.F90" "_cap.F90.o" capobj "./${CMAKE_FILES_DIRECTORY}/ccppphys.dir/${cap}")
    list(APPEND CAPOBJS ${capobj})
endforeach(cap)

add_custom_command(TARGET ccppphys
                   PRE_LINK
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/pgifix.py --cmake ${CAPOBJS}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   COMMENT "Running pgifix_wrapper.py over all scheme caps")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
# *DH end hack for PGI compiler

if (PROJECT STREQUAL "CCPP-FV3")
  # Define where to install the library
  install(TARGETS ccppphys
          EXPORT ccppphys-targets
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib
          RUNTIME DESTINATION lib
  )
  # Export our configuration
  install(EXPORT ccppphys-targets
          FILE ccppphys-config.cmake
          DESTINATION lib/cmake
  )
  if(STATIC)
    # Define where to install the C headers and Fortran modules
    #install(FILES ${HEADERS_C} DESTINATION include)
    install(FILES ${MODULES_F90} DESTINATION include)
  endif(STATIC)
endif (PROJECT STREQUAL "CCPP-FV3")
