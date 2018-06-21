cmake_minimum_required(VERSION 2.8.11)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif(POLICY CMP0042)

PROJECT(CCPP-FV3)
set(PROJECT "CCPP-FV3")

ENABLE_LANGUAGE(Fortran)

include(CMakeForceCompiler)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# Set include directories for make, also set CCPP_... for external projects
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/framework/src)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/physics)
SET(CCPP_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/framework/src;${CMAKE_BINARY_DIR}/physics")
SET(CCPP_LIB_DIRS "${CMAKE_BINARY_DIR}/framework/src;${CMAKE_BINARY_DIR}/physics")

#------------------------------------------------------------------------------
# Generate Compiler flags for C/CXX/Fortran - set to match NEMSfv3gfs flags
# for DEBUG and PROD mode for all three compilers

# Set basic compiler optimization levels, depending on build type
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  ADD_DEFINITIONS(-DDEBUG)
  set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS_DEBUG}")
  set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS_DEBUG}")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS_DEBUG}")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
  if (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    set (CMAKE_C_FLAGS_RELEASE       "-O1 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_CXX_FLAGS_RELEASE     "-O1 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_Fortran_FLAGS_RELEASE "-O1 -fPIC" CACHE STRING "" FORCE)
  else (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
      set (CMAKE_C_FLAGS_RELEASE       "-O2 -fPIC" CACHE STRING "" FORCE)
      set (CMAKE_CXX_FLAGS_RELEASE     "-O2 -fPIC" CACHE STRING "" FORCE)
      set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
  endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS_RELEASE}")
  set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS_RELEASE}")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS_RELEASE}")
endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")

# Set compiler-specific options, depending on build type
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fcray-pointer -ffree-line-length-none -fno-range-check")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fbacktrace -fdefault-real-8 -fdefault-double-8")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -cpp")
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ffpe-trap=invalid,zero,overflow -fcheck=bounds -fbacktrace -fno-range-check")
  endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mcray=pointer -Kieee -byteswapio -i4 -r8 -Mfprelaxed=div -Mfprelaxed=sqrt")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -tp=haswell -traceback")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mpreprocess")
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Ktrap=fp -Mbounds")
  endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -i4 -real-size 64 -no-prec-div -no-prec-sqrt -xCORE-AVX2 -qno-opt-dynamic-align")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fpp")
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fstack-protector-all -fpe0 -traceback -debug extended -ftrapuv")
  elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3")
  endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message (FATAL_ERROR "CCPP has only been compiled with gfortran, pgf90 and ifort.")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

#------------------------------------------------------------------------------
# The Fortran compiler/linker flag inserted by cmake to create shared libraries
# with the Intel compiler is deprecated (-i_dynamic), correct here.
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  string(REPLACE "-i_dynamic" "-shared-intel"
         CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
         "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")
  string(REPLACE "-i_dynamic" "-shared-intel"
         CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS
         "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")
endif()

#------------------------------------------------------------------------------
# CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/framework/cmake")

#------------------------------------------------------------------------------
# Add -DCCPP preprocessor flag (needed to preprocess GFS_typedefs.F90 from FV3)
ADD_DEFINITIONS(-DCCPP)

#------------------------------------------------------------------------------
# Set MPI flags for C/C++/Fortran preprocessor
if (MPI)
  ADD_DEFINITIONS(-DMPI)
endif (MPI)

#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran compiler and preprocessor
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Enable OpenMP support for C/C++/Fortran compiler")
  ADD_DEFINITIONS(-DOPENMP)
else (OPENMP)
  message (STATUS "Disable OpenMP support for C/C++/Fortran compiler")
endif (OPENMP)

#------------------------------------------------------------------------------
# Set SIONlib flags for Fortran compiler, preprocessor and linker
if (SIONLIB)
  execute_process(COMMAND ${SIONLIB}/bin/sionconfig --mpi --cflags --f90
                  RESULT_VARIABLE exit_code
                  OUTPUT_VARIABLE SIONlib_Fortran_FLAGS
                  ERROR_VARIABLE stderr
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if (NOT exit_code EQUAL "0")
     message (FATAL_ERROR "Execution of '${SIONLIB}/bin/sionconfig --mpi --cflags --f90' failed")
  endif (NOT exit_code EQUAL "0")
  execute_process(COMMAND ${SIONLIB}/bin/sionconfig --mpi --libs --f90
                  RESULT_VARIABLE exit_code
                  OUTPUT_VARIABLE SIONlib_Fortran_LIBS
                  ERROR_VARIABLE stderr
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if (NOT exit_code EQUAL "0")
     message (FATAL_ERROR "Execution of '${SIONLIB}/bin/sionconfig --mpi --libs --f90' failed")
  endif (NOT exit_code EQUAL "0")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${SIONlib_Fortran_FLAGS} ${SIONlib_Fortran_LIBS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SIONlib_Fortran_LIBS}")
  set (CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS} ${SIONlib_Fortran_LIBS}")
  ADD_DEFINITIONS(-DSION)
  message (STATUS "Enable SIONlib support for Fortran compiler")
  message (STATUS "CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS}")
  message (STATUS "CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS}")
  message (STATUS "CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS ${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")
else (SIONLIB)
  message (STATUS "Disable SIONlib support for Fortran compiler")
endif (SIONLIB)

#------------------------------------------------------------------------------
# Build CCPP framework and physics

SET(CCPP_FRAMEWORK_SRC ${CMAKE_SOURCE_DIR}/framework)
ADD_SUBDIRECTORY(${CCPP_FRAMEWORK_SRC} ${CMAKE_BINARY_DIR}/framework)

# Add directory containing compiled CCPP Fortran modules to included directories for physics
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${CMAKE_BINARY_DIR}/framework/src")

SET(CCPP_PHYSICS_SRC ${CMAKE_SOURCE_DIR}/physics)
ADD_SUBDIRECTORY(${CCPP_PHYSICS_SRC} ${CMAKE_BINARY_DIR}/physics)
