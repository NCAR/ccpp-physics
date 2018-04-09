#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.11)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
    project(ccppphys VERSION 0.0.1)
else(POLICY CMP0048)
    project(ccppphys)
    set(PROJECT_VERSION 0.0.1)
    set(PROJECT_VERSION_MAJOR 0)
    set(PROJECT_VERSION_MINOR 0)
    set(PROJECT_VERSION_PATCH 1)
endif(POLICY CMP0048)

#------------------------------------------------------------------------------
set(PACKAGE "ccppphys")
set(AUTHORS  "Grant J. Firl" "Dom Heinzeller")

#------------------------------------------------------------------------------
# Enable Fortran
enable_language(Fortran)

#------------------------------------------------------------------------------
# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../ccpp-framework/cmake")

#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Enable OpenMP support for C/C++/Fortran compiler")
else(OPENMP)
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
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                 "MinSizeRel" "RelWithDebInfo")
endif()

#------------------------------------------------------------------------------
# By default we want a shared library
option(BUILD_SHARED_LIBS "Build a shared library" ON)

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
# Set the sources
set(SOURCES
    ./physics/machine.F
    ./physics/physcons.f90
    ./physics/funcphys.f90
    ./GFS_layer/GFS_typedefs.F90
    ./physics/radlw_param.f
    ./physics/physparam.f
    ./physics/radsw_param.f
    ./physics/radiation_aerosols.f
    ./physics/iounitdef.f
    ./physics/gfs_phy_tracer_config.f
    ./physics/radiation_clouds.f
    ./physics/module_bfmicrophysics.f
    ./physics/radiation_gases.f
    ./physics/ozne_def.f
    ./physics/radiation_surface.f
    ./physics//h2o_def.f
    ./physics/sfcsub.F
    ./physics/mersenne_twister.f
    ./physics/rascnvv2.f
    ./physics/GFDL_parse_tracers.F90
    ./physics/wam_f107_kp_mod.f90
    ./physics/num_parthds.F
    ./physics/gocart_tracer_config_stub.f
    ./GFS_layer/GFS_initialize_scm.F90
    ./physics/cldwat2m_micro.F
    ./physics/wv_saturation.F
    ./physics/aer_cloud.F
    ./physics/rad_initialize.f
    ./GFS_layer/GFS_radiation_driver.F90
    ./physics/radcons.f90
    ./physics/radiation_astronomy.f
    ./physics/radsw_main.f
    ./physics/radsw_datatb.f
    ./physics/GFS_rrtmg_pre.F90
    ./physics/GFS_rrtmg_post.F90
    ./physics/rrtmg_sw_pre.F90
    ./physics/rrtmg_sw_post.F90
    ./physics/rrtmg_lw_pre.F90
    ./physics/rrtmg_lw_post.F90
    ./physics/radlw_main.f
    ./physics/radlw_datatb.f
    ./physics/set_soilveg.f
    ./physics/namelist_soilveg.f
    ./physics/GFS_phys_time_vary.f90
    ./physics/ozinterp.f90
    ./physics/h2ointerp.f90
    ./physics/gcycle.f90
    ./physics/GFS_rad_time_vary.f90
    ./physics/GFS_radupdate.f90
    ./physics/GFS_suite_interstitial.ccpp.f90
    ./physics/get_prs_fv3.f90
    ./physics/sfc_sice.f
    ./physics/dcyc2.f
    ./physics/GFS_surface_generic.f90
    ./physics/GFS_PBL_generic.f90
    ./physics/sfc_drv.f
    ./physics/sflx.f
    ./physics/sfc_diff.f
    ./physics/GFS_surface_loop_control.f
    ./physics/sfc_nst.f
    ./physics/date_def.f
    ./physics/module_nst_water_prop.f90
    ./physics/module_nst_parameters.f90
    ./physics/module_nst_model.f90
    ./physics/sfc_diag.f
    ./physics/moninedmf.f
    ./physics/mfpbl.f
    ./physics/tridi.f
    ./physics/gwdps.f
    ./physics/rayleigh_damp.f
    ./physics/ozphys.f
    ./physics/GFS_DCNV_generic.f90
    ./physics/GFS_zhao_carr_pre.f90
    ./physics/mfdeepcnv.f
    ./physics/gwdc.f
    ./physics/GFS_SCNV_generic.f90
    ./physics/mfshalcnv.f
    ./physics/cnvc90.f
    ./physics/GFS_MP_generic_pre.f90
    ./physics/gscond.f
    ./physics/precpd.f
    ./physics/GFS_calpreciptype.f90
    ./physics/GFS_MP_generic_post.f90
    ./physics/cnvc90_cap.F90
    ./physics/lsm_noah_pre_cap.F90
    ./physics/GFS_DCNV_generic_post_cap.F90
    ./physics/ozphys_pre_cap.F90
    ./physics/gwdps_pre_cap.F90
    ./physics/ozphys_cap.F90
    ./physics/GFS_suite_interstitial_rad_reset_cap.F90
    ./physics/GFS_MP_generic_post_cap.F90
    ./physics/rrtmg_sw_cap.F90
    ./physics/dcyc2t3_post_cap.F90
    ./physics/sfc_nst_cap.F90
    ./physics/rrtmg_lw_cap.F90
    ./physics/GFS_zhao_carr_pre_cap.F90
    ./physics/GFS_PBL_generic_pre_cap.F90
    ./physics/gwdps_post_cap.F90
    ./physics/GFS_suite_interstitial_1_cap.F90
    ./physics/sfc_sice_pre_cap.F90
    ./physics/GFS_suite_interstitial_2_cap.F90
    ./physics/dcyc2t3_pre_cap.F90
    ./physics/rrtmg_lw_pre_cap.F90
    ./physics/GFS_rad_time_vary_cap.F90
    ./physics/sasas_shal_post_cap.F90
    ./physics/sfc_sice_cap.F90
    ./physics/GFS_SCNV_generic_pre_cap.F90
    ./physics/GFS_PBL_generic_post_cap.F90
    ./physics/sasas_shal_cap.F90
    ./physics/GFS_surface_loop_control_part2_cap.F90
    ./physics/GFS_rrtmg_pre_cap.F90
    ./physics/GFS_phys_time_vary_2_cap.F90
    ./physics/GFS_phys_time_vary_1_cap.F90
    ./physics/dcyc2t3_cap.F90
    ./physics/gwdc_cap.F90
    ./physics/get_prs_fv3_cap.F90
    ./physics/GFS_surface_generic_post_cap.F90
    ./physics/rrtmg_sw_post_cap.F90
    ./physics/edmf_cap.F90
    ./physics/sfc_ex_coef_cap.F90
    ./physics/GFS_surface_loop_control_part1_cap.F90
    ./physics/GFS_surface_loop_control_part0_cap.F90
    ./physics/rrtmg_lw_post_cap.F90
    ./physics/rayleigh_damp_cap.F90
    ./physics/lsm_noah_post_cap.F90
    ./physics/lsm_noah_cap.F90
    ./physics/GFS_suite_interstitial_phys_reset_cap.F90
    ./physics/gwdps_cap.F90
    ./physics/sfc_nst_pre_cap.F90
    ./physics/GFS_SCNV_generic_post_cap.F90
    ./physics/sfc_diag_cap.F90
    ./physics/zhaocarr_precpd_cap.F90
    ./physics/sasas_deep_cap.F90
    ./physics/GFS_initialize_scm_cap.F90
    ./physics/GFS_DCNV_generic_pre_cap.F90
    ./physics/gwdc_pre_cap.F90
    ./physics/zhaocarr_gscond_cap.F90
    ./physics/GFS_suite_interstitial_3_cap.F90
    ./physics/GFS_calpreciptype_cap.F90
    ./physics/GFS_MP_generic_pre_cap.F90
    ./physics/sfc_sice_post_cap.F90
    ./physics/gwdc_post_cap.F90
    ./physics/GFS_suite_update_stateout_cap.F90
    ./physics/sfc_nst_post_cap.F90
    ./physics/GFS_surface_generic_pre_cap.F90
    ./physics/ozphys_post_cap.F90
    ./physics/GFS_rrtmg_post_cap.F90
    ./physics/get_phi_fv3_cap.F90
    ./physics/rrtmg_sw_pre_cap.F90
)

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(f_flags -ffree-line-length-none)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

  SET_SOURCE_FILES_PROPERTIES(./physics/module_bfmicrophysics.f ./physics/rascnvv2.f ./physics/sflx.f ./physics/sfc_diff.f ./physics/sfc_diag.f PROPERTIES COMPILE_FLAGS -fdefault-real-8)
  SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_model.f90 ./physics/GFS_calpreciptype.f90 PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -ffree-form")
  SET_SOURCE_FILES_PROPERTIES(./physics/mersenne_twister.f PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fno-range-check")
  SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_water_prop.f90 PROPERTIES COMPILE_FLAGS "-ffree-line-length-none -fdefault-real-8 -ffree-form")
  SET_SOURCE_FILES_PROPERTIES(./physics/aer_cloud.F ./physics/wv_saturation.F ./physics/cldwat2m_micro.F PROPERTIES COMPILE_FLAGS "-DNEMS_GSM -fdefault-real-8 -fdefault-double-8")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")

  SET_SOURCE_FILES_PROPERTIES(./physics/module_bfmicrophysics.f ./physics/rascnvv2.f ./physics/sflx.f ./physics/sfc_diff.f ./physics/sfc_diag.f PROPERTIES COMPILE_FLAGS -r8)
  SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_model.f90 ./physics/GFS_calpreciptype.f90 PROPERTIES COMPILE_FLAGS "-r8 -free")
  SET_SOURCE_FILES_PROPERTIES(./physics/mersenne_twister.f PROPERTIES COMPILE_FLAGS "-r8 -ftz")
  SET_SOURCE_FILES_PROPERTIES(./physics/module_nst_water_prop.f90 PROPERTIES COMPILE_FLAGS "-extend-source 132 -r8 -free")
  SET_SOURCE_FILES_PROPERTIES(./physics/aer_cloud.F ./physics/wv_saturation.F ./physics/cldwat2m_micro.F PROPERTIES COMPILE_FLAGS "-DNEMS_GSM -r8")
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message (FATAL_ERROR "This program has only been compiled with gfortran and ifort. If another compiler is needed, the appropriate flags must be added in ${GFS_PHYS_SRC}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

#apply general fortran tags to all fortran source files
if(${CMAKE_VERSION} LESS 3.3)
  string (REPLACE ";" " " f_flags_str "${f_flags}")
  SET_PROPERTY(SOURCE ${GFS_phys_source_code} APPEND_STRING PROPERTY COMPILE_FLAGS " ${f_flags_str}")
else(${CMAKE_VERSION} LESS 3.3)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${f_flags}>")
endif (${CMAKE_VERSION} LESS 3.3)

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

#------------------------------------------------------------------------------
# Add the auto-generated caps
#add_custom_command(
#    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/scm_cap.f90
#    DEPENDS ${CCPP_MKCAP}
#    COMMAND ${CCPP_MKCAP} -o ${CMAKE_CURRENT_BINARY_DIR}/scm_cap.f90
#                             ${CMAKE_CURRENT_SOURCE_DIR}/scheme.xml
#)
#list(APPEND SOURCES ${CMAKE_CURRENT_BINARY_DIR}/scm_test1_cap.f90)


add_library(ccppphys ${SOURCES})
target_link_libraries(ccppphys LINK_PUBLIC ${LIBS} w3 sp bacio)
set_target_properties(ccppphys PROPERTIES VERSION ${PROJECT_VERSION}
                                     SOVERSION ${PROJECT_VERSION_MAJOR}
                                     COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                                     LINK_FLAGS "${CMAKE_Fortran_FLAGS}")
