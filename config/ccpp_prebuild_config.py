#!/usr/bin/env python

# CCPP prebuild config for GFDL Finite-Volume Cubed-Sphere Model (FV3)


###############################################################################
# Definitions                                                                 #
###############################################################################

HOST_MODEL_IDENTIFIER = "FV3"

# Add all files with metadata tables on the host model side and in CCPP,
# relative to basedir = top-level directory of host model. This includes
# kind and type definitions used in CCPP physics. Also add any internal
# dependencies of these files to the list.
VARIABLE_DEFINITION_FILES = [
    # actual variable definition files
    'framework/src/ccpp_types.F90',
    'physics/physics/machine.F',
    'physics/physics/radsw_param.f',
    'physics/physics/radlw_param.f',
    'physics/physics/h2o_def.f',
    'physics/physics/radiation_surface.f',
    'physics/physics/module_ozphys.F90',
    'data/CCPP_typedefs.F90',
    'data/GFS_typedefs.F90',
    'data/CCPP_data.F90',
    ]

TYPEDEFS_NEW_METADATA = {
    'ccpp_types' : {
        'ccpp_t' : 'cdata',
        'ccpp_types' : '',
        },
    'machine' : {
        'machine' : '',
        },
    'module_radlw_parameters' : {
        'module_radsw_parameters' : '',
        },
    'module_radlw_parameters' : {
        'module_radlw_parameters' : '',
        },
    'module_ozphys' : {
        'module_ozphys' : '',
        'ty_ozphys'     : '',
        },
    'CCPP_typedefs' : {
        'GFS_interstitial_type' : 'GFS_Interstitial(cdata%thrd_no)',
        'GFDL_interstitial_type' : 'GFDL_interstitial',
        'CCPP_typedefs' : '',
        },
    'CCPP_data' : {
        'CCPP_data' : '',
        },
    'GFS_typedefs' : {
        'GFS_control_type'      : 'GFS_Control',
        'GFS_data_type'         : 'GFS_Data(cdata%blk_no)',
        'GFS_diag_type'         : 'GFS_Data(cdata%blk_no)%Intdiag',
        'GFS_tbd_type'          : 'GFS_Data(cdata%blk_no)%Tbd',
        'GFS_sfcprop_type'      : 'GFS_Data(cdata%blk_no)%Sfcprop',
        'GFS_coupling_type'     : 'GFS_Data(cdata%blk_no)%Coupling',
        'GFS_statein_type'      : 'GFS_Data(cdata%blk_no)%Statein',
        'GFS_cldprop_type'      : 'GFS_Data(cdata%blk_no)%Cldprop',
        'GFS_radtend_type'      : 'GFS_Data(cdata%blk_no)%Radtend',
        'GFS_grid_type'         : 'GFS_Data(cdata%blk_no)%Grid',
        'GFS_stateout_type'     : 'GFS_Data(cdata%blk_no)%Stateout',
        'GFS_typedefs' : '',
        },
    }

# Add all physics scheme files relative to basedir
SCHEME_FILES = [
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    'physics/physics/GFS_DCNV_generic_pre.F90',
    'physics/physics/GFS_DCNV_generic_post.F90',
    'physics/physics/GFS_GWD_generic_pre.F90',
    'physics/physics/GFS_GWD_generic_post.F90',
    'physics/physics/GFS_MP_generic_pre.F90',
    'physics/physics/GFS_MP_generic_post.F90',
    'physics/physics/GFS_PBL_generic_pre.F90',
    'physics/physics/GFS_PBL_generic_post.F90',
    'physics/physics/GFS_SCNV_generic_pre.F90',
    'physics/physics/GFS_SCNV_generic_post.F90',
    'physics/physics/GFS_debug.F90',
    'physics/physics/GFS_phys_time_vary.fv3.F90',
    'physics/physics/GFS_rad_time_vary.fv3.F90',
    'physics/physics/GFS_radiation_surface.F90',
    'physics/physics/GFS_rrtmg_post.F90',
    'physics/physics/GFS_rrtmg_pre.F90',
    'physics/physics/GFS_rrtmg_setup.F90',
    'physics/physics/GFS_stochastics.F90',
    'physics/physics/GFS_suite_interstitial_rad_reset.F90',
    'physics/physics/GFS_suite_interstitial_phys_reset.F90',
    'physics/physics/GFS_suite_interstitial_1.F90',
    'physics/physics/GFS_suite_interstitial_2.F90',
    'physics/physics/GFS_suite_stateout_reset.F90',
    'physics/physics/GFS_suite_stateout_update.F90',
    'physics/physics/GFS_suite_interstitial_3.F90',
    'physics/physics/GFS_suite_interstitial_4.F90',
    'physics/physics/GFS_suite_interstitial_5.F90',
    'physics/physics/GFS_surface_generic_pre.F90',
    'physics/physics/GFS_surface_generic_post.F90',
    'physics/physics/GFS_surface_composites_pre.F90',
    'physics/physics/GFS_surface_composites_inter.F90',
    'physics/physics/GFS_surface_composites_post.F90',
    'physics/physics/GFS_surface_loop_control_part1.F90',
    'physics/physics/GFS_surface_loop_control_part2.F90',
    'physics/physics/GFS_time_vary_pre.fv3.F90',
    'physics/physics/GFS_physics_post.F90',
    'physics/physics/cires_ugwp.F90',
    'physics/physics/cires_ugwp_post.F90',
    'physics/physics/unified_ugwp.F90',
    'physics/physics/unified_ugwp_post.F90',
    'physics/physics/ugwpv1_gsldrag.F90',
    'physics/physics/ugwpv1_gsldrag_post.F90',
    'physics/physics/cnvc90.f',
    'physics/physics/cs_conv_pre.F90',
    'physics/physics/cs_conv.F90',
    'physics/physics/cs_conv_post.F90',
    'physics/physics/cs_conv_aw_adj.F90',
    'physics/physics/cu_ntiedtke_pre.F90',
    'physics/physics/cu_ntiedtke.F90',
    'physics/physics/cu_ntiedtke_post.F90',
    'physics/physics/dcyc2t3.f',
    'physics/physics/drag_suite.F90',
    'physics/physics/shoc.F90',
    'physics/physics/get_prs_fv3.F90',
    'physics/physics/get_phi_fv3.F90',
    'physics/physics/gfdl_cloud_microphys.F90',
    'physics/physics/fv_sat_adj.F90',
    'physics/physics/gfdl_sfc_layer.F90',
    'physics/physics/zhaocarr_gscond.f',
    'physics/physics/gwdc_pre.f',
    'physics/physics/gwdc.f',
    'physics/physics/gwdc_post.f',
    'physics/physics/gwdps.f',
    'physics/physics/h2ophys.f',
    'physics/physics/samfdeepcnv.f',
    'physics/physics/samfshalcnv.f',
    'physics/physics/sascnvn.F',
    'physics/physics/shalcnv.F',
    'physics/physics/maximum_hourly_diagnostics.F90',
    'physics/physics/m_micro.F90',
    'physics/physics/m_micro_pre.F90',
    'physics/physics/m_micro_post.F90',
    'physics/physics/cu_gf_driver_pre.F90',
    'physics/physics/cu_gf_driver.F90',
    'physics/physics/cu_gf_driver_post.F90',
    'physics/physics/cu_c3_driver_pre.F90',
    'physics/physics/cu_c3_driver.F90',
    'physics/physics/cu_c3_driver_post.F90',
    'physics/physics/hedmf.f',
    'physics/physics/moninshoc.f',
    'physics/physics/satmedmfvdif.F',
    'physics/physics/satmedmfvdifq.F',
    'physics/physics/shinhongvdif.F90',
    'physics/physics/ysuvdif.F90',
    'physics/physics/mynnedmf_wrapper.F90',
    'physics/physics/mynnsfc_wrapper.F90',
    'physics/physics/sgscloud_radpre.F90',
    'physics/physics/sgscloud_radpost.F90',
    'physics/physics/myjsfc_wrapper.F90',
    'physics/physics/myjpbl_wrapper.F90',
    'physics/physics/mp_thompson_pre.F90',
    'physics/physics/mp_thompson.F90',
    'physics/physics/mp_thompson_post.F90',
    'physics/physics/mp_nssl.F90',
    'physics/physics/zhaocarr_precpd.f',
    'physics/physics/radlw_main.F90',
    'physics/physics/radsw_main.F90',
    'physics/physics/rascnv.F90',
    'physics/physics/rayleigh_damp.f',
    'physics/physics/rrtmg_lw_post.F90',
    'physics/physics/rrtmg_lw_pre.F90',
    'physics/physics/rrtmg_sw_post.F90',
    'physics/physics/rad_sw_pre.F90',
    'physics/physics/sfc_diag.f',
    'physics/physics/sfc_diag_post.F90',
    'physics/physics/lsm_ruc.F90',
    'physics/physics/sfc_cice.f',
    'physics/physics/sfc_diff.f',
    'physics/physics/lsm_noah.f',
    'physics/physics/noahmpdrv.F90',
    'physics/physics/noahmpdrv_time_vary.F90',
    'physics/physics/flake_driver.F90',
    'physics/physics/clm_lake.f90',
    'physics/physics/sfc_nst_pre.f90',
    'physics/physics/sfc_nst.f90',
    'physics/physics/sfc_nst_post.f90',
    'physics/physics/sfc_ocean.F',
    'physics/physics/sfc_sice.f',
    # HAFS FER_HIRES
    'physics/physics/mp_fer_hires.F90',
    # SMOKE
    'physics/physics/smoke_dust/rrfs_smoke_wrapper.F90',
    'physics/physics/smoke_dust/rrfs_smoke_postpbl.F90',
    # RRTMGP
    'physics/physics/rrtmgp_aerosol_optics.F90',
    'physics/physics/rrtmgp_lw_main.F90',
    'physics/physics/rrtmgp_sw_main.F90',
    'physics/physics/GFS_rrtmgp_setup.F90',
    'physics/physics/GFS_rrtmgp_pre.F90',
    'physics/physics/GFS_cloud_diagnostics.F90',
    'physics/physics/GFS_rrtmgp_cloud_mp.F90',
    'physics/physics/GFS_rrtmgp_cloud_overlap.F90',
    'physics/physics/GFS_rrtmgp_post.F90'
    ]

# Default build dir, relative to current working directory,
# if not specified as command-line argument
DEFAULT_BUILD_DIR = 'build'

# Auto-generated makefile/cmakefile snippets that contain all type definitions
TYPEDEFS_MAKEFILE   = '{build_dir}/physics/CCPP_TYPEDEFS.mk'
TYPEDEFS_CMAKEFILE  = '{build_dir}/physics/CCPP_TYPEDEFS.cmake'
TYPEDEFS_SOURCEFILE = '{build_dir}/physics/CCPP_TYPEDEFS.sh'

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE   = '{build_dir}/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE  = '{build_dir}/physics/CCPP_SCHEMES.cmake'
SCHEMES_SOURCEFILE = '{build_dir}/physics/CCPP_SCHEMES.sh'

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE   = '{build_dir}/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE  = '{build_dir}/physics/CCPP_CAPS.cmake'
CAPS_SOURCEFILE = '{build_dir}/physics/CCPP_CAPS.sh'

# Directory where to put all auto-generated physics caps
CAPS_DIR = '{build_dir}/physics'

# Directory where the suite definition files are stored
SUITES_DIR = 'suites'

# Directory where to write static API to
STATIC_API_DIR = '{build_dir}/physics'
STATIC_API_CMAKEFILE = '{build_dir}/physics/CCPP_STATIC_API.cmake'
STATIC_API_SOURCEFILE = '{build_dir}/physics/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
# used by metadata2html.py for generating scientific documentation
METADATA_HTML_OUTPUT_DIR = '{build_dir}/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = '{build_dir}/physics/CCPP_VARIABLES_FV3.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = '{build_dir}/framework/doc/DevelopersGuide/CCPP_VARIABLES_FV3.tex'
