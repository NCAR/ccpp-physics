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
    'physics/physics/hooks/machine.F',
    'physics/physics/Radiation/RRTMG/radsw_param.f',
    'physics/physics/Radiation/RRTMG/radlw_param.f',
    'physics/physics/photochem/module_ozphys.F90',
    'physics/physics/photochem/module_h2ophys.F90',
    'physics/physics/SFC_Models/Land/Noahmp/lnd_iau_mod.F90',
    'data/CCPP_typedefs.F90',
    'data/GFS_typedefs.F90',
    'data/CCPP_data.F90',
    ]

TYPEDEFS_NEW_METADATA = {
    'ccpp_types' : {
        'ccpp_t' : 'cdata',
        'MPI_Comm' : '',
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
    'module_h2ophys' : {
        'module_h2ophys' : '',
        'ty_h2ophys'     : '',
        },
    'land_iau_mod' : {
        'land_iau_mod' : '',
        'land_iau_external_data_type' : '',
        'land_iau_state_type' : '',
        'land_iau_control_type' : '',
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
        'GFS_statein_type'      : 'GFS_Statein',
        'GFS_stateout_type'     : 'GFS_Stateout',
        'GFS_grid_type'         : 'GFS_Grid',
        'GFS_tbd_type'          : 'GFS_Tbd',
        'GFS_cldprop_type'      : 'GFS_Cldprop',
        'GFS_sfcprop_type'      : 'GFS_Sfcprop',
        'GFS_radtend_type'      : 'GFS_Radtend',
        'GFS_coupling_type'     : 'GFS_Coupling',
        'GFS_diag_type'         : 'GFS_Intdiag',
        'GFS_typedefs' : '',
        },
    }

# Add all physics scheme files relative to basedir
SCHEME_FILES = [
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_DCNV_generic_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_DCNV_generic_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_GWD_generic_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_GWD_generic_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_MP_generic_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_MP_generic_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_PBL_generic_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_PBL_generic_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_SCNV_generic_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_SCNV_generic_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_debug.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_phys_time_vary.fv3.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_photochemistry.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rad_time_vary.fv3.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_radiation_surface.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_setup.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_setup.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_cloud_diagnostics.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_cloud_mp.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_cloud_overlap.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_stochastics.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_rad_reset.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_phys_reset.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_1.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_2.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_stateout_reset.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_stateout_update.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_3.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_4.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_5.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_generic_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_generic_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_inter.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part1.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part2.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_time_vary_pre.fv3.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/cnvc90.f',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/dcyc2t3.f',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/maximum_hourly_diagnostics.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_physics_post.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/sgscloud_radpre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/sgscloud_radpost.F90',
    'physics/physics/CONV/Chikira_Sugiyama/cs_conv_pre.F90',
    'physics/physics/CONV/Chikira_Sugiyama/cs_conv.F90',
    'physics/physics/CONV/Chikira_Sugiyama/cs_conv_post.F90',
    'physics/physics/CONV/Chikira_Sugiyama/cs_conv_aw_adj.F90',
    'physics/physics/CONV/nTiedtke/cu_ntiedtke_pre.F90',
    'physics/physics/CONV/nTiedtke/cu_ntiedtke.F90',
    'physics/physics/CONV/nTiedtke/cu_ntiedtke_post.F90',
    'physics/physics/CONV/SAMF/samfdeepcnv.f',
    'physics/physics/CONV/SAMF/samfshalcnv.f',
    'physics/physics/CONV/SAS/sascnvn.F',
    'physics/physics/CONV/SAS/shalcnv.F',
    'physics/physics/CONV/Grell_Freitas/cu_gf_driver_pre.F90',
    'physics/physics/CONV/Grell_Freitas/cu_gf_driver.F90',
    'physics/physics/CONV/Grell_Freitas/cu_gf_driver_post.F90',
    'physics/physics/CONV/C3/cu_c3_driver_pre.F90',
    'physics/physics/CONV/C3/cu_c3_driver.F90',
    'physics/physics/CONV/C3/cu_c3_driver_post.F90',
    'physics/physics/CONV/RAS/rascnv.F90',
    'physics/physics/GWD/cires_ugwp.F90',
    'physics/physics/GWD/cires_ugwp_post.F90',
    'physics/physics/GWD/unified_ugwp.F90',
    'physics/physics/GWD/unified_ugwp_post.F90',
    'physics/physics/GWD/ugwpv1_gsldrag.F90',
    'physics/physics/GWD/ugwpv1_gsldrag_post.F90',
    'physics/physics/GWD/drag_suite.F90',
    'physics/physics/GWD/gwdc_pre.f',
    'physics/physics/GWD/gwdc.f',
    'physics/physics/GWD/gwdc_post.f',
    'physics/physics/GWD/gwdps.f',
    'physics/physics/GWD/rayleigh_damp.f',
    'physics/physics/photochem/module_h2ophys.F90',
    'physics/physics/photochem/module_ozphys.F90',
    'physics/physics/MP/Ferrier_Aligo/mp_fer_hires.F90',
    'physics/physics/MP/GFDL/gfdl_cloud_microphys.F90',
    'physics/physics/MP/GFDL/fv_sat_adj.F90',
    'physics/physics/MP/Morrison_Gettelman/m_micro.F90',
    'physics/physics/MP/Morrison_Gettelman/m_micro_pre.F90',
    'physics/physics/MP/Morrison_Gettelman/m_micro_post.F90',
    'physics/physics/MP/NSSL/mp_nssl.F90',
    'physics/physics/MP/Thompson/mp_thompson_pre.F90',
    'physics/physics/MP/Thompson/mp_thompson.F90',
    'physics/physics/MP/Thompson/mp_thompson_post.F90',
    'physics/physics/MP/Zhao_Carr/zhaocarr_gscond.f',
    'physics/physics/MP/Zhao_Carr/zhaocarr_precpd.f',
    'physics/physics/PBL/HEDMF/hedmf.f',
    'physics/physics/PBL/SHOC/moninshoc.f',
    'physics/physics/PBL/SHOC/shoc.F90',
    'physics/physics/PBL/MYJ/myjpbl_wrapper.F90',
    'physics/physics/PBL/MYNN_EDMF/mynnedmf_wrapper.F90',
    'physics/physics/PBL/SATMEDMF/satmedmfvdif.F',
    'physics/physics/PBL/SATMEDMF/satmedmfvdifq.F',
    'physics/physics/PBL/YSU/ysuvdif.F90',
    'physics/physics/PBL/saYSU/shinhongvdif.F90',
    'physics/physics/Radiation/RRTMG/radsw_main.F90',
    'physics/physics/Radiation/RRTMG/radlw_main.F90',
    'physics/physics/Radiation/RRTMG/rrtmg_lw_post.F90',
    'physics/physics/Radiation/RRTMG/rrtmg_sw_post.F90',
    'physics/physics/Radiation/RRTMG/rad_sw_pre.F90',
    'physics/physics/Radiation/RRTMGP/rrtmgp_aerosol_optics.F90',
    'physics/physics/Radiation/RRTMGP/rrtmgp_lw_main.F90',
    'physics/physics/Radiation/RRTMGP/rrtmgp_sw_main.F90',
    'physics/physics/SFC_Layer/GFDL/gfdl_sfc_layer.F90',
    'physics/physics/SFC_Layer/MYNN/mynnsfc_wrapper.F90',
    'physics/physics/SFC_Layer/MYJ/myjsfc_wrapper.F90',
    'physics/physics/SFC_Layer/UFS/sfc_diag.f',
    'physics/physics/SFC_Layer/UFS/sfc_diag_post.F90',
    'physics/physics/SFC_Layer/UFS/sfc_diff.f',
    'physics/physics/SFC_Layer/UFS/sfc_nst_pre.f90',
    'physics/physics/SFC_Layer/UFS/sfc_nst.f90',
    'physics/physics/SFC_Layer/UFS/sfc_nst_post.f90',
    'physics/physics/SFC_Models/Land/RUC/lsm_ruc.F90',
    'physics/physics/SFC_Models/SeaIce/CICE/sfc_cice.f',
    'physics/physics/SFC_Models/Land/sfc_land.F90',
    'physics/physics/SFC_Models/Land/Noah/lsm_noah.f',
    'physics/physics/SFC_Models/Land/Noahmp/noahmpdrv.F90',
    'physics/physics/SFC_Models/Lake/Flake/flake_driver.F90',
    'physics/physics/SFC_Models/Lake/CLM/clm_lake.f90',
    'physics/physics/SFC_Models/Ocean/UFS/sfc_ocean.F',
    'physics/physics/SFC_Models/SeaIce/CICE/sfc_sice.f',
    'physics/physics/smoke_dust/rrfs_smoke_wrapper.F90',
    'physics/physics/smoke_dust/rrfs_smoke_postpbl.F90',
    'physics/physics/tools/get_prs_fv3.F90',
    'physics/physics/tools/get_phi_fv3.F90'
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
