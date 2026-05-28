#!/usr/bin/env python

# CCPP prebuild config for CCPP Single Column Model (SCM)


###############################################################################
# Definitions                                                                 #
###############################################################################

HOST_MODEL_IDENTIFIER = "SCM"

# Add all files with metadata tables on the host model side and in CCPP,
# relative to basedir = top-level directory of host model. This includes
# kind and type definitions used in CCPP physics. Also add any internal
# dependencies of these files to the list.
VARIABLE_DEFINITION_FILES = [
    # actual variable definition files
    'ccpp/framework/src/ccpp_types.F90',
    'ccpp/physics/physics/hooks/machine.F',
    'ccpp/physics/physics/Radiation/RRTMG/radsw_param.f',
    'ccpp/physics/physics/Radiation/RRTMG/radlw_param.f',
    'ccpp/physics/physics/photochem/h2o_def.f',
    'ccpp/physics/physics/photochem/module_ozphys.F90',
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/module_ccpp_suite_simulator.F90',
    'scm/src/CCPP_typedefs.F90',
    'scm/src/GFS_typedefs.F90',
    'scm/src/scm_kinds.F90',
    'scm/src/scm_type_defs.F90',
    'scm/src/scm_physical_constants.F90',
    'scm/src/scm_utils.F90', #no definitions, but scm_type_defs.F90 uses a module from this file
    ]

TYPEDEFS_NEW_METADATA = {
    'ccpp_types' : {
        'ccpp_types' : '',
        'MPI_Comm' : '',
        'ccpp_t' : 'cdata',
        },
    'machine' : {
        'machine' : '',
        },
    'module_radsw_parameters' : {
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
        'GFS_interstitial_type' : 'physics%Interstitial',
        'CCPP_typedefs' : '',
        },
    'GFS_typedefs' : {
        'GFS_diag_type' : 'physics%Diag',
        'GFS_control_type' : 'physics%Model',
        'GFS_cldprop_type' : 'physics%Cldprop',
        'GFS_tbd_type' : 'physics%Tbd',
        'GFS_sfcprop_type' : 'physics%Sfcprop',
        'GFS_coupling_type' : 'physics%Coupling',
        'GFS_statein_type' : 'physics%Statein',
        'GFS_radtend_type' : 'physics%Radtend',
        'GFS_grid_type' : 'physics%Grid',
        'GFS_stateout_type' : 'physics%Stateout',
        'GFS_typedefs' : '',
        },
    'scm_physical_constants' : {
        'scm_physical_constants' : '',
        },
    'scm_type_defs' : {
        'scm_type_defs' : '',
        'physics_type' : 'physics',
        },
    'module_ccpp_suite_simulator' : {
        'base_physics_process' : '',
        'module_ccpp_suite_simulator' : '',
        },
    }

# Add all physics scheme files relative to basedir
SCHEME_FILES = [
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_DCNV_generic_pre.F90'         ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_DCNV_generic_post.F90'        ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_GWD_generic_pre.F90'          ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_GWD_generic_post.F90'         ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_MP_generic_pre.F90'           ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_MP_generic_post.F90'          ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_PBL_generic_pre.F90'          ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_PBL_generic_post.F90'         ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_SCNV_generic_pre.F90'         ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_SCNV_generic_post.F90'        ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_phys_time_vary.scm.F90'       ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rad_time_vary.scm.F90'        ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_radiation_surface.F90'        ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_post.F90'               ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_pre.F90'                ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_setup.F90'              ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_setup.F90'             ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_pre.F90'               ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_cloud_diagnostics.F90'        ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_cloud_mp.F90'          ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_cloud_overlap.F90'     ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_post.F90'              ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_rad_reset.F90',
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_phys_reset.F90',
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_1.F90'     ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_2.F90'     ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_stateout_reset.F90'     ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_stateout_update.F90'    ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_3.F90'     ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_4.F90'     ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_5.F90'     ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_generic_pre.F90'      ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_generic_post.F90'     ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_pre.F90'   ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_inter.F90' ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_post.F90'  ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part1.F90' ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part2.F90' ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_time_vary_pre.scm.F90'        ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/cnvc90.f'                         ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/dcyc2t3.f'                        ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/maximum_hourly_diagnostics.F90'   ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_physics_post.F90'             ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/sgscloud_radpre.F90'              ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/sgscloud_radpost.F90'             ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/scm_sfc_flux_spec.F90'            ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_ccpp_suite_sim_pre.F90'       ,
    'ccpp/physics/physics/Interstitials/UFS_SCM_NEPTUNE/ccpp_suite_simulator.F90'         ,
#    'ccpp/physics/physics/bl_acm.F90'                       ,
    'ccpp/physics/physics/CONV/Chikira_Sugiyama/cs_conv_pre.F90',
    'ccpp/physics/physics/CONV/Chikira_Sugiyama/cs_conv.F90',
    'ccpp/physics/physics/CONV/Chikira_Sugiyama/cs_conv_post.F90',
    'ccpp/physics/physics/CONV/Chikira_Sugiyama/cs_conv_aw_adj.F90',
    'ccpp/physics/physics/CONV/nTiedtke/cu_ntiedtke_pre.F90',
    'ccpp/physics/physics/CONV/nTiedtke/cu_ntiedtke.F90',
    'ccpp/physics/physics/CONV/nTiedtke/cu_ntiedtke_post.F90',
    'ccpp/physics/physics/CONV/SAMF/samfdeepcnv.f',
    'ccpp/physics/physics/CONV/SAMF/samfshalcnv.f',
    'ccpp/physics/physics/CONV/SAS/sascnvn.F',
    'ccpp/physics/physics/CONV/SAS/shalcnv.F',
    'ccpp/physics/physics/CONV/Grell_Freitas/cu_gf_driver_pre.F90',
    'ccpp/physics/physics/CONV/Grell_Freitas/cu_gf_driver.F90',
    'ccpp/physics/physics/CONV/Grell_Freitas/cu_gf_driver_post.F90',
    'ccpp/physics/physics/CONV/C3/cu_c3_driver_pre.F90',
    'ccpp/physics/physics/CONV/C3/cu_c3_driver.F90',
    'ccpp/physics/physics/CONV/C3/cu_c3_driver_post.F90',
    'ccpp/physics/physics/CONV/RAS/rascnv.F90',
    'ccpp/physics/physics/GWD/cires_ugwp.F90',
    'ccpp/physics/physics/GWD/cires_ugwp_post.F90',
    'ccpp/physics/physics/GWD/unified_ugwp.F90',
    'ccpp/physics/physics/GWD/unified_ugwp_post.F90',
    'ccpp/physics/physics/GWD/ugwpv1_gsldrag.F90',
    'ccpp/physics/physics/GWD/ugwpv1_gsldrag_post.F90',
    'ccpp/physics/physics/GWD/drag_suite.F90',
    'ccpp/physics/physics/GWD/gwdc_pre.f',
    'ccpp/physics/physics/GWD/gwdc.f',
    'ccpp/physics/physics/GWD/gwdc_post.f',
    'ccpp/physics/physics/GWD/gwdps.f',
    'ccpp/physics/physics/GWD/rayleigh_damp.f',
    'ccpp/physics/physics/photochem/h2ophys.f',
    'ccpp/physics/physics/photochem/module_ozphys.F90',
    'ccpp/physics/physics/MP/Ferrier_Aligo/mp_fer_hires.F90',
    'ccpp/physics/physics/MP/GFDL/gfdl_cloud_microphys.F90',
    'ccpp/physics/physics/MP/GFDL/fv_sat_adj.F90',
    'ccpp/physics/physics/MP/Morrison_Gettelman/m_micro.F90',
    'ccpp/physics/physics/MP/Morrison_Gettelman/m_micro_pre.F90',
    'ccpp/physics/physics/MP/Morrison_Gettelman/m_micro_post.F90',
    'ccpp/physics/physics/MP/NSSL/mp_nssl.F90',
    'ccpp/physics/physics/MP/Thompson/mp_thompson_pre.F90',
    'ccpp/physics/physics/MP/Thompson/mp_thompson.F90',
    'ccpp/physics/physics/MP/Thompson/mp_thompson_post.F90',
    'ccpp/physics/physics/MP/Zhao_Carr/zhaocarr_gscond.f',
    'ccpp/physics/physics/MP/Zhao_Carr/zhaocarr_precpd.f',
    'ccpp/physics/physics/PBL/HEDMF/hedmf.f',
    'ccpp/physics/physics/PBL/SHOC/moninshoc.f',
    'ccpp/physics/physics/PBL/SHOC/shoc.F90',
    'ccpp/physics/physics/PBL/MYJ/myjpbl_wrapper.F90',
    'ccpp/physics/physics/PBL/MYNN_EDMF/mynnedmf_wrapper.F90',
    'ccpp/physics/physics/PBL/SATMEDMF/satmedmfvdif.F',
    'ccpp/physics/physics/PBL/SATMEDMF/satmedmfvdifq.F',
    'ccpp/physics/physics/PBL/YSU/ysuvdif.F90',
    'ccpp/physics/physics/PBL/saYSU/shinhongvdif.F90',
    'ccpp/physics/physics/Radiation/RRTMG/radsw_main.F90',
    'ccpp/physics/physics/Radiation/RRTMG/radlw_main.F90',
    'ccpp/physics/physics/Radiation/RRTMG/rrtmg_lw_post.F90',
    'ccpp/physics/physics/Radiation/RRTMG/rrtmg_sw_post.F90',
    'ccpp/physics/physics/Radiation/RRTMG/rad_sw_pre.F90',
    'ccpp/physics/physics/Radiation/RRTMGP/rrtmgp_aerosol_optics.F90',
    'ccpp/physics/physics/Radiation/RRTMGP/rrtmgp_lw_main.F90',
    'ccpp/physics/physics/Radiation/RRTMGP/rrtmgp_sw_main.F90',
    'ccpp/physics/physics/SFC_Layer/GFDL/gfdl_sfc_layer.F90',
    'ccpp/physics/physics/SFC_Layer/MYNN/mynnsfc_wrapper.F90',
    'ccpp/physics/physics/SFC_Layer/MYJ/myjsfc_wrapper.F90',
    'ccpp/physics/physics/SFC_Layer/UFS/sfc_diag.f',
    'ccpp/physics/physics/SFC_Layer/UFS/sfc_diag_post.F90',
    'ccpp/physics/physics/SFC_Layer/UFS/sfc_diff.f',
    'ccpp/physics/physics/SFC_Layer/UFS/sfc_nst_pre.f90',
    'ccpp/physics/physics/SFC_Layer/UFS/sfc_nst.f90',
    'ccpp/physics/physics/SFC_Layer/UFS/sfc_nst_post.f90',
    'ccpp/physics/physics/SFC_Models/Land/RUC/lsm_ruc.F90',
    'ccpp/physics/physics/SFC_Models/SeaIce/CICE/sfc_cice.f',
    'ccpp/physics/physics/SFC_Models/Land/sfc_land.F90',
    'ccpp/physics/physics/SFC_Models/Land/Noah/lsm_noah.f',
    'ccpp/physics/physics/SFC_Models/Land/Noahmp/noahmpdrv.F90',
    'ccpp/physics/physics/SFC_Models/Lake/Flake/flake_driver.F90',
    'ccpp/physics/physics/SFC_Models/Lake/CLM/clm_lake.f90',
    'ccpp/physics/physics/SFC_Models/Ocean/UFS/sfc_ocean.F',
    'ccpp/physics/physics/SFC_Models/SeaIce/CICE/sfc_sice.f',
    'ccpp/physics/physics/smoke_dust/rrfs_smoke_wrapper.F90',
    'ccpp/physics/physics/smoke_dust/rrfs_smoke_postpbl.F90',
    'ccpp/physics/physics/tools/get_prs_fv3.F90',
    'ccpp/physics/physics/tools/get_phi_fv3.F90'
    ]

# Default build dir, relative to current working directory,
# if not specified as command-line argument
DEFAULT_BUILD_DIR = 'scm/bin'

# Auto-generated makefile/cmakefile snippets that contain all type definitions
TYPEDEFS_MAKEFILE   = '{build_dir}/ccpp/physics/CCPP_TYPEDEFS.mk'
TYPEDEFS_CMAKEFILE  = '{build_dir}/ccpp/physics/CCPP_TYPEDEFS.cmake'
TYPEDEFS_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_TYPEDEFS.sh'

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE = '{build_dir}/ccpp/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE = '{build_dir}/ccpp/physics/CCPP_SCHEMES.cmake'
SCHEMES_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_SCHEMES.sh'

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE = '{build_dir}/ccpp/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE = '{build_dir}/ccpp/physics/CCPP_CAPS.cmake'
CAPS_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_CAPS.sh'

# Directory where to put all auto-generated physics caps
CAPS_DIR = '{build_dir}/ccpp/physics/physics'

# Directory where the suite definition files are stored
SUITES_DIR = 'ccpp/suites'

# Directory where to write static API to
STATIC_API_DIR = 'scm/src/'
STATIC_API_CMAKEFILE = 'scm/src/CCPP_STATIC_API.cmake'
STATIC_API_SOURCEFILE = 'scm/src/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
METADATA_HTML_OUTPUT_DIR = 'ccpp/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = 'ccpp/physics/CCPP_VARIABLES_SCM.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = 'ccpp/framework/doc/DevelopersGuide/CCPP_VARIABLES_SCM.tex'
