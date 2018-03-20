
# All CCPP schemes are defined here.
#
# This file is auto-generated using ccpp_prebuild.py
# at compile time, do not edit manually.
#
SCHEMES_F =

SCHEMES_F90 = \
	   ./GFS_layer/GFS_initialize_scm.F90 \
	   ./physics/GFS_rrtmg_pre.F90 \
	   ./physics/rrtmg_sw_pre.F90 \
	   ./physics/rrtmg_sw_post.F90 \
	   ./physics/rrtmg_lw_pre.F90 \
	   ./physics/rrtmg_lw_post.F90 \
	   ./physics/GFS_rrtmg_post.F90

SCHEMES_f = \
	   ./physics/radsw_main.f \
	   ./physics/radlw_main.f \
	   ./physics/sfc_sice.f \
	   ./physics/dcyc2.f \
	   ./physics/sfc_drv.f \
	   ./physics/sfc_diff.f \
	   ./physics/GFS_surface_loop_control.f \
	   ./physics/sfc_nst.f \
	   ./physics/sfc_diag.f \
	   ./physics/moninedmf.f \
	   ./physics/gwdps.f \
	   ./physics/rayleigh_damp.f \
	   ./physics/ozphys.f \
	   ./physics/mfdeepcnv.f \
	   ./physics/gwdc.f \
	   ./physics/mfshalcnv.f \
	   ./physics/cnvc90.f \
	   ./physics/gscond.f \
	   ./physics/precpd.f

SCHEMES_f90 = \
	   ./physics/GFS_phys_time_vary.f90 \
	   ./physics/GFS_rad_time_vary.f90 \
	   ./physics/GFS_suite_interstitial.ccpp.f90 \
	   ./physics/get_prs_fv3.f90 \
	   ./physics/GFS_surface_generic.f90 \
	   ./physics/GFS_PBL_generic.f90 \
	   ./physics/GFS_DCNV_generic.f90 \
	   ./physics/GFS_zhao_carr_pre.f90 \
	   ./physics/GFS_SCNV_generic.f90 \
	   ./physics/GFS_MP_generic_pre.f90 \
	   ./physics/GFS_calpreciptype.f90 \
	   ./physics/GFS_MP_generic_post.f90