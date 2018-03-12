SHELL = /bin/sh

inside_nems := $(wildcard ../../../conf/configure.nems)
ifneq ($(strip $(inside_nems)),)
    include ../../../conf/configure.nems
else
    exist_configure_fv3 := $(wildcard ../conf/configure.fv3)
    ifneq ($(strip $(exist_configure_fv3)),)
        include ../conf/configure.fv3
    else
        $(error "../conf/configure.fv3 file is missing. Run ./configure")
    endif
    $(info )
    $(info Build standalone FV3 gfsphysics ...)
    $(info )
endif

ifneq (,$(findstring MACOSX,$(CPPDEFS)))
LIBRARY = libgfsphys.dylib
else
LIBRARY = libgfsphys.so
endif
VER_MAJOR = 1
VER_MINOR = 0
VER_PATCH = 0

FFLAGS   += -I../fms -I../fms/include -fPIC

CPPDEFS += -DNEW_TAUCTMAX -DSMALL_PE -DNEMS_GSM

ifneq (,$(findstring CCPP,$(CPPDEFS)))
    # IPD with full CCPP - not making use of IPD steps
    IPD_DRIVER_CAP = ./IPD_layer/IPD_driver_cap.F90
    IPD_CCPP_DRIVER = ./IPD_layer/IPD_CCPP_driver.F90
else
    # IPD without CCPP
    IPD_DRIVER_CAP =
    IPD_CCPP_DRIVER =
endif

ifneq (,$(findstring CCPP,$(CPPDEFS)))
    GFS_SUITE_INTERSTITIAL = ./physics/GFS_suite_interstitial.ccpp.f90
    GFS_PHYSICS_DRIVER =
    GFS_RADIATION_DRIVER = ./GFS_layer/GFS_radiation_driver.F90
else
    GFS_SUITE_INTERSTITIAL = ./physics/GFS_suite_interstitial.ipd.f90
    GFS_PHYSICS_DRIVER = ./GFS_layer/GFS_physics_driver.F90
    GFS_RADIATION_DRIVER = ./GFS_layer/GFS_radiation_driver.F90
endif

SRCS_f   =  \
	   ./physics/cnvc90.f                                                        \
	   ./physics/co2hc.f                                                         \
	   ./physics/date_def.f                                                      \
	   ./physics/dcyc2.f                                                         \
	   ./physics/dcyc2.pre.rad.f                                                 \
	   ./physics/efield.f                                                        \
	   ./physics/get_prs.f                                                       \
	   ./physics/gfs_phy_tracer_config.f                                         \
	   ./physics/gocart_tracer_config_stub.f                                     \
	   ./physics/gscond.f                                                        \
	   ./physics/gscondp.f                                                       \
	   ./physics/gwdc.f                                                          \
	   ./physics/gwdps.f                                                         \
	   ./physics/h2o_def.f                                                       \
	   ./physics/h2oc.f                                                          \
	   ./physics/h2ohdc.f                                                        \
	   ./physics/h2ophys.f                                                       \
	   ./physics/ideaca.f                                                        \
	   ./physics/idea_co2.f                                                      \
	   ./physics/idea_composition.f                                              \
	   ./physics/idea_dissipation.f                                              \
	   ./physics/idea_h2o.f                                                      \
	   ./physics/idea_ion.f                                                      \
	   ./physics/idea_o2_o3.f                                                    \
	   ./physics/idea_phys.f                                                     \
	   ./physics/idea_solar_heating.f                                            \
	   ./physics/idea_tracer.f                                                   \
	   ./physics/iounitdef.f                                                     \
	   ./physics/lrgsclr.f                                                       \
	   ./physics/mersenne_twister.f                                              \
	   ./physics/mfdeepcnv.f                                                     \
	   ./physics/mfpbl.f                                                         \
	   ./physics/mfshalcnv.f                                                     \
	   ./physics/module_bfmicrophysics.f                                         \
	   ./physics/moninedmf.f                                                     \
	   ./physics/moninp.f                                                        \
	   ./physics/moninp1.f                                                       \
	   ./physics/moninq.f                                                        \
	   ./physics/moninq1.f                                                       \
	   ./physics/moninshoc.f                                                     \
	   ./physics/mstadb.f                                                        \
	   ./physics/mstadbtn.f                                                      \
	   ./physics/mstcnv.f                                                        \
	   ./physics/namelist_soilveg.f                                              \
	   ./physics/ozne_def.f                                                      \
	   ./physics/ozphys.f                                                        \
	   ./physics/ozphys_2015.f                                                   \
	   ./physics/physparam.f                                                     \
	   ./physics/precpd.f                                                        \
	   ./physics/precpd_shoc.f                                                   \
	   ./physics/precpdp.f                                                       \
	   ./physics/progt2.f                                                        \
	   ./physics/progtm_module.f                                                 \
	   ./physics/rad_initialize.f                                                \
	   ./physics/radiation_aerosols.f                                            \
	   ./physics/radiation_astronomy.f                                           \
	   ./physics/radiation_clouds.f                                              \
	   ./physics/radiation_gases.f                                               \
	   ./physics/radiation_surface.f                                             \
	   ./physics/radlw_datatb.f                                                  \
	   ./physics/radlw_main.f                                                    \
	   ./physics/radlw_param.f                                                   \
	   ./physics/radsw_datatb.f                                                  \
	   ./physics/radsw_main.f                                                    \
	   ./physics/radsw_param.f                                                   \
	   ./physics/rascnvv2.f                                                      \
	   ./physics/rayleigh_damp.f                                                 \
	   ./physics/rayleigh_damp_mesopause.f                                       \
	   ./physics/sascnv.f                                                        \
	   ./physics/sascnvn.f                                                       \
	   ./physics/set_soilveg.f                                                   \
	   ./physics/GFS_surface_loop_control.f                                      \
	   ./physics/sfc_cice.f                                                      \
	   ./physics/sfc_diag.f                                                      \
	   ./physics/sfc_diff.f                                                      \
	   ./physics/sfc_drv.f                                                       \
	   ./physics/sfc_land.f                                                      \
	   ./physics/sfc_nst.f                                                       \
	   ./physics/sfc_ocean.f                                                     \
	   ./physics/sfc_sice.f                                                      \
	   ./physics/sflx.f                                                          \
	   ./physics/shalcnv.f                                                       \
	   ./physics/shalcv.f                                                        \
	   ./physics/shalcv_opr.f                                                    \
	   ./physics/tracer_const_h.f                                                \
	   ./physics/tridi.f                                                         \
	   ./physics/tridi2t3.f

SRCS_f90 = \
	   ./physics/GFS_calpreciptype.f90                                           \
	   ./physics/GFS_MP_generic_post.f90                                         \
	   ./physics/GFS_MP_generic_pre.f90                                          \
	   ./physics/GFS_zhao_carr_pre.f90                                           \
	   ./physics/GFS_rad_time_vary.f90                                           \
	   ./physics/GFS_radupdate.f90                                               \
	   ./physics/cs_conv.f90                                                     \
	   ./physics/funcphys.f90                                                    \
	   ./physics/gcm_shoc.f90                                                    \
	   ./physics/gcycle.f90                                                      \
	   ./physics/get_prs_fv3.f90                                                 \
	   ./physics/GFS_DCNV_generic.f90                                            \
	   ./physics/GFS_SCNV_generic.f90                                            \
	   ./physics/GFS_PBL_generic.f90                                             \
	   $(GFS_SUITE_INTERSTITIAL)                                                 \
	   ./physics/GFS_phys_time_vary.f90                                          \
	   ./physics/GFS_stochastics.f90                                             \
	   ./physics/GFS_surface_generic.f90                                         \
	   ./physics/h2ointerp.f90                                                   \
	   ./physics/m_micro_driver.f90                                              \
	   ./physics/module_nst_model.f90                                            \
	   ./physics/module_nst_parameters.f90                                       \
	   ./physics/module_nst_water_prop.f90                                       \
	   ./physics/ozinterp.f90                                                    \
	   ./physics/physcons.f90                                                    \
	   ./physics/radcons.f90                                                     \
	   ./physics/wam_f107_kp_mod.f90                                             \
	   ./physics/GFS_debug.f90

SRCS_F   = \
	   ./physics/aer_cloud.F                                                     \
	   ./physics/cldmacro.F                                                      \
	   ./physics/cldwat2m_micro.F                                                \
	   ./physics/machine.F                                                       \
	   ./physics/sfcsub.F                                                        \
	   ./physics/num_parthds.F                                                   \
	   ./physics/wv_saturation.F

SRCS_F90 = \
	   ./physics/GFDL_parse_tracers.F90                                          \
	   ./GFS_layer/GFS_abstraction_layer.F90                                     \
	   ./GFS_layer/GFS_diagnostics.F90                                           \
	   ./GFS_layer/GFS_driver.F90                                                \
	   ./physics/GFS_rrtmg_pre.F90                                               \
	   ./physics/GFS_rrtmg_post.F90                                              \
	   ./physics/rrtmg_sw_pre.F90                                                \
	   ./physics/rrtmg_sw_post.F90                                               \
	   ./physics/rrtmg_lw_pre.F90                                                \
	   ./physics/rrtmg_lw_post.F90                                               \
	   $(GFS_PHYSICS_DRIVER)                                                     \
	   $(GFS_RADIATION_DRIVER)                                                   \
	   ./GFS_layer/GFS_restart.F90                                               \
	   ./GFS_layer/GFS_typedefs.F90                                              \
	   ./IPD_layer/IPD_driver.F90                                                \
	   ./IPD_layer/IPD_typedefs.F90                                              \
	   $(IPD_CCPP_DRIVER)

SRCS_c   =

ifneq (,$(findstring CCPP,$(CPPDEFS)))
include ./CCPP_CAPS.mk
CAPS_F90 += \
	   $(IPD_DRIVER_CAP)
#include ./CCPP_SCHEMES.mk
# Schemes not yet used, all hardcoded above
else
CAPS_F90 = \
	   $(IPD_DRIVER_CAP)
endif

DEPEND_FILES = $(SRCS_f) $(SRCS_f90) $(SRCS_F) $(SRCS_F90) $(CAPS_F90)

OBJS_f   = $(SRCS_f:.f=.o)
OBJS_f90 = $(SRCS_f90:.f90=.o)
OBJS_F   = $(SRCS_F:.F=.o)
OBJS_F90 = $(SRCS_F90:.F90=.o)
OBJS_c   = $(SRCS_c:.c=.o)

OBJS = $(OBJS_f) $(OBJS_f90) $(OBJS_F) $(OBJS_F90) $(OBJS_c)

CAPS = $(CAPS_F90:.F90=.o)

all default: depend $(LIBRARY)

ifneq (,$(findstring MACOSX,$(CPPDEFS)))
LIBRARY_FULL_NAME = $(subst .dylib,.$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH).dylib,$(LIBRARY))
$(LIBRARY): $(OBJS) $(CAPS)
	$(FC) -shared $(OBJS) $(CAPS) $(LDFLAGS) $(NCEPLIBS) -o $(LIBRARY_FULL_NAME)
	ln -sf $(LIBRARY_FULL_NAME) $(LIBRARY)
	ln -sf $(LIBRARY_FULL_NAME) $(subst .dylib,.$(VER_MAJOR).dylib,$(LIBRARY))
else
$(LIBRARY): $(OBJS) $(CAPS)
	$(FC) -shared -Wl,-soname,$(LIBRARY).$(VER_MAJOR) $(OBJS) $(CAPS) $(LDFLAGS) $(NCEPLIBS) -o $(LIBRARY).$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH)
	ln -sf $(LIBRARY).$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH) $(LIBRARY)
	ln -sf $(LIBRARY).$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH) $(LIBRARY).$(VER_MAJOR)
endif

# this is the place to override default (implicit) compilation rules
# and create specific (explicit) rules

# this has no effect, because radiation_aerosols.o is in physics, not in gfsphys
./radiation_aerosols.o : ./gfsphys/radiation_aerosols.f
	$(FC) $(CPPDEFS) $(FFLAGS) $(OTHER_FFLAGS) -xCORE-AVX-I -c $< -o $@

# Reduce optimization for sfc_sice for bit-for-bit reproducibility
FFLAGS_REDUCED_OPT=$(subst -O2,-O1,$(subst -xCORE-AVX2,-xCORE-AVX-I,$(FFLAGS)))
./physics/sfc_sice.o : ./physics/sfc_sice.f
	$(FC) $(CPPDEFS) $(FFLAGS_REDUCED_OPT) $(OTHER_FFLAGS) -c $< -o $@

./GFS_layer/GFS_diagnostics.o : ./GFS_layer/GFS_diagnostics.F90
	$(FC) $(CPPDEFS) $(FFLAGS) $(OTHER_FFLAGS) -O0 -c $< -o $@

ifneq (,$(findstring PGIFIX,$(CPPDEFS)))
$(CAPS):
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHER_FFLAGS) -c $< -o $@
	# Apply a fix specific to the PGI compiler (rename objects in cap object files)
	./pgifix.py $@
else
$(CAPS):
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHER_FFLAGS) -c $< -o $@
endif

# Do preprocessing of the IPD-CCPP driver in two steps to be
# able to look at the actual .f90 file that gets compiled
ifneq (,$(findstring CCPP,$(CPPDEFS)))
./IPD_layer/IPD_CCPP_driver.o: ./IPD_layer/IPD_CCPP_driver.F90
	$(CPP) $(CPPDEFS) $(CPPFLAGS) $< > $*.tmp.f90
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -c $*.tmp.f90 -o $@
endif

.PHONY: clean
clean:
	@echo "Cleaning gfsphysics  ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o */*.o *.tmp.f90 */*.tmp.f90 *.mod *.lst *.i depend
	$(RM) -f $(LIBRARY).$(VER_MAJOR)
	$(RM) -f $(LIBRARY).$(VER_MAJOR).$(VER_MINOR).$(VER_PATCH)

MKDEPENDS = ../mkDepends.pl
include ../conf/make.rules

include ./depend

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
   -include depend
endif
