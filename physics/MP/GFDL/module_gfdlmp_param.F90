! #########################################################################################
! #########################################################################################
module module_gfdlmp_param
  use machine, only: kind_phys
  implicit none

  ! #######################################################################################
  ! Data container for GFDL MP runtime configurations information (i.e. Namelist)
  ! #######################################################################################
  type ty_gfdlmp_config
     ! V3 configuration
     integer :: ntimes    = 1  !< cloud microphysics sub cycles (v3)
     integer :: nconds    = 1  !< condensation sub cycles (v3)
     integer :: inflag    = 1  !< ice nucleation scheme (v3)
                               !< 1: Hong et al. (2004)
                               !< 2: Meyers et al. (1992)
                               !< 3: Meyers et al. (1992)
                               !< 4: Cooper (1986)
                               !< 5: Fletcher (1962)
     integer :: igflag    = 3  !< ice generation scheme (v3)
                               !< 1: WSM6
                               !< 2: WSM6 with 0 at 0 C
                               !< 3: WSM6 with 0 at 0 C and fixed value at - 10 C
                               !< 4: combination of 1 and 3
     integer :: ifflag    = 1  !< ice fall scheme (v3)
                               !< 1: Deng and Mace (2008)
                               !< 2: Heymsfield and Donner (1990)
     integer :: rewflag   = 1  !< cloud water effective radius scheme (v3)
                               !< 1: Martin et al. (1994)
                               !< 2: Martin et al. (1994), GFDL revision
                               !< 3: Kiehl et al. (1994)
                               !< 4: effective radius
     integer :: rerflag   = 1  !< rain effective radius scheme (v3)
                               !< 1: effective radius
     integer :: resflag   = 1  !< snow effective radius scheme (v3)
                               !< 1: effective radius
     integer :: regflag   = 1  !< graupel effective radius scheme (v3)
                               !< 1: effective radius
     integer :: radr_flag = 1  !< radar reflectivity for rain (v3)
                               !< 1: Mark Stoelinga (2005)
                               !< 2: Smith et al. (1975), Tong and Xue (2005)
                               !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
     integer :: rads_flag = 1  !< radar reflectivity for snow (v3)
                               !< 1: Mark Stoelinga (2005)
                               !< 2: Smith et al. (1975), Tong and Xue (2005)
                               !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
     integer :: radg_flag = 1  !< radar reflectivity for graupel (v3)
                               !< 1: Mark Stoelinga (2005)
                               !< 2: Smith et al. (1975), Tong and Xue (2005)
                               !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
     integer :: sedflag   = 1  !< sedimentation scheme (v3)
                               !< 1: implicit scheme
                               !< 2: explicit scheme
                               !< 3: lagrangian scheme
                               !< 4: combined implicit and lagrangian scheme
     integer :: vdiffflag = 1  !< wind difference scheme in accretion (v3)
                               !< 1: Wisner et al. (1972)
                               !< 2: Mizuno (1990)
                               !< 3: Murakami (1990)     
     !
     real :: tice_mlt    = 273.16  !< can set ice melting temperature to 268 based on observation (Kay et al. 2016) (K) (v3)
     real :: tau_gmlt    = 600.0   !< graupel melting time scale (s) (v3)
     real :: tau_wbf     = 300.0   !< graupel melting time scale (s) (v3)
     real :: c_pracw     = 0.8     !< cloud water to rain accretion efficiency (v3)
     real :: c_praci     = 1.0     !< cloud ice to rain accretion efficiency (v3)
     real :: c_pracs     = 1.0     !< snow to rain accretion efficiency (v3)
     real :: c_psacw     = 1.0     !< cloud water to snow accretion efficiency (v3)
     real :: c_pgacw     = 1.0     !< cloud water to graupel accretion efficiency  (v3)
     real :: c_psacr     = 1.0     !< rain to snow accretion efficiency  (v3)
     real :: c_pgaci     = 0.05    !< cloud ice to graupel accretion efficiency (was 0.1 in ZETAC) (v3)
     real :: c_pgacr     = 1.0     !< rain to graupel accretion efficiency  (v3)
     real :: is_fac      = 0.2     !< cloud ice sublimation temperature factor (v3)
     real :: ss_fac      = 0.2     !< snow sublimation temperature factor (v3)
     real :: gs_fac      = 0.2     !< graupel sublimation temperature factor (v3)
     real :: rh_fac_evap = 10.0    !< cloud water evaporation relative humidity factor (v3)
     real :: rh_fac_cond = 10.0    !< cloud water condensation relative humidity factor (v3)
     real :: sed_fac     = 1.0     !< coefficient for sedimentation fall, scale from 1.0 (implicit) to 0.0 (lagrangian)  (v3)
     real :: vw_fac      = 1.0     !< (v3)
     real :: vw_max      = 0.01    !< maximum fall speed for cloud water (m/s) (V3)
     real :: xr_a        = 0.25    !< p value in Xu and Randall (1996) (v3)
     real :: xr_b        = 100.0   !< alpha_0 value in Xu and Randall (1996) (v3)
     real :: xr_c        = 0.49    !< gamma value in Xu and Randall (1996) (v3)
     real :: te_err      = 1.e-5   !< 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time (v3)
     real :: tw_err      = 1.e-8   !< 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time (v3)
     real :: rh_thres    = 0.75    !< minimum relative humidity for cloud fraction (v3)
     real :: rhc_cevap   = 0.85    !< maximum relative humidity for cloud water evaporation (v3)
     real :: rhc_revap   = 0.85    !< maximum relative humidity for rain evaporation (v3)
     real :: f_dq_p      = 1.0     !< cloud fraction adjustment for supersaturation (v3)
     real :: f_dq_m      = 1.0     !< cloud fraction adjustment for undersaturation (v3)
     real :: fi2s_fac    = 1.0     !< maximum sink of cloud ice to form snow: 0-1 (v3)
     real :: fi2g_fac    = 1.0     !< maximum sink of cloud ice to form graupel: 0-1 (v3)
     real :: fs2g_fac    = 1.0     !< maximum sink of snow to form graupel: 0-1 (v3)
     real :: beta        = 1.22    !< defined in Heymsfield and Mcfarquhar (1996) (v3)
     real :: n0w_exp     = 41      !< intercept parameter (exponent) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994) (v3)
     real :: n0i_exp     = 18      !< intercept parameter (exponent) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015) (v3)
     real :: n0r_exp     = 6       !< intercept parameter (exponent) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948) (v3)
     real :: n0s_exp     = 6       !< intercept parameter (exponent) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958) (v3)
     real :: n0g_exp     = 6       !< intercept parameter (exponent) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979) (v3)
     real :: n0h_exp     = 4       !< intercept parameter (exponent) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975) (v3)
     real :: muw         = 6.0     !< shape parameter of cloud water in Gamma distribution (Martin et al. 1994) (v3)
     real :: mui         = 3.35    !< shape parameter of cloud ice in Gamma distribution (McFarquhar et al. 2015) (v3)
     real :: mur         = 1.0     !< shape parameter of rain in Gamma distribution (Marshall and Palmer 1948) (v3)
     real :: mus         = 1.0     !< shape parameter of snow in Gamma distribution (Gunn and Marshall 1958) (v3)
     real :: mug         = 1.0     !< shape parameter of graupel in Gamma distribution (Houze et al. 1979) (v3)
     real :: muh         = 1.0     !< shape parameter of hail in Gamma distribution (Federer and Waldvogel 1975) (v3)
     real :: alinw       = 3.e7    !< "a" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990) (v3)
     real :: alini       = 7.e2    !< "a" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990) (v3)
     real :: alinr       = 842.0   !< "a" in Lin et al. (1983) for rain (Liu and Orville 1969) (v3)
     real :: alins       = 4.8     !< "a" in Lin et al. (1983) for snow (straka 2009) (v3)
     real :: aling       = 1.0     !< "a" in Lin et al. (1983), similar to a, but for graupel (Pruppacher and Klett 2010) (v3)
     real :: alinh       = 1.0     !< "a" in Lin et al. (1983), similar to a, but for hail (Pruppacher and Klett 2010) (v3)
     real :: blinw       = 2.0     !< "b" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990) (v3)
     real :: blini       = 1.0     !< "b" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990) (v3)
     real :: blinr       = 0.8     !< "b" in Lin et al. (1983) for rain (Liu and Orville 1969) (v3)
     real :: blins       = 0.25    !< "b" in Lin et al. (1983) for snow (straka 2009) (v3)
     real :: bling       = 0.5     !< "b" in Lin et al. (1983), similar to b, but for graupel (Pruppacher and Klett 2010) (v3)
     real :: blinh       = 0.5     !< "b" in Lin et al. (1983), similar to b, but for hail (Pruppacher and Klett 2010) (v3)
     real :: rewfac      = 1.0     !< this is a tuning parameter to compromise the inconsistency between (v3)
                                   !< GFDL MP's PSD and cloud water radiative property's PSD assumption.
                                   !< after the cloud water radiative property's PSD is rebuilt,
                                   !< this parameter should be 1.0.
     real :: reifac      = 1.0     !< this is a tuning parameter to compromise the inconsistency between (v3)
                                   !< GFDL MP's PSD and cloud ice radiative property's PSD assumption.
                                   !< after the cloud ice radiative property's PSD is rebuilt,
                                   !< this parameter should be 1.0.
     real :: n0w_sig     = 1.1     !< intercept parameter (significand) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994) (v3)
     real :: n0i_sig     = 1.3     !< intercept parameter (significand) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015) (v3)
     real :: n0r_sig     = 8.0     !< intercept parameter (significand) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948) (v3)
     real :: n0s_sig     = 3.0     !< intercept parameter (significand) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958) (v3)
     real :: n0g_sig     = 4.0     !< intercept parameter (significand) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)(v3)
     real :: n0h_sig     = 4.0     !< intercept parameter (significand) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975) (v3)
     
     ! v3 cloud microphysics switches
     logical :: const_vw             = .false. !< if .true., the constants are specified by v * _fac (v3)
     logical :: liq_ice_combine      = .false. !< combine all liquid water, combine all solid water (v3)
     logical :: snow_grauple_combine = .true.  !< combine snow and graupel (v3)
     logical :: use_rhc_cevap        = .false. !< cap of rh for cloud water evaporation (V3)
     logical :: use_rhc_revap        = .false. !< cap of rh for rain evaporation (V3)
     logical :: do_cld_adj           = .false. !< do cloud fraction adjustment  (v3)
     logical :: do_evap_timescale    = .true.  !< whether to apply a timescale to evaporation (v3)
     logical :: do_cond_timescale    = .false. !< whether to apply a timescale to condensation (v3)
     logical :: consv_checker        = .false. !< turn on energy and water conservation checker (v3)
     logical :: do_warm_rain_mp      = .false. !< do warm rain cloud microphysics only (v3)
     logical :: do_wbf               = .false. !< do Wegener Bergeron Findeisen process (v3)
     logical :: do_psd_water_fall    = .false. !< calculate cloud water terminal velocity based on PSD (v3)
     logical :: do_psd_ice_fall      = .false. !< calculate cloud ice terminal velocity based on PSD (v3)
     logical :: do_psd_water_num     = .false. !< calculate cloud water number concentration based on PSD (v3)
     logical :: do_psd_ice_num       = .false. !< calculate cloud ice number concentration based on PSD (v3)
     logical :: do_new_acc_water     = .false. !< perform the new accretion for cloud water (v3)
     logical :: do_new_acc_ice       = .false. !< perform the new accretion for cloud ice (v3)
     logical :: cp_heating           = .false. !< update temperature based on constant pressure (v3)
     logical :: delay_cond_evap      = .false. !< do condensation evaporation only at the last time step (v3)
     logical :: do_subgrid_proc      = .true.  !< do temperature sentive high vertical resolution processes (v3)
     logical :: fast_fr_mlt          = .true.  !< do freezing and melting in fast microphysics (v3)
     logical :: fast_dep_sub         = .true.  !< do deposition and sublimation in fast microphysics (v3)
     logical :: do_sedi_uv           = .true.  !< transport of horizontal momentum in sedimentation (v3)
     logical :: do_sedi_melt         = .true.  !< melt cloud ice, snow, and graupel during sedimentation (v3)
     
     ! v1/v3 configuration
     real :: cld_min     = 0.05    !< minimum cloud fraction (v1/v3)
     real :: tice        = 273.16  !< set tice = 165. to trun off ice - phase phys (kessler emulator) (v3=273.15)
     real :: t_min       = 178.    !< min temp to freeze - dry all water vapor (v1/v3)
     real :: t_sub       = 184.    !< min temp for sublimation of cloud ice (v1/v3)
     real :: mp_time     = 150.0   !< maximum micro - physics time step (sec) (v1/v3)
     real :: rh_inc      = 0.25    !< rh increment for complete evaporation of cloud water and cloud ice (v1/v3)
     real :: rh_inr      = 0.25    !< rh increment for minimum evaporation of rain (v1/v3)
     real :: rh_ins      = 0.25    !< rh increment for sublimation of snow (v1/v3)
     real :: tau_r2g     = 900.    !< rain freezing during fast_sat (v1/v3)
     real :: tau_smlt    = 900.    !< snow melting (v1/v3)
     real :: tau_g2r     = 600.    !< graupel melting to rain (v1/v3)
     real :: tau_imlt    = 600.    !< cloud ice melting (v1/v3)
     real :: tau_i2s     = 1000.   !< cloud ice to snow auto-conversion (v1/v3)
     real :: tau_l2r     = 900.    !< cloud water to rain auto-conversion (v1/v3)
     real :: tau_v2l     = 150.    !< water vapor to cloud water (condensation) (v1/v3)
     real :: tau_l2v     = 300.    !< cloud water to water vapor (evaporation) (v1/v3)
     real :: tau_g2v     = 900.    !< graupel sublimation  (v1/v3)
     real :: tau_v2g     = 21600.  !< graupel deposition -- make it a slow process (v1)
     real :: dw_land     = 0.20    !< base value for subgrid deviation / variability over land (v1/v3)
     real :: dw_ocean    = 0.10    !< base value for subgrid deviation / variability over ocean (v1/v3)
     real :: ccn_o       = 90.     !< prescribed ccn over ocean (cm^ - 3) (v1/v3)
     real :: ccn_l       = 270.    !< prescribed ccn over land (cm^ - 3) (v1/v3)
     real :: rthresh     = 10.0e-6 !< critical cloud drop radius (micro m) (v3=20.0e-6)
     real :: sat_adj0    = 0.90    !< adjustment factor (0: no, 1: full) during fast_sat_adj (v1/v3)
     real :: qc_crt      = 5.0e-8  !< mini condensate mixing ratio to allow partial cloudiness (v1)
     real :: qi_lim      = 1.      !< cloud ice limiter to prevent large ice build up (v1/v3)
     real :: ql_mlt      = 2.0e-3  !< max value of cloud water allowed from melted cloud ice (v1/v3)
     real :: qs_mlt      = 1.0e-6  !< max cloud water due to snow melt (v1/v3)
     real :: ql_gen      = 1.0e-3  !< max cloud water generation during remapping step if fast_sat_adj = .t.  (v1/v3)
     real :: qi_gen      = 1.82e-6 !< max cloud ice generation during remapping step (V1. Computed internally in V3) (v1/v3)
     real :: ql0_max     = 2.0e-3  !< max cloud water condensate value (auto converted to rain) (v1/v3)
     real :: qi0_max     = 1.0e-4  !< max cloud ice condensatevalue (by other sources) (v1/v3)
     real :: qi0_crt     = 1.0e-4  !< cloud ice to snow autoconversion threshold (was 1.e-4);
                                   !< qi0_crt is highly dependent on horizontal resolution (v1/v3)
     real :: qr0_crt     = 1.0e-4  !< rain to snow or graupel/hail threshold (v1)
                                   !< lfo used * mixing ratio * = 1.e-4 (hail in lfo)
     real :: qs0_crt     = 1.0e-3  !< snow to graupel density threshold (0.6e-3 in purdue lin scheme) (v1/v3)
     real :: c_paut      = 0.55    !< autoconversion cloud water to rain (use 0.5 to reduce autoconversion) (v1/v3)
     real :: c_psaci     = 0.02    !< accretion: cloud ice to snow (was 0.1 in zetac) (v3 = 0.05)
     real :: c_piacr     = 5.0     !< accretion: rain to ice: (v1)
     real :: c_cracw     = 0.9     !< rain accretion efficiency (v1)
     real :: c_pgacs     = 2.0e-3  !< snow to graupel "accretion" eff. (was 0.1 in zetac) (v3 = 0.01)
     real :: alin        = 842.0   !< "a" in lin1983 (v1)
     real :: clin        = 4.8     !< "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs) (v1)
     real :: vi_fac      = 1.0     !< if const_vi: 1 / 3 (v1/v3)
     real :: vs_fac      = 1.0     !< if const_vs: 1. (v1/v3)
     real :: vg_fac      = 1.0     !< if const_vg: 2. (v1/v3)
     real :: vr_fac      = 1.0     !< if const_vr: 4. (v1/v3)
     real :: vi_max      = 0.5     !< max fall speed for ice (v1/v3)
     real :: vs_max      = 5.0     !< max fall speed for snow (v1/v3)
     real :: vg_max      = 8.0     !< max fall speed for graupel (v1/v3)
     real :: vr_max      = 12.     !< max fall speed for rain (v1/v3)
     real :: rewmin      = 5.0     !< minimum effective radii (liquid) (v1/v3)
     real :: rewmax      = 10.0    !< maximum effective radii (liquid) (v3= = 15.0)
     real :: reimin      = 10.0    !< minimum effective radii (ice) (v1/v3)
     real :: reimax      = 150.0   !< maximum effective radii (ice) (v1/v3)
     real :: rermin      = 10.0    !< minimum effective radii (rain) (v3=15.0)
     real :: rermax      = 10000.0 !< maximum effective radii (rain) (v1/v3)
     real :: resmin      = 150.0   !< minimum effective radii (snow) (v1.v3)
     real :: resmax      = 10000.0 !< maximum effective radii (snow) (v1/v3)
     real :: regmin      = 300.0   !< minimum effective radii (graupel) (v3=150.)
     real :: regmax      = 10000.0 !< maximum effective radii (graupel) (v1/v3)
     
     ! v1/v3 cloud microphysics switches
     logical :: const_vi             = .false. !< if .true., the constants are specified by v * _fac (v1/v3)
     logical :: const_vs             = .false. !< if .true., the constants are specified by v * _fac (v1/v3)
     logical :: const_vg             = .false. !< if .true., the constants are specified by v * _fac (v1/v3)
     logical :: const_vr             = .false. !< if .true., the constants are specified by v * _fac (v1/v3)
     logical :: fast_sat_adj         = .false. !< has fast saturation adjustments (v1)
     logical :: z_slope_liq          = .true.  !< use linear mono slope for autocconversions (v1/v3)
     logical :: z_slope_ice          = .false. !< use linear mono slope for autocconversions (v3 = .true.)
     logical :: use_ccn              = .false. !< must be true when prog_ccn is false (v1)
     logical :: use_ppm              = .false. !< use ppm fall scheme (v1)
     logical :: mono_prof            = .true.  !< perform terminal fall with mono ppm scheme (v1)
     logical :: mp_print             = .false. !< cloud microphysics debugging printout (v1)
     logical :: do_hail              = .false. !< use hail parameters instead of graupel (v1/v3)
     logical :: do_qa                = .true.  !< do inline cloud fraction (v1/v3)
     logical :: rad_snow             = .true.  !< consider snow in cloud fraciton calculation (v1/v3)
     logical :: rad_graupel          = .true.  !< consider graupel in cloud fraction calculation (v1/v3)
     logical :: rad_rain             = .true.  !< consider rain in cloud fraction calculation (v1/v3)
     logical :: de_ice               = .false. !< to prevent excessive build - up of cloud ice from external sources (v1)
     logical :: sedi_transport       = .true.  !< transport of momentum in sedimentation (v1)
     logical :: do_sedi_w            = .false. !< transport of vertical motion in sedimentation (v3 = .true.)
     logical :: do_sedi_heat         = .true.  !< transport of heat in sedimentation (v1/v3)
     logical :: prog_ccn             = .false. !< do prognostic ccn (yi ming's method) (v1/v3)
     logical :: fix_negative         = .false. !< fix negative water species (v3 = .true.)
     logical :: tintqs               = .false. !< (v1)

     ! v1/v3 integers
     integer :: icloud_f  = 0  !< GFDL cloud scheme (v1/v3)
                               !< 0: subgrid variability based scheme
                               !< 1: same as 0, but for old fvgfs implementation
                               !< 2: binary cloud scheme
                               !< 3: extension of 0
     integer :: irain_f   = 0  !< cloud water to rain auto conversion scheme (v1/v3)
                               !< 0: subgrid variability based scheme
                               !< 1: no subgrid varaibility
     integer :: reiflag   = 1  !< cloud ice effective radius scheme
                               !< 1: Heymsfield and Mcfarquhar (1996) (v3 = 5)
                               !< 2: Donner et al. (1997)
                               !< 3: Fu (2007)
                               !< 4: Kristjansson et al. (2000)
                               !< 5: Wyser (1998)
                               !< 6: Sun and Rikus (1999), Sun (2001)
                               !< 7: effective radius
   contains
     generic,   public :: setup => setup_v1, setup_v3
     procedure, private :: setup_v1
     procedure, private :: setup_v3
     
  end type ty_gfdlmp_config
  !
  type(ty_gfdlmp_config) :: cfg

  public :: cfg
contains
  ! #######################################################################################
  !
  ! #######################################################################################
  function setup_v1(cfg, mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land,      &
                    dw_ocean,vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, &
                    vi_max, vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max,     &
                    qi0_max, qi0_crt, qr0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr,      &
                    const_vi, const_vs, const_vg, const_vr, use_ccn, rthresh, ccn_l,      &
                    ccn_o, qc_crt, tau_g2v, tau_v2g, sat_adj0, c_piacr, tau_imlt, tau_v2l,&
                    tau_l2v, tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs,  &
                    z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin, tice,        &
                    rad_snow, rad_graupel, rad_rain, cld_min, use_ppm, mono_prof,         &
                    do_sedi_heat, sedi_transport, do_sedi_w, de_ice, icloud_f, irain_f,   &
                    mp_print, reiflag, rewmin, rewmax, reimin, reimax, rermin, rermax,    &
                    resmin, resmax, regmin, regmax, tintqs, do_hail) result(err_message)
    !
    class(ty_gfdlmp_config), intent(inout)  :: cfg
    character(len=128) :: err_message
    logical, intent(in) :: do_qa, fix_negative, fast_sat_adj, const_vi, const_vs,         &
         const_vg, const_vr, use_ccn, z_slope_liq, z_slope_ice, prog_ccn, rad_snow,       &
         rad_graupel, rad_rain, use_ppm, mono_prof, mp_print, do_hail, tintqs,            &
         do_sedi_heat, do_sedi_w, sedi_transport, de_ice
    real, intent(in) :: mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r,                &
         tau_g2v, tau_v2g, tau_imlt, tau_v2l, tau_l2v, tau_i2s, tau_l2r, dw_land,         &
         dw_ocean, vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, vi_max, vs_max, vg_max,        &
         vr_max, ql0_max, qi0_max, qs0_crt, qi0_crt, qr0_crt,qc_crt, qs_mlt, qi_gen,      &
         rh_inc, rh_ins, rh_inr, rthresh, ccn_l, ccn_o, sat_adj0, c_piacr, qi_lim,        &
         ql_gen, c_paut, c_psaci, c_pgacs, c_cracw, alin, clin, tice, cld_min, rewmin,    &
         rewmax, reimin, reimax, rermin, rermax, resmin, resmax, regmin, regmax
    integer, intent(in) :: icloud_f, irain_f, reiflag

    ! initialize error message
    err_message = ""
  
    cfg%mp_time        = mp_time
    cfg%t_min          = t_min
    cfg%t_sub          = t_sub
    cfg%tau_r2g        = tau_r2g
    cfg%tau_smlt       = tau_smlt
    cfg%tau_g2r        = tau_g2r
    cfg%dw_land        = dw_land
    cfg%dw_ocean       = dw_ocean
    cfg%vi_fac         = vi_fac
    cfg%vr_fac         = vr_fac
    cfg%vs_fac         = vs_fac
    cfg%vg_fac         = vg_fac
    cfg%ql_mlt         = ql_mlt
    cfg%do_qa          = do_qa
    cfg%fix_negative   = fix_negative
    cfg%vi_max         = vi_max
    cfg%vs_max         = vs_max
    cfg%vg_max         = vg_max
    cfg%vr_max         = vr_max
    cfg%qs_mlt         = qs_mlt
    cfg%qs0_crt        = qs0_crt
    cfg%qi_gen         = qi_gen
    cfg%ql0_max        = ql0_max
    cfg%qi0_max        = qi0_max
    cfg%qi0_crt        = qi0_crt
    cfg%qr0_crt        = qr0_crt
    cfg%fast_sat_adj   = fast_sat_adj
    cfg%rh_inc         = rh_inc
    cfg%rh_ins         = rh_ins
    cfg%rh_inr         = rh_inr
    cfg%const_vi       = const_vi
    cfg%const_vs       = const_vs
    cfg%const_vg       = const_vg
    cfg%const_vr       = const_vr
    cfg%use_ccn        = use_ccn
    cfg%rthresh        = rthresh
    cfg%ccn_l          = ccn_l
    cfg%ccn_o          = ccn_o
    cfg%qc_crt         = qc_crt
    cfg%tau_g2v        = tau_g2v
    cfg%tau_v2g        = tau_v2g
    cfg%sat_adj0       = sat_adj0
    cfg%c_piacr        = c_piacr
    cfg%tau_imlt       = tau_imlt
    cfg%tau_v2l        = tau_v2l
    cfg%tau_l2v        = tau_l2v
    cfg%tau_i2s        = tau_i2s
    cfg%tau_l2r        = tau_l2r
    cfg%qi_lim         = qi_lim
    cfg%ql_gen         = ql_gen
    cfg%c_paut         = c_paut
    cfg%c_psaci        = c_psaci
    cfg%c_pgacs        = c_pgacs
    cfg%z_slope_liq    = z_slope_liq
    cfg%z_slope_ice    = z_slope_ice
    cfg%prog_ccn       = prog_ccn
    cfg%c_cracw        = c_cracw
    cfg%alin           = alin
    cfg%clin           = clin
    cfg%tice           = tice
    cfg%rad_snow       = rad_snow
    cfg%rad_graupel    = rad_graupel
    cfg%rad_rain       = rad_rain
    cfg%cld_min        = cld_min
    cfg%use_ppm        = use_ppm
    cfg%mono_prof      = mono_prof
    cfg%do_sedi_heat   = do_sedi_heat
    cfg%sedi_transport = sedi_transport
    cfg%do_sedi_w      = do_sedi_w
    cfg%de_ice         = de_ice
    cfg%icloud_f       = icloud_f
    cfg%irain_f        = irain_f
    cfg%mp_print       = mp_print
    cfg%reiflag        = reiflag
    cfg%rewmin         = rewmin
    cfg%rewmax         = rewmax
    cfg%reimin         = reimin
    cfg%reimax         = reimax
    cfg%rermin         = rermin
    cfg%rermax         = rermax
    cfg%resmin         = resmin
    cfg%resmax         = resmax
    cfg%regmin         = regmin
    cfg%regmax         = regmax
    cfg%tintqs         = tintqs
    cfg%do_hail        = do_hail

  end function setup_v1
  
  ! #######################################################################################
  !
  ! #######################################################################################
  function setup_v3(cfg) result(err_message)
    class(ty_gfdlmp_config), intent(inout)  :: cfg
    character(len=128) :: err_message

    ! initialize error message
    err_message = ""
    
  end function setup_v3

end module module_gfdlmp_param
