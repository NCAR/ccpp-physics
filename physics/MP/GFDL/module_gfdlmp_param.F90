! #########################################################################################
! #########################################################################################
module module_gfdlmp_param
  use machine, only: kind_phys
  implicit none

  ! #######################################################################################
  ! Data container for GFDL MP runtime configurations information (i.e. Namelist)
  ! #######################################################################################
  type ty_gfdlmp_v1_config
     real :: cld_min  = 0.05       !< minimum cloud fraction
     real :: tice     = 273.16     !< set tice = 165. to trun off ice - phase phys (kessler emulator)
     real :: t_min    = 178.       !< min temp to freeze - dry all water vapor
     real :: t_sub    = 184.       !< min temp for sublimation of cloud ice
     real :: mp_time  = 150.       !< maximum micro - physics time step (sec)
     real :: rh_inc = 0.25         !< rh increment for complete evaporation of cloud water and cloud ice
     real :: rh_inr = 0.25         !< rh increment for minimum evaporation of rain
     real :: rh_ins = 0.25         !< rh increment for sublimation of snow
     real :: tau_r2g  = 900.       !< rain freezing during fast_sat
     real :: tau_smlt = 900.       !< snow melting
     real :: tau_g2r  = 600.       !< graupel melting to rain
     real :: tau_imlt = 600.       !< cloud ice melting
     real :: tau_i2s  = 1000.      !< cloud ice to snow auto-conversion
     real :: tau_l2r  = 900.       !< cloud water to rain auto-conversion
     real :: tau_v2l  = 150.       !< water vapor to cloud water (condensation)
     real :: tau_l2v  = 300.       !< cloud water to water vapor (evaporation)
     real :: tau_g2v  = 900.       !< graupel sublimation
     real :: tau_v2g  = 21600.     !< graupel deposition -- make it a slow process
     real :: dw_land  = 0.20       !< value for subgrid deviation / variability over land
     real :: dw_ocean = 0.10       !< base value for ocean
     real :: ccn_o = 90.           !< ccn over ocean (cm^ - 3)
     real :: ccn_l = 270.          !< ccn over land (cm^ - 3)
     real :: rthresh  = 10.0e-6    !< critical cloud drop radius (micro m)
     real :: sat_adj0 = 0.90       !< adjustment factor (0: no, 1: full) during fast_sat_adj
     real :: qc_crt   = 5.0e-8     !< mini condensate mixing ratio to allow partial cloudiness
     real :: qi_lim   = 1.         !< cloud ice limiter to prevent large ice build up
     real :: ql_mlt   = 2.0e-3     !< max value of cloud water allowed from melted cloud ice
     real :: qs_mlt   = 1.0e-6     !< max cloud water due to snow melt
     real :: ql_gen   = 1.0e-3     !< max cloud water generation during remapping step if fast_sat_adj = .t.
     real :: qi_gen   = 1.82e-6    !< max cloud ice generation during remapping step (V1 ONLY. Computed internally in V3)
     real :: ql0_max = 2.0e-3      !< max cloud water value (auto converted to rain)
     real :: qi0_max = 1.0e-4      !< max cloud ice value (by other sources)
     real :: qi0_crt = 1.0e-4      !< cloud ice to snow autoconversion threshold (was 1.e-4);
                                   !< qi0_crt is highly dependent on horizontal resolution
     real :: qr0_crt = 1.0e-4      !< rain to snow or graupel/hail threshold
                                   !< lfo used * mixing ratio * = 1.e-4 (hail in lfo)
     real :: qs0_crt = 1.0e-3      !< snow to graupel density threshold (0.6e-3 in purdue lin scheme)
     real :: c_paut  = 0.55        !< autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
     real :: c_psaci = 0.02        !< accretion: cloud ice to snow (was 0.1 in zetac)
     real :: c_piacr = 5.0         !< accretion: rain to ice:
     real :: c_cracw = 0.9         !< rain accretion efficiency
     real :: c_pgacs = 2.0e-3      !< snow to graupel "accretion" eff. (was 0.1 in zetac)
     real :: alin = 842.0          !< "a" in lin1983
     real :: clin = 4.8            !< "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)
     real :: vi_fac = 1.           !< if const_vi: 1 / 3
     real :: vs_fac = 1.           !< if const_vs: 1.
     real :: vg_fac = 1.           !< if const_vg: 2.
     real :: vr_fac = 1.           !< if const_vr: 4.
     real :: vi_max = 0.5          !< max fall speed for ice
     real :: vs_max = 5.0          !< max fall speed for snow
     real :: vg_max = 8.0          !< max fall speed for graupel
     real :: vr_max = 12.          !< max fall speed for rain
     real :: rewmin = 5.0          !< minimum effective radii (liquid)
     real :: rewmax = 10.0         !< maximum effective radii (liquid)
     real :: reimin = 10.0         !< minimum effective radii (ice)
     real :: reimax = 150.0        !< maximum effective radii (ice)
     real :: rermin = 10.0         !< minimum effective radii (rain)
     real :: rermax = 10000.0      !< maximum effective radii (rain)
     real :: resmin = 150.0        !< minimum effective radii (snow)
     real :: resmax = 10000.0      !< maximum effective radii (snow)
     real :: regmin = 300.0        !< minimum effective radii (graupel)
     real :: regmax = 10000.0      !< maximum effective radii (graupel)
     !
     logical :: const_vi       = .false. !< if .t. the constants are specified by v * _fac
     logical :: const_vs       = .false. !< if .t. the constants are specified by v * _fac
     logical :: const_vg       = .false. !< if .t. the constants are specified by v * _fac
     logical :: const_vr       = .false. !< if .t. the constants are specified by v * _fac
     logical :: fast_sat_adj   = .false. !< has fast saturation adjustments
     logical :: z_slope_liq    = .true.  !< use linear mono slope for autocconversions
     logical :: z_slope_ice    = .false. !< use linear mono slope for autocconversions
     logical :: use_ccn        = .false. !< must be true when prog_ccn is false
     logical :: use_ppm        = .false. !< use ppm fall scheme
     logical :: mono_prof      = .true.  !< perform terminal fall with mono ppm scheme
     logical :: mp_print       = .false. !< cloud microphysics debugging printout
     logical :: do_hail        = .false. !< use hail parameters instead of graupel
     logical :: do_qa          = .true.  !< do inline cloud fraction
     logical :: rad_snow       = .true.  !< consider snow in cloud fraciton calculation
     logical :: rad_graupel    = .true.  !< consider graupel in cloud fraction calculation
     logical :: rad_rain       = .true.  !< consider rain in cloud fraction calculation
     logical :: de_ice         = .false. !< to prevent excessive build - up of cloud ice from external sources
     logical :: sedi_transport = .true.  !< transport of momentum in sedimentation
     logical :: do_sedi_w      = .false. !< transport of vertical motion in sedimentation
     logical :: do_sedi_heat   = .true.  !< transport of heat in sedimentation
     logical :: prog_ccn       = .false. !< do prognostic ccn (yi ming's method)
     logical :: fix_negative   = .false. !< fix negative water species
     logical :: tintqs         = .false. !<
     !
     integer :: icloud_f = 0          !< GFDL cloud scheme
                                      !< 0: subgrid variability based scheme
                                      !< 1: same as 0, but for old fvgfs implementation
                                      !< 2: binary cloud scheme
                                      !< 3: extension of 0
     integer :: irain_f = 0           !< cloud water to rain auto conversion scheme
                                      !< 0: subgrid variability based scheme
                                      !< 1: no subgrid varaibility
     integer :: reiflag = 1           !< cloud ice effective radius scheme
                                      !< 1: Heymsfield and Mcfarquhar (1996)
                                      !< 2: Donner et al. (1997)
                                      !< 3: Fu (2007)
                                      !< 4: Kristjansson et al. (2000)
                                      !< 5: Wyser (1998)
                                      !< 6: Sun and Rikus (1999), Sun (2001)
                                      !< 7: effective radius
   contains
     generic,   public  :: update => update_cfg_v1
     procedure, private :: update_cfg_v1
  end type ty_gfdlmp_v1_config

  ! #######################################################################################
  type ty_gfdlmp_v3_config
     !
     real :: cld_min  = 0.05       !< minimum cloud fraction
     real :: tice     = 273.15     !< freezing temperature (K): ref: GFDL, GFS
     real :: t_min    = 178.       !< min temp to freeze - dry all water vapor
     real :: t_sub    = 184.       !< min temp for sublimation of cloud ice
     real :: mp_time  = 150.       !< maximum micro - physics time step (sec)
     real :: rh_inc = 0.25         !< rh increment for complete evaporation of cloud water and cloud ice
     real :: rh_inr = 0.25         !< rh increment for minimum evaporation of rain
     real :: rh_ins = 0.25         !< rh increment for sublimation of snow
     real :: tau_r2g = 900.0 ! rain freezing to graupel time scale (s)
     real :: tau_smlt = 900.0 ! snow melting time scale (s)
     real :: tau_imlt = 1200.0 ! cloud ice melting time scale (s) 
     real :: tau_i2s = 1000.0 ! cloud ice to snow autoconversion time scale (s)
     real :: tau_l2r = 900.0 ! cloud water to rain autoconversion time scale (s)
     real :: tau_v2l = 150.0 ! water vapor to cloud water condensation time scale (s)
     real :: tau_l2v = 300.0 ! cloud water to water vapor evaporation time scale (s)
     real :: dw_land = 0.20 ! base value for subgrid deviation / variability over land
     real :: dw_ocean = 0.10 ! base value for subgrid deviation / variability over ocean
     real :: ccn_o = 90.0 ! ccn over ocean (1/cm^3)
     real :: ccn_l = 270.0 ! ccn over land (1/cm^3)
     real :: rthresh = 20.0e-6 ! critical cloud drop radius (micron) for autoconversion
     real :: sat_adj0 = 0.90  ! adjustment factor (0: no, 1: full) during fast_sat_adj
     real :: qi_lim = 1.0 ! cloud ice limiter (0: no, 1: full, >1: extra) to prevent large ice build up
     real :: ql_mlt = 2.0e-3 ! maximum cloud water allowed from melted cloud ice (kg/kg)
     real :: qs_mlt = 1.0e-6 ! maximum cloud water allowed from melted snow (kg/kg)
     real :: ql_gen = 1.0e-3 ! maximum cloud water generation during remapping step (kg/kg)
     real :: qi_gen =  1.82e-6 ! max cloud ice generation during remapping step
     real :: ql0_max = 2.0e-3 ! maximum cloud water value (autoconverted to rain) (kg/kg)
     real :: qi0_max = 1.0e-4 ! maximum cloud ice value (autoconverted to snow) (kg/m^3)
     real :: qi0_crt = 1.0e-4 ! cloud ice to snow autoconversion threshold (kg/m^3)
     real :: qs0_crt = 1.0e-3 ! snow to graupel autoconversion threshold (0.6e-3 in Purdue Lin scheme) (kg/m^3)
     real :: c_paut = 0.55 ! cloud water to rain autoconversion efficiency
     real :: c_psacw = 1.0 ! cloud water to snow accretion efficiency
     real :: c_psaci = 0.05 ! cloud ice to snow accretion efficiency (was 0.1 in ZETAC)
     real :: c_pracw = 0.8 ! cloud water to rain accretion efficiency
     real :: c_praci = 1.0 ! cloud ice to rain accretion efficiency
     real :: c_pgacw = 1.0 ! cloud water to graupel accretion efficiency
     real :: c_pgaci = 0.05 ! cloud ice to graupel accretion efficiency (was 0.1 in ZETAC)
     real :: c_pracs = 1.0 ! snow to rain accretion efficiency
     real :: c_psacr = 1.0 ! rain to snow accretion efficiency
     real :: c_pgacr = 1.0 ! rain to graupel accretion efficiency
     real :: c_pgacs = 0.01 ! snow to graupel accretion efficiency (was 0.1 in ZETAC)

     real :: alinw = 3.e7 ! "a" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
     real :: alini = 7.e2 ! "a" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
     real :: alinr = 842.0 ! "a" in Lin et al. (1983) for rain (Liu and Orville 1969)
     real :: alins = 4.8 ! "a" in Lin et al. (1983) for snow (straka 2009)
     real :: aling = 1.0 ! "a" in Lin et al. (1983), similar to a, but for graupel (Pruppacher and Klett 2010)
     real :: alinh = 1.0 ! "a" in Lin et al. (1983), similar to a, but for hail (Pruppacher and Klett 2010)

     real :: blinw = 2.0 ! "b" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
     real :: blini = 1.0 ! "b" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
     real :: blinr = 0.8 ! "b" in Lin et al. (1983) for rain (Liu and Orville 1969)
     real :: blins = 0.25 ! "b" in Lin et al. (1983) for snow (straka 2009)
     real :: bling = 0.5 ! "b" in Lin et al. (1983), similar to b, but for graupel (Pruppacher and Klett 2010)
     real :: blinh = 0.5 ! "b" in Lin et al. (1983), similar to b, but for hail (Pruppacher and Klett 2010)
     real :: vw_fac = 1.0
     real :: vi_fac = 1.0 ! IFS: if const_vi: 1 / 3
     real :: vs_fac = 1.0 ! IFS: if const_vs: 1.
     real :: vg_fac = 1.0 ! IFS: if const_vg: 2.
     real :: vr_fac = 1.0 ! IFS: if const_vr: 4.
     real :: vw_max = 0.01 ! maximum fall speed for cloud water (m/s)
     real :: vi_max = 0.5 ! maximum fall speed for cloud ice (m/s)
     real :: vs_max = 5.0 ! maximum fall speed for snow (m/s)
     real :: vg_max = 8.0 ! maximum fall speed for graupel (m/s)
     real :: vr_max = 12.0 ! maximum fall speed for rain (m/s)
     real :: rewmin = 5.0, rewmax = 15.0 ! minimum and maximum effective radius for cloud water (micron)
     real :: reimin = 10.0, reimax = 150.0 ! minimum and maximum effective radius for cloud ice (micron)
     real :: rermin = 15.0, rermax = 10000.0 ! minimum and maximum effective radius for rain (micron)
     real :: resmin = 150.0, resmax = 10000.0 ! minimum and maximum effective radius for snow (micron)
     real :: regmin = 150.0, regmax = 10000.0 ! minimum and maximum effective radius for graupel
     real :: tice_mlt = 273.16 ! can set ice melting temperature to 268 based on observation (Kay et al. 2016) (K)
     real :: tau_gmlt = 600.0 ! graupel melting time scale (s)
     real :: tau_wbf = 300.0 ! graupel melting time scale (s)
     real :: tau_revp = 0.0 ! rain evaporation time scale (s)
     real :: is_fac = 0.2 ! cloud ice sublimation temperature factor
     real :: ss_fac = 0.2 ! snow sublimation temperature factor
     real :: gs_fac = 0.2 ! graupel sublimation temperature factor
     real :: rh_fac_evap = 10.0 ! cloud water evaporation relative humidity factor
     real :: rh_fac_cond = 10.0 ! cloud water condensation relative humidity factor
     real :: sed_fac = 1.0 ! coefficient for sedimentation fall, scale from 1.0 (implicit) to 0.0 (lagrangian)
     real :: xr_a = 0.25 ! p value in Xu and Randall (1996)
     real :: xr_b = 100.0 ! alpha_0 value in Xu and Randall (1996)
     real :: xr_c = 0.49 ! gamma value in Xu and Randall (1996)
     real :: te_err = 1.e-5 ! 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time
     real :: tw_err = 1.e-8 ! 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time
     real :: rh_thres = 0.75 ! minimum relative humidity for cloud fraction
     real :: rhc_cevap = 0.85 ! maximum relative humidity for cloud water evaporation
     real :: rhc_revap = 0.85 ! maximum relative humidity for rain evaporation
     real :: f_dq_p = 1.0 ! cloud fraction adjustment for supersaturation
     real :: f_dq_m = 1.0 ! cloud fraction adjustment for undersaturation
     real :: fi2s_fac = 1.0 ! maximum sink of cloud ice to form snow: 0-1
     real :: fi2g_fac = 1.0 ! maximum sink of cloud ice to form graupel: 0-1
     real :: fs2g_fac = 1.0 ! maximum sink of snow to form graupel: 0-1

     logical :: const_vw = .false. ! if .ture., the constants are specified by v * _fac
     logical :: const_vi = .false. ! if .ture., the constants are specified by v * _fac
     logical :: const_vs = .false. ! if .ture., the constants are specified by v * _fac
     logical :: const_vg = .false. ! if .ture., the constants are specified by v * _fac
     logical :: const_vr = .false. ! if .ture., the constants are specified by v * _fac
     logical :: z_slope_liq = .true. ! use linear mono slope for autocconversions
     logical :: z_slope_ice = .true. ! use linear mono slope for autocconversions
     logical :: do_hail = .false. ! use hail parameters instead of graupel
     logical :: do_qa = .true. ! do inline cloud fraction
     logical :: rad_snow = .true. ! include snow in cloud fraciton calculation
     logical :: rad_graupel = .true. ! include graupel in cloud fraction calculation
     logical :: rad_rain = .true. ! include rain in cloud fraction calculation
     logical :: do_sedi_uv = .true. ! transport of horizontal momentum in sedimentation
     logical :: do_sedi_w = .true. ! transport of vertical momentum in sedimentation
     logical :: do_sedi_heat = .true. ! transport of heat in sedimentation
     logical :: do_sedi_melt = .true. ! melt cloud ice, snow, and graupel during sedimentation
     logical :: prog_ccn = .false. ! do prognostic ccn (Yi Ming's method)
     logical :: fix_negative = .true. ! fix negative water species
     logical :: tintqs         = .false. !<
     logical :: liq_ice_combine = .false. ! combine all liquid water, combine all solid water
     logical :: snow_grauple_combine = .true. ! combine snow and graupel
     logical :: use_rhc_cevap = .false. ! cap of rh for cloud water evaporation
     logical :: use_rhc_revap = .false. ! cap of rh for rain evaporation
     logical :: do_cld_adj = .false. ! do cloud fraction adjustment
    logical :: do_evap_timescale = .true. ! whether to apply a timescale to evaporation
    logical :: do_cond_timescale = .false. ! whether to apply a timescale to condensation    
    logical :: consv_checker = .false. ! turn on energy and water conservation checker
    logical :: do_warm_rain_mp = .false. ! do warm rain cloud microphysics only
    logical :: do_wbf = .false. ! do Wegener Bergeron Findeisen process
    logical :: do_psd_water_fall = .false. ! calculate cloud water terminal velocity based on PSD
    logical :: do_psd_ice_fall = .false. ! calculate cloud ice terminal velocity based on PSD
    logical :: do_psd_water_num = .false. ! calculate cloud water number concentration based on PSD
    logical :: do_psd_ice_num = .false. ! calculate cloud ice number concentration based on PSD
    logical :: do_new_acc_water = .false. ! perform the new accretion for cloud water
    logical :: do_new_acc_ice = .false. ! perform the new accretion for cloud ice
    logical :: cp_heating = .false. ! update temperature based on constant pressure
    logical :: delay_cond_evap = .false. ! do condensation evaporation only at the last time step
    logical :: do_subgrid_proc = .true. ! do temperature sentive high vertical resolution processes
    logical :: fast_fr_mlt = .true. ! do freezing and melting in fast microphysics
    logical :: fast_dep_sub = .true. ! do deposition and sublimation in fast microphysics
    integer :: ntimes = 1 ! cloud microphysics sub cycles
    integer :: nconds = 1 ! condensation sub cycles
    integer :: icloud_f = 0 ! GFDL cloud scheme
                            ! 0: subgrid variability based scheme
                            ! 1: same as 0, but for old fvgfs implementation
                            ! 2: binary cloud scheme
                            ! 3: extension of 0
    integer :: irain_f = 0 ! cloud water to rain auto conversion scheme
                           ! 0: subgrid variability based scheme
                           ! 1: no subgrid varaibility
    integer :: reiflag = 5 ! cloud ice effective radius scheme
                           ! 1: Heymsfield and Mcfarquhar (1996)
                           ! 2: Donner et al. (1997)
                           ! 3: Fu (2007)
                           ! 4: Kristjansson et al. (2000)
                           ! 5: Wyser (1998)
                           ! 6: Sun and Rikus (1999), Sun (2001)
                           ! 7: effective radius
     !
     integer :: inflag = 1             !< ice nucleation scheme
                                      !< 1: Hong et al. (2004)
                                      !< 2: Meyers et al. (1992)
                                      !< 3: Meyers et al. (1992)
                                      !< 4: Cooper (1986)
                                      !, 5: Fletcher (1962)
     !
     integer :: igflag = 3             !< ice generation scheme
                                      !< 1: WSM6
                                      !< 2: WSM6 with 0 at 0 C
                                      !< 3: WSM6 with 0 at 0 C and fixed value at - 10 C
                                      !< 4: combination of 1 and 3
     !
     integer :: ifflag = 1            !< ice fall scheme
                                      !< 1: Deng and Mace (2008)
                                      !< 2: Heymsfield and Donner (1990)
     !
     integer :: rewflag = 1           !< cloud water effective radius scheme
                                      !< 1: Martin et al. (1994)
                                      !< 2: Martin et al. (1994), GFDL revision
                                      !< 3: Kiehl et al. (1994)
                                      !< 4: effective radius
     !
     integer :: rerflag = 1           !< rain effective radius scheme
                                      !< 1: effective radius
     !
     integer :: resflag = 1           !< snow effective radius scheme
                                      !< 1: effective radius
     !
     integer :: regflag = 1           !< graupel effective radius scheme
                                      !< 1: effective radius
     !
     integer :: radr_flag = 1         !< radar reflectivity for rain
                                      !< 1: Mark Stoelinga (2005)
                                      !< 2: Smith et al. (1975), Tong and Xue (2005)
                                      !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
     !
     integer :: rads_flag = 1         !< radar reflectivity for snow
                                      !< 1: Mark Stoelinga (2005)
                                      !< 2: Smith et al. (1975), Tong and Xue (2005)
                                      !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
     !
     integer :: radg_flag = 1         !< radar reflectivity for graupel
                                      !< 1: Mark Stoelinga (2005)
                                      !< 2: Smith et al. (1975), Tong and Xue (2005)
                                      !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
     !
     integer :: sedflag   = 1         !< sedimentation scheme
                                      !< 1: implicit scheme
                                      !< 2: explicit scheme
                                      !< 3: lagrangian scheme
                                      !< 4: combined implicit and lagrangian scheme
     !
     integer :: vdiffflag = 1         !< wind difference scheme in accretion
                                      !< 1: Wisner et al. (1972)
                                      !< 2: Mizuno (1990)
                                      !< 3: Murakami (1990)

     real :: n0w_sig = 1.1 ! intercept parameter (significand) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
     real :: n0i_sig = 1.3 ! intercept parameter (significand) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
     real :: n0r_sig = 8.0 ! intercept parameter (significand) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
     real :: n0s_sig = 3.0 ! intercept parameter (significand) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
     real :: n0g_sig = 4.0 ! intercept parameter (significand) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
     real :: n0h_sig = 4.0 ! intercept parameter (significand) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)
     
     real :: n0w_exp = 41 ! intercept parameter (exponent) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
     real :: n0i_exp = 18 ! intercept parameter (exponent) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
     real :: n0r_exp = 6 ! intercept parameter (exponent) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
     real :: n0s_exp = 6 ! intercept parameter (exponent) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
     real :: n0g_exp = 6 ! intercept parameter (exponent) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
     real :: n0h_exp = 4 ! intercept parameter (exponent) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)
     
     real :: muw = 6.0 ! shape parameter of cloud water in Gamma distribution (Martin et al. 1994)
     real :: mui = 3.35 ! shape parameter of cloud ice in Gamma distribution (McFarquhar et al. 2015)
     real :: mur = 1.0 ! shape parameter of rain in Gamma distribution (Marshall and Palmer 1948)
     real :: mus = 1.0 ! shape parameter of snow in Gamma distribution (Gunn and Marshall 1958)
     real :: mug = 1.0 ! shape parameter of graupel in Gamma distribution (Houze et al. 1979)
     real :: muh = 1.0 ! shape parameter of hail in Gamma distribution (Federer and Waldvogel 1975)
     real :: beta = 1.22 ! defined in Heymsfield and Mcfarquhar (1996)
     
     real :: rewfac = 1.0  !< this is a tuning parameter to compromise the inconsistency between
                           !< GFDL MP's PSD and cloud water radiative property's PSD assumption.
                           !< after the cloud water radiative property's PSD is rebuilt,
                           !< this parameter should be 1.0.
     real :: reifac = 1.0  !< this is a tuning parameter to compromise the inconsistency between
                           !< GFDL MP's PSD and cloud ice radiative property's PSD assumption.
                           !< after the cloud ice radiative property's PSD is rebuilt,
                           !< this parameter should be 1.0.
   contains
     generic,   public  :: update => update_cfg_v3
     procedure, private :: update_cfg_v3
  end type ty_gfdlmp_v3_config

  type(ty_gfdlmp_v3_config) :: cfg_v3
  type(ty_gfdlmp_v1_config) :: cfg_v1
  
contains
  ! #######################################################################################
  !
  ! #######################################################################################
  function update_cfg_v1(cfg, mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, &
       dw_ocean,vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max,      &
       vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max, qi0_crt,        &
       qr0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr, const_vi, const_vs, const_vg,       &
       const_vr, use_ccn, rthresh, ccn_l, ccn_o, qc_crt, tau_g2v, tau_v2g, sat_adj0,      &
       c_piacr, tau_imlt, tau_v2l, tau_l2v, tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut,     &
       c_psaci, c_pgacs, z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin, tice,   &
       rad_snow, rad_graupel, rad_rain, cld_min, use_ppm, mono_prof, do_sedi_heat,        &
       sedi_transport, do_sedi_w, de_ice, icloud_f, irain_f, mp_print, reiflag, rewmin,   &
       rewmax, reimin, reimax, rermin, rermax, resmin, resmax, regmin, regmax, tintqs,    &
       do_hail) result(err_message)
    
    class(ty_gfdlmp_v1_config), intent(inout)  :: cfg
    character(len=128) :: err_message
    logical, intent(in) :: do_qa, fix_negative, fast_sat_adj, const_vi, const_vs,         &
         const_vg, const_vr, use_ccn, z_slope_liq, z_slope_ice, prog_ccn, rad_snow,       &
         rad_graupel, rad_rain, use_ppm, mono_prof, mp_print, do_hail, tintqs,            &
         do_sedi_heat, do_sedi_w, sedi_transport, de_ice
    real, intent(in) :: mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r,                &
         tau_g2v, tau_v2g, tau_imlt, tau_v2l, tau_l2v, tau_i2s, tau_l2r, dw_land,         &
         dw_ocean,  vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, vi_max, vs_max, vg_max,       &
         vr_max, ql0_max, qi0_max, qs0_crt, qi0_crt, qr0_crt,qc_crt, qs_mlt, qi_gen,      &
         rh_inc, rh_ins, rh_inr, rthresh, ccn_l, ccn_o, sat_adj0, c_piacr, qi_lim, ql_gen,&
         c_paut, c_psaci, c_pgacs, c_cracw, alin, clin, tice, cld_min, rewmin, rewmax,    &
         reimin, reimax, rermin, rermax, resmin, resmax, regmin, regmax
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

  end function update_cfg_v1
  
  ! #######################################################################################
  !
  ! #######################################################################################
  function update_cfg_v3(cfg, t_min, t_sub, tau_r2g, tau_smlt, tau_gmlt, dw_land, dw_ocean,&
       vw_fac, vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vw_max,       &
       vi_max, vs_max, vg_max, vr_max, qs_mlt, qs0_crt, ql0_max, qi0_max, qi0_crt, ifflag,&
       rh_inc, rh_ins, rh_inr, const_vw, const_vi, const_vs, const_vg, const_vr, rthresh, &
       ccn_l, ccn_o, igflag, c_paut, tau_imlt, tau_v2l, tau_l2v, tau_i2s, tau_l2r, qi_lim,&
       ql_gen, do_hail, inflag, c_psacw, c_psaci, c_pracs, c_psacr, c_pgacr, c_pgacs,     &
       c_pgacw, c_pgaci, z_slope_liq, z_slope_ice, prog_ccn, c_pracw, c_praci, rad_snow,  &
       rad_graupel, rad_rain, cld_min, sedflag, sed_fac, do_sedi_uv, do_sedi_w,           &
       do_sedi_heat, icloud_f, irain_f, xr_a, xr_b, xr_c, ntimes, tau_revp, tice_mlt,     &
       do_cond_timescale, mp_time, consv_checker, te_err, tw_err, use_rhc_cevap,          &
       use_rhc_revap, tau_wbf, do_warm_rain_mp, rh_thres, f_dq_p, f_dq_m, do_cld_adj,     &
       rhc_cevap, rhc_revap, beta, liq_ice_combine, rewflag, reiflag, rerflag, resflag,   &
       regflag, rewmin, rewmax, reimin, reimax, rermin, rermax, resmin, resmax, regmin,   &
       regmax, fs2g_fac, fi2s_fac, fi2g_fac, do_sedi_melt, radr_flag, rads_flag,          &
       radg_flag, do_wbf, do_psd_water_fall, do_psd_ice_fall, n0w_sig, n0i_sig, n0r_sig,  &
       n0s_sig, n0g_sig, n0h_sig, n0w_exp, n0i_exp, n0r_exp, n0s_exp, n0g_exp, n0h_exp,   &
       muw, mui, mur, mus, mug, muh, alinw, alini, alinr, alins, aling, alinh, blinw,     &
       blini, blinr, blins, bling, blinh, do_new_acc_water, do_new_acc_ice, is_fac,       &
       ss_fac, gs_fac, rh_fac_evap, rh_fac_cond, snow_grauple_combine, do_psd_water_num,  &
       do_psd_ice_num, vdiffflag, rewfac, reifac, cp_heating, nconds, do_evap_timescale,  &
       delay_cond_evap, do_subgrid_proc, fast_fr_mlt, fast_dep_sub, qi_gen, sat_adj0)     &
       result(err_message)
    
    class(ty_gfdlmp_v3_config), intent(inout)  :: cfg
    character(len=128) :: err_message
    logical, intent(in) :: do_qa, fix_negative, const_vw, const_vi, const_vs, const_vg,   &
         const_vr,  z_slope_liq, z_slope_ice, prog_ccn, rad_snow, rad_graupel, rad_rain,  &
         do_hail, do_sedi_heat, do_sedi_w, do_sedi_uv, do_sedi_melt,                      &
         do_cond_timescale, consv_checker, use_rhc_cevap, use_rhc_revap, do_warm_rain_mp, &
         do_cld_adj, liq_ice_combine, do_wbf, do_psd_water_fall, do_psd_ice_fall,         &
         do_psd_water_num, do_psd_ice_num, do_new_acc_water, do_new_acc_ice,              &
         snow_grauple_combine, cp_heating, do_evap_timescale, delay_cond_evap,            &
         do_subgrid_proc, fast_fr_mlt, fast_dep_sub
    real, intent(in) :: mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_gmlt, tau_imlt,     &
         tau_v2l, tau_l2v, tau_i2s, tau_l2r, dw_land, dw_ocean, vw_fac,  vi_fac, vr_fac,  &
         vs_fac, vg_fac, ql_mlt, vw_max, vi_max, vs_max, vg_max, vr_max, ql0_max, qi0_max,&
         qs0_crt, qi0_crt, qs_mlt, qi_gen, rh_inc, rh_ins, rh_inr, rthresh, ccn_l, ccn_o, &
         sat_adj0, qi_lim, ql_gen, sed_fac, c_paut, c_psaci, c_psacw, c_psacr, c_pracs,   &
         c_pgacs, c_pgacr, c_pgacw, c_pgaci, c_pracw, c_praci, cld_min,                   &
         rewmin, rewmax, reimin, reimax, rermin, rermax, resmin, resmax, regmin, regmax,  &
         xr_a, xr_b, xr_c, tau_revp, tice_mlt, te_err, tw_err, tau_wbf, rh_thres, f_dq_p, &
         f_dq_m, rhc_cevap, rhc_revap, beta, fs2g_fac, fi2s_fac, fi2g_fac, n0w_sig,       &
         n0i_sig, n0r_sig, n0s_sig, n0g_sig, n0h_sig, n0w_exp, n0i_exp, n0r_exp, n0s_exp, &
         n0g_exp, n0h_exp, muw, mui, mur, mus, mug, muh, alinw, alini, alinr, alins,      &
         aling, alinh, blinw, blini, blinr, blins, bling, blinh, is_fac, ss_fac, gs_fac,  &
         rh_fac_evap, rh_fac_cond, rewfac, reifac
    integer, intent(in) :: icloud_f, irain_f, reiflag, ifflag, igflag, inflag, sedflag,   &
         ntimes, rewflag, rerflag, resflag, regflag, radr_flag, rads_flag, radg_flag,     &
         vdiffflag, nconds
    
    ! initialize error message
    err_message = ""

    cfg%t_min                = t_min
    cfg%t_sub                = t_sub
    cfg%tau_r2g              = tau_r2g
    cfg%tau_smlt             = tau_smlt
    cfg%tau_gmlt             = tau_gmlt
    cfg%dw_land              = dw_land
    cfg%dw_ocean             = dw_ocean
    cfg%vw_fac               = vw_fac
    cfg%vi_fac               = vi_fac
    cfg%vr_fac               = vr_fac
    cfg%vs_fac               = vs_fac
    cfg%vg_fac               = vg_fac
    cfg%ql_mlt               = ql_mlt
    cfg%do_qa                = do_qa
    cfg%fix_negative         = fix_negative
    cfg%vw_max               = vw_max
    cfg%vi_max               = vi_max
    cfg%vs_max               = vs_max
    cfg%vg_max               = vg_max
    cfg%vr_max               = vr_max
    cfg%qs_mlt               = qs_mlt
    cfg%qs0_crt              = qs0_crt
    cfg%ql0_max              = ql0_max
    cfg%qi0_max              = qi0_max
    cfg%qi0_crt              = qi0_crt
    cfg%ifflag               = ifflag
    cfg%rh_inc               = rh_inc
    cfg%rh_ins               = rh_ins
    cfg%rh_inr               = rh_inr
    cfg%const_vw             = const_vw
    cfg%const_vi             = const_vi
    cfg%const_vs             = const_vs
    cfg%const_vg             = const_vg
    cfg%const_vr             = const_vr
    cfg%rthresh              = rthresh
    cfg%ccn_l                = ccn_l
    cfg%ccn_o                = ccn_o
    cfg%igflag               = igflag
    cfg%c_paut               = c_paut
    cfg%tau_imlt             = tau_imlt
    cfg%tau_v2l              = tau_v2l
    cfg%tau_l2v              = tau_l2v
    cfg%tau_i2s              = tau_i2s
    cfg%tau_l2r              = tau_l2r
    cfg%qi_lim               = qi_lim
    cfg%ql_gen               = ql_gen
    cfg%do_hail              = do_hail
    cfg%inflag               = inflag
    cfg%c_psacw              = c_psacw
    cfg%c_psaci              = c_psaci
    cfg%c_pracs              = c_pracs
    cfg%c_psacr              = c_psacr
    cfg%c_pgacr              = c_pgacr
    cfg%c_pgacs              = c_pgacs
    cfg%c_pgacw              = c_pgacw
    cfg%c_pgaci              = c_pgaci
    cfg%z_slope_liq          = z_slope_liq
    cfg%z_slope_ice          = z_slope_ice
    cfg%prog_ccn             = prog_ccn
    cfg%c_pracw              = c_pracw
    cfg%c_praci              = c_praci
    cfg%rad_snow             = rad_snow
    cfg%rad_graupel          = rad_graupel
    cfg%rad_rain             = rad_rain
    cfg%cld_min              = cld_min
    cfg%sedflag              = sedflag
    cfg%sed_fac              = sed_fac
    cfg%do_sedi_uv           = do_sedi_uv
    cfg%do_sedi_w            = do_sedi_w
    cfg%do_sedi_heat         = do_sedi_heat
    cfg%icloud_f             = icloud_f
    cfg%irain_f              = irain_f
    cfg%xr_a                 = xr_a
    cfg%xr_b                 = xr_b
    cfg%xr_c                 = xr_c
    cfg%ntimes               = ntimes
    cfg%tau_revp             = tau_revp
    cfg%tice_mlt             = tice_mlt
    cfg%do_cond_timescale    = do_cond_timescale
    cfg%mp_time              = mp_time
    cfg%consv_checker        = consv_checker
    cfg%te_err               = te_err
    cfg%tw_err               = tw_err
    cfg%use_rhc_cevap        = use_rhc_cevap
    cfg%use_rhc_revap        = use_rhc_revap
    cfg%tau_wbf              = tau_wbf
    cfg%do_warm_rain_mp      = do_warm_rain_mp
    cfg%rh_thres             = rh_thres
    cfg%f_dq_p               = f_dq_p
    cfg%f_dq_m               = f_dq_m
    cfg%do_cld_adj           = do_cld_adj
    cfg%rhc_cevap            = rhc_cevap
    cfg%rhc_revap            = rhc_revap
    cfg%beta                 = beta
    cfg%liq_ice_combine      = liq_ice_combine
    cfg%rewflag              = rewflag
    cfg%reiflag              = reiflag
    cfg%rerflag              = rerflag
    cfg%resflag              = resflag
    cfg%regflag              = regflag
    cfg%rewmin               = rewmin
    cfg%rewmax               = rewmax
    cfg%reimin               = reimin
    cfg%reimax               = reimax
    cfg%rermin               = rermin
    cfg%rermax               = rermax
    cfg%resmin               = resmin
    cfg%resmax               = resmax
    cfg%regmin               = regmin
    cfg%regmax               = regmax
    cfg%fs2g_fac             = fs2g_fac
    cfg%fi2s_fac             = fi2s_fac
    cfg%fi2g_fac             = fi2g_fac
    cfg%do_sedi_melt         = do_sedi_melt
    cfg%radr_flag            = radr_flag
    cfg%rads_flag            = rads_flag
    cfg%radg_flag            = radg_flag
    cfg%do_wbf               = do_wbf
    cfg%do_psd_water_fall    = do_psd_water_fall
    cfg%do_psd_ice_fall      = do_psd_ice_fall
    cfg%n0w_sig              = n0w_sig
    cfg%n0i_sig              = n0i_sig
    cfg%n0r_sig              = n0r_sig
    cfg%n0s_sig              = n0s_sig
    cfg%n0g_sig              = n0g_sig
    cfg%n0h_sig              = n0h_sig
    cfg%n0w_exp              = n0w_exp
    cfg%n0i_exp              = n0i_exp
    cfg%n0r_exp              = n0r_exp
    cfg%n0s_exp              = n0s_exp
    cfg%n0g_exp              = n0g_exp
    cfg%n0h_exp              = n0h_exp
    cfg%muw                  = muw
    cfg%mui                  = mui
    cfg%mur                  = mur
    cfg%mus                  = mus
    cfg%mug                  = mug
    cfg%muh                  = muh
    cfg%alinw                = alinw
    cfg%alini                = alini
    cfg%alinr                = alinr
    cfg%alins                = alins
    cfg%aling                = aling
    cfg%alinh                = alinh
    cfg%blinw                = blinw
    cfg%blini                = blini
    cfg%blinr                = blinr
    cfg%blins                = blins
    cfg%bling                = bling
    cfg%blinh                = blinh
    cfg%do_new_acc_water     = do_new_acc_water
    cfg%do_new_acc_ice       = do_new_acc_ice
    cfg%is_fac               = is_fac
    cfg%ss_fac               = ss_fac
    cfg%gs_fac               = gs_fac
    cfg%rh_fac_evap          = rh_fac_evap
    cfg%rh_fac_cond          = rh_fac_cond
    cfg%snow_grauple_combine = snow_grauple_combine
    cfg%do_psd_water_num     = do_psd_water_num
    cfg%do_psd_ice_num       = do_psd_ice_num
    cfg%vdiffflag            = vdiffflag
    cfg%rewfac               = rewfac
    cfg%reifac               = reifac
    cfg%cp_heating           = cp_heating
    cfg%nconds               = nconds
    cfg%do_evap_timescale    = do_evap_timescale
    cfg%delay_cond_evap      = delay_cond_evap
    cfg%do_subgrid_proc      = do_subgrid_proc
    cfg%fast_fr_mlt          = fast_fr_mlt
    cfg%fast_dep_sub         = fast_dep_sub
    cfg%qi_gen               = qi_gen
    cfg%sat_adj0             = sat_adj0
    
  end function update_cfg_v3

end module module_gfdlmp_param
