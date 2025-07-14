! #########################################################################################
! #########################################################################################
module module_gfdlmp_param
  use machine, only: kind_phys
  implicit none
  public :: read_gfdlmp_nml
  private

  ! #####################################################################################
  ! GFDL MP Version 1 parameters.
  ! #####################################################################################
  real(kind_phys) :: tau_g2r        = 600.    !< graupel melting to rain time scale (s)
  real(kind_phys) :: tau_g2v        = 900.    !< graupel sublimation time scale (s)
  real(kind_phys) :: tau_v2g        = 21600.  !< graupel deposition -- make it a slow process time scale (s)
  real(kind_phys) :: qc_crt         = 5.0e-8  !< mini condensate mixing ratio to allow partial cloudiness
  real(kind_phys) :: qr0_crt        = 1.0e-4  !< rain to snow or graupel/hail threshold
                                              !< lfo used * mixing ratio * = 1.e-4 (hail in lfo)
  real(kind_phys) :: c_piacr        = 5.0     !< accretion: rain to ice:
  real(kind_phys) :: c_cracw        = 0.9     !< rain accretion efficiency
  real(kind_phys) :: alin           = 842.0   !< "a" in lin1983
  real(kind_phys) :: clin           = 4.8     !< "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)
  logical         :: fast_sat_adj   = .false. !< has fast saturation adjustments 
  logical         :: use_ccn        = .false. !< must be true when prog_ccn is false
  logical         :: use_ppm        = .false. !< use ppm fall scheme
  logical         :: mono_prof      = .true.  !< perform terminal fall with mono ppm scheme
  logical         :: mp_print       = .false. !< cloud microphysics debugging printout
  logical         :: de_ice         = .false. !< to prevent excessive build - up of cloud ice from external sources
  logical         :: sedi_transport = .true.  !< transport of momentum in sedimentation

  ! #####################################################################################
  ! GFDL MP common (v1/v3) parameters
  ! #####################################################################################
  real(kind_phys) :: cld_min        = 0.05    !< (v1/v3) minimum cloud fraction
  real(kind_phys) :: t_min          = 178.    !< (v1/v3) min temp to freeze - dry all water vapor
  real(kind_phys) :: t_sub          = 184.    !< (v1/v3) min temp for sublimation of cloud ice
  real(kind_phys) :: mp_time        = 150.    !< (v1/v3) maximum micro - physics time step (sec)
  real(kind_phys) :: rh_inc         = 0.25    !< (v1/v3) rh increment for complete evaporation of cloud water and cloud ice
  real(kind_phys) :: rh_inr         = 0.25    !< (v1/v3) rh increment for minimum evaporation of rain
  real(kind_phys) :: rh_ins         = 0.25    !< (v1/v3) rh increment for sublimation of snow
  real(kind_phys) :: tau_r2g        = 900.    !< (v1/v3) rain freezing during fast_sat time scale (s)
  real(kind_phys) :: tau_smlt       = 900.    !< (v1/v3) snow melting time scale (s)
  real(kind_phys) :: tau_i2s        = 1000.   !< (v1/v3) cloud ice to snow auto-conversion time scale (s)
  real(kind_phys) :: tau_l2r        = 900.    !< (v1/v3) cloud water to rain auto-conversion time scale (s)
  real(kind_phys) :: tau_v2l        = 150.    !< (v1/v3) water vapor to cloud water (condensation) time scale (s)
  real(kind_phys) :: tau_l2v        = 300.    !< (v1/v3) cloud water to water vapor (evaporation) time scale (s)
  real(kind_phys) :: dw_land        = 0.20    !< (v1/v3) value for subgrid deviation / variability over land
  real(kind_phys) :: dw_ocean       = 0.10    !< (v1/v3) base value for ocean
  real(kind_phys) :: ccn_o          = 90.     !< (v1/v3) ccn over ocean (cm^ - 3)
  real(kind_phys) :: ccn_l          = 270.    !< (v1/v3) ccn over land (cm^ - 3)
  real(kind_phys) :: sat_adj0       = 0.90    !< (v1/v3) adjustment factor (0: no, 1: full) during fast_sat_adj
  real(kind_phys) :: qi_lim         = 1.      !< (v1/v3) cloud ice limiter (0: no, 1: full, >1: extra) to prevent large ice build up
  real(kind_phys) :: ql_mlt         = 2.0e-3  !< (v1/v3) max value of cloud water allowed from melted cloud ice
  real(kind_phys) :: qs_mlt         = 1.0e-6  !< (v1/v3) max cloud water due to snow melt
  real(kind_phys) :: ql_gen         = 1.0e-3  !< (v1/v3) max cloud water generation during remapping step if fast_sat_adj = .t.
  real(kind_phys) :: qi_gen         = 1.82e-6 !< (v1/v3) max cloud ice generation during remapping step (V1 ONLY. Computed internally in V3)
  real(kind_phys) :: ql0_max        = 2.0e-3  !< (v1/v3) max cloud water value (auto converted to rain)
  real(kind_phys) :: qi0_max        = 1.0e-4  !< (v1/v3) max cloud ice value (by other sources)
  real(kind_phys) :: qi0_crt        = 1.0e-4  !< (v1/v3) cloud ice to snow autoconversion threshold (was 1.e-4);
                                              !< qi0_crt is highly dependent on horizontal resolution
  real(kind_phys) :: qs0_crt        = 1.0e-3  !< (v1/v3) snow to graupel density threshold (0.6e-3 in purdue lin scheme)
  real(kind_phys) :: c_paut         = 0.55    !< (v1/v3) autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
  real(kind_phys) :: vi_fac         = 1.      !< (v1/v3) if const_vi: 1 / 3
  real(kind_phys) :: vs_fac         = 1.      !< (v1/v3) if const_vs: 1.
  real(kind_phys) :: vg_fac         = 1.      !< (v1/v3) if const_vg: 2.
  real(kind_phys) :: vr_fac         = 1.      !< (v1/v3) if const_vr: 4.
  real(kind_phys) :: vr_max         = 12.     !< (v1/v3) max fall speed for rain
  real(kind_phys) :: rewmin         = 5.0     !< (v1/v3) minimum effective radii (liquid)
  real(kind_phys) :: reimin         = 10.0    !< (v1/v3) minimum effective radii (ice)
  real(kind_phys) :: reimax         = 150.0   !< (v1/v3) maximum effective radii (ice)
  real(kind_phys) :: rermax         = 10000.0 !< (v1/v3) maximum effective radii (rain)
  real(kind_phys) :: resmin         = 150.0   !< (v1/v3) minimum effective radii (snow)
  real(kind_phys) :: resmax         = 10000.0 !< (v1/v3) maximum effective radii (snow)
  real(kind_phys) :: regmax         = 10000.0 !< (v1/v3) maximum effective radii (graupel)
  !
  logical :: const_vi               = .false. !< (v1/v3) if .t. the constants are specified by v * _fac
  logical :: const_vs               = .false. !< (v1/v3) if .t. the constants are specified by v * _fac
  logical :: const_vg               = .false. !< (v1/v3) if .t. the constants are specified by v * _fac
  logical :: const_vr               = .false. !< (v1/v3) if .t. the constants are specified by v * _fac
  logical :: z_slope_liq            = .true.  !< (v1/v3) use linear mono slope for autocconversions
  logical :: do_hail                = .false. !< (v1/v3) use hail parameters instead of graupel
  logical :: do_qa                  = .true.  !< (v1/v3) do inline cloud fraction
  logical :: rad_snow               = .true.  !< (v1/v3) consider snow in cloud fraciton calculation
  logical :: rad_graupel            = .true.  !< (v1/v3) consider graupel in cloud fraction calculation
  logical :: rad_rain               = .true.  !< (v1/v3) consider rain in cloud fraction calculation
  logical :: do_sedi_heat           = .true.  !< (v1/v3) transport of heat in sedimentation
  logical :: prog_ccn               = .false. !< (v1/v3) do prognostic ccn (yi ming's method)
  logical :: tintqs                 = .false. !< (v1/v3)
  !
  integer :: icloud_f               = 0       !< (v1/v3) GFDL cloud scheme
                                              !< 0: subgrid variability based scheme
                                              !< 1: same as 0, but for old fvgfs implementation
                                              !< 2: binary cloud scheme
                                              !< 3: extension of 0
  integer :: irain_f                = 0       !< (v1/v3) cloud water to rain auto conversion scheme
                                              !< 0: subgrid variability based scheme
                                              !< 1: no subgrid varaibility

  ! #####################################################################################
  ! GFDL MP common (v1/v3) parameters, with different default values
  ! #####################################################################################
  real(kind_phys) :: tice           = 273.16  !< freezing temperature (K): ref: GFDL, GFS (DJS: V3=273.15)
  real(kind_phys) :: tau_imlt       = 600.    !< cloud ice melting time scale (s) (DJS: V3=1200.)
  real(kind_phys) :: rthresh        = 10.0e-6 !< critical cloud drop radius (micro m) (DJS: v3=20.0e-6)
  real(kind_phys) :: c_psaci        = 0.02    !< accretion: cloud ice to snow (was 0.1 in zetac) (DJS: v3=0.05)
  real(kind_phys) :: c_pgacs        = 2.0e-3  !< snow to graupel "accretion" eff. (was 0.1 in zetac) (DJS: v3=0.01)
  real(kind_phys) :: vi_max         = 0.5     !< max fall speed for ice (DJS: v3=1.0)
  real(kind_phys) :: vs_max         = 5.0     !< max fall speed for snow (DJS: v3=2.0)
  real(kind_phys) :: vg_max         = 8.0     !< max fall speed for graupel (DJS: v3=12.0)
  real(kind_phys) :: rewmax         = 10.0    !< maximum effective radii (liquid) (DJS: v3=15.0)
  real(kind_phys) :: rermin         = 10.0    !< minimum effective radii (rain) (DJS: v3=15.0)
  real(kind_phys) :: regmin         = 300.0   !< minimum effective radii (graupel) (DJS: v3=150.0)
  logical :: z_slope_ice            = .false. !< use linear mono slope for autocconversions (DJS: v3=.true.)
  logical :: do_sedi_w              = .false. !< transport of vertical motion in sedimentation (DJS: v3=.true.)
  logical :: fix_negative           = .false. !< fix negative water species (DJS: v3=.true.)
  integer :: reiflag                = 1       !< cloud ice effective radius scheme (DJS: v3=5)
                                              !< 1: Heymsfield and Mcfarquhar (1996)
                                              !< 2: Donner et al. (1997)
                                              !< 3: Fu (2007)
                                              !< 4: Kristjansson et al. (2000)
                                              !< 5: Wyser (1998)
                                              !< 6: Sun and Rikus (1999), Sun (2001)
                                              !< 7: effective radius

  ! #####################################################################################
  ! GFDL MP Version 3 parameters
  ! #####################################################################################
  logical :: const_vw               = .false. !< if .ture., the constants are specified by v * _fac
  logical :: do_sedi_uv             = .true.  !< transport of horizontal momentum in sedimentation
  logical :: do_sedi_melt           = .true.  !< melt cloud ice, snow, and graupel during sedimentation
  logical :: liq_ice_combine        = .false. !< combine all liquid water, combine all solid water
  logical :: snow_grauple_combine   = .true.  !< combine snow and graupel
  logical :: use_rhc_cevap          = .false. !< cap of rh for cloud water evaporation
  logical :: use_rhc_revap          = .false. !< cap of rh for rain evaporation
  logical :: do_cld_adj             = .false. !< do cloud fraction adjustment
  logical :: do_evap_timescale      = .true.  !< whether to apply a timescale to evaporation
  logical :: do_cond_timescale      = .false. !< whether to apply a timescale to condensation    
  logical :: consv_checker          = .false. !< turn on energy and water conservation checker
  logical :: do_warm_rain_mp        = .false. !< do warm rain cloud microphysics only
  logical :: do_wbf                 = .false. !< do Wegener Bergeron Findeisen process
  logical :: do_psd_water_fall      = .false. !< calculate cloud water terminal velocity based on PSD
  logical :: do_psd_ice_fall        = .false. !< calculate cloud ice terminal velocity based on PSD
  logical :: do_psd_water_num       = .false. !< calculate cloud water number concentration based on PSD
  logical :: do_psd_ice_num         = .false. !< calculate cloud ice number concentration based on PSD
  logical :: do_new_acc_water       = .false. !< perform the new accretion for cloud water
  logical :: do_new_acc_ice         = .false. !< perform the new accretion for cloud ice
  logical :: cp_heating             = .false. !< update temperature based on constant pressure
  logical :: delay_cond_evap        = .false. !< do condensation evaporation only at the last time step
  logical :: do_subgrid_proc        = .true.  !< do temperature sentive high vertical resolution processes
  logical :: fast_fr_mlt            = .true.  !< do freezing and melting in fast microphysics
  logical :: fast_dep_sub           = .true.  !< do deposition and sublimation in fast microphysics
  integer :: ntimes                 = 1       !< cloud microphysics sub cycles
  integer :: nconds                 = 1       !< condensation sub cycles
  integer :: inflag                 = 1       !< ice nucleation scheme
                                              !< 1: Hong et al. (2004)
                                              !< 2: Meyers et al. (1992)
                                              !< 3: Meyers et al. (1992)
                                              !< 4: Cooper (1986)
                                              !< 5: Fletcher (1962)
  integer :: igflag                 = 3       !< ice generation scheme
                                              !< 1: WSM6
                                              !< 2: WSM6 with 0 at 0 C
                                              !< 3: WSM6 with 0 at 0 C and fixed value at - 10 C
                                              !< 4: combination of 1 and 3
  integer :: ifflag                 = 1       !< ice fall scheme
                                              !< 1: Deng and Mace (2008)
                                              !< 2: Heymsfield and Donner (1990)
  integer :: rewflag                = 1       !< cloud water effective radius scheme
                                              !< 1: Martin et al. (1994)
                                              !< 2: Martin et al. (1994), GFDL revision
                                              !< 3: Kiehl et al. (1994)
                                              !< 4: effective radius
  integer :: rerflag                = 1       !< rain effective radius scheme
                                              !< 1: effective radius
  integer :: resflag                = 1       !< snow effective radius scheme
                                              !< 1: effective radius
  integer :: regflag                = 1       !< graupel effective radius scheme
                                              !< 1: effective radius
  integer :: radr_flag              = 1       !< radar reflectivity for rain
                                              !< 1: Mark Stoelinga (2005)
                                              !< 2: Smith et al. (1975), Tong and Xue (2005)
                                              !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
  integer :: rads_flag              = 1       !< radar reflectivity for snow
                                              !< 1: Mark Stoelinga (2005)
                                              !< 2: Smith et al. (1975), Tong and Xue (2005)
                                              !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
  integer :: radg_flag              = 1       !< radar reflectivity for graupel
                                              !< 1: Mark Stoelinga (2005)
                                              !< 2: Smith et al. (1975), Tong and Xue (2005)
                                              !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
  integer :: sedflag                = 1       !< sedimentation scheme
                                              !< 1: implicit scheme
                                              !< 2: explicit scheme
                                              !< 3: lagrangian scheme
                                              !< 4: combined implicit and lagrangian scheme
  integer :: vdiffflag              = 1       !< wind difference scheme in accretion
                                              !< 1: Wisner et al. (1972)
                                              !< 2: Mizuno (1990)
                                              !< 3: Murakami (1990)
  real(kind_phys) :: c_psacw        = 1.0     !< cloud water to snow accretion efficiency
  real(kind_phys) :: c_pracw        = 0.8     !< cloud water to rain accretion efficiency
  real(kind_phys) :: c_praci        = 1.0     !< cloud ice to rain accretion efficiency
  real(kind_phys) :: c_pgacw        = 1.0     !< cloud water to graupel accretion efficiency
  real(kind_phys) :: c_pgaci        = 0.05    !< cloud ice to graupel accretion efficiency (was 0.1 in ZETAC)
  real(kind_phys) :: c_pracs        = 1.0     !< snow to rain accretion efficiency
  real(kind_phys) :: c_psacr        = 1.0     !< rain to snow accretion efficiency
  real(kind_phys) :: c_pgacr        = 1.0     !< rain to graupel accretion efficiency
  real(kind_phys) :: alinw          = 3.e7    !< "a" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
  real(kind_phys) :: alini          = 7.e2    !< "a" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
  real(kind_phys) :: alinr          = 842.0   !< "a" in Lin et al. (1983) for rain (Liu and Orville 1969)
  real(kind_phys) :: alins          = 4.8     !< "a" in Lin et al. (1983) for snow (straka 2009)
  real(kind_phys) :: aling          = 1.0     !< "a" in Lin et al. (1983), similar to a, but for graupel (Pruppacher and Klett 2010)
  real(kind_phys) :: alinh          = 1.0     !< "a" in Lin et al. (1983), similar to a, but for hail (Pruppacher and Klett 2010)
  real(kind_phys) :: blinw          = 2.0     !< "b" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
  real(kind_phys) :: blini          = 1.0     !< "b" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
  real(kind_phys) :: blinr          = 0.8     !< "b" in Lin et al. (1983) for rain (Liu and Orville 1969)
  real(kind_phys) :: blins          = 0.25    !< "b" in Lin et al. (1983) for snow (straka 2009)
  real(kind_phys) :: bling          = 0.5     !< "b" in Lin et al. (1983), similar to b, but for graupel (Pruppacher and Klett 2010)
  real(kind_phys) :: blinh          = 0.5     !< "b" in Lin et al. (1983), similar to b, but for hail (Pruppacher and Klett 2010)
  real(kind_phys) :: vw_fac         = 1.0     !<
  real(kind_phys) :: vw_max         = 0.01    !< maximum fall speed for cloud water (m/s)
  real(kind_phys) :: tice_mlt       = 273.16  !< can set ice melting temperature to 268 based on observation (Kay et al. 2016) (K)
  real(kind_phys) :: tau_gmlt       = 600.0   !< graupel melting time scale (s)
  real(kind_phys) :: tau_wbf        = 300.0   !< graupel melting time scale (s)
  real(kind_phys) :: tau_revp       = 0.0     !< rain evaporation time scale (s)
  real(kind_phys) :: is_fac         = 0.2     !< cloud ice sublimation temperature factor
  real(kind_phys) :: ss_fac         = 0.2     !< snow sublimation temperature factor
  real(kind_phys) :: gs_fac         = 0.2     !< graupel sublimation temperature factor
  real(kind_phys) :: rh_fac_evap    = 10.0    !< cloud water evaporation relative humidity factor
  real(kind_phys) :: rh_fac_cond    = 10.0    !< cloud water condensation relative humidity factor
  real(kind_phys) :: sed_fac        = 1.0     !< coefficient for sedimentation fall, scale from 1.0 (implicit) to 0.0 (lagrangian)
  real(kind_phys) :: xr_a           = 0.25    !< p value in Xu and Randall (1996)
  real(kind_phys) :: xr_b           = 100.0   !< alpha_0 value in Xu and Randall (1996)
  real(kind_phys) :: xr_c           = 0.49    !< gamma value in Xu and Randall (1996)
  real(kind_phys) :: te_err         = 1.e-5   !< 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time
  real(kind_phys) :: tw_err         = 1.e-8   !< 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time
  real(kind_phys) :: rh_thres       = 0.75    !< minimum relative humidity for cloud fraction
  real(kind_phys) :: rhc_cevap      = 0.85    !< maximum relative humidity for cloud water evaporation
  real(kind_phys) :: rhc_revap      = 0.85    !< maximum relative humidity for rain evaporation
  real(kind_phys) :: f_dq_p         = 1.0     !< cloud fraction adjustment for supersaturation
  real(kind_phys) :: f_dq_m         = 1.0     !< cloud fraction adjustment for undersaturation
  real(kind_phys) :: fi2s_fac       = 1.0     !< maximum sink of cloud ice to form snow: 0-1
  real(kind_phys) :: fi2g_fac       = 1.0     !< maximum sink of cloud ice to form graupel: 0-1
  real(kind_phys) :: fs2g_fac       = 1.0     !< maximum sink of snow to form graupel: 0-1
  real(kind_phys) :: n0w_sig        = 1.1     !< intercept parameter (significand) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
  real(kind_phys) :: n0i_sig        = 1.3     !< intercept parameter (significand) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
  real(kind_phys) :: n0r_sig        = 8.0     !< intercept parameter (significand) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
  real(kind_phys) :: n0s_sig        = 3.0     !< intercept parameter (significand) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
  real(kind_phys) :: n0g_sig        = 4.0     !< intercept parameter (significand) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
  real(kind_phys) :: n0h_sig        = 4.0     !< intercept parameter (significand) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)
  real(kind_phys) :: n0w_exp        = 41      !< intercept parameter (exponent) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
  real(kind_phys) :: n0i_exp        = 18      !< intercept parameter (exponent) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
  real(kind_phys) :: n0r_exp        = 6       !< intercept parameter (exponent) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
  real(kind_phys) :: n0s_exp        = 6       !< intercept parameter (exponent) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
  real(kind_phys) :: n0g_exp        = 6       !< intercept parameter (exponent) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
  real(kind_phys) :: n0h_exp        = 4       !< intercept parameter (exponent) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)    
  real(kind_phys) :: muw            = 6.0     !< shape parameter of cloud water in Gamma distribution (Martin et al. 1994)
  real(kind_phys) :: mui            = 3.35    !< shape parameter of cloud ice in Gamma distribution (McFarquhar et al. 2015)
  real(kind_phys) :: mur            = 1.0     !< shape parameter of rain in Gamma distribution (Marshall and Palmer 1948)
  real(kind_phys) :: mus            = 1.0     !< shape parameter of snow in Gamma distribution (Gunn and Marshall 1958)
  real(kind_phys) :: mug            = 1.0     !< shape parameter of graupel in Gamma distribution (Houze et al. 1979)
  real(kind_phys) :: muh            = 1.0     !< shape parameter of hail in Gamma distribution (Federer and Waldvogel 1975)
  real(kind_phys) :: beta           = 1.22    !< defined in Heymsfield and Mcfarquhar (1996)
  real(kind_phys) :: rewfac         = 1.0     !< this is a tuning parameter to compromise the inconsistency between
                                              !< GFDL MP's PSD and cloud water radiative property's PSD assumption.
                                              !< after the cloud water radiative property's PSD is rebuilt,
                                              !< this parameter should be 1.0.
  real(kind_phys) :: reifac         = 1.0     !< this is a tuning parameter to compromise the inconsistency between
                                              !< GFDL MP's PSD and cloud ice radiative property's PSD assumption.
                                              !< after the cloud ice radiative property's PSD is rebuilt,
                                              !< this parameter should be 1.0.

  ! #######################################################################################
  ! NAMELISTS
  ! #######################################################################################
  
  ! V1 namelist
  namelist / gfdl_cloud_microphysics_nml /                                                &
       mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, vi_fac,      &
       vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max, vs_max, vg_max,       &
       vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max, qi0_crt, qr0_crt, fast_sat_adj, &
       rh_inc, rh_ins, rh_inr, const_vi, const_vs, const_vg, const_vr, use_ccn, rthresh,  &
       ccn_l, ccn_o, qc_crt, tau_g2v, tau_v2g, sat_adj0, c_piacr, tau_imlt, tau_v2l,      &
       tau_l2v, tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, z_slope_liq,  &
       z_slope_ice, prog_ccn, c_cracw, alin, clin, tice, rad_snow, rad_graupel, rad_rain, &
       cld_min, use_ppm, mono_prof, do_sedi_heat, sedi_transport, do_sedi_w, de_ice,      &
       icloud_f, irain_f, mp_print, reiflag, rewmin, rewmax, reimin, reimax, rermin,      &
       rermax, resmin, resmax, regmin, regmax, tintqs, do_hail

  ! V3 Namelist
  namelist / gfdl_cloud_microphysics_v3_nml /                                             &
       t_min, t_sub, tau_r2g, tau_smlt, tau_gmlt, dw_land, dw_ocean, vw_fac, vi_fac,      &
       vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vw_max, vi_max, vs_max,       &
       vg_max, vr_max, qs_mlt, qs0_crt, ql0_max, qi0_max, qi0_crt, ifflag, rh_inc, rh_ins,&
       rh_inr, const_vw, const_vi, const_vs, const_vg, const_vr, rthresh, ccn_l, ccn_o,   &
       igflag, c_paut, tau_imlt, tau_v2l, tau_l2v, tau_i2s, tau_l2r, qi_lim, ql_gen,      &
       do_hail, inflag, c_psacw, c_psaci, c_pracs, c_psacr, c_pgacr, c_pgacs, c_pgacw,    &
       c_pgaci, z_slope_liq, z_slope_ice, prog_ccn, c_pracw, c_praci, rad_snow,           &
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
       delay_cond_evap, do_subgrid_proc, fast_fr_mlt, fast_dep_sub
  !
  public                                                                                  &
       tau_g2r, tau_g2v, tau_v2g, qc_crt, qr0_crt, c_piacr, c_cracw, alin, clin,          &
       fast_sat_adj, use_ccn, use_ppm, mono_prof, mp_print, de_ice, sedi_transport,       &
       t_min, t_sub, tau_r2g, tau_smlt, tau_gmlt, dw_land, dw_ocean, vw_fac, vi_fac,      &
       vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vw_max, vi_max, vs_max,       &
       vg_max, vr_max, qs_mlt, qs0_crt, ql0_max, qi0_max, qi0_crt, ifflag, rh_inc, rh_ins,&
       rh_inr, const_vw, const_vi, const_vs, const_vg, const_vr, rthresh, ccn_l, ccn_o,   &
       igflag, c_paut, tau_imlt, tau_v2l, tau_l2v, tau_i2s, tau_l2r, qi_lim, ql_gen,      &
       do_hail, inflag, c_psacw, c_psaci, c_pracs, c_psacr, c_pgacr, c_pgacs, c_pgacw,    &
       c_pgaci, z_slope_liq, z_slope_ice, prog_ccn, c_pracw, c_praci, rad_snow,           &
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
       delay_cond_evap, do_subgrid_proc, fast_fr_mlt, fast_dep_sub, qi_gen, sat_adj0,     &
       tice, tintqs
contains

  ! #######################################################################################
  ! Procedure to read GFDLMP namelists
  ! #######################################################################################
  subroutine read_gfdlmp_nml(errmsg, errflg, unit, input_nml_file, fn_nml, version, iostat)
    
    character(len = *), intent(in ), optional  :: input_nml_file(:)
    character(len = *), intent(in ), optional  :: fn_nml
    integer,            intent(in ), optional  :: unit
    integer,            intent(in ), optional  :: version
    integer,            intent(out), optional  :: iostat
    character(len=*),   intent(out), optional  :: errmsg
    integer,            intent(out), optional  :: errflg
    logical :: exists
    ! Make sure that all inputs to read appropriate NML are provided, if not use default
    ! parameters
    if (present(unit)           .and. present(iostat) .and. &
        present(input_nml_file) .and. present(fn_nml) .and. &
        present(version)        .and. present(errflg) .and. &
        present(errmsg)) then

       if ((version .ne. 1) .and. (version .ne. 3)) then
          write (6, *) 'gfdl - mp :: invalid scheme version number'
          errflg = 1
          errmsg = 'ERROR(module_gfdlmp_param): invalid scheme version number'
          return
       endif
    
#ifdef INTERNAL_FILE_NML
       if (version==1) read (input_nml_file, nml = gfdl_cloud_microphysics_nml)
       if (version==3) read (input_nml_file, nml = gfdl_cloud_microphysics_v3_nml)
#else
       inquire (file = trim (fn_nml), exist = exists)
       if (.not. exists) then
          write (6, *) 'gfdl - mp :: namelist file: ', trim (fn_nml), ' does not exist'
          errflg = 1
          errmsg = 'ERROR(module_gfdlmp_param): namelist file '//trim (fn_nml)//' does not exist'
          return
       else
          open (unit = unit, file = fn_nml, action = 'read' , status = 'old', iostat = iostat)
       endif
       rewind (unit)
       if (version==1) read (unit, nml = gfdl_cloud_microphysics_nml)
       if (version==3) read (unit, nml = gfdl_cloud_microphysics_v3_nml)
       close (unit)
#endif
    endif
  end subroutine read_gfdlmp_nml  
  !
end module module_gfdlmp_param
