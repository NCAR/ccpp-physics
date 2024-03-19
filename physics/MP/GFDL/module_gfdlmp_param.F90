module module_gfdlmp_param
  use machine, only: kind_phys
  implicit none

  type ty_gfdlmp_config
     !
     real(kind_phys) :: cld_min       !< minimum cloud fraction
     real(kind_phys) :: tice          !< set tice = 165. to trun off ice - phase phys (kessler emulator)
     real(kind_phys) :: t_min         !< min temp to freeze - dry all water vapor
     real(kind_phys) :: t_sub         !< min temp for sublimation of cloud ice
     real(kind_phys) :: tice_mlt      !< can set ice melting temperature to 268 based on observation (Kay et al. 2016) (K) 
     real(kind_phys) :: mp_time       !< maximum micro - physics time step (sec)

     !
     real(kind_phys) :: rh_inc        !< rh increment for complete evaporation of cloud water and cloud ice
     real(kind_phys) :: rh_inr        !< rh increment for minimum evaporation of rain
     real(kind_phys) :: rh_ins        !< rh increment for sublimation of snow

     !
     real(kind_phys) :: tau_r2g       !< rain freezing during fast_sat
     real(kind_phys) :: tau_smlt      !< snow melting
     real(kind_phys) :: tau_g2r       !< graupel melting to rain
     real(kind_phys) :: tau_imlt      !< cloud ice melting
     real(kind_phys) :: tau_i2s       !< cloud ice to snow auto-conversion
     real(kind_phys) :: tau_l2r       !< cloud water to rain auto-conversion
     real(kind_phys) :: tau_v2l       !< water vapor to cloud water (condensation)
     real(kind_phys) :: tau_l2v       !< cloud water to water vapor (evaporation)
     real(kind_phys) :: tau_g2v       !< graupel sublimation
     real(kind_phys) :: tau_v2g       !< graupel deposition -- make it a slow process
     real(kind_phys) :: tau_gmlt      !< graupel melting time scale (s)
     real(kind_phys) :: tau_wbf       !< graupel melting time scale (s)
     real(kind_phys) :: tau_revp      !< rain evaporation time scale (s)

     ! horizontal subgrid variability
     real(kind_phys) :: dw_land       !< value for subgrid deviation / variability over land
     real(kind_phys) :: dw_ocean      !< base value for ocean

     ! prescribed ccn
     real(kind_phys) :: ccn_o         !< ccn over ocean (cm^ - 3)
     real(kind_phys) :: ccn_l         !< ccn over land (cm^ - 3)

     !
     real(kind_phys) :: rthresh       !< critical cloud drop radius (micro m)
     real(kind_phys) :: sat_adj0      !< adjustment factor (0: no, 1: full) during fast_sat_adj
     real(kind_phys) :: qc_crt        !< mini condensate mixing ratio to allow partial cloudiness
     real(kind_phys) :: qi_lim        !< cloud ice limiter to prevent large ice build up
     real(kind_phys) :: ql_mlt        !< max value of cloud water allowed from melted cloud ice
     real(kind_phys) :: qs_mlt        !< max cloud water due to snow melt
     real(kind_phys) :: ql_gen        !< max cloud water generation during remapping step if fast_sat_adj = .t.
     real(kind_phys) :: qi_gen        !< max cloud ice generation during remapping step (V1 ONLY. Computed internally in V3)

     ! cloud condensate upper bounds: "safety valves" for ql & qi
     real(kind_phys) :: ql0_max       !< max cloud water value (auto converted to rain)
     real(kind_phys) :: qi0_max       !< max cloud ice value (by other sources)

     !
     real(kind_phys) :: qi0_crt       !< cloud ice to snow autoconversion threshold (was 1.e-4);
                                      !< qi0_crt is highly dependent on horizontal resolution
     real(kind_phys) :: qr0_crt       !< rain to snow or graupel/hail threshold
                                      !< lfo used * mixing ratio * = 1.e-4 (hail in lfo)
     real(kind_phys) :: qs0_crt       !< snow to graupel density threshold (0.6e-3 in purdue lin scheme)

     !
     real(kind_phys) :: c_paut        !< autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
     real(kind_phys) :: c_psacw       !< cloud water to snow accretion efficiency 
     real(kind_phys) :: c_psaci       !< accretion: cloud ice to snow (was 0.1 in zetac)
     real(kind_phys) :: c_pracw       !< cloud water to rain accretion efficiency
     real(kind_phys) :: c_praci       !< cloud ice to rain accretion efficiency 
     real(kind_phys) :: c_pracs       !< snow to rain accretion efficiency
     real(kind_phys) :: c_piacr       !< accretion: rain to ice:
     real(kind_phys) :: c_cracw       !< rain accretion efficiency
     real(kind_phys) :: c_pgacs       !< snow to graupel "accretion" eff. (was 0.1 in zetac)
     real(kind_phys) :: c_pgacw       !< cloud water to graupel accretion efficiency 
     real(kind_phys) :: c_psacr       !< rain to snow accretion efficiency
     real(kind_phys) :: c_pgaci       !< cloud ice to graupel accretion efficiency (was 0.1 in ZETAC)
     real(kind_phys) :: c_pgacr       !< rain to graupel accretion efficiency 

     !
     real(kind_phys) :: is_fac        !< cloud ice sublimation temperature factor
     real(kind_phys) :: ss_fac        !< snow sublimation temperature factor
     real(kind_phys) :: gs_fac        !< graupel sublimation temperature factor
     
     !
     real(kind_phys) :: rh_fac_evap   !< cloud water evaporation relative humidity factor
     real(kind_phys) :: rh_fac_cond   !< cloud water condensation relative humidity factor

     ! decreasing clin to reduce csacw (so as to reduce cloud water --- > snow)
     real(kind_phys) :: alin          !< "a" in lin1983
     real(kind_phys) :: clin          !< "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)

     ! fall velocity tuning constants:
     logical :: const_vi              !< if .t. the constants are specified by v * _fac
     logical :: const_vs              !< if .t. the constants are specified by v * _fac
     logical :: const_vg              !< if .t. the constants are specified by v * _fac
     logical :: const_vr              !< if .t. the constants are specified by v * _fac

     !
     logical :: liq_ice_combine       !< combine all liquid water, combine all solid water
     logical :: snow_grauple_combine  !< combine snow and graupel
     logical :: use_rhc_cevap         !< cap of rh for cloud water evaporation (V3)
     logical :: use_rhc_revap         !< cap of rh for rain evaporation (V3)
     !
     real(kind_phys) :: sed_fac       !< coefficient for sedimentation fall, scale from 1.0 (implicit) to 0.0 (lagrangian) 
     real(kind_phys) :: vw_fac        !<
     real(kind_phys) :: vi_fac        !< if const_vi: 1 / 3
     real(kind_phys) :: vs_fac        !< if const_vs: 1.
     real(kind_phys) :: vg_fac        !< if const_vg: 2.
     real(kind_phys) :: vr_fac        !< if const_vr: 4.

     ! upper bounds of fall speed (with variable speed option)
     real(kind_phys) :: vw_max        !< maximum fall speed for cloud water (m/s) (V3)
     real(kind_phys) :: vi_max        !< max fall speed for ice
     real(kind_phys) :: vs_max        !< max fall speed for snow
     real(kind_phys) :: vg_max        !< max fall speed for graupel
     real(kind_phys) :: vr_max        !< max fall speed for rain
     !
     real(kind_phys) :: xr_a          !< p value in Xu and Randall (1996)
     real(kind_phys) :: xr_b          !< alpha_0 value in Xu and Randall (1996)
     real(kind_phys) :: xr_c          !< gamma value in Xu and Randall (1996)
     !
     real(kind_phys) :: te_err        !< 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time
     real(kind_phys) :: tw_err        !< 64bit: 1.e-14, 32bit: 1.e-7; turn off to save computer time
     real(kind_phys) :: rh_thres      !< minimum relative humidity for cloud fraction
     real(kind_phys) :: rhc_cevap     !< maximum relative humidity for cloud water evaporation
     real(kind_phys) :: rhc_revap     !< maximum relative humidity for rain evaporation 
     real(kind_phys) :: f_dq_p        !< cloud fraction adjustment for supersaturation
     real(kind_phys) :: f_dq_m        !< cloud fraction adjustment for undersaturation
     real(kind_phys) :: fi2s_fac      !< maximum sink of cloud ice to form snow: 0-1
     real(kind_phys) :: fi2g_fac      !< maximum sink of cloud ice to form graupel: 0-1
     real(kind_phys) :: fs2g_fac      !< maximum sink of snow to form graupel: 0-1

     ! cloud microphysics switchers
     logical :: fast_sat_adj          !< has fast saturation adjustments
     logical :: z_slope_liq           !< use linear mono slope for autocconversions
     logical :: z_slope_ice           !< use linear mono slope for autocconversions
     logical :: use_ccn               !< must be true when prog_ccn is false
     logical :: use_ppm               !< use ppm fall scheme
     logical :: mono_prof             !< perform terminal fall with mono ppm scheme
     logical :: mp_print              !< cloud microphysics debugging printout
     logical :: do_hail               !< use hail parameters instead of graupel
     logical :: do_qa                 !< do inline cloud fraction
     logical :: rad_snow              !< consider snow in cloud fraciton calculation
     logical :: rad_graupel           !< consider graupel in cloud fraction calculation
     logical :: rad_rain              !< consider rain in cloud fraction calculation
     logical :: do_cld_adj            !< do cloud fraction adjustment 
     logical :: do_evap_timescale     !< whether to apply a timescale to evaporation 
     logical :: do_cond_timescale     !< whether to apply a timescale to condensation 
     logical :: consv_checker         !< turn on energy and water conservation checker  
     logical :: do_warm_rain_mp       !< do warm rain cloud microphysics only
     logical :: do_wbf                !< do Wegener Bergeron Findeisen process 
     logical :: do_psd_water_fall     !< calculate cloud water terminal velocity based on PSD
     logical :: do_psd_ice_fall       !< calculate cloud ice terminal velocity based on PSD
     logical :: do_psd_water_num      !< calculate cloud water number concentration based on PSD
     logical :: do_psd_ice_num        !< calculate cloud ice number concentration based on PSD
     logical :: do_new_acc_water      !< perform the new accretion for cloud water
     logical :: do_new_acc_ice        !< perform the new accretion for cloud ice
     logical :: cp_heating            !< update temperature based on constant pressure
     logical :: delay_cond_evap       !< do condensation evaporation only at the last time step
     logical :: do_subgrid_proc       !< do temperature sentive high vertical resolution processes
     logical :: fast_fr_mlt           !< do freezing and melting in fast microphysics
     logical :: fast_dep_sub          !< do deposition and sublimation in fast microphysics
     integer :: ntimes                !< cloud microphysics sub cycles
     integer :: nconds                !< condensation sub cycles
     !
     integer :: icloud_f              !< GFDL cloud scheme
                                      !< 0: subgrid variability based scheme
                                      !< 1: same as 0, but for old fvgfs implementation
                                      !< 2: binary cloud scheme
                                      !< 3: extension of 0
     !
     integer :: irain_f               !< cloud water to rain auto conversion scheme
                                      !< 0: subgrid variability based scheme
                                      !< 1: no subgrid varaibility
     !
     integer :: inflag                !< ice nucleation scheme
                                      !< 1: Hong et al. (2004)
                                      !< 2: Meyers et al. (1992)
                                      !< 3: Meyers et al. (1992)
                                      !< 4: Cooper (1986)
                                      !, 5: Fletcher (1962)
     !
     integer :: igflag                !< ice generation scheme
                                      !< 1: WSM6
                                      !< 2: WSM6 with 0 at 0 C
                                      !< 3: WSM6 with 0 at 0 C and fixed value at - 10 C
                                      !< 4: combination of 1 and 3
     !
     integer :: ifflag                !< ice fall scheme
                                      !< 1: Deng and Mace (2008)
                                      !< 2: Heymsfield and Donner (1990)
     !
     integer :: rewflag               !< cloud water effective radius scheme
                                      !< 1: Martin et al. (1994)
                                      !< 2: Martin et al. (1994), GFDL revision
                                      !< 3: Kiehl et al. (1994)
                                      !< 4: effective radius
     !
     integer :: reiflag               !< cloud ice effective radius scheme
                                      !< 1: Heymsfield and Mcfarquhar (1996)
                                      !< 2: Donner et al. (1997)
                                      !< 3: Fu (2007)
                                      !< 4: Kristjansson et al. (2000)
                                      !< 5: Wyser (1998)
                                      !< 6: Sun and Rikus (1999), Sun (2001)
                                      !< 7: effective radius
     !
     integer :: rerflag               !< rain effective radius scheme
                                      !< 1: effective radius
     !
     integer :: resflag               !< snow effective radius scheme
                                      !< 1: effective radius
     !
     integer :: regflag               !< graupel effective radius scheme
                                      !< 1: effective radius
     !
     integer :: radr_flag             !< radar reflectivity for rain
                                      !< 1: Mark Stoelinga (2005)
                                      !< 2: Smith et al. (1975), Tong and Xue (2005)
                                      !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
     !
     integer :: rads_flag             !< radar reflectivity for snow
                                      !< 1: Mark Stoelinga (2005)
                                      !< 2: Smith et al. (1975), Tong and Xue (2005)
                                      !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
     !
     integer :: radg_flag             !< radar reflectivity for graupel
                                      !< 1: Mark Stoelinga (2005)
                                      !< 2: Smith et al. (1975), Tong and Xue (2005)
                                      !< 3: Marshall-Palmer formula (https://en.wikipedia.org/wiki/DBZ_(meteorology))
     !
     integer :: sedflag               !< sedimentation scheme
                                      !< 1: implicit scheme
                                      !< 2: explicit scheme
                                      !< 3: lagrangian scheme
                                      !< 4: combined implicit and lagrangian scheme
     !
     integer :: vdiffflag             !< wind difference scheme in accretion
                                      !< 1: Wisner et al. (1972)
                                      !< 2: Mizuno (1990)
                                      !< 3: Murakami (1990)
     !
     real(kind_phys) :: n0w_sig       !< intercept parameter (significand) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
     real(kind_phys) :: n0i_sig       !< intercept parameter (significand) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
     real(kind_phys) :: n0r_sig       !< intercept parameter (significand) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
     real(kind_phys) :: n0s_sig       !< intercept parameter (significand) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
     real(kind_phys) :: n0g_sig       !< intercept parameter (significand) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
     real(kind_phys) :: n0h_sig       !< intercept parameter (significand) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)
     !
     real(kind_phys) :: n0w_exp       !< intercept parameter (exponent) of cloud water (Lin et al. 1983) (1/m^4) (Martin et al. 1994)
     real(kind_phys) :: n0i_exp       !< intercept parameter (exponent) of cloud ice (Lin et al. 1983) (1/m^4) (McFarquhar et al. 2015)
     real(kind_phys) :: n0r_exp       !< intercept parameter (exponent) of rain (Lin et al. 1983) (1/m^4) (Marshall and Palmer 1948)
     real(kind_phys) :: n0s_exp       !< intercept parameter (exponent) of snow (Lin et al. 1983) (1/m^4) (Gunn and Marshall 1958)
     real(kind_phys) :: n0g_exp       !< intercept parameter (exponent) of graupel (Rutledge and Hobbs 1984) (1/m^4) (Houze et al. 1979)
     real(kind_phys) :: n0h_exp       !< intercept parameter (exponent) of hail (Lin et al. 1983) (1/m^4) (Federer and Waldvogel 1975)
     !
     real(kind_phys) :: muw           !< shape parameter of cloud water in Gamma distribution (Martin et al. 1994)
     real(kind_phys) :: mui           !< shape parameter of cloud ice in Gamma distribution (McFarquhar et al. 2015)
     real(kind_phys) :: mur           !< shape parameter of rain in Gamma distribution (Marshall and Palmer 1948)
     real(kind_phys) :: mus           !< shape parameter of snow in Gamma distribution (Gunn and Marshall 1958)
     real(kind_phys) :: mug           !< shape parameter of graupel in Gamma distribution (Houze et al. 1979)
     real(kind_phys) :: muh           !< shape parameter of hail in Gamma distribution (Federer and Waldvogel 1975)
     !
     real(kind_phys) :: alinw         !< "a" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
     real(kind_phys) :: alini         !< "a" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
     real(kind_phys) :: alinr         !< "a" in Lin et al. (1983) for rain (Liu and Orville 1969)
     real(kind_phys) :: alins         !< "a" in Lin et al. (1983) for snow (straka 2009)
     real(kind_phys) :: aling         !< "a" in Lin et al. (1983), similar to a, but for graupel (Pruppacher and Klett 2010)
     real(kind_phys) :: alinh         !< "a" in Lin et al. (1983), similar to a, but for hail (Pruppacher and Klett 2010)
     !
     real(kind_phys) :: blinw         !< "b" in Lin et al. (1983) for cloud water (Ikawa and Saito 1990)
     real(kind_phys) :: blini         !< "b" in Lin et al. (1983) for cloud ice (Ikawa and Saita 1990)
     real(kind_phys) :: blinr         !< "b" in Lin et al. (1983) for rain (Liu and Orville 1969)
     real(kind_phys) :: blins         !< "b" in Lin et al. (1983) for snow (straka 2009)
     real(kind_phys) :: bling         !< "b" in Lin et al. (1983), similar to b, but for graupel (Pruppacher and Klett 2010)
     real(kind_phys) :: blinh         !< "b" in Lin et al. (1983), similar to b, but for hail (Pruppacher and Klett 2010)
     !
     logical :: de_ice                !< to prevent excessive build - up of cloud ice from external sources
     logical :: sedi_transport        !< transport of momentum in sedimentation
     logical :: do_sedi_uv            !< transport of horizontal momentum in sedimentation
     logical :: do_sedi_w             !< transport of vertical motion in sedimentation
     logical :: do_sedi_heat          !< transport of heat in sedimentation
     logical :: do_sedi_melt          !< melt cloud ice, snow, and graupel during sedimentation
     logical :: prog_ccn              !< do prognostic ccn (yi ming's method)
     logical :: fix_negative          !< fix negative water species
     logical :: tintqs                !<
     !
     real(kind_phys) :: beta          !< defined in Heymsfield and Mcfarquhar (1996)

     real(kind_phys) :: rewmin        !< minimum effective radii (liquid)
     real(kind_phys) :: rewmax        !< maximum effective radii (liquid)
     real(kind_phys) :: reimin        !< minimum effective radii (ice)
     real(kind_phys) :: reimax        !< maximum effective radii (ice)
     real(kind_phys) :: rermin        !< minimum effective radii (rain)
     real(kind_phys) :: rermax        !< maximum effective radii (rain)
     real(kind_phys) :: resmin        !< minimum effective radii (snow)
     real(kind_phys) :: resmax        !< maximum effective radii (snow)
     real(kind_phys) :: regmin        !< minimum effective radii (graupel)
     real(kind_phys) :: regmax        !< maximum effective radii (graupel)
     !
     real(kind_phys) :: rewfac        !< this is a tuning parameter to compromise the inconsistency between
                                      !< GFDL MP's PSD and cloud water radiative property's PSD assumption.
                                      !< after the cloud water radiative property's PSD is rebuilt,
                                      !< this parameter should be 1.0.
     real(kind_phys) :: reifac        !< this is a tuning parameter to compromise the inconsistency between
                                      !< GFDL MP's PSD and cloud ice radiative property's PSD assumption.
                                      !< after the cloud ice radiative property's PSD is rebuilt,
                                      !< this parameter should be 1.0.
  end type ty_gfdlmp_config
  !
  type(ty_gfdlmp_config) :: cfg
  
contains
  
  subroutine setup_v3(cfg, mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
                      vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max,  &
                      vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max,    &
                      qi0_crt, qr0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr, const_vi,     &
                      const_vs, const_vg, const_vr, use_ccn, rthresh, ccn_l, ccn_o, qc_crt, &
                      tau_g2v, tau_v2g, sat_adj0, c_piacr, tau_imlt, tau_v2l, tau_l2v,      &
                      tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs,           &
                      z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin, tice,        &
                      rad_snow, rad_graupel, rad_rain, cld_min, use_ppm, mono_prof,         &
                      do_sedi_heat, sedi_transport, do_sedi_w, de_ice, icloud_f, irain_f,   &
                      mp_print, reiflag, rewmin, rewmax, reimin, reimax, rermin, rermax,    &
                      resmin, resmax, regmin, regmax, tintqs, do_hail)
    !
    class(ty_gfdlmp_config), intent(in)  :: cfg
    logical, intent(in) :: do_qa, fix_negative, fast_sat_adj, const_vi,      &
         const_vs, const_vg, const_vr, use_ccn, z_slope_liq, z_slope_ice,    &
         prog_ccn, rad_snow, rad_graupel, rad_rain, use_ppm, mono_prof,      &
         mp_print, do_hail, tintqs
    real(kind_phys), intent(in) :: mp_time, t_min, t_sub, tau_r2g, tau_smlt, &
         tau_g2r, tau_g2v, tau_v2g, tau_imlt, tau_v2l, tau_l2v, tau_i2s,     &
         tau_l2r, dw_land, dw_ocean,  vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt,&
         vi_max, vs_max, vg_max, vr_max, ql0_max, qi0_max, qs0_crt, qi0_crt, &
         qr0_crt,qc_crt, qs_mlt, qi_gen, rh_inc, rh_ins, rh_inr, rthresh,    &
         ccn_l, ccn_o, sat_adj0, c_piacr, qi_lim, ql_gen, c_paut, c_psaci,   &
         c_pgacs, c_cracw, alin, clin, tice, cld_min, do_sedi_heat,          &
         do_sedi_w, sedi_transport, de_ice, rewmin, rewmax, reimin, reimax,  &
         rermin, rermax, resmin, resmax, regmin, regmax
    integer, intent(in) :: icloud_f, irain_f, reiflag
  
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

  end subroutine setup_v3

end module module_gfdlmp_param
