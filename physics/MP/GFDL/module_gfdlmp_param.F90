! #########################################################################################
! #########################################################################################
module module_gfdlmp_param
  use machine, only: kind_phys
  implicit none

  ! #######################################################################################
  ! Data container for GFDL MP runtime configurations information (i.e. Namelist)
  ! #######################################################################################
  type ty_gfdlmp_config

     ! ####################################################################################
     ! GFDL MP Version 1 parameters.
     ! ####################################################################################
     real :: tau_g2r  = 600.       !< graupel melting to rain time scale (s)
     real :: tau_g2v  = 900.       !< graupel sublimation time scale (s)
     real :: tau_v2g  = 21600.     !< graupel deposition -- make it a slow process time scale (s)
     real :: qc_crt   = 5.0e-8     !< mini condensate mixing ratio to allow partial cloudiness
     real :: qr0_crt = 1.0e-4      !< rain to snow or graupel/hail threshold
                                   !< lfo used * mixing ratio * = 1.e-4 (hail in lfo)
     real :: c_piacr = 5.0         !< accretion: rain to ice:
     real :: c_cracw = 0.9         !< rain accretion efficiency
     real :: alin = 842.0          !< "a" in lin1983
     real :: clin = 4.8            !< "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)
     logical :: fast_sat_adj   = .false. !< has fast saturation adjustments 
     logical :: use_ccn        = .false. !< must be true when prog_ccn is false
     logical :: use_ppm        = .false. !< use ppm fall scheme
     logical :: mono_prof      = .true.  !< perform terminal fall with mono ppm scheme
     logical :: mp_print       = .false. !< cloud microphysics debugging printout
     logical :: de_ice         = .false. !< to prevent excessive build - up of cloud ice from external sources
     logical :: sedi_transport = .true.  !< transport of momentum in sedimentation
     
     ! ####################################################################################
     ! GFDL MP common (v1/v3) parameters
     ! ####################################################################################
     real :: cld_min  = 0.05       !< (v1/v3) minimum cloud fraction
     real :: tice     = 273.16     !< (DIF for v3) freezing temperature (K): ref: GFDL, GFS (DJS: V3=273.15)
     real :: t_min    = 178.       !< (v1/v3) min temp to freeze - dry all water vapor
     real :: t_sub    = 184.       !< (v1/v3) min temp for sublimation of cloud ice
     real :: mp_time  = 150.       !< (v1/v3) maximum micro - physics time step (sec)
     real :: rh_inc = 0.25         !< (v1/v3) rh increment for complete evaporation of cloud water and cloud ice
     real :: rh_inr = 0.25         !< (v1/v3) rh increment for minimum evaporation of rain
     real :: rh_ins = 0.25         !< (v1/v3) rh increment for sublimation of snow
     real :: tau_r2g  = 900.       !< (v1/v3) rain freezing during fast_sat time scale (s)
     real :: tau_smlt = 900.       !< (v1/v3) snow melting time scale (s)
     real :: tau_imlt = 600.       !< (DIF for v3) cloud ice melting time scale (s) (DJS: V3=1200.)
     real :: tau_i2s  = 1000.      !< (v1/v3) cloud ice to snow auto-conversion time scale (s)
     real :: tau_l2r  = 900.       !< (v1/v3) cloud water to rain auto-conversion time scale (s)
     real :: tau_v2l  = 150.       !< (v1/v3) water vapor to cloud water (condensation) time scale (s)
     real :: tau_l2v  = 300.       !< (v1/v3) cloud water to water vapor (evaporation) time scale (s)
     real :: dw_land  = 0.20       !< (v1/v3) value for subgrid deviation / variability over land
     real :: dw_ocean = 0.10       !< (v1/v3) base value for ocean
     real :: ccn_o = 90.           !< (v1/v3) ccn over ocean (cm^ - 3)
     real :: ccn_l = 270.          !< (v1/v3) ccn over land (cm^ - 3)
     real :: rthresh  = 10.0e-6    !< (DIF for v3) critical cloud drop radius (micro m) (DJS: v3=20.0e-6)
     real :: sat_adj0 = 0.90       !< (v1/v3) adjustment factor (0: no, 1: full) during fast_sat_adj
     real :: qi_lim   = 1.         !< (v1/v3) cloud ice limiter (0: no, 1: full, >1: extra) to prevent large ice build up
     real :: ql_mlt   = 2.0e-3     !< (v1/v3) max value of cloud water allowed from melted cloud ice
     real :: qs_mlt   = 1.0e-6     !< (v1/v3) max cloud water due to snow melt
     real :: ql_gen   = 1.0e-3     !< (v1/v3) max cloud water generation during remapping step if fast_sat_adj = .t.
     real :: qi_gen   = 1.82e-6    !< (v1/v3) max cloud ice generation during remapping step (V1 ONLY. Computed internally in V3)
     real :: ql0_max = 2.0e-3      !< (v1/v3) max cloud water value (auto converted to rain)
     real :: qi0_max = 1.0e-4      !< (v1/v3) max cloud ice value (by other sources)
     real :: qi0_crt = 1.0e-4      !< (v1/v3) cloud ice to snow autoconversion threshold (was 1.e-4);
                                   !< qi0_crt is highly dependent on horizontal resolution
     real :: qs0_crt = 1.0e-3      !< (v1/v3) snow to graupel density threshold (0.6e-3 in purdue lin scheme)
     real :: c_paut  = 0.55        !< (v1/v3) autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
     real :: c_psaci = 0.02        !< (DIF for v3) accretion: cloud ice to snow (was 0.1 in zetac) (DJS: v3=0.05)
     real :: c_pgacs = 2.0e-3      !< (DIF for v3) snow to graupel "accretion" eff. (was 0.1 in zetac) (DJS: v3=0.01)
     real :: vi_fac = 1.           !< (v1/v3) if const_vi: 1 / 3
     real :: vs_fac = 1.           !< (v1/v3) if const_vs: 1.
     real :: vg_fac = 1.           !< (v1/v3) if const_vg: 2.
     real :: vr_fac = 1.           !< (v1/v3) if const_vr: 4.
     real :: vi_max = 0.5          !< (v1/v3) max fall speed for ice
     real :: vs_max = 5.0          !< (v1/v3) max fall speed for snow
     real :: vg_max = 8.0          !< (v1/v3) max fall speed for graupel
     real :: vr_max = 12.          !< (v1/v3) max fall speed for rain
     real :: rewmin = 5.0          !< (v1/v3) minimum effective radii (liquid)
     real :: rewmax = 10.0         !< (DIF for v3) maximum effective radii (liquid) (DJS: v3=15.0)
     real :: reimin = 10.0         !< (v1/v3) minimum effective radii (ice)
     real :: reimax = 150.0        !< (v1/v3) maximum effective radii (ice)
     real :: rermin = 10.0         !< (DIF for v3) minimum effective radii (rain) (DJS: v3=15.0)
     real :: rermax = 10000.0      !< (v1/v3) maximum effective radii (rain)
     real :: resmin = 150.0        !< (v1/v3) minimum effective radii (snow)
     real :: resmax = 10000.0      !< (v1/v3) maximum effective radii (snow)
     real :: regmin = 300.0        !< (DIF for v3) minimum effective radii (graupel) (DJS: v3=150.0)
     real :: regmax = 10000.0      !< (v1/v3) maximum effective radii (graupel)
     !
     logical :: const_vi       = .false. !< (v1/v3) if .t. the constants are specified by v * _fac
     logical :: const_vs       = .false. !< (v1/v3) if .t. the constants are specified by v * _fac
     logical :: const_vg       = .false. !< (v1/v3) if .t. the constants are specified by v * _fac
     logical :: const_vr       = .false. !< (v1/v3) if .t. the constants are specified by v * _fac
     logical :: z_slope_liq    = .true.  !< (v1/v3) use linear mono slope for autocconversions
     logical :: z_slope_ice    = .false. !< (DIF for v3) use linear mono slope for autocconversions (DJS: v3=.true.)
     logical :: do_hail        = .false. !< (v1/v3) use hail parameters instead of graupel
     logical :: do_qa          = .true.  !< (v1/v3) do inline cloud fraction
     logical :: rad_snow       = .true.  !< (v1/v3) consider snow in cloud fraciton calculation
     logical :: rad_graupel    = .true.  !< (v1/v3) consider graupel in cloud fraction calculation
     logical :: rad_rain       = .true.  !< (v1/v3) consider rain in cloud fraction calculation
     logical :: do_sedi_w      = .false. !< (DIF for v3) transport of vertical motion in sedimentation (DJS: v3=.true.)
     logical :: do_sedi_heat   = .true.  !< (v1/v3) transport of heat in sedimentation
     logical :: prog_ccn       = .false. !< (v1/v3) do prognostic ccn (yi ming's method)
     logical :: fix_negative   = .false. !< (DIF for v3) fix negative water species (DJS: v3=.true.)
     logical :: tintqs         = .false. !< (v1/v3)
     !
     integer :: icloud_f = 0          !< (v1/v3) GFDL cloud scheme
                                      !< 0: subgrid variability based scheme
                                      !< 1: same as 0, but for old fvgfs implementation
                                      !< 2: binary cloud scheme
                                      !< 3: extension of 0
     integer :: irain_f = 0           !< (v1/v3) cloud water to rain auto conversion scheme
                                      !< 0: subgrid variability based scheme
                                      !< 1: no subgrid varaibility
     integer :: reiflag = 1           !< (DIF for v3) cloud ice effective radius scheme (DJS: v3=5)
                                      !< 1: Heymsfield and Mcfarquhar (1996)
                                      !< 2: Donner et al. (1997)
                                      !< 3: Fu (2007)
                                      !< 4: Kristjansson et al. (2000)
                                      !< 5: Wyser (1998)
                                      !< 6: Sun and Rikus (1999), Sun (2001)
                                      !< 7: effective radius

     ! ####################################################################################
     ! GFDL MP Version 3 parameters
     ! ####################################################################################
     real :: c_psacw = 1.0 ! cloud water to snow accretion efficiency
     real :: c_pracw = 0.8 ! cloud water to rain accretion efficiency
     real :: c_praci = 1.0 ! cloud ice to rain accretion efficiency
     real :: c_pgacw = 1.0 ! cloud water to graupel accretion efficiency
     real :: c_pgaci = 0.05 ! cloud ice to graupel accretion efficiency (was 0.1 in ZETAC)
     real :: c_pracs = 1.0 ! snow to rain accretion efficiency
     real :: c_psacr = 1.0 ! rain to snow accretion efficiency
     real :: c_pgacr = 1.0 ! rain to graupel accretion efficiency
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
     real :: vw_max = 0.01 ! maximum fall speed for cloud water (m/s)
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
     logical :: do_sedi_uv = .true. ! transport of horizontal momentum in sedimentation
     logical :: do_sedi_melt = .true. ! melt cloud ice, snow, and graupel during sedimentation
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
     !
     integer :: inflag = 1            !< ice nucleation scheme
                                      !< 1: Hong et al. (2004)
                                      !< 2: Meyers et al. (1992)
                                      !< 3: Meyers et al. (1992)
                                      !< 4: Cooper (1986)
                                      !, 5: Fletcher (1962)
     !
     integer :: igflag = 3            !< ice generation scheme
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
     procedure, pass   :: read => read_namelist_v1,read_namelist_v3
     generic,   public :: read(formatted) => read
  end type ty_gfdlmp_config

  ! Instance of DDT (Used by fv_sat_adj.F90 in FV3 dycore prior to physics initialization)
  type(ty_gfdlmp_config) :: cfg

  public :: cfg

contains
  ! #######################################################################################
  !
  ! #######################################################################################
  subroutine read_namelist_v1(self, unit, iotype, v_list, iostat, iomsg)
    class(ty_gfdlmp_config), intent(inout) :: self
    integer,      intent(in)    :: unit
    character(*), intent(in)    :: iotype
    integer,      intent(in)    :: v_list(:)
    integer,      intent(out)   :: iostat
    character(*), intent(inout) :: iomsg
 
    character(len=200) :: err_msg
    
    read (unit, *, iostat=iostat, iomsg=err_msg) self%mp_time
    if (iostat/=0 ) write(*,*) "Error 1: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%t_min
    if (iostat/=0 ) write(*,*) "Error 2: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%t_sub
    if (iostat/=0 ) write(*,*) "Error 3: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_r2g
    if (iostat/=0 ) write(*,*) "Error 4: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_smlt
    if (iostat/=0 ) write(*,*) "Error 5: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_g2r
    if (iostat/=0 ) write(*,*) "Error 6: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%dw_land
    if (iostat/=0 ) write(*,*) "Error 7: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%dw_ocean
    if (iostat/=0 ) write(*,*) "Error 8: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vi_fac
    if (iostat/=0 ) write(*,*) "Error 9: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vr_fac
    if (iostat/=0 ) write(*,*) "Error 10: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vs_fac
    if (iostat/=0 ) write(*,*) "Error 11: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vg_fac
    if (iostat/=0 ) write(*,*) "Error 12: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ql_mlt
    if (iostat/=0 ) write(*,*) "Error 13: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_qa
    if (iostat/=0 ) write(*,*) "Error 14: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%fix_negative
    if (iostat/=0 ) write(*,*) "Error 15: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vi_max
    if (iostat/=0 ) write(*,*) "Error 16: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vs_max
    if (iostat/=0 ) write(*,*) "Error 17: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vg_max
    if (iostat/=0 ) write(*,*) "Error 18: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vr_max
    if (iostat/=0 ) write(*,*) "Error 19: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qs_mlt
    if (iostat/=0 ) write(*,*) "Error 20: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qs0_crt
    if (iostat/=0 ) write(*,*) "Error 21: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qi_gen
    if (iostat/=0 ) write(*,*) "Error 22: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ql0_max
    if (iostat/=0 ) write(*,*) "Error 23: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qi0_max
    if (iostat/=0 ) write(*,*) "Error 24: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qi0_crt
    if (iostat/=0 ) write(*,*) "Error 25: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qr0_crt
    if (iostat/=0 ) write(*,*) "Error 26: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%fast_sat_adj
    if (iostat/=0 ) write(*,*) "Error 27: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rh_inc
    if (iostat/=0 ) write(*,*) "Error 28: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rh_ins
    if (iostat/=0 ) write(*,*) "Error 29: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rh_inr
    if (iostat/=0 ) write(*,*) "Error 30: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%const_vi
    if (iostat/=0 ) write(*,*) "Error 31: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%const_vs
    if (iostat/=0 ) write(*,*) "Error 32: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%const_vg
    if (iostat/=0 ) write(*,*) "Error 33: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%const_vr
    if (iostat/=0 ) write(*,*) "Error 34: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%use_ccn
    if (iostat/=0 ) write(*,*) "Error 35: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rthresh
    if (iostat/=0 ) write(*,*) "Error 36: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ccn_l
    if (iostat/=0 ) write(*,*) "Error 37: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ccn_o
    if (iostat/=0 ) write(*,*) "Error 38: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qc_crt
    if (iostat/=0 ) write(*,*) "Error 39: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_g2v
    if (iostat/=0 ) write(*,*) "Error 40: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_v2g
    if (iostat/=0 ) write(*,*) "Error 41: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%sat_adj0
    if (iostat/=0 ) write(*,*) "Error 42: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_piacr
    if (iostat/=0 ) write(*,*) "Error 43: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_imlt
    if (iostat/=0 ) write(*,*) "Error 44: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_v2l
    if (iostat/=0 ) write(*,*) "Error 45: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_l2v
    if (iostat/=0 ) write(*,*) "Error 46: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_i2s
    if (iostat/=0 ) write(*,*) "Error 47: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_l2r
    if (iostat/=0 ) write(*,*) "Error 48: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qi_lim
    if (iostat/=0 ) write(*,*) "Error 49: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ql_gen
    if (iostat/=0 ) write(*,*) "Error 50: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_paut
    if (iostat/=0 ) write(*,*) "Error 51: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_psaci
    if (iostat/=0 ) write(*,*) "Error 52: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_pgacs
    if (iostat/=0 ) write(*,*) "Error 53: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%z_slope_liq
    if (iostat/=0 ) write(*,*) "Error 54: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%z_slope_ice
    if (iostat/=0 ) write(*,*) "Error 55: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%prog_ccn
    if (iostat/=0 ) write(*,*) "Error 56: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_cracw
    if (iostat/=0 ) write(*,*) "Error 57: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%alin
    if (iostat/=0 ) write(*,*) "Error 58: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%clin
    if (iostat/=0 ) write(*,*) "Error 59: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tice
    if (iostat/=0 ) write(*,*) "Error 60: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rad_snow
    if (iostat/=0 ) write(*,*) "Error 61: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rad_graupel
    if (iostat/=0 ) write(*,*) "Error 62: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rad_rain
    if (iostat/=0 ) write(*,*) "Error 63: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%cld_min
    if (iostat/=0 ) write(*,*) "Error 64: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%use_ppm
    if (iostat/=0 ) write(*,*) "Error 65: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%mono_prof
    if (iostat/=0 ) write(*,*) "Error 66: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_sedi_heat
    if (iostat/=0 ) write(*,*) "Error 67: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%sedi_transport
    if (iostat/=0 ) write(*,*) "Error 68: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_sedi_w
    if (iostat/=0 ) write(*,*) "Error 69: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%de_ice
    if (iostat/=0 ) write(*,*) "Error 70: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%icloud_f
    if (iostat/=0 ) write(*,*) "Error 71: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%irain_f
    if (iostat/=0 ) write(*,*) "Error 72: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%mp_print
    if (iostat/=0 ) write(*,*) "Error 73: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%reiflag
    if (iostat/=0 ) write(*,*) "Error 74: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rewmin
    if (iostat/=0 ) write(*,*) "Error 75: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rewmax
    if (iostat/=0 ) write(*,*) "Error 76: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%reimin
    if (iostat/=0 ) write(*,*) "Error 77: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%reimax
    if (iostat/=0 ) write(*,*) "Error 78: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rermin
    if (iostat/=0 ) write(*,*) "Error 79: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rermax
    if (iostat/=0 ) write(*,*) "Error 80: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%resmin
    if (iostat/=0 ) write(*,*) "Error 81: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%resmax
    if (iostat/=0 ) write(*,*) "Error 82: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%regmin
    if (iostat/=0 ) write(*,*) "Error 83: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%regmax
    if (iostat/=0 ) write(*,*) "Error 84: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tintqs
    if (iostat/=0 ) write(*,*) "Error 85: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_hail
    if (iostat/=0 ) write(*,*) "Error 86: " // err_msg
    
  end subroutine read_namelist_v1

  ! #######################################################################################
  !
  ! #######################################################################################
  subroutine read_namelist_v3(self, unit, iotype, v_list, iostat, iomsg)
    class(ty_gfdlmp_config), intent(inout)  :: self
    integer,      intent(in)    :: unit
    character(*), intent(in)    :: iotype
    integer,      intent(in)    :: v_list(:)
    integer,      intent(out)   :: iostat
    character(*), intent(inout) :: iomsg
    character(len=200) :: err_msg

    read (unit, *, iostat=iostat, iomsg=err_msg) self%t_min
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%t_sub
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_r2g
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_smlt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_gmlt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%dw_land
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%dw_ocean
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vw_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vi_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vr_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vs_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vg_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ql_mlt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_qa
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%fix_negative
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vw_max
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vi_max
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vs_max
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vg_max
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vr_max
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qs_mlt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qs0_crt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ql0_max
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qi0_max
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qi0_crt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ifflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rh_inc
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rh_ins
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rh_inr
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%const_vw
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%const_vi
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%const_vs
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%const_vg
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%const_vr
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rthresh
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ccn_l
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ccn_o
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%igflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_paut
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_imlt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_v2l
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_l2v
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_i2s
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_l2r
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%qi_lim
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ql_gen
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_hail
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%inflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_psacw
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_psaci
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_pracs
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_psacr
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_pgacr
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_pgacs
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_pgacw
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_pgaci
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%z_slope_liq
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%z_slope_ice
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%prog_ccn
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_pracw
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%c_praci
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rad_snow
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rad_graupel
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rad_rain
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%cld_min
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%sedflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%sed_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_sedi_uv
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_sedi_w
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_sedi_heat
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%icloud_f
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%irain_f
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%xr_a
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%xr_b
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%xr_c
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ntimes
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_revp
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tice_mlt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_cond_timescale
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%mp_time
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%consv_checker
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%te_err
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tw_err
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%use_rhc_cevap
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%use_rhc_revap
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%tau_wbf
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_warm_rain_mp
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rh_thres
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%f_dq_p
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%f_dq_m
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_cld_adj
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rhc_cevap
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rhc_revap
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%beta
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%liq_ice_combine
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rewflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%reiflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rerflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%resflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%regflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rewmin
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rewmax
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%reimin
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%reimax
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rermin
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rermax
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%resmin
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%resmax
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%regmin
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%regmax
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%fs2g_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%fi2s_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%fi2g_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_sedi_melt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%radr_flag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rads_flag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%radg_flag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_wbf
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_psd_water_fall
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_psd_ice_fall
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0w_sig
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0i_sig
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0r_sig
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0s_sig
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0g_sig
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0h_sig
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0w_exp
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0i_exp
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0r_exp
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0g_exp
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%n0h_exp
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%muw
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%mui
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%mur
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%mus
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%mug
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%muh
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%alinw
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%alini
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%alinr
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%alins
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%aling
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%alinh
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%blinw
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%blini
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%blinr
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%blins
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%bling
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%blinh
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_new_acc_water
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_new_acc_ice
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%is_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%ss_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%gs_fac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rh_fac_evap
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rh_fac_cond
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%snow_grauple_combine
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_psd_water_num
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_psd_ice_num
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%vdiffflag
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%rewfac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%reifac
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%cp_heating
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%nconds
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_evap_timescale
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%delay_cond_evap
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%do_subgrid_proc
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%fast_fr_mlt
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    read (unit, *, iostat=iostat, iomsg=err_msg) self%fast_dep_sub
    if (iostat/=0 ) write(*,*) "ERROR Reading GFDLMP V3 Namelist:: " // err_msg
    !
  end subroutine read_namelist_v3
  !
end module module_gfdlmp_param
