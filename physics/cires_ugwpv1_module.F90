
module  cires_ugwpv1_module

!
!   driver is called after pbl & before chem-parameterizations
!    it uses ugwp_common (like phys_cons) and some module-param od solvers/sources init-modules
!....................................................................................
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
!...................................................................................
!
!
    use machine,            only : kind_phys
    use  ugwp_common,       only : arad, pi, pi2, hpscale, rhp, rhp2, rh4, rhp4, khp, hpskm 
    use  ugwp_wmsdis_init,  only : ilaunch, nslope, lhmet, lzmax, lzmin, lzstar   
    use  ugwp_wmsdis_init,  only : tau_min, tamp_mpa   

    implicit none
    logical            :: module_is_initialized

    character(len=8)   :: strsolver='pss-1986'
    logical            :: do_physb_gwsrcs = .false.        ! control for physics-based GW-sources
    logical            :: do_rfdamp       = .false.        ! control for Rayleigh friction inside ugwp_driver
    integer, parameter :: idebug_gwrms=0                   ! control for diag computaions pw wind-temp GW-rms and MF fluxs   
    logical, parameter :: do_adjoro = .false.
    
    real(kind=kind_phys), parameter    ::  max_kdis = 450.                 ! 400 m2/s
    real(kind=kind_phys), parameter    ::  max_axyz = 450.e-5              ! 400 m/s/day
    real(kind=kind_phys), parameter    ::  max_eps =  max_kdis*4.e-4       ! max_kdis*BN2 
    real(kind=kind_phys), parameter    ::  maxdudt =  max_axyz
    real(kind=kind_phys), parameter    ::  maxdtdt =  max_eps*1.e-3        ! max_kdis*BN2/cp
    real(kind=kind_phys), parameter    ::  dked_min = 0.01
    real(kind=kind_phys), parameter    ::  dked_max = max_kdis    
!    
!
! Pr = Kv/Kt < 1  for upper layers; Pr_mol = 1./1.95 check it
! 
    real(kind=kind_phys), parameter    :: Pr_kvkt = 1./1.                             ! kv/kt = 1./3.    
    real(kind=kind_phys), parameter    :: Pr_kdis = Pr_kvkt/(1.+Pr_kvkt)
    
    real(kind=kind_phys), parameter    :: iPr_ktgw =1./3., iPr_spgw=iPr_ktgw 
    real(kind=kind_phys), parameter    :: iPr_turb =1./3., iPr_mol =1.95
 
    real(kind=kind_phys), parameter    :: cd_ulim = 1.0                 ! critical level precision or Lz ~ 0 ~dz of model 
    real(kind=kind_phys), parameter    :: linsat  = 1.00
    real(kind=kind_phys), parameter    :: linsat2 = linsat*linsat
        
    real(kind=kind_phys), parameter    :: ricrit = 0.25
    real(kind=kind_phys), parameter    :: frcrit = 0.50


    integer               :: knob_ugwp_version = 1    
    integer               :: knob_ugwp_solver=1              ! 1, 2, 3, 4 - (linsat, ifs_2010, ad_gfdl, dsp_dis)
    integer, dimension(4) :: knob_ugwp_source=(/1,0,1,0/)    ! [1,0,1,1]  - (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_wvspec=(/1,32,32,32/) !  number of waves for- (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_azdir=(/2,4,4,4/)     ! number of wave azimuths for-(oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_stoch=(/0,0,0,0/)     !  0 - deterministic ; 1 - stochastic
    real(kind=kind_phys),    dimension(4) :: knob_ugwp_effac=(/1.,1.,1.,1./) !  efficiency factors for- (oro, fronts, conv, imbf-owp]

    integer               :: knob_ugwp_doaxyz=1            ! 1 -gwdrag
    integer               :: knob_ugwp_doheat=1            ! 1 -gwheat
    integer               :: knob_ugwp_dokdis=0            ! 1 -gwmixing
    integer               :: knob_ugwp_ndx4lh = 2          ! n-number  of  "unresolved" "n*dx" for lh_gw
    integer               :: knob_ugwp_nslope = 1          ! spectral"growth" S-slope of GW-energy spectra mkz^S    
 
    real(kind=kind_phys)                  :: knob_ugwp_palaunch = 500.e2   ! fixed pressure layer in Pa for "launch" of NGWs
    real(kind=kind_phys)                  :: knob_ugwp_lzmax = 12.5e3      ! 12.5 km max-VERT-WL of GW-spectra
    real(kind=kind_phys)                  :: knob_ugwp_lzstar = 2.0e3      ! UTLS  mstar = 6.28/lzstar  2-2.5 km 
    real(kind=kind_phys)                  :: knob_ugwp_lzmin = 1.5e3       ! 1.5 km  min-VERT-WL of GW-spectra       
    real(kind=kind_phys)                  :: knob_ugwp_taumin = 0.25e-3
    real(kind=kind_phys)                  :: knob_ugwp_tauamp = 7.75e-3    ! range from 30.e-3 to 3.e-3 ( space-borne values)
    real(kind=kind_phys)                  :: knob_ugwp_lhmet  = 200.e3     ! 200 km
    
    logical               :: knob_ugwp_tlimb  = .true.    
    character(len=8)      :: knob_ugwp_orosolv='pss-1986'  
    
    real(kind=kind_phys)  :: kxw = 6.28/200.e3              ! single horizontal wavenumber of ugwp schemes
!
    integer  :: ugwp_azdir
    integer  :: ugwp_stoch

    integer  :: ugwp_src
    integer  :: ugwp_nws
    
    real(kind=kind_phys)     :: ugwp_effac
!
    integer  :: launch_level = 55
!
    namelist /cires_ugwp_nml/ knob_ugwp_solver, knob_ugwp_source,knob_ugwp_wvspec, knob_ugwp_azdir,      &
            knob_ugwp_stoch,  knob_ugwp_effac,knob_ugwp_doaxyz,  knob_ugwp_doheat, knob_ugwp_dokdis,     &
            knob_ugwp_ndx4lh, knob_ugwp_version, knob_ugwp_palaunch, knob_ugwp_nslope,  knob_ugwp_lzmax, &
	    knob_ugwp_lzmin,  knob_ugwp_lzstar,  knob_ugwp_lhmet, knob_ugwp_tauamp, knob_ugwp_taumin,    &
	    knob_ugwp_tlimb,  knob_ugwp_orosolv  

!
! allocatable arrays, initilized during "cires_ugwp_init" &
!                     released   during "cires_ugwp_finalize"
!
   real(kind=kind_phys), allocatable :: kvg(:), ktg(:), krad(:), kion(:)
   real(kind=kind_phys), allocatable :: zkm(:), pmb(:)
   real(kind=kind_phys), allocatable :: rfdis(:), rfdist(:)
!
! RF-not active now
!   
   integer                           :: levs_rf
   real(kind=kind_phys)              :: pa_rf, tau_rf
!
! simple modulation of tau_ngw by the total rain/precip  strength
!   
   real(kind=kind_phys), parameter    :: rain_max=8.e-5, rain_lat=41.0, rain_lim=1.e-5  
   real(kind=kind_phys), parameter    :: w_merra = 1.0,  w_nomerra = 1.-w_merra, w_rain =1.
   real(kind=kind_phys), parameter    :: mtau_rain = 1.e-3, ft_min =0.5, ft_max=2 
   real(kind=kind_phys), parameter    :: tau_ngw_max = 20.e-3                          ! 20 mPa 
   real(kind=kind_phys), parameter    :: tau_ngw_min = .20e-3                          ! .2 mPa    
!   
! Bushell et al. (2015) tau = tau_rainum (~3.8 km) x sqrt(Precip/base_rainum) 
!    
   real(kind=kind_phys), parameter    :: tau_rainum  = 0.7488e-3                       ! 0.74 mPa
   real(kind=kind_phys), parameter    :: base_rainum = 0.1e-5                          ! ~0.1 mm/day  
   real(kind=kind_phys), parameter    :: pbase_um =1./sqrt(base_rainum) * tau_rainum   !  
   integer, parameter :: metoum_rain = 0
!=================================================================
! switches that can ba activated for NGW physics include/omit
!
! rotational, non-hydrostatic and eddy-dissipative 
!   F_coriol   F_nonhyd             F_kds
!===================================================
   real(kind=kind_phys), parameter :: F_coriol=1.0                    ! Coriolis effects
   real(kind=kind_phys), parameter :: F_nonhyd=1.0                    ! Nonhydrostatic waves
   real(kind=kind_phys), parameter :: F_kds   =0.0                    ! Eddy mixing due to GW-unstable below

  
   contains
!
!-----------------------------------------------------------------------
!
! init  of cires_ugwp   (_init)  called from CCPP cap file
!
! ---------------------------------------------------------------------------------
! non-ccpp ....  
!
!  subroutine cires_ugwp_init_v1 (me, master, nlunit, logunit, jdat_gfs, fn_nml2,   &
!              lonr, latr, levs, ak, bk, pref, dtp)      
!-----------------------------------------------------------------------------------

  subroutine cires_ugwpv1_init (me, master, nlunit, logunit, jdat_gfs, con_pi,  &
              con_rerth, fn_nml2, input_nml_file, lonr, latr, levs, ak, bk,     &
              pref, dtp, errmsg, errflg)
!
!  input_nml_file ='input.nml'=fn_nml   ..... OLD_namelist and cdmvgwd(4) Corrected Bug Oct 4
!
    use  netcdf
    use  ugwp_oro_init,     only :  init_oro_gws
    use  ugwp_conv_init,    only :  init_conv_gws
    use  ugwp_fjet_init,    only :  init_fjet_gws
    use  ugwp_okw_init,     only :  init_okw_gws
    use  ugwp_lsatdis_init, only :  initsolv_lsatdis    
    
    use  ugwp_wmsdis_init,  only :  initsolv_wmsdis           
    use  ugwp_wmsdis_init,  only :  ilaunch, nslope, lhmet, lzmax, lzmin, lzstar   
    use  ugwp_wmsdis_init,  only :  tau_min, tamp_mpa 
   
    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master
    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit
    integer, intent (in) :: lonr
    integer, intent (in) :: levs
    integer, intent (in) :: latr
    integer, intent (in) :: jdat_gfs(8)
    real(kind=kind_phys),    intent (in) :: ak(levs+1), bk(levs+1), pref
    real(kind=kind_phys),    intent (in) :: dtp
!
! consider to retire them
!
    real(kind=kind_phys),    intent (in) :: con_pi, con_rerth
 
    character(len=64), intent (in) :: fn_nml2
    character(len=*),  intent (in) :: input_nml_file(:)

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

!    character,  intent (in) :: input_nml_file
!   
    integer :: ios
    logical :: exists
    
    integer :: ncid,  iernc, vid, dimid, status
    integer :: k
    integer :: ddd_ugwp,    curday_ugwp 
!    integer :: version

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml = cires_ugwp_nml)
#else
    if (me == master) print *, trim (fn_nml2), ' GW-namelist file '
    inquire (file =trim (fn_nml2) , exist = exists)
!
    if (.not. exists) then
        write(errmsg,'(3a)') 'cires_ugwpv1_init: namelist file: ', trim (fn_nml2), ' does not exist'
        errflg = 1
        return
    else
        open (unit = nlunit, file = trim(fn_nml2), action = 'read', status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read   (nlunit, nml = cires_ugwp_nml)
    close  (nlunit)
#endif

    strsolver= knob_ugwp_orosolv     
    
    curday_ugwp = jdat_gfs(1)*10000 + jdat_gfs(2)*100 +jdat_gfs(3)    
    call calendar_ugwp(jdat_gfs(1), jdat_gfs(2), jdat_gfs(3), ddd_ugwp)
        
! write version number and namelist to log file
    if (me == master) then
        write (logunit, *) " ================================================================== "
        write (logunit, *) "CCPP cires_ugwp_namelist_extended_v1"
        write (logunit, nml = cires_ugwp_nml)
        write (logunit, *) " ================================================================== "
	
        write (6, *) " ================================================================== "
        write (6, *) "CCPP cires_ugwp_namelist_extended_v1"
        write (6, nml = cires_ugwp_nml)		
        write (6, *) " ================================================================== "	
        write (6, *) "calendar_ugwp ddd_ugwp=", ddd_ugwp
        write (6, *) "calendar_ugwp curday_ugwp=", curday_ugwp
        write (6, *) " ================================================================== "	
        write (6, *) 	 ddd_ugwp, ' jdat_gfs ddd of year '	
    endif
!
! effective kxw - resolution-aware
!
! 
    kxw  = pi2/knob_ugwp_lhmet
!
!
! init global background dissipation for ugwp -> 4d-variable for fv3wam linked with pbl-vert_diff
!
!
!   allocate(fcor(latr), fcor2(latr)  )
!
    allocate( kvg(levs+1),   ktg(levs+1)  )
    allocate( krad(levs+1),  kion(levs+1) )
    allocate( zkm(levs),   pmb(levs) )
    
!
! ak -pa  bk-dimensionless  from surf => tol_lid_pressure =0
!
    
    do k=1, levs
       pmb(k) = ak(k) + pref*bk(k)    ! Pa -unit  Pref = 1.e5, pmb = Pa
       zkm(k) = -hpskm*alog(pmb(k)/pref)
    enddo   
   
!
! find ilaunch
!
   if (knob_ugwp_palaunch > 900.e2) then
     write(errmsg,'(a,e16.7)') 'cires_ugwpv1_init: unrealistic value of knob_ugwp_palaunch', knob_ugwp_palaunch
     errflg = 1
     return
   endif

   do k=levs, 1, -1
     if (pmb(k) .gt. knob_ugwp_palaunch ) exit
   enddo

      launch_level = max(k-1, 5)    ! above 5-layers from the surface
   if (me == master) then
      print *, 'cires_ugwpv1 klev_ngw =', launch_level, nint(pmb(launch_level))
   endif   
!
! Part-1 :init_global_gwdis               again "damn"-con_p
!
    call init_global_gwdis(levs, zkm, pmb, kvg, ktg, krad, kion,  me,  master)
			       
!
! Part-2 :init_SOURCES_gws
!
    
!    
! call init-solver for "stationary" multi-wave spectra and sub-grid oro
!
    call init_oro_gws( knob_ugwp_wvspec(1), knob_ugwp_azdir(1), &
         knob_ugwp_stoch(1), knob_ugwp_effac(1), lonr, kxw       )
!
! call init-sources for "non-sationary" multi-wave spectra
!
    do_physb_gwsrcs=.true.

    IF (do_physb_gwsrcs) THEN

!      if (me == master) print *, ' do_physb_gwsrcs ',  do_physb_gwsrcs, ' in cires_ugwp_init_modv1 '
      if (knob_ugwp_wvspec(4) > 0) then
! okw
        call init_okw_gws(knob_ugwp_wvspec(4), knob_ugwp_azdir(4), &
                          knob_ugwp_stoch(4), knob_ugwp_effac(4),  &
                          lonr, kxw )
!        if (me == master) print *, ' init_okw_gws '
      endif

      if (knob_ugwp_wvspec(3) > 0) then
! fronts
        call init_fjet_gws(knob_ugwp_wvspec(3), knob_ugwp_azdir(3), &
                           knob_ugwp_stoch(3), knob_ugwp_effac(3),  &
                           lonr, kxw )
!        if (me == master) print *, ' init_fjet_gws '
      endif

      if (knob_ugwp_wvspec(2) > 0) then
! conv :   con_pi, con_rerth,
        call init_conv_gws(knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                           knob_ugwp_stoch(2), knob_ugwp_effac(2),  &
                           lonr, kxw              )
!        if (me == master)   &
!           print *, ' init_convective GWs ', knob_ugwp_wvspec(2), knob_ugwp_azdir(2)

      endif

     ENDIF   !IF (do_physb_gwsrcs)
    
!======================
! Part-3 :init_SOLVERS
! =====================
!
! call init-solvers for "broad" non-stationary multi-wave spectra
!
    if   (knob_ugwp_solver==1) then
!
      kxw = pi2/lhmet
      call initsolv_lsatdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                            knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw )
    endif
     if   (knob_ugwp_solver == 2) then 
!
! re-assign from namelists 
!     
       nslope = knob_ugwp_nslope             ! the GW sprctral slope at small-m             
       lzstar = knob_ugwp_lzstar
       lzmax  = knob_ugwp_lzmax 
       lzmin  = knob_ugwp_lzmin      
       lhmet  = knob_ugwp_lhmet
       tamp_mpa =knob_ugwp_tauamp             !amplitude for GEOS-5/MERRA-2
       tau_min  =knob_ugwp_taumin             ! min of GW MF 0.25 mPa                        
       ilaunch = launch_level

       kxw = pi2/lhmet
       
       call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2),           &
           knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw, knob_ugwp_version)
	   
     endif


!======================
    module_is_initialized = .true.
    if (me == master) print *, ' CIRES_ugwpV1 is initialized ', module_is_initialized

    end subroutine cires_ugwpv1_init

    
!=============================================

     subroutine cires_ugwp_advance
!-----------------------------------------------------------------------
!   FV3-dycore and CCPP-physics has limited options to
!   add "horizontal" gradients of winds and temp-re to
! compute GW-triggers:    reserved option if  it will be funded ......
!
! the day-to-day variable sources/spectra and diagnostics for stochastic "triggers"
!     
! diagnose GW-source functions * FGF + OKWP + SGO/CONV from IAU-fields
!     and use for stochastic GWP-sources "memory"
!
!   this option is not active  due to "weak" flexibility
!     in communication between "ccpp/gfsphysics" and FV3-dycore
!        extension of State%in is needed to pass horizontal gradients
!        winds and temperature to compute "spontatneous" GW triggers
!-----------------------------------------------------------------------
      implicit none
!
! update GW sources and dissipation
!  a) physics-based GW triggers eliminated from cires_ugwpv1_triggers.F90
!  b) stochastic-based spectra and amplitudes is not considered 
!  c) use "memory" on GW-spectra from previous time-step is not considered 
!  d) update "background" dissipation of GWs as needed (option for FV3WAM)
!
     end subroutine cires_ugwp_advance
 
!      
! -----------------------------------------------------------------------
! finalize  of cires_ugwp_dealloc
! -----------------------------------------------------------------------


  subroutine cires_ugwp_dealloc
!
! deallocate sources/spectra & some diagnostics need to find where "deaalocate them"
! before "end" of the FV3GFS
!
    implicit none
!
!   deallocate arrays employed in:
!     cires_ugwp_advance / cires_ugwp_driver / cires_ugwp_init
!
   if (allocated (kvg)) deallocate (kvg)
   if (allocated (ktg)) deallocate (ktg)
   if (allocated (krad)) deallocate (krad)   
   if (allocated (kion)) deallocate (kion)   
   if (allocated (zkm)) deallocate (zkm) 
   if (allocated (pmb)) deallocate (pmb) 
!   if (allocated (ugwp_taulat))  deallocate(ugwp_taulat) 
!   if (allocated (tau_limb)) deallocate (tau_limb)   
!   if (allocated (days_limb)) deallocate(days_limb)                    
 
 
   end subroutine cires_ugwp_dealloc

!
!
    subroutine calendar_ugwp(yr, mm, dd, ddd_ugwp)
!
! computes day of year to get tau_limb forcing written with 1-day precision
!    
    implicit none
    integer, intent(in) :: yr, mm, dd
    integer :: ddd_ugwp
    
    integer ::  iw3jdn
    integer :: jd1, jddd
    jd1     = iw3jdn(yr,1,1)
    jddd    = iw3jdn(yr,mm,dd)
    ddd_ugwp = jddd-jd1+1
    
    end subroutine calendar_ugwp
    
        
    subroutine ngwflux_update(me, master, im, levs, kdt, ddd, curdate, &
	   tau_ddd, xlatd, sinlat,coslat, rain, tau_ngw)
	    
       use machine, only: kind_phys	           
       implicit none
!input
       	 
       integer, intent(in) :: me,  master                                         !, jdat(8) 
       integer, intent(in) :: im,   levs, kdt     
       integer, intent(in) :: ddd, curdate  
       
!       integer, intent(in), dimension(im)                :: j1_tau,    j2_tau
!       real(kind=kind_phys), intent(in), dimension(im)   :: ddy_j2tau, ddy_j1tau  
                            
       real(kind=kind_phys), intent(in), dimension(im)   :: xlatd, sinlat,coslat      
       real(kind=kind_phys), intent(in), dimension(im)   :: rain, tau_ddd
        
       real(kind=kind_phys), intent(inout), dimension(im) ::  tau_ngw           
!
! locals
!  

       integer                                :: i, j1, j2, k, it1, it2, iday
       real(kind=kind_phys)                   :: tem,  tx1, tx2, w1, w2, wlat, rw1, rw2
       real(kind=kind_phys)                   :: tau_rain, flat_rain, tau_3dt      
      
! 

! code below inside cires_tauamf_data.F90	 
!            it1 = 2
!         do iday=1, ntau_d2t
!	    if (float(ddd) .lt. days_limb(iday) ) then
!	    it2 = iday
!	    exit
!	    endif
!	 enddo
!	 it2 = min(it2,ntau_d2t)	 
!	 it1 = max(it2-1,1)
!	 if (it2 > ntau_d2t ) then
!	  print *, ' it1, it2, ntau_d2t ', it1, it2, ntau_d2t
!	  stop
!	 endif
!	 w2 = (float(ddd)-days_limb(it1))/(days_limb(it2)-days_limb(it1))
!	 w1 = 1.0-w2
!      do i=1, im	 
!	 j1 = j1_tau(i)
!	 j2 = j2_tau(i)
!	 tx1 = tau_limb(j1, it1)*ddy_j1tau(i)+tau_limb(j2, it1)*ddy_j2tau(i)
!	 tx2 = tau_limb(j1, it2)*ddy_j1tau(i)+tau_limb(j2, it2)*ddy_j2tau(i)	 
!	 tau_ddd(i) =  tx1*w1 + w2*tx2
!
! add modulattion by the total "rain"-strength Yudin et al.(2020-FV3GFS) and Bushell et al. (2015-UM/METO)
!
       do i=1, im        
	 tau_3dt = tau_ngw(i) * w_merra + w_nomerra *tau_ddd(i)
         	 
	 if (w_rain > 0. .and. rain(i) > 0.) then 
	 
             wlat =  abs(xlatd(i))
	     	 	 	 
	  if (wlat <= rain_lat .and. rain(i) > rain_lim) then 	 
              flat_rain = wlat/rain_lat
	      rw1 = 0.75 * flat_rain  ; rw2 = 1.-rw1
	      
	      tau_rain = tau_3dt * rw1 + rw2 * mtau_rain*min(rain_max, rain(i))/rain_lim
	      tau_rain = tau_3dt*(1.-w_rain) + w_rain* tau_rain
!
! restict variations from the "tau_ngw" without precip-impact
!	      
!      real, parameter    :: ft_min =0.5*tau_g5 < tau_rain < ft_max =2. *tau_g5	     
! 	      
	      if (tau_rain < ft_min *tau_3dt) tau_rain = ft_min *tau_3dt
	      if (tau_rain > ft_max *tau_3dt) tau_rain = ft_max *tau_3dt  
    
	      tau_3dt  =  tau_rain  
	          
	  endif 
	  if (metoum_rain == 1) then	     
	     tau_rain  =  min( sqrt(rain(i))*pbase_um, tau_ngw_max)
	     tau_3dt = max(tau_ngw_min, tau_rain)
	  endif
	 endif 
	  tau_ngw(i) = tau_3dt    
      enddo   
                        
      end subroutine ngwflux_update        
!    
 end module cires_ugwpv1_module

