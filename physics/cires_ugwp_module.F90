!
module  cires_ugwpv0_module

!
!   driver is called after pbl & before chem-parameterizations
!

    implicit none
    logical            :: module_is_initialized

    logical            :: do_physb_gwsrcs = .false.        ! control for physics-based GW-sources
    logical            :: do_rfdamp       = .false.        ! control for Rayleigh friction inside ugwp_driver

    real, parameter    :: arad=6370.e3
    real, parameter    :: pi = atan(1.0)
    real, parameter    :: pi2 = 2.*pi
    real, parameter    :: hps   = 7000.
    real, parameter    :: hpskm = hps/1000.
!
    real               :: kxw = 6.28e-3/100.               ! single horizontal wavenumber of ugwp schemes
    real, parameter    :: ricrit = 0.25
    real, parameter    :: frcrit = 0.50
    real, parameter    :: linsat = 1.00
    real, parameter    :: linsat2 = linsat*linsat
!

    integer               :: knob_ugwp_solver=1            ! 1, 2, 3, 4 - (linsat, ifs_2010, ad_gfdl, dsp_dis)
    integer, dimension(4) :: knob_ugwp_source              ! [1,1,1,0]  - (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_wvspec              !  number of waves for- (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_azdir               !   number of wave azimuths for- (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_stoch               !  1 - deterministic ; 0 - stochastic
    real,    dimension(4) :: knob_ugwp_effac               !  efficiency factors for- (oro, fronts, conv, imbf-owp]

    integer               :: knob_ugwp_doaxyz=1            ! 1 -gwdrag
    integer               :: knob_ugwp_doheat=1            ! 1 -gwheat
    integer               :: knob_ugwp_dokdis=0            ! 1 -gwmixing
    integer               :: knob_ugwp_ndx4lh = 2          ! n-number  of  "unresolved" "n*dx" for lh_gw
!
    integer  :: ugwp_azdir
    integer  :: ugwp_stoch

    integer  :: ugwp_src
    integer  :: ugwp_nws
    real     :: ugwp_effac

!
    data knob_ugwp_source / 1,0, 1, 0 /                    !  oro-conv-fjet-okw-taub_lat:      1-active 0-off
    data knob_ugwp_wvspec /1,32,32,32/                     !  number of waves for- (oro, fronts, conv, imbf-owp, taulat]
    data knob_ugwp_azdir  /2, 4, 4,4/                      !  number of wave azimuths for- (oro, fronts, conv, imbf-okwp]
    data knob_ugwp_stoch  /0, 0, 0,0/                      !  0 - deterministic ; 1 - stochastic, non-activated option
    data knob_ugwp_effac  /1.,1.,1.,1./                    !  efficiency factors for- (oro, fronts, conv, imbf-owp]
    integer  :: knob_ugwp_version = 0                      !  version control had sense under IPD in CCPP=> to SUITES
    integer  :: launch_level = 55
!
    namelist /cires_ugwp_nml/ knob_ugwp_solver, knob_ugwp_source,knob_ugwp_wvspec, knob_ugwp_azdir,  &
            knob_ugwp_stoch,  knob_ugwp_effac,knob_ugwp_doaxyz,  knob_ugwp_doheat, knob_ugwp_dokdis, &
            knob_ugwp_ndx4lh, knob_ugwp_version, launch_level

!&cires_ugwp_nml
! knob_ugwp_solver=2
! knob_ugwp_source=1,1,1,0
! knob_ugwp_wvspec=1,32,32,32
! knob_ugwp_azdir =2, 4, 4,4
! knob_ugwp_stoch =0, 0, 0,0
! knob_ugwp_effac=1, 1, 1,1
! knob_ugwp_doaxyz=1
! knob_ugwp_doheat=1
! knob_ugwp_dokdis=0
! knob_ugwp_ndx4lh=4
!/
!
! allocatable arrays, initilized during "cires_ugwp_init" &
!                     released   during "cires_ugwp_finalize"
!
   real, allocatable :: kvg(:), ktg(:), krad(:), kion(:)
   real, allocatable :: zkm(:), pmb(:)
   real, allocatable :: rfdis(:), rfdist(:)
   integer           :: levs_rf
   real              :: pa_rf, tau_rf
!
! limiters
!
   real, parameter ::  max_kdis = 400.              ! 400 m2/s
   real, parameter ::  max_axyz = 400.e-5           ! 400 m/s/day
   real, parameter ::  max_eps =  max_kdis*4.e-7    ! ~16   K/day
!
!======================================================================
   real, parameter :: F_coriol=1                    ! Coriolis effects
   real, parameter :: F_nonhyd=1                    ! Nonhydrostatic waves
   real, parameter :: F_kds   =0                    ! Eddy mixing due to GW-unstable below
   real, parameter :: iPr_ktgw =1./3., iPr_spgw=iPr_ktgw 
   real, parameter :: iPr_turb =1./3., iPr_mol =1.95
   real, parameter :: rhp1=1./hps, rhp2=0.5*rhp1, rhp4 = rhp2*rhp2
   real, parameter :: khp =  0.287*rhp1             ! R/Cp/Hp
   real, parameter :: cd_ulim = 1.0                 ! critical level precision or Lz ~ 0 ~dz of model

   contains
!
! -----------------------------------------------------------------------
!
! init  of cires_ugwp   (_init)  called from GFS_driver.F90
!
! -----------------------------------------------------------------------
   subroutine cires_ugwpv0_mod_init (me, master, nlunit, input_nml_file, logunit, &
              fn_nml, lonr, latr, levs, ak, bk, pref, dtp, cdmvgwd, cgwf,    &
              pa_rf_in, tau_rf_in)

    use  ugwpv0_oro_init,     only :  init_oro_gws_v0
    use  ugwpv0_wmsdis_init,  only :  initsolv_wmsdis_v0, ilaunch
    use  ugwpv0_lsatdis_init, only :  initsolv_lsatdis_v0
    
    implicit none

    integer,              intent (in) :: me
    integer,              intent (in) :: master
    integer,              intent (in) :: nlunit
    character (len = *),  intent (in) :: input_nml_file(:)
    integer,              intent (in) :: logunit
    character(len=64),    intent (in) :: fn_nml
    integer,              intent (in) :: lonr
    integer,              intent (in) :: levs
    integer,              intent (in) :: latr
    real,                 intent (in) :: ak(levs+1), bk(levs+1), pref
    real,                 intent (in) :: dtp
    real,                 intent (in) :: cdmvgwd(2), cgwf(2)             ! "scaling" controls for "old" GFS-GW schemes
    real,                 intent (in) :: pa_rf_in, tau_rf_in

    integer :: ios
    logical :: exists
    real    :: dxsg
    integer :: k

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml = cires_ugwp_nml)
#else
    if (me == master) print *, trim (fn_nml), ' GW-namelist file '
    
    inquire (file =trim (fn_nml) , exist = exists)

    if (.not. exists) then
       if (me == master) &
        write (6, *) 'separate ugwp :: namelist file: ', trim (fn_nml), ' does not exist'
    else
        open (unit = nlunit, file = trim(fn_nml), action = 'read', status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read   (nlunit, nml = cires_ugwp_nml)
    close  (nlunit)
#endif
    
!
    ilaunch = launch_level
    pa_rf   = pa_rf_in
    tau_rf  = tau_rf_in

! write version number and namelist to log file
    if (me == master) then
        write (logunit, *) " ================================================================== "
        write (logunit, *) "cires_ugwp_cires"
        write (logunit, nml = cires_ugwp_nml)
    endif
!
! effective kxw - resolution-aware
!
    dxsg =  pi2*arad/float(lonr) * knob_ugwp_ndx4lh
!
    allocate( kvg(levs+1),   ktg(levs+1)  )
    allocate( krad(levs+1),  kion(levs+1) )        
    allocate( zkm(levs),   pmb(levs) )
    allocate( rfdis(levs), rfdist(levs) )
!
! ak -pa  bk-dimensionless  from surf => tol_lid_pressure =0
!
    do k=1, levs
       pmb(k) = 1.e0*(ak(k) + pref*bk(k))    ! Pa -unit  Pref = 1.e5
       zkm(k) = -hpskm*alog(pmb(k)/pref)
    enddo
!
! Part-1 :init_global_gwdis
!
    call init_global_gwdis_v0(levs, zkm, pmb, kvg, ktg, krad, kion)

!
! Part-2 :init_SOURCES_gws -- only orowaves, but ugwp-v0 is based on gwdps.f of EMC
!
    
!    
! call init-solver for "stationary" multi-wave spectra and sub-grid oro
!
    call init_oro_gws_v0( knob_ugwp_wvspec(1), knob_ugwp_azdir(1),    &
         knob_ugwp_stoch(1), knob_ugwp_effac(1), lonr, kxw, cdmvgwd )
!
! call init-sources for "non-sationary" multi-wave spectra
!
    do_physb_gwsrcs=.true.

!======================
! Part-3 :init_SOLVERS
! =====================
!
! call init-solvers for "broad" non-stationary multi-wave spectra
!
    if   (knob_ugwp_solver==1) then
!
      call initsolv_lsatdis_v0(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                            knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw )
    endif
     if   (knob_ugwp_solver==2) then 

       call initsolv_wmsdis_v0(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                            knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw)
     endif


!======================
    module_is_initialized = .true.

    end subroutine cires_ugwpv0_mod_init
!      
! -----------------------------------------------------------------------
! finalize  of cires_ugwp   (_finalize)
! -----------------------------------------------------------------------

  subroutine cires_ugwpv0_mod_finalize
!
! deallocate sources/spectra & some diagnostics need to find where "deaalocate them"
! before "end" of the FV3GFS
!
    implicit none
!
!   deallocate arrays employed in V0
!
    deallocate( kvg,   ktg  )
    deallocate( krad,  kion )
    deallocate( zkm,   pmb  )
    deallocate( rfdis, rfdist)

   end subroutine cires_ugwpv0_mod_finalize
!
 end module cires_ugwpv0_module

