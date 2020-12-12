
module  cires_ugwp_module_v1

!
!   driver is called after pbl & before chem-parameterizations
!    it uses ugwp_common (like phys_cons) and some module-param od solvers/sources init-modules
!....................................................................................
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
!...................................................................................
!
!
    use ugwp_common_v1, only : arad, pi, pi2, hpscale, rhp, rhp2, rh4 
    implicit none
    logical            :: module_is_initialized
!logical               :: do_ugwp   = .false.              ! control => true - ugwp false old gws + rayeleigh friction
    character(len=8)   :: strsolver='pss-1986'
    logical            :: do_physb_gwsrcs = .false.        ! control for physics-based GW-sources
    logical            :: do_rfdamp       = .false.        ! control for Rayleigh friction inside ugwp_driver
    integer, parameter :: idebug_gwrms=1                   ! control for diag computaions pw wind-temp GW-rms and MF fluxs   
    logical, parameter :: do_adjoro = .false.
    real, parameter    ::  max_kdis = 250.              ! 400 m2/s
    real, parameter    ::  max_axyz = 250.e-5           ! 400 m/s/day
    real, parameter    ::  max_eps =  max_kdis*4.e-7    ! ~16   K/day  max_kdis*BN2/cp 
    real, parameter    ::  maxdudt =  max_axyz
    real, parameter    ::  maxdtdt =  max_eps    
    real, parameter    ::  dked_min = 0.01
    real, parameter    ::  dked_max = max_kdis  
  
 
    real, parameter    :: hps   = hpscale
    real, parameter    :: hpskm = hps/1000.
!

    real, parameter    :: ricrit = 0.25
    real, parameter    :: frcrit = 0.50
    real, parameter    :: linsat = 1.00
    real, parameter    :: linsat2 = linsat*linsat
!
!    integer               ::    curday_ugwp                !  yyyymmdd    20150101 
!    integer               ::    ddd_ugwp                   !  ddd of year from 1-366
    
    integer               :: knob_ugwp_solver=1              ! 1, 2, 3, 4 - (linsat, ifs_2010, ad_gfdl, dsp_dis)
    integer, dimension(4) :: knob_ugwp_source=(/1,0,1,0/)    ! [1,0,1,1]  - (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_wvspec=(/1,32,32,32/) !  number of waves for- (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_azdir=(/2,4,4,4/)     !   number of wave azimuths for- (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_stoch=(/0,0,0,0/)     !  0 - deterministic ; 1 - stochastic
    real,    dimension(4) :: knob_ugwp_effac=(/1.,1.,1.,1./) !  efficiency factors for- (oro, fronts, conv, imbf-owp]

    integer               :: knob_ugwp_doaxyz=1            ! 1 -gwdrag
    integer               :: knob_ugwp_doheat=1            ! 1 -gwheat
    integer               :: knob_ugwp_dokdis=0            ! 1 -gwmixing
    integer               :: knob_ugwp_ndx4lh = 2          ! n-number  of  "unresolved" "n*dx" for lh_gw
    integer               :: knob_ugwp_nslope = 1          ! spectral"growth" S-slope of GW-energy spectra mkz^S    
 
    real                  :: knob_ugwp_palaunch = 500.e2   ! fixed pressure layer in Pa for "launch" of NGWs
    real                  :: knob_ugwp_lzmax = 12.5e3      ! 12.5 km max-VERT-WL of GW-spectra
    real                  :: knob_ugwp_lzstar = 2.0e3      ! UTLS  mstar = 6.28/lzstar  2-2.5 km 
    real                  :: knob_ugwp_lzmin = 1.5e3       ! 1.5 km  min-VERT-WL of GW-spectra       
    real                  :: knob_ugwp_taumin = 0.25e-3
    real                  :: knob_ugwp_tauamp = 7.75e-3    ! range from 30.e-3 to 3.e-3 ( space-borne values)
    real                  :: knob_ugwp_lhmet  = 200.e3     ! 200 km
!
    real                  :: kxw = pi2/200.e3              ! single horizontal wavenumber of ugwp schemes
!
! tune-ups for qbo
!	    
    real                  :: knob_ugwp_qbolev = 500.e2   ! fixed pressure layer in Pa for "launch" of conv-GWs
    real                  :: knob_ugwp_qbosin = 1.86     ! semiannual cycle of tau_qbo_src in radians 
    real                  :: knob_ugwp_qbotav = 2.285e-3 ! additional to "climate" for QBO-sg forcing 
    real                  :: knob_ugwp_qboamp = 1.191e-3 ! additional to "climate" QBO 
    real                  :: knob_ugwp_qbotau = 10.      !  relaxation time scale in days
    real                  :: knob_ugwp_qbolat = 15.      !  qbo-domain for extra-forcing
    real                  :: knob_ugwp_qbowid = 7.5      !  qbo-attenuation for extra-forcing
    character(len=8)      :: knob_ugwp_orosolv='pss-1986'

    character(len=255)    :: ugwp_qbofile =  'qbo_zmf_2009_2018.nc'
    character(len=255)    :: ugwp_taufile =  'ugwp_limb_tau.nc' 
    
!    character(len=250)     :: knob_ugwp_qbofile='qbo_zmf_2009_2018.nc'!
!    character(len=250)     :: knob_ugwp_amffile='mern_zmf_amf_12month.nc'
!     character(len=255)    :: file_limb_tab='ugwp_limb_tau.nc'

!     integer, parameter    :: ny_tab=73, nt_tab=14
!     real,    parameter    :: rdy_tab = 1./2.5,  rdd_tab = 1./30. 
!     real                  :: days_tab(nt_tab), lat_tab(ny_tab)
!     real                  :: abmf_tab(ny_tab,nt_tab)   

    integer  :: ugwp_azdir
    integer  :: ugwp_stoch

    integer  :: ugwp_src
    integer  :: ugwp_nws
    real     :: ugwp_effac

!
    integer  :: knob_ugwp_version = 0
    integer  :: launch_level = 55
!
    namelist /cires_ugwp_nml/ knob_ugwp_solver, knob_ugwp_source,knob_ugwp_wvspec, knob_ugwp_azdir,  &
            knob_ugwp_stoch,  knob_ugwp_effac,knob_ugwp_doaxyz,  knob_ugwp_doheat, knob_ugwp_dokdis, &
            knob_ugwp_ndx4lh, knob_ugwp_version, knob_ugwp_palaunch, knob_ugwp_nslope,  knob_ugwp_lzmax, &
	    knob_ugwp_lzmin,  knob_ugwp_lzstar, knob_ugwp_lhmet, knob_ugwp_tauamp, knob_ugwp_taumin,  &
	    knob_ugwp_qbolev, knob_ugwp_qbosin, knob_ugwp_qbotav, knob_ugwp_qboamp, knob_ugwp_qbotau, &
	    knob_ugwp_qbolat, knob_ugwp_qbowid, knob_ugwp_orosolv

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
! tabulated GW-sources
!
   integer           :: ntau_d1y, ntau_d2t, nqbo_d1y, nqbo_d2z, nqbo_d3t  
   real, allocatable :: ugwp_taulat(:),  ugwp_qbolat(:)
   real, allocatable :: tau_limb(:,:), days_limb(:)
   real, allocatable :: uzmf_merra(:,:,:),  days_merra(:),  pmb127(:)
   real, allocatable :: uqboe(:,:)  
   real, allocatable ::  days_y4ddd(:),  zkm127(:)
   real, allocatable ::  tau_qbo(:),   stau_qbo(:) 
   integer,allocatable ::  days_y4md(:) 
   real, allocatable   ::   vert_qbo(:)
   
!
! limiters
!
   real, parameter    :: latqbo =20., widqbo=15., taurel = 21600.        
   integer, parameter :: kz2 = 127-7, kz1= 127-49, kz5=5    !  64km - 18km   
!   

!======================================================================
   real, parameter :: F_coriol=1                    ! Coriolis effects
   real, parameter :: F_nonhyd=1                    ! Nonhydrostatic waves
   real, parameter :: F_kds   =0                    ! Eddy mixing due to GW-unstable below
   real, parameter :: iPr_ktgw =1./3., iPr_spgw=iPr_ktgw 
   real, parameter :: iPr_turb =1./3., iPr_mol =1.95
   real, parameter :: rhp1=1./hps, rh2=0.5*rhp1, rhp4 = rh2*rh2
   real, parameter :: khp =  0.287*rhp1             ! R/Cp/Hp
   real, parameter :: cd_ulim = 1.0                 ! critical level precision or Lz ~ 0 ~dz of model

   contains
!
! -----------------------------------------------------------------------
!
! init  of cires_ugwp   (_init)  called from CCPP cap file
!
! -----------------------------------------------------------------------



  subroutine cires_ugwp_init_v1 (me, master, nlunit, logunit, jdat_gfs, con_pi,  &
              con_rerth, fn_nml2, lonr, latr, levs, ak, bk, pref, dtp, cdmvgwd,  &
              cgwf, pa_rf_in, tau_rf_in, errmsg, errflg)
!
!  input_nml_file ='input.nml'=fn_nml   ..... OLD_namelist and cdmvgwd(4) Corrected Bug Oct 4
!
    use  netcdf
    use  ugwp_oro_init_v1,     only :  init_oro_gws
    use  ugwp_conv_init_v1,    only :  init_conv_gws
    use  ugwp_fjet_init_v1,    only :  init_fjet_gws
    use  ugwp_okw_init_v1,     only :  init_okw_gws
    use  ugwp_wmsdis_init_v1,  only :  initsolv_wmsdis
    
    use  ugwp_lsatdis_init_v1, only :  initsolv_lsatdis
       
    
    use  ugwp_wmsdis_init_v1,  only : ilaunch, nslope, lhmet, lzmax, lzmin, lzstar   
    use  ugwp_wmsdis_init_v1,  only : tau_min, tamp_mpa    
    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master
    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit
    integer, intent (in) :: lonr
    integer, intent (in) :: levs
    integer, intent (in) :: latr
    integer, intent (in) :: jdat_gfs(8)
    real,    intent (in) :: ak(levs+1), bk(levs+1), pref
    real,    intent (in) :: dtp
    real,    intent (in) :: cdmvgwd(2), cgwf(2)             ! "scaling" controls for "old" GFS-GW  dims(2) !!!
    real,    intent (in) :: pa_rf_in, tau_rf_in, con_pi, con_rerth
 
    character(len=64), intent (in) :: fn_nml2
    character(len=64), parameter   :: fn_nml='input.nml'

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

!    character,  intent (in) :: input_nml_file
!    integer, parameter :: logunit =  6
    integer :: ios
    logical :: exists
    real    :: dxsg
    
    integer :: ncid,  iernc, vid, dimid, status         
    integer :: k
    integer :: ddd_ugwp,    curday_ugwp 
    real, dimension(6) :: avqbo = (/0.05, 0.1, 0.25, 0.5, 0.75, 0.95/)
!
    if (me == master) print *, trim (fn_nml), ' GW-namelist file '
    inquire (file =trim (fn_nml) , exist = exists)
!
    if (.not. exists) then
       if (me == master) &
        write (6, *) 'separate ugwp :: namelist file: ', trim (fn_nml), ' does not exist'
    else
        open (unit = nlunit, file = trim(fn_nml), action = 'read', status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read   (nlunit, nml = cires_ugwp_nml)
    close  (nlunit)
!

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0


    strsolver= knob_ugwp_orosolv     
    pa_rf   = pa_rf_in
    tau_rf  = tau_rf_in
    
    curday_ugwp = jdat_gfs(1)*10000 + jdat_gfs(2)*100 +jdat_gfs(3)    
    call calendar_ugwp(jdat_gfs(1), jdat_gfs(2), jdat_gfs(3), ddd_ugwp)
        
! write version number and namelist to log file
    if (me == master) then
        write (logunit, *) " ================================================================== "
        write (logunit, *) "cires_ugwp_namelist_extended_v1"
        write (logunit, nml = cires_ugwp_nml)
        write (logunit, *) " ================================================================== "
	
        write (6, *) " ================================================================== "
        write (6, *) "cires_ugwp_namelist_extended_v1"
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
    dxsg =  pi2*arad/float(lonr) * knob_ugwp_ndx4lh
    kxw  = pi2/knob_ugwp_lhmet
!
!   kxw  =  pi2/dxsg
!
! init global background dissipation for ugwp -> 4d-variable for fv3wam linked with pbl-vert_diff
!

!   allocate(fcor(latr), fcor2(latr)  )
!
    allocate( kvg(levs+1),   ktg(levs+1)  )
    allocate( krad(levs+1),  kion(levs+1) )        
    allocate( zkm(levs),   pmb(levs) )
    allocate( rfdis(levs), rfdist(levs) )
    
    allocate (vert_qbo(levs))    
    
!
! ak -pa  bk-dimensionless  from surf => tol_lid_pressure =0
!
    
    do k=1, levs
       pmb(k) = 1.e0*(ak(k) + pref*bk(k))    ! Pa -unit  Pref = 1.e5, pmb = Pa
       zkm(k) = -hpskm*alog(pmb(k)/pref)
    enddo   
    vert_qbo(1:levs) = 0.   
     
   do k=kz1, kz2        		  
      vert_qbo(k)=1.
      if (k.le.(kz1+kz5))  vert_qbo(k)  = avqbo(k+1-kz1)
      if (k.ge.(kz2-kz5))  vert_qbo(k) =  avqbo(kz2+1-k)
      if (me == master) print *, 'vertqbo', vert_qbo(k), zkm(k) 		   
    enddo
    
!
! find ilaunch
!
    
   do k=levs, 1, -1
     if (pmb(k) .gt. knob_ugwp_palaunch ) exit
   enddo

      launch_level = max(k-1, 5)    ! above 5-layers from the surface
    
!
! Part-1 :init_global_gwdis_v1
!
    call init_global_gwdis_v1(levs, zkm, pmb, kvg, ktg, krad, kion, con_pi,     &
                              pa_rf, tau_rf,  me,  master)
    call rf_damp_init_v1  (levs, pa_rf, tau_rf, dtp, pmb, rfdis, rfdist, levs_rf)
!
! Part-2 :init_SOURCES_gws
!
    
!    
! call init-solver for "stationary" multi-wave spectra and sub-grid oro
!
    call init_oro_gws( knob_ugwp_wvspec(1), knob_ugwp_azdir(1), &
         knob_ugwp_stoch(1), knob_ugwp_effac(1), lonr, kxw, cdmvgwd )
!
! call init-sources for "non-sationary" multi-wave spectra
!
    do_physb_gwsrcs=.true.

    IF (do_physb_gwsrcs) THEN

      if (me == master) print *, ' do_physb_gwsrcs ',  do_physb_gwsrcs, ' in cires_ugwp_init_v1 '
      if (knob_ugwp_wvspec(4) > 0) then
! okw
        call init_okw_gws(knob_ugwp_wvspec(4), knob_ugwp_azdir(4), &
                          knob_ugwp_stoch(4), knob_ugwp_effac(4),  &
                          con_pi, lonr, kxw )
        if (me == master) print *, ' init_okw_gws '
      endif

      if (knob_ugwp_wvspec(3) > 0) then
! fronts
        call init_fjet_gws(knob_ugwp_wvspec(3), knob_ugwp_azdir(3), &
                           knob_ugwp_stoch(3), knob_ugwp_effac(3),  &
                           con_pi, lonr, kxw )
        if (me == master) print *, ' init_fjet_gws '
      endif

      if (knob_ugwp_wvspec(2) > 0) then
! conv
        call init_conv_gws(knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                           knob_ugwp_stoch(2), knob_ugwp_effac(2),  &
                           con_pi, con_rerth, lonr, kxw, cgwf )
        if (me == master)   &
           print *, ' init_convective GWs cgwf', knob_ugwp_wvspec(2), knob_ugwp_azdir(2)

      endif

     ENDIF   !IF (do_physb_gwsrcs)
!
!
! Tabulated sources 
!
!      goto 121
      
      iernc=NF90_OPEN(trim(ugwp_taufile), nf90_nowrite, ncid)
     
       if(iernc.ne.0) then         
          write(errmsg,'(*(a))') "Cannot open file_limb_tab data-file ",  &
                                    trim(ugwp_taufile)
          errflg = 1
          return
        else


       status = nf90_inq_dimid(ncid, "lat", DimID)
!      if (status /= nf90_noerr) call handle_err(status)
!
       status = nf90_inquire_dimension(ncid, DimID,  len =ntau_d1y )
       
       status = nf90_inq_dimid(ncid, "days", DimID)
       status = nf90_inquire_dimension(ncid, DimID,  len =ntau_d2t )
           if (me == master)  print *, ntau_d1y, ntau_d2t, ' dimd-tlimb '
        allocate (ugwp_taulat(ntau_d1y ), days_limb(ntau_d2t))
        allocate ( tau_limb (ntau_d1y, ntau_d2t ))         
                  
	iernc=nf90_inq_varid( ncid, 'DAYS', vid )
        iernc= nf90_get_var( ncid, vid, days_limb)
	iernc=nf90_inq_varid( ncid, 'LATS', vid )
        iernc= nf90_get_var( ncid, vid, ugwp_taulat)
	iernc=nf90_inq_varid( ncid, 'ABSMF', vid )
        iernc= nf90_get_var( ncid, vid, tau_limb)
			
	iernc=nf90_close(ncid)
	
	endif
!	
       iernc=NF90_OPEN(trim(ugwp_qbofile), nf90_nowrite, ncid)
     
       if(iernc.ne.0) then         
          write(errmsg,'(*(a))') "Cannot open qbofile data-file ",  &
                                    trim(ugwp_qbofile)
          errflg = 1 
          return
        else
	
       status = nf90_inq_dimid(ncid, "lat", DimID)
       status = nf90_inquire_dimension(ncid, DimID, len =nqbo_d1y )
       status = nf90_inq_dimid(ncid, "lev", DimID)
       status = nf90_inquire_dimension(ncid, DimID, len =nqbo_d2z)       
       status = nf90_inq_dimid(ncid, "days", DimID)
       status = nf90_inquire_dimension(ncid, DimID, len =nqbo_d3t )
       if (me == master) print *, nqbo_d1y, nqbo_d2z, nqbo_d3t, ' dims tauqbo '
        allocate (ugwp_qbolat(nqbo_d1y ), days_merra(nqbo_d3t) )
        allocate (zkm127(nqbo_d2z), pmb127(nqbo_d2z))
        allocate ( uzmf_merra (nqbo_d1y, nqbo_d2z, nqbo_d3t )) 
        allocate ( uqboe (nqbo_d2z, nqbo_d3t )) 	        
        allocate (days_y4ddd(nqbo_d3t), days_y4md(nqbo_d3t) )
	allocate (tau_qbo(nqbo_d3t), stau_qbo(nqbo_d3t) )  
	               
	iernc=nf90_inq_varid( ncid, 'DAYS', vid )	
        iernc= nf90_get_var( ncid, vid, days_merra)
	
	iernc=nf90_inq_varid( ncid, 'Y4MD', vid )	
        iernc= nf90_get_var( ncid, vid, days_y4md)
	
	iernc=nf90_inq_varid( ncid, 'Y4DDD', vid )	
        iernc= nf90_get_var( ncid, vid, days_y4ddd)		

	iernc=nf90_inq_varid( ncid, 'LATS', vid )
        iernc= nf90_get_var( ncid, vid, ugwp_qbolat)

	iernc=nf90_inq_varid( ncid, 'LEVS', vid )
        iernc= nf90_get_var( ncid, vid, zkm127)	

	
	iernc=nf90_inq_varid( ncid, 'UQBO', vid )
        iernc= nf90_get_var( ncid, vid, uzmf_merra)
	
	iernc=nf90_inq_varid( ncid, 'TAUQBO', vid )
        iernc= nf90_get_var( ncid, vid, tau_qbo)	
	
	iernc=nf90_inq_varid( ncid, 'STAUQBO', vid )
        iernc= nf90_get_var( ncid, vid, stau_qbo)
	iernc=nf90_inq_varid( ncid, 'UQBOE', vid )
        iernc= nf90_get_var( ncid, vid, uqboe)				
	iernc=nf90_close(ncid)	
      endif	
      
     if (me == master)	then
        print *
        print *, ' ugwp_tabulated files input '
	print *, ' ugwp_taulat ', ugwp_taulat
	print *, ' days ', days_limb
        print *, ' TAU-limb ', maxval(tau_limb)*1.e3, minval(tau_limb)*1.e3
        print *, ' TAU-qbo ', maxval(stau_qbo)*1.e3, minval(stau_qbo)*1.e3	
	print *, ' YMD-qbo ', maxval(days_y4md), minval(days_y4md)
	print *, ' YDDD-qbo ', maxval(days_y4ddd), minval(days_y4ddd)
	print *, ' uzmf_merra ',maxval(uzmf_merra), minval(uzmf_merra)
	print *, ' uEq_merra ',maxval(uqboe), minval(uqboe)
	print *	
     endif
     
!
121  continue
!    endif   ! tabulated sources SABER/HIRDLS/QBO
    
!======================
! Part-3 :init_SOLVERS
! =====================
!
! call init-solvers for "broad" non-stationary multi-wave spectra
!
    if   (knob_ugwp_solver==1) then
!
      call initsolv_lsatdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                            knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw )
    endif
     if   (knob_ugwp_solver==2) then 
!
! re-assign from namelists 
!     
       nslope = knob_ugwp_nslope             ! the GW sprctral slope at small-m             
       lzstar = knob_ugwp_lzstar
       lzmax  = knob_ugwp_lzmax 
       lzmin  = knob_ugwp_lzmin      
       lhmet  = knob_ugwp_lhmet
       tamp_mpa =knob_ugwp_tauamp           !amplitude for GEOS-5/MERRA-2
       tau_min  =knob_ugwp_taumin             ! min of GW MF 0.25 mPa                        
       ilaunch = launch_level
       kxw = pi2/lhmet
       call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
                   knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw)
     endif
!
! other solvers not yet tested for fv3gfs
!
!<    if   (knob_ugwp_solver==3) call init_dspdis
!<    if   (knob_ugwp_solver==4) call init_adodis
! 

!======================
    module_is_initialized = .true.
    if (me == master)   print *, ' CIRES-ugwp-V1 is initialized ', module_is_initialized

    end subroutine cires_ugwp_init_v1

    
!=============================================


     subroutine cires_ugwp_advance
!-----------------------------------------------------------------------
!
!   options for the day-to-day variable sources/spectra + diagnostics
!           for stochastic "triggers"
!   diagnose GW-source functions * FGF + OKWP + SGO/CONV from IAU-fields
!     or use for stochastic GWP-sources "memory"
!-----------------------------------------------------------------------
      implicit none
!
! update sources
!  a) physics-based triggers for multi-wave
!  b) stochastic-based spectra and amplitudes
!  c) use "memory" on GW-spectra from previous time-step
!  d) update "background" GW  dissipation as needed 
!
     end subroutine cires_ugwp_advance
 
!      
! -----------------------------------------------------------------------
! finalize  of cires_ugwp   (_finalize)
! -----------------------------------------------------------------------


  subroutine cires_ugwp_finalize
!
! deallocate sources/spectra & some diagnostics need to find where "deaalocate them"
! before "end" of the FV3GFS
!
    implicit none
!
!   deallocate arrays employed in:
!     cires_ugwp_advance / cires_ugwp_driver / cires_ugwp_init
!
    deallocate( kvg,   ktg  )
    deallocate( krad,  kion )
    deallocate( zkm,   pmb  )
    deallocate( rfdis, rfdist)
    deallocate(ugwp_taulat,  ugwp_qbolat)
    deallocate(tau_limb, uzmf_merra)
    deallocate(days_limb,  days_merra,  pmb127)
 
   end subroutine cires_ugwp_finalize

!
!
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
    
    
    subroutine cires_indx_ugwp (npts, me, master, dlat,j1_tau,j2_tau, w1_j1tau, w2_j2tau, &
                                j1_qbo,j2_qbo, w1_j1qbo, w2_j2qbo, dexp_latqbo  )   
	 
    implicit none
!    
! ntau_d1y, ntau_d2t, nqbo_d1y, nqbo_d2z, nqbo_d3t
! ugwp_taulat(:),  ugwp_qbolat(:), ugwp_merlat(:)
!
      integer    ::   npts, me, master
      integer, dimension(npts) ::  j1_tau,j2_tau, j1_qbo, j2_qbo
      real ,   dimension(npts) ::  dlat, w1_j1tau, w2_j2tau, w1_j1qbo, w2_j2qbo
      real ,   dimension(npts) ::  dexp_latqbo  
      real    :: widqbo2, xabs          
!
      integer i,j, j1, j2
!      
! weights for tau_limb  w1_j1tau, w2_j2tau
!
      do j=1,npts
        j2_qbo(j) = nqbo_d1y
        do i=1, nqbo_d1y
          if (dlat(j) < ugwp_qbolat(i)) then
            j2_qbo(j) = i
            exit
          endif
        enddo
	
        
        j2_qbo(j) = min(j2_qbo(j),nqbo_d1y)
        j1_qbo(j) = max(j2_qbo(j)-1,1)	
	
        if (j1_qbo(j) /= j2_qbo(j) ) then
          w2_j2qbo(j) = (dlat(j)  - ugwp_qbolat(j1_qbo(j))) &
                 / (ugwp_qbolat(j2_qbo(j))-ugwp_qbolat(j1_qbo(j)))
       	 
        else
          w2_j2qbo(j) = 1.0
        endif
          w1_j1qbo(j) = 1.0 - w2_j2qbo(j)	
	  
!
      enddo	 
!	   
! weights for tau_limb  w1_j1tau, w2_j2tau
!
      do j=1,npts
        j2_tau(j) = ntau_d1y
        do i=1,ntau_d1y
          if (dlat(j) < ugwp_taulat(i)) then
            j2_tau(j) = i
            exit
          endif
        enddo
	
      
        j2_tau(j) = min(j2_tau(j),ntau_d1y)
        j1_tau(j) = max(j2_tau(j)-1,1)	
	
        if (j1_tau(j) /= j2_tau(j) ) then
          w2_j2tau(j) = (dlat(j)  - ugwp_taulat(j1_tau(j))) &
                 / (ugwp_taulat(j2_tau(j))-ugwp_taulat(j1_tau(j)))
       	 
        else
          w2_j2tau(j) = 1.0
        endif
          w1_j1tau(j) = 1.0 -	w2_j2tau(j)	

      enddo
 		widqbo2 =1./widqbo/widqbo     
       do j=1,npts         
       	    dexp_latqbo(j) =0.		 		 
	    xabs =abs(dlat(j))
	if (xabs .le. latqbo) then
	    dexp_latqbo(j) = exp(-xabs*xabs*widqbo2)
	    if (xabs .le. 4.0 ) dexp_latqbo(j)  =1.
!	    print *, ' indx_ugwp dexp=', dexp_latqbo(j), nint(dlat(j))
        endif		 
       enddo     
      
      if (me == master ) then
222   format( 2x, 'vay-wqbo',  I4, 5(2x, F10.3))     
223   format( 2x, 'vay-limb',  I4, 5(2x, F10.3))
      print *, 'vay_indx_ugwp ',  size(dlat),  ' npts ', npts
        do j=1,npts
	  j1 =  j1_tau(j)
	  j2 =  j2_tau(j)	  
	 write(6,223) j,  ugwp_taulat(j1),  dlat(j),  ugwp_taulat(j2), w2_j2tau(j), w1_j1tau(j)
	enddo 
	print *
        do j=1,npts
	  j1 =  j1_qbo(j)
	  j2 =  j2_qbo(j)	  
	 write(6,222) j,  ugwp_qbolat(j1),  dlat(j),  ugwp_qbolat(j2), w2_j2qbo(j), w1_j1qbo(j)
	enddo          
      endif 	 
    end subroutine cires_indx_ugwp 

!    
 end module cires_ugwp_module_v1

