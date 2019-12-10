!>\file sfc_drv_ruc.F90 
!!  This file contains the RUC land surface scheme driver.

module lsm_ruc

        use machine,           only: kind_phys

        use namelist_soilveg_ruc
        use set_soilveg_ruc_mod,  only: set_soilveg_ruc
        use module_soil_pre
        use module_sf_ruclsm

        implicit none

        private

        public :: lsm_ruc_init, lsm_ruc_run, lsm_ruc_finalize

      contains

!> This subroutine calls set_soilveg_ruc() to specify vegetation and soil parameters for 
!! a given soil and land-use classification.
!! \section arg_table_lsm_ruc_init Argument Table
!! \htmlinclude lsm_ruc_init.html
!!
      subroutine lsm_ruc_init (me, isot, ivegsrc, nlunit,  &
     &                         errmsg, errflg)

      implicit none

      integer,              intent(in)  :: me, isot, ivegsrc, nlunit
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0 

      !--- initialize soil vegetation
      call set_soilveg_ruc(me, isot, ivegsrc, nlunit)

      end subroutine lsm_ruc_init

!! \section arg_table_lsm_ruc_finalize Argument Table
!! \htmlinclude lsm_ruc_finalize.html
!!
      subroutine lsm_ruc_finalize (errmsg, errflg)

      implicit none

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

  end subroutine lsm_ruc_finalize

! ===================================================================== !
!  lsm_ruc_run:                                                         !
!  RUC Surface Model - WRF4.0 version                                   ! 
!  program history log:                                                 !
!    may  2018  -- tanya smirnova                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horiz dimention and num of used pts      1    !
!     km       - integer, vertical soil layer dimension            9    !
!     ps       - real, surface pressure (pa)                       im   !
!     t1       - real, surface layer mean temperature (k)          im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     soiltyp  - integer, soil type (integer index)                im   !
!     vegtype  - integer, vegetation type (integer index)          im   !
!     sigmaf   - real, areal fractional cover of green vegetation  im   !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     dswflx   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     snet     - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     delt     - real, time interval (second)                      1    !
!     tg3      - real, deep soil temperature (k)                   im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, sfc layer 1 mean pressure (pa)              im   !
!     prslki   - real, dimensionless exner function at layer 1     im   !
!     zf       - real, height of bottom layer (m)                  im   !
!     wind       real, surface layer wind speed (m/s)              im   !
!     slopetyp - integer, class of sfc slope (integer index)       im   !
!     shdmin   - real, min fractional coverage of green veg        im   !
!     shdmax   - real, max fractnl cover of green veg (not used)   im   !
!     snoalb   - real, upper bound on max albedo over deep snow    im   !
!     sfalb    - real, mean sfc diffused sw albedo with effect          !
!                      of snow (fractional)                        im   !
!     flag_iter- logical,                                          im   !
!     flag_guess-logical,                                          im   !
!     isot     - integer, sfc soil type data source zobler or statsgo   !
!     ivegsrc  - integer, sfc veg type data source umd or igbp          !
!     smois    - real, total soil moisture content (fractional)   im,km !
!                                                                       !
!  input/outputs:                                                       !
!     weasd    - real, water equivalent accumulated snow depth (mm) im  !
!     snwdph   - real, snow depth (water equiv) over land          im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     srflag   - real, snow/rain flag for precipitation or mixed-phase
!                      precipitation fraction (depends on MP)      im   !
!     tslb     - real, soil temp (k)                              im,km !
!     sh2o     - real, liquid soil moisture                       im,km !
!     canopy   - real, canopy moisture content (mm)                im   !
!     trans    - real, total plant transpiration (m/s)             im   !
!     tsurf    - real, surface skin temperature (after iteration)  im   !
!                                                                       !
!  outputs:                                                             !
!     sncovr1  - real, snow cover over land (fractional)           im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     drain    - real, subsurface runoff (m/s)                     im   !
!     evap     - real, latent heat flux in kg kg-1 m s-1           im   !
!     runof    - real, surface runoff (m/s)                        im   !
!     evbs     - real, direct soil evaporation (m/s)               im   !
!     evcw     - real, canopy water evaporation (m/s)              im   !
!     sbsno    - real, sublimation/deposit from snopack (m/s)      im   !
!     stm      - real, total soil column moisture content (m)      im   !
!     zorl     - real, surface roughness                           im   !
!     wetness  - real, normalized soil wetness                     im   !
!                                                                       !
!  ====================    end of description    =====================  !

!> \defgroup lsm_ruc_group GSD RUC LSM Model
!! This module contains the RUC Land Surface Model developed by NOAA/GSD
!! (Smirnova et al. 2016 \cite Smirnova_2016).
#if 0
!> \section arg_table_lsm_ruc_run Argument Table
!! \htmlinclude lsm_ruc_run.html
!!
#endif
!>\section gen_lsmruc GSD RUC LSM General Algorithm
! DH* TODO - make order of arguments the same as in the metadata table
      subroutine lsm_ruc_run                                            & ! inputs
     &     ( iter, me, master, kdt, im, nlev, lsoil_ruc, lsoil, zs,     &
     &       t1, q1, qc, soiltyp, vegtype, sigmaf,                      &
     &       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
     &       prsl1, zf, wind, shdmin, shdmax, alvwf, alnwf,             &
     &       snoalb, sfalb, flag_iter, flag_guess, isot, ivegsrc, fice, &
     &       smc, stc, slc, lsm_ruc, lsm, land, islimsk,                &
     &       imp_physics, imp_physics_gfdl, imp_physics_thompson,       &
     &       smcwlt2, smcref2, do_mynnsfclay,                           &
     &       con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, con_fvirt,& !  constants
     &       weasd, snwdph, tskin, tskin_ocn,                           & !  in/outs
     &       rainnc, rainc, ice, snow, graupel,                         & ! in
     &       srflag, smois, tslb, sh2o, keepfr, smfrkeep,               & ! in/outs, on RUC levels
     &       canopy, trans, tsurf, tsnow, zorl,                         &
     &       sfcqc, sfcdew, tice, sfcqv,                                &
     &       sncovr1, qsurf, gflux, drain, evap, hflx,                  & ! outputs
     &       rhosnf, runof, runoff, srunoff,                            &
     &       chh, cmm, evbs, evcw, sbsno, stm, wetness,                 &
     &       acsnow, snowfallac,                                        &
     &       flag_init, flag_restart, errmsg, errflg                    &
     &     )

      implicit none

!  ---  constant parameters:
      real(kind=kind_phys), parameter :: rhoh2o  = 1000.0
      real(kind=kind_phys), parameter :: stbolt  = 5.670400e-8

!  ---  input:
      integer, intent(in) :: me, master
      integer, intent(in) :: im, nlev, iter, lsoil_ruc, lsoil, kdt, isot, ivegsrc
      integer, intent(in) :: lsm_ruc, lsm
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson

      real (kind=kind_phys), dimension(im,lsoil), intent(inout) :: smc,stc,slc

      real (kind=kind_phys), dimension(im), intent(in) ::               &
     &       t1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,        &
     &       ch, prsl1, wind, shdmin, shdmax,                           &
     &       snoalb, alvwf, alnwf, zf, qc, q1

      real (kind=kind_phys),  intent(in) :: delt
      real (kind=kind_phys),  intent(in) :: con_cp, con_rv, con_g,      &
                                            con_pi, con_rd,             &
                                            con_hvap, con_fvirt

      logical, dimension(im), intent(in) :: flag_iter, flag_guess, land
      integer, dimension(im), intent(in) :: islimsk ! sea/land/ice mask (=0/1/2)
      logical,                intent(in) :: do_mynnsfclay

!  ---  in/out:
      integer, dimension(im), intent(inout) :: soiltyp, vegtype
      real (kind=kind_phys), dimension(lsoil_ruc) :: dzs
      real (kind=kind_phys), dimension(lsoil_ruc), intent(inout   ) :: zs
      real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
     &       snwdph, tskin, tskin_ocn,                                  &
     &       srflag, canopy, trans, tsurf, zorl, tsnow,                 &
     &       sfcqc, sfcqv, sfcdew, fice, tice, sfalb, smcwlt2, smcref2
!  ---  in
      real (kind=kind_phys), dimension(im), intent(in) ::               &
     &       rainnc, rainc, ice, snow, graupel
!  ---  in/out:
!  --- on RUC levels
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) ::         &
     &       smois, tslb, sh2o, keepfr, smfrkeep

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(inout) :: sncovr1,   &
     &       qsurf , gflux , evap , runof , drain ,                     &
     &       runoff, srunoff, hflx, cmm, chh,                           &
     &       rhosnf, evbs, evcw, sbsno, stm, wetness,                   &
     &       acsnow, snowfallac

      logical,          intent(in)  :: flag_init, flag_restart
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rch, rho,                 &
     &       q0, qs1, weasd_old, snwdph_old,                            &
     &       tprcp_old, srflag_old, tskin_old, canopy_old,              &
     &       tsnow_old, snowfallac_old, acsnow_old, sfalb_old,          &
     &       sfcqv_old, sfcqc_old, wetness_old, zorl_old, sncovr1_old

      real (kind=kind_phys), dimension(lsoil_ruc) :: et

      real (kind=kind_phys), dimension(im,lsoil_ruc,1) :: smsoil,       &
           slsoil, stsoil, smfrsoil, keepfrsoil

      real (kind=kind_phys), dimension(im,lsoil_ruc) :: smois_old,      &
     &       tslb_old, sh2o_old, keepfr_old, smfrkeep_old

      real (kind=kind_phys),dimension (im,1,1)      ::                  &
     &     conflx2, sfcprs, sfctmp, q2, qcatm, rho2 
      real (kind=kind_phys),dimension (im,1)        ::                  &
     &     albbck, alb, chs, flhc, flqc, wet, smmax, cmc,               &
     &     dew, drip,  ec, edir, ett, lh, esnow, etp, qfx,              &
     &     acceta, ffrozp, lwdn, prcp, xland, xice,                     &
     &     graupelncv, snowncv, rainncv, raincv,                        &
     &     solnet, sfcexc,                                              &
     &     runoff1, runoff2, acrunoff,                                  &
     &     sfcems, hfx, shdfac, shdmin1d, shdmax1d,                     &
     &     sneqv, snoalb1d, snowh, snoh, tsnav,                         &
     &     snomlt, sncovr, soilw, soilm, ssoil, soilt, tbot,            &
     &     xlai, swdn, z0, znt, rhosnfr, infiltr,                       &
     &     precipfr, snfallac, acsn,                                    &
     &     qsfc, qsg, qvg, qcg, soilt1, chklowq

      real (kind=kind_phys) :: xice_threshold

      character(len=256) :: llanduse  !< Land-use dataset.  Valid values are :
                                      !! "USGS" (USGS 24/27 category dataset) and
                                      !! "MODIFIED_IGBP_MODIS_NOAH" (MODIS 20-category dataset)

      integer :: nscat, nlcat
      real (kind=kind_phys), dimension(:,:,:), allocatable :: landusef !< fractional landuse
      real (kind=kind_phys), dimension(:,:,:), allocatable :: soilctop !< fractional soil type

      integer :: nsoil, iswater, isice
      integer, dimension (1:im,1:1) :: stype, vtype
      integer :: ipr

! local
      integer :: ims,ime, its,ite, jms,jme, jts,jte, kms,kme, kts,kte
      integer :: l, k, i, j,  fractional_seaice

      logical :: flag(im)
      logical :: rdlai2d, myj, frpcpn
      logical :: debug_print
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ipr = 10

      debug_print=.false.

      chklowq = 1.

      if (isot == 1) then
        nscat = 19 ! stasgo
      else
        nscat = 9  ! zobler
      endif
      allocate(soilctop(im,nscat,1))

      if(ivegsrc == 1) then
        nlcat = 20  ! IGBP - "MODI-RUC"
      else
        nlcat = 13
      endif
      allocate(landusef(im,nlcat,1))

      if(debug_print) then
        write (0,*)'RUC LSM run'
        write (0,*)'noah soil temp',ipr,stc(ipr,:)
        write (0,*)'noah soil mois',ipr,smc(ipr,:)
        write (0,*)'soiltyp=',ipr,soiltyp(ipr)
        write (0,*)'vegtype=',ipr,vegtype(ipr)
        write (0,*)'kdt, iter =',kdt,iter
        write (0,*)'flag_init =',flag_init
        write (0,*)'flag_restart =',flag_restart
      endif
 
!> - Call rucinit() at the first time step and the first interation
!! for RUC initialization,then overwrite Noah soil fields
!! with initialized RUC soil fields for output.
      if(flag_init .and. iter==1) then
        if (debug_print) write (0,'(a,i0,a,l)') 'RUC LSM initialization, kdt = ', kdt, ', flag_restart = ', flag_restart

        call rucinit          (flag_restart, im, lsoil_ruc, lsoil, nlev, & ! in
                               isot, soiltyp, vegtype, fice,             & ! in
                               land, tskin, tskin_ocn, tg3,              & ! in
                               smc, slc, stc,                            & ! in
                               smcref2, smcwlt2,                         & ! inout
                               lsm_ruc, lsm,                             & ! in
                               zs, sh2o, smfrkeep, tslb, smois, wetness, & ! out
                               me, master, errmsg, errflg)

      endif ! flag_init=.true.,iter=1
!-- end of initialization

      ims = 1
      its = 1
      ime = 1
      ite = 1
      jms = 1
      jts = 1
      jme = 1
      jte = 1
      kms = 1
      kts = 1
      kme = 1
      kte = 1

      ! mosaic_lu=mosaic_soil=0, set in set_soilveg_ruc.F90
      ! set mosaic_lu=mosaic_soil=1 when fractional land and soil 
      ! categories available
      ! for now set fractions of differnet landuse and soil types 
      ! in the grid cell to zero

      landusef (:,:,:) = 0.0
      soilctop (:,:,:) = 0.0

      ! -- number of soil categories          
      !if(isot == 1) then
      !nscat = 19 ! stasgo
      !else
      !nscat = 9  ! zobler
      !endif
      !> - Set parameters for IGBP land-use data.
      if(ivegsrc == 1) then
        llanduse = 'MODI-RUC'  ! IGBP
        iswater = 17
        isice = 15
      else
        write(errmsg, '(a,i0)') 'Logic error in sfc_drv_ruc_run: iswater/isice not configured for ivegsrc=', ivegsrc
        errflg = 1
        return
      endif

      fractional_seaice = 1
      if ( fractional_seaice == 0 ) then
        xice_threshold = 0.5
      else if ( fractional_seaice == 1 ) then
        xice_threshold = 0.02
      endif

      nsoil = lsoil_ruc

      do i  = 1, im ! i - horizontal loop
        ! reassign smcref2 and smcwlt2 to RUC values
        if(.not. land(i)) then
          !water and sea ice
          smcref2 (i) = 1.
          smcwlt2 (i) = 0.
        else
          !land 
          smcref2 (i) = REFSMC(soiltyp(i))
          smcwlt2 (i) = WLTSMC(soiltyp(i))
        endif
      enddo

      do i  = 1, im ! i - horizontal loop
        !> - Set flag for land and ice points.
        !- 10may19 - ice points are turned off.
        flag(i) = land(i)
        if (land(i) .and. (vegtype(i)==iswater .or. (vegtype(i)==isice.and.islimsk(i)==2))) then
            !write(errmsg,'(a,i0,a,i0)') 'Logic error in sfc_drv_ruc_run: for i=', i, &
            !           ', land(i) is true but vegtype(i) is water or ice: ', vegtype(i)
            !errflg = 1
            !return
            if(flag_init .and. iter==1) then
              write(0,'(a,i0,a,i0)') 'Warning: in sfc_drv_ruc_run: for i=', i, &
                ', land(i) is true but vegtype(i) is water or ice: ', vegtype(i)
            end if
        end if
      enddo

      do i  = 1, im ! i - horizontal loop
        if (flag(i) .and. flag_guess(i)) then
          !> - Save land-related prognostic fields for guess run.
          !if(me==0 .and. i==ipr) write (0,*)'before call to RUC guess run', i
          weasd_old(i)       = weasd(i)
          snwdph_old(i)      = snwdph(i)
          tskin_old(i)       = tskin(i)
          canopy_old(i)      = canopy(i)
          !tprcp_old(i)       = tprcp(i)
          srflag_old(i)      = srflag(i)
          tsnow_old(i)       = tsnow(i)
          snowfallac_old(i)  = snowfallac(i)
          acsnow_old(i)      = acsnow(i)
          sfalb_old(i)       = sfalb(i)
          sfcqv_old(i)       = sfcqv(i)
          sfcqc_old(i)       = sfcqc(i)
          wetness_old(i)     = wetness(i)
          zorl_old(i)        = zorl(i)
          sncovr1_old(i)     = sncovr1(i)
          do k = 1, lsoil_ruc
            smois_old(i,k)  = smois(i,k)
            tslb_old(i,k)   = tslb(i,k)
            sh2o_old(i,k)   = sh2o(i,k)
            keepfr_old(i,k) = keepfr(i,k)
            smfrkeep_old(i,k) = smfrkeep(i,k)
          enddo
        endif
      enddo

!  --- ...  initialization block

      do j  = 1, 1
      do i  = 1, im ! i - horizontal loop
        if (flag_iter(i) .and. flag(i)) then
          !if(me==0 .and. i==ipr) write (0,*)'iter run', iter, i, flag_iter(i),flag_guess(i)
          evap (i)  = 0.0
          hflx (i)  = 0.0
          gflux(i)  = 0.0
          drain(i)  = 0.0
          canopy(i) = max(canopy(i), 0.0)
          sfcdew(i) = 0.0

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno(i)  = 0.0

          !local i,j arrays
          dew(i,j)      = 0.0
          soilm(i,j)    = 0.0
          smmax(i,j)    = 0.0
          hfx(i,j)      = 0.0
          qfx(i,j)      = 0.0
          lh(i,j)       = 0.0
          acsn(i,j)     = 0.0
          sfcexc(i,j)   = 0.0
          acceta(i,j)   = 0.0
          ssoil(i,j)    = 0.0
          snomlt(i,j)   = 0.0
          infiltr(i,j)  = 0.0
          runoff1(i,j)  = 0.0
          runoff2(i,j)  = 0.0
          acrunoff(i,j) = 0.0
          snfallac(i,j) = 0.0
          rhosnfr(i,j)  = 0.0
          precipfr(i,j) = 0.0

        endif
      enddo ! i=1,im
      enddo

!  --- ...  initialize atm. forcing variables

      do i  = 1, im
        if (flag_iter(i) .and. flag(i)) then
          q0(i)   = max(q1(i)/(1.-q1(i)), 1.e-8)   !* q1=specific humidity at level 1 (kg/kg)
          rho(i) = prsl1(i) / (con_rd*t1(i)*(1.0+con_fvirt*q0(i)))
          qs1(i) = rslf(prsl1(i),t1(i))  !* qs1=sat. mixing ratio at level 1 (kg/kg)
          q0 (i) = min(qs1(i), q0(i))
        endif
      enddo ! i

!> - Prepare variables to run RUC LSM: 
!!  -   1. configuration information (c):
!!\n  \a ffrozp  - fraction of frozen precipitation
!!\n  \a frpcpn  - .true. if mixed phase precipitation available
!!\n  \a 1:im - horizontal_loop_extent
!!\n  \a fice    - fraction of sea-ice in the grid cell
!!\n  \a delt    - timestep (sec) (dt should not exceed 3600 secs) 
!!\n  \a conflx2 - height (\f$m\f$) above ground of atmospheric forcing variables
!!\n  \a lsoil_ruc - number of soil layers (= 6 or 9)
!!\n  \a zs      - the depth of each soil level (\f$m\f$)

      ! Set flag for mixed phase precipitation depending on microphysics scheme.
      ! For GFDL and Thompson, srflag is fraction of frozen precip for convective+explicit precip.
      if (imp_physics==imp_physics_gfdl .or. imp_physics==imp_physics_thompson) then
        frpcpn = .true.
      else
        frpcpn = .false.
      endif

      do j  = 1, 1    ! 1:1
      do i  = 1, im   ! i - horizontal loop
        if (flag_iter(i) .and. flag(i)) then

        if (frpcpn) then
          ffrozp(i,j) = srflag(i)
        else
          ffrozp(i,j) = real(nint(srflag(i)),kind_phys)
        endif

        !tgs - for now set rdlai2d to .false., WRF has LAI maps, and RUC LSM
        !      uses rdlai2d = .true.
        rdlai2d = .false.
        !if( .not. rdlai2d) xlai = lai_data(vtype)

        conflx2(i,1,j)  = zf(i) * 2. ! factor 2. is needed to get the height of
                                     ! atm. forcing inside RUC LSM (inherited
                                     ! from WRF)

!>  -   2. forcing data (f):
!!\n  \a sfcprs  - pressure at height zf above ground (pascals)
!!\n  \a sfctmp  - air temperature (\f$K\f$) at height zf above ground
!!\n  \a q2      - pressure at height zf above ground (pascals)
!!\n  \a qcatm   - cloud water mising ration at height zf above ground (\f$kg !kg^{-1}\f$)
!!\n  \a rho2    - air density at height zf above ground (pascals)

        sfcprs(i,1,j)  = prsl1(i)
        sfctmp(i,1,j)  = t1(i)
        q2(i,1,j)      = q0(i)
        qcatm(i,1,j)   = max(0., qc(i))
        rho2(i,1,j)    = rho(i)

!!\n  \a lwdn    - lw dw radiation flux at surface (\f$W m^{-2}\f$)
!!\n  \a swdn    - sw dw radiation flux at surface (\f$W m^{-2}\f$)
!!\n  \a solnet  - net sw radiation flux (dn-up) (\f$W m^{-2}\f$)
!!\n  \a prcp    - time-step total precip (\f$kg m^{-2} \f$)
!!\n  \a raincv  - time-step convective precip (\f$kg m^{-2} \f$)
!!\n  \a rainncv - time-step non-convective precip (\f$kg m^{-2} \f$)
!!\n  \a graupelncv - time-step graupel (\f$kg m^{-2} \f$)
!!\n  \a snowncv - time-step snow (\f$kg m^{-2} \f$)
!!\n  \a precipfr - time-step precipitation in solod form (\f$kg m^{-2} \f$)
!!\n  \a qsfc    - specific humidity at surface (\f$kg kg^{-1}\f$)
!!\n  \a qvg     - water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n  \a qsg     - saturated water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n  \a qcg     - cloud water mixing ratio at surface (\f$kg kg^{-1}\f$)

        lwdn(i,j)   = dlwflx(i)         !..downward lw flux at sfc in w/m2
        swdn(i,j)   = dswsfc(i)         !..downward sw flux at sfc in w/m2
        solnet(i,j) = dswsfc(i)*(1.-sfalb(i)) !snet(i) !..net sw rad flx (dn-up) at sfc in w/m2

        ! all precip input to RUC LSM is in [mm]
        !prcp(i,j)       = rhoh2o * tprcp(i)                   ! tprcp in [m] - convective plus explicit
        !raincv(i,j)     = rhoh2o * rainc(i)                   ! total time-step convective precip
        !rainncv(i,j)    = rhoh2o * max(rain(i)-rainc(i),0.0)  ! total time-step explicit precip 
        prcp(i,j)       = rhoh2o * (rainc(i)+rainnc(i))        ! [mm] - convective plus explicit
        raincv(i,j)     = rhoh2o * rainc(i)                    ! [mm] - total time-step convective precip
        rainncv(i,j)    = rhoh2o * rainnc(i)                   ! [mm] - total time-step explicit precip 
        graupelncv(i,j) = rhoh2o * graupel(i)
        snowncv(i,j)    = rhoh2o * snow(i)
        ! ice not used
        ! precipfr(i,j)   = rainncv(i,j) * ffrozp(i,j)

        qvg(i,j)    = sfcqv(i)
        qsfc(i,j)   = sfcqv(i)/(1.+sfcqv(i))
        qsg(i,j)    = rslf(prsl1(i),tsurf(i))
        qcg(i,j)    = sfcqc(i) 

!>  -   3. canopy/soil characteristics (s):
!!\n \a vegtyp  - vegetation type (integer index)                   -> vtype
!!\n \a soiltyp - soil type (integer index)                         -> stype
!!\n \a shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
!!\n \a shdmin  - minimum areal fractional coverage of green vegetation -> shdmin1d
!!\n \a shdmax  - maximum areal fractional coverage of green vegetation -> shdmax1d
!!\n \a sfcems  -  surface emmisivity                                   -> sfcemis
!!\n \a 0.5*(alvwf + alnwf) - backround snow-free surface albedo (fraction)         -> albbck
!!\n \a snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d
!!\n \a sfalb  - surface albedo including snow effect (unitless fraction) -> alb
!!\n \a tbot    - bottom soil temperature (local yearly-mean sfc air temp)

        if(ivegsrc == 1) then   ! IGBP - MODIS
        !> - Prepare land/ice/water masks for RUC LSM
        !> - for land only
            vtype(i,j) = vegtype(i)
            stype(i,j) = soiltyp(i)
            xland(i,j) = 1.
            xice(i,j)  = 0.
        else
          write (0,*)'MODIS landuse is not available'
        endif

        ! --- units %
        shdfac(i,j) = sigmaf(i)*100.
        shdmin1d(i,j) = shdmin(i)*100.
        shdmax1d(i,j) = shdmax(i)*100.

        sfcems(i,j) = sfcemis(i)

        snoalb1d(i,j) = snoalb(i)
        albbck(i,j)   = max(0.01, 0.5 * (alvwf(i) + alnwf(i))) 
        alb(i,j)      = sfalb(i)

        tbot(i,j) = tg3(i)

!>  -   4. history (state) variables (h):
!!\n \a cmc        - canopy moisture content (\f$mm\f$)
!!\n \a soilt = tskin - ground/canopy/snowpack effective skin temperature (\f$K\f$)
!!\n \a soilt1 = snowpack temperature at the bottom of the 1st layer (\f$K\f$)
!!\n \a tslb(lsoil_ruc) - soil temp (\f$K\f$)                                    -> stsoil
!!\n \a smois(lsoil_ruc) - total soil moisture content (volumetric fraction)     -> smsoil
!!\n \a sh2o(lsoil_ruc) - unfrozen soil moisture content (volumetric fraction)   -> slsoil
!!\n \a smfrsoil(lsoil_ruc) - frozen soil moisture content (volumetric fraction) -> smfrsoil
!!\n \a keepfrflag(lsoil_ruc) - flag for frozen soil physics: 0. or 1.
!!\n \a wet        - soil moisture availability at surface
!!\n \a snowh      - actual snow depth (\f$m\f$)
!!\n \a sneqv      - liquid water-equivalent snow depth (\f$m\f$)
!!\n \a sncovr     - fraction of snow in the grid cell
!!\n \a ch         - surface exchange coefficient for heat (\f$m s^{-1}\f$)      -> chs
!!\n \a z0         - surface roughness (\f$m\f$)     -> zorl(\f$cm\f$)

        cmc(i,j) = canopy(i)            !  [mm] 
        soilt(i,j) = tsurf(i)            ! clu_q2m_iter
        ! sanity check for snow temperature tsnow
        if (tsnow(i) > 0. .and. tsnow(i) < 273.15) then
          soilt1(i,j) = tsnow(i)
        else
          soilt1(i,j) = tsurf(i)
        endif

          tsnav(i,j) = 0.5*(soilt(i,j) + soilt1(i,j)) - 273.15

        do k = 1, lsoil_ruc
          smsoil  (i,k,j) = smois(i,k)
          slsoil  (i,k,j) = sh2o(i,k)
          stsoil  (i,k,j) = tslb(i,k)
          smfrsoil(i,k,j) = smfrkeep(i,k)
          keepfrsoil(i,k,j) = keepfr(i,k)
        enddo

        if(stype(i,j) .ne. 14) then
           ! land
           if (wetness(i) > 0.) then
            wet(i,j) = wetness(i)
           else
            wet(i,j) = max(0.0001,smsoil(i,1,j)/0.3)
           endif
        else
          ! water
          wet(i,j) = 1.
        endif

        snowh(i,j) = snwdph(i) * 0.001         ! convert from mm to m
        sneqv(i,j) = weasd(i)                  ! [mm]

        snfallac(i,j) = snowfallac(i)
        acsn(i,j)     = acsnow(i)

        ! -- sanity checks on sneqv and snowh
        if (sneqv(i,j) /= 0.0 .and. snowh(i,j) == 0.0) then
          snowh(i,j) = 0.003 * sneqv(i,j) ! snow density ~300 kg m-3 
        endif

        if (snowh(i,j) /= 0.0 .and. sneqv(i,j) == 0.0) then
          sneqv(i,j) = 300. * snowh(i,j) ! snow density ~300 kg m-3 
        endif

        if (sneqv(i,j) > 0. .and. snowh(i,j) > 0.) then
          if(sneqv(i,j)/snowh(i,j) > 950.) then
            sneqv(i,j) = 300. * snowh(i,j)
          endif
        endif

        sncovr(i,j) = sncovr1(i)

        chs(i,j)    = ch(i) * wind(i) ! compute conductance 
        flhc(i,j)   = chs(i,j) * rho(i) * con_cp * (1. + 0.84*q2(i,1,j))
        flqc(i,j)   = chs(i,j) * rho(i) * wet(i,j)
        ! for output
        cmm(i)      = cm(i) * wind(i)
        chh(i)      = chs(i,j) * rho(i)
        !

        !  ---- ... outside sflx, roughness uses cm as unit
        z0(i,j)  = zorl(i)/100.
        znt(i,j) = zorl(i)/100.

        if(debug_print) then
          !if(me==0 .and. i==ipr) then
            write (0,*)'before RUC smsoil = ',smsoil(i,:,j), i,j
            write (0,*)'stsoil = ',stsoil(i,:,j), i,j
            write (0,*)'soilt = ',soilt(i,j), i,j
            write (0,*)'wet = ',wet(i,j), i,j
            write (0,*)'soilt1 = ',soilt1(i,j), i,j
            write (0,*)'delt =',delt
            write (0,*)'kdt =',kdt
            write (0,*)'flag_init =',flag_init
            write (0,*)'flag_restart =',flag_restart
            write (0,*)'nsoil =',nsoil
            write (0,*)'frpcpn =',frpcpn
            write (0,*)'zs =',zs
            write (0,*)'graupelncv(i,j) =',i,j,graupelncv(i,j)
            write (0,*)'snowncv(i,j) =',i,j,snowncv(i,j)
            write (0,*)'rainncv(i,j) =',i,j,rainncv(i,j)
            write (0,*)'raincv(i,j) =',i,j,raincv(i,j)
            write (0,*)'prcp(i,j) =',i,j,prcp(i,j)
            write (0,*)'sneqv(i,j) =',i,j,sneqv(i,j)
            write (0,*)'snowh(i,j) =',i,j,snowh(i,j)
            write (0,*)'sncovr(i,j) =',i,j,sncovr(i,j)
            write (0,*)'ffrozp(i,j) =',i,j,ffrozp(i,j)
            write (0,*)'conflx2(i,1,j) =',i,j,conflx2(i,1,j)
            write (0,*)'sfcprs(i,1,j) =',i,j,sfcprs(i,1,j)
            write (0,*)'sfctmp(i,1,j) =',i,j,sfctmp(i,1,j)
            write (0,*)'q2(i,1,j) =',i,j,q2(i,1,j)
            write (0,*)'qcatm(i,1,j) =',i,j,qcatm(i,1,j)
            write (0,*)'rho2(i,1,j) =',i,j,rho2(i,1,j)
            write (0,*)'lwdn(i,j) =',i,j,lwdn(i,j)
            write (0,*)'solnet(i,j) =',i,j,solnet(i,j)
            write (0,*)'sfcems(i,j) =',i,j,sfcems(i,j)
            write (0,*)'chklowq(i,j) =',i,j,chklowq(i,j)
            write (0,*)'chs(i,j) =',i,j,chs(i,j)
            write (0,*)'flqc(i,j) =',i,j,flqc(i,j)
            write (0,*)'flhc(i,j) =',i,j,flhc(i,j)
            write (0,*)'wet(i,j) =',i,j,wet(i,j)
            write (0,*)'cmc(i,j) =',i,j,cmc(i,j)
            write (0,*)'shdfac(i,j) =',i,j,shdfac(i,j)
            write (0,*)'alb(i,j) =',i,j,alb(i,j)
            write (0,*)'znt(i,j) =',i,j,znt(i,j)
            write (0,*)'z0(i,j) =',i,j,z0(i,j)
            write (0,*)'snoalb1d(i,j) =',i,j,snoalb1d(i,j)
            write (0,*)'alb(i,j) =',i,j,alb(i,j)
            write (0,*)'landusef(i,:,j) =',i,j,landusef(i,:,j)
            write (0,*)'soilctop(i,:,j) =',i,j,soilctop(i,:,j)
            write (0,*)'nlcat=',nlcat
            write (0,*)'nscat=',nscat
            write (0,*)'qsfc(i,j) =',i,j,qsfc(i,j)
            write (0,*)'qvg(i,j) =',i,j,qvg(i,j)
            write (0,*)'qsg(i,j) =',i,j,qsg(i,j)
            write (0,*)'qcg(i,j) =',i,j,qcg(i,j)
            write (0,*)'dew(i,j) =',i,j,dew(i,j)
            write (0,*)'soilt(i,j) =',i,j,soilt(i,j)
            write (0,*)'tskin(i) =',i,j,tskin(i)
            write (0,*)'soilt1(i,j) =',i,j,soilt1(i,j)
            write (0,*)'tsnav(i,j) =',i,j,tsnav(i,j)
            write (0,*)'tbot(i,j) =',i,j,tbot(i,j)
            write (0,*)'vtype(i,j) =',i,j,vtype(i,j)
            write (0,*)'stype(i,j) =',i,j,stype(i,j)
            write (0,*)'xland(i,j) =',i,j,xland(i,j)
            write (0,*)'xice(i,j) =',i,j,xice(i,j)
            write (0,*)'iswater=',iswater
            write (0,*)'isice=',isice
            write (0,*)'xice_threshold=',xice_threshold
            write (0,*)'con_cp=',con_cp
            write (0,*)'con_rv=',con_rv
            write (0,*)'con_rd=',con_rd
            write (0,*)'con_g=',con_g
            write (0,*)'con_pi=',con_pi
            write (0,*)'con_hvap=',con_hvap
            write (0,*)'stbolt=',stbolt
            write (0,*)'smsoil(i,:,j)=',i,j,smsoil(i,:,j)
            write (0,*)'slsoil(i,:,j)=',i,j,slsoil(i,:,j)
            write (0,*)'stsoil(i,:,j)=',i,j,stsoil(i,:,j)
            write (0,*)'smfrsoil(i,:,j)=',i,j,smfrsoil(i,:,j)
            write (0,*)'keepfrsoil(i,:,j)=',i,j,keepfrsoil(i,:,j)
            write (0,*)'soilm(i,j) =',i,j,soilm(i,j)
            write (0,*)'smmax(i,j) =',i,j,smmax(i,j)
            write (0,*)'hfx(i,j) =',i,j,hfx(i,j)
            write (0,*)'qfx(i,j) =',i,j,qfx(i,j)
            write (0,*)'lh(i,j) =',i,j,lh(i,j)
            write (0,*)'infiltr(i,j) =',i,j,infiltr(i,j)
            write (0,*)'runoff1(i,j) =',i,j,runoff1(i,j)
            write (0,*)'runoff2(i,j) =',i,j,runoff2(i,j)
            write (0,*)'acrunoff(i,j) =',i,j,acrunoff(i,j)
            write (0,*)'sfcexc(i,j) =',i,j,sfcexc(i,j)
            write (0,*)'acceta(i,j) =',i,j,acceta(i,j)
            write (0,*)'ssoil(i,j) =',i,j,ssoil(i,j)
            write (0,*)'snfallac(i,j) =',i,j,snfallac(i,j)
            write (0,*)'acsn(i,j) =',i,j,acsn(i,j)
            write (0,*)'snomlt(i,j) =',i,j,snomlt(i,j)
            write (0,*)'shdmin1d(i,j) =',i,j,shdmin1d(i,j)
            write (0,*)'shdmax1d(i,j) =',i,j,shdmax1d(i,j)
            write (0,*)'rdlai2d =',rdlai2d
          !endif
        endif

!> - Call RUC LSM lsmruc(). 
      call lsmruc( delt, flag_init, flag_restart, kdt, iter, nsoil,          &
     &          graupelncv(i,j), snowncv(i,j), rainncv(i,j), raincv(i,j),    &
     &          zs, prcp(i,j), sneqv(i,j), snowh(i,j), sncovr(i,j),          &
     &          ffrozp(i,j), frpcpn,                                         &
     &          rhosnfr(i,j), precipfr(i,j),                                 &
!  ---  inputs:
     &          conflx2(i,1,j), sfcprs(i,1,j), sfctmp(i,1,j), q2(i,1,j),     &
     &          qcatm(i,1,j), rho2(i,1,j),                                   &
     &          lwdn(i,j), solnet(i,j), sfcems(i,j), chklowq(i,j),           &
     &          chs(i,j), flqc(i,j), flhc(i,j),                              &
!  ---  input/outputs:
     &          wet(i,j), cmc(i,j), shdfac(i,j), alb(i,j), znt(i,j),         &
     &          z0(i,j), snoalb1d(i,j), albbck(i,j),                         &
!     &          z0, snoalb1d, alb, xlai,                                     &
     &          landusef(i,:,j), nlcat,                                      &
!  --- mosaic_lu and mosaic_soil are moved to the namelist
!     &          mosaic_lu, mosaic_soil,                                      &
     &          soilctop(i,:,j), nscat,                                      &
     &          qsfc(i,j), qsg(i,j), qvg(i,j), qcg(i,j), dew(i,j),           &
     &          soilt1(i,j),                                                 &
     &          tsnav(i,j), tbot(i,j), vtype(i,j), stype(i,j), xland(i,j),   &
     &          iswater, isice, xice(i,j), xice_threshold,                   &
!  ---  constants
     &          con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, stbolt,     &
!  ---  input/outputs:
     &          smsoil(i,:,j), slsoil(i,:,j), soilm(i,j), smmax(i,j),        &
     &          stsoil(i,:,j), soilt(i,j), hfx(i,j), qfx(i,j), lh(i,j),      &
     &          infiltr(i,j), runoff1(i,j), runoff2(i,j), acrunoff(i,j),     &
     &          sfcexc(i,j), acceta(i,j), ssoil(i,j),                        &
     &          snfallac(i,j), acsn(i,j), snomlt(i,j),                       &
     &          smfrsoil(i,:,j),keepfrsoil(i,:,j), .false.,                  &
     &          shdmin1d(i,j), shdmax1d(i,j), rdlai2d,                       &
     &          ims,ime, jms,jme, kms,kme,                                   &
     &          its,ite, jts,jte, kts,kte                                    )

        if(debug_print) then
          write (0,*)'after sneqv(i,j) =',i,j,sneqv(i,j)
          write (0,*)'after snowh(i,j) =',i,j,snowh(i,j)
          write (0,*)'after sncovr(i,j) =',i,j,sncovr(i,j)
          write (0,*)'after vtype(i,j) =',i,j,vtype(i,j)
          write (0,*)'after stype(i,j) =',i,j,stype(i,j)
          write (0,*)'after wet(i,j) =',i,j,wet(i,j)
          write (0,*)'after cmc(i,j) =',i,j,cmc(i,j)
          write (0,*)'after qsfc(i,j) =',i,j,qsfc(i,j)
          write (0,*)'after qvg(i,j) =',i,j,qvg(i,j)
          write (0,*)'after qsg(i,j) =',i,j,qsg(i,j)
          write (0,*)'after qcg(i,j) =',i,j,qcg(i,j)
          write (0,*)'after dew(i,j) =',i,j,dew(i,j)
          write (0,*)'after soilt(i,j) =',i,j,soilt(i,j)
          write (0,*)'after tskin(i) =',i,j,tskin(i)
          write (0,*)'after soilt1(i,j) =',i,j,soilt1(i,j)
          write (0,*)'after tsnav(i,j) =',i,j,tsnav(i,j)
          write (0,*)'after smsoil(i,:,j)=',i,j,smsoil(i,:,j)
          write (0,*)'after slsoil(i,:,j)=',i,j,slsoil(i,:,j)
          write (0,*)'after stsoil(i,:,j)=',i,j,stsoil(i,:,j)
          write (0,*)'after smfrsoil(i,:,j)=',i,j,smfrsoil(i,:,j)
          write (0,*)'after keepfrsoil(i,:,j)=',i,j,keepfrsoil(i,:,j)
          write (0,*)'after soilm(i,j) =',i,j,soilm(i,j)
          write (0,*)'after smmax(i,j) =',i,j,smmax(i,j)
          write (0,*)'after hfx(i,j) =',i,j,hfx(i,j)
          write (0,*)'after qfx(i,j) =',i,j,qfx(i,j)
          write (0,*)'after lh(i,j) =',i,j,lh(i,j)
          write (0,*)'after infiltr(i,j) =',i,j,infiltr(i,j)
          write (0,*)'after runoff1(i,j) =',i,j,runoff1(i,j)
          write (0,*)'after runoff2(i,j) =',i,j,runoff2(i,j)
          write (0,*)'after ssoil(i,j) =',i,j,ssoil(i,j)
          write (0,*)'after snfallac(i,j) =',i,j,snfallac(i,j)
          write (0,*)'after acsn(i,j) =',i,j,acsn(i,j)
          write (0,*)'after snomlt(i,j) =',i,j,snomlt(i,j)
        endif


!> - RUC LSM: prepare variables for return to parent model and unit conversion.
!>  -   6. output (o):
!!\n \a lh     - actual latent heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a hfx    - sensible heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a ssoil   - soil heat flux (\f$W m^{-2}\f$: negative if downward from surface)
!!\n \a runoff1 - surface runoff (\f$m s^{-1}\f$), not infiltrating the surface
!!\n \a runoff2 - subsurface runoff (\f$m s^{-1}\f$), drainage out bottom
!!\n \a snoh    - phase-change heat flux from snowmelt (w m-2)
!
        if(debug_print) then
          !if(me==0.and.i==ipr) then
            write (0,*)'after  RUC smsoil = ',smsoil(i,:,j), i, j
            write (0,*)'stsoil = ',stsoil(i,:,j), i,j
            write (0,*)'soilt = ',soilt(i,j), i,j
            write (0,*)'wet = ',wet(i,j), i,j
            write (0,*)'soilt1 = ',soilt1(i,j), i,j
            write (0,*)'rhosnfr = ',rhosnfr(i,j), i,j
          !endif
        endif

        ! Interstitial
        evap(i)   = qfx(i,j) / rho(i)           ! kinematic
        hflx(i)   = hfx(i,j) / (con_cp*rho(i))  ! kinematic
        gflux(i)  = ssoil(i,j)

        !evbs(i)  = edir(i,j)
        !evcw(i)  = ec(i,j)
        !trans(i) = ett(i,j)
        !sbsno(i) = esnow(i,j)
        !snohf(i) = snoh(i,j)

        sfcdew(i)  = dew(i,j)
        qsurf(i)   = qsfc(i,j)
        sncovr1(i) = sncovr(i,j)
        stm(i)     = soilm(i,j)
        tsurf(i)   = soilt(i,j)
        tice(i)    = tsurf(i)
        
        runof (i)  = runoff1(i,j) * 1000.0 ! unit conversion (from m s-1 to mm s-1 and kg m-2 s-1)
        drain (i)  = runoff2(i,j) * 1000.0 ! unit conversion (from m s-1 to mm s-1 and kg m-2 s-1)

        wetness(i) = wet(i,j)

        ! State variables
        tsnow(i)   = soilt1(i,j)
        sfcqc(i)  = qcg(i,j)
        sfcqv(i)  = qvg(i,j)
        rhosnf(i) = rhosnfr(i,j)

        ! --- ... accumulated total runoff and surface runoff
        runoff(i)  = runoff(i)  + (drain(i)+runof(i)) * delt          ! kg m-2
        srunoff(i) = srunoff(i) + runof(i) * delt                     ! kg m-2

        ! --- ... accumulated frozen precipitation (accumulation in lsmruc)
        snowfallac(i) = snfallac(i,j) ! kg m-2
        acsnow(i)     = acsn(i,j)     ! kg m-2

        !  --- ...  unit conversion (from m to mm)
        snwdph(i)  = snowh(i,j) * 1000.0

        canopy(i)  = cmc(i,j)   ! mm
        weasd(i)   = sneqv(i,j) ! mm
        sncovr1(i) = sncovr(i,j)
        !  ---- ... outside RUC LSM, roughness uses cm as unit 
        !  (update after snow's effect)
        zorl(i) = znt(i,j)*100.
        sfalb(i)= alb(i,j)

        do k = 1, lsoil_ruc
          smois(i,k)  = smsoil(i,k,j)
          sh2o(i,k)   = slsoil(i,k,j)
          tslb(i,k)   = stsoil(i,k,j)
          keepfr(i,k)   = keepfrsoil(i,k,j)
          smfrkeep(i,k) = smfrsoil(i,k,j)
        enddo

!  --- ...  do not return the following output fields to parent model
!    ec      - canopy water evaporation (m s-1)
!    edir    - direct soil evaporation (m s-1)
!    et(nsoil)-plant transpiration from a particular root layer (m s-1)
!    ett     - total plant transpiration (m s-1)
!    esnow   - sublimation from (or deposition to if <0) snowpack (m s-1)
!    drip    - through-fall of precip and/or dew in excess of canopy
!              water-holding capacity (m)
!    dew     - dewfall (or frostfall for t<273.15) (m)
!    snomlt  - snow melt (m) (water equivalent)
!    sncovr  - fractional snow cover (unitless fraction, 0-1)
!              for a given soil layer at the end of a time step
!    xlai    - leaf area index (dimensionless)
!    soilw   - available soil moisture in root zone (unitless fraction
!              between smcwlt and smcmax)
!    soilm   - total soil column moisture content (frozen+unfrozen) (m)
!    nroot   - number of root layers, a function of veg type, determined
!              in subroutine redprm.

        endif   ! end if_flag_iter_and_flag_block
      enddo   ! j
      enddo   ! i

!> - Restore land-related prognostic fields for guess run.
      do j  = 1, 1
      do i  = 1, im
        if (flag(i)) then
          if(debug_print) write (0,*)'end ',i,flag_guess(i),flag_iter(i)
          if (flag_guess(i)) then
            if(debug_print) write (0,*)'guess run'
            weasd(i)       = weasd_old(i)
            snwdph(i)      = snwdph_old(i)
            tskin(i)       = tskin_old(i)
            canopy(i)      = canopy_old(i)
            !tprcp(i)       = tprcp_old(i)
            srflag(i)      = srflag_old(i)
            tsnow(i)       = tsnow_old(i)
            snowfallac(i)  = snowfallac_old(i)
            acsnow(i)      = acsnow_old(i)
            sfalb(i)       = sfalb_old(i)
            sfcqv(i)       = sfcqv_old(i)
            sfcqc(i)       = sfcqc_old(i)
            wetness(i)     = wetness_old(i)
            zorl(i)        = zorl_old(i)
            sncovr1(i)     = sncovr1_old(i)
            do k = 1, lsoil_ruc
              smois(i,k) = smois_old(i,k)
              tslb(i,k) = tslb_old(i,k)
              sh2o(i,k) = sh2o_old(i,k)
              keepfr(i,k) = keepfr_old(i,k)
              smfrkeep(i,k) = smfrkeep_old(i,k)
            enddo
          else
            if(debug_print) write (0,*)'iter run', i,j, tskin(i),tsurf(i)
            tskin(i) = tsurf(i)
            tice (i) = tsurf(i)
          endif
        endif
      enddo  ! i
      enddo  ! j
!
      deallocate(soilctop)
      deallocate(landusef)
!
      !! Update standard (Noah LSM) soil variables for physics
      !! that require these variables (e.g. sfc_sice), independent
      !! of whether it is a land point or not
      !do i  = 1, im
      !  if (land(i)) then
      !    do k = 1, lsoil
      !      smc(i,k) = smois(i,k)
      !      slc(i,k) = sh2o(i,k)
      !      stc(i,k) = tslb(i,k)
      !    enddo
      !  endif
      !enddo
      !
      !write(0,*) "DH DEBUG: i, k, land(i), smc(i,k), slc(i,k), stc(i,k):"
      !do i  = 1, im
      !  do k = 1, lsoil
      !    write(0,'(2i5,1x,l,1x,3e20.10)'), i, k, land(i), smc(i,k), slc(i,k), stc(i,k)
      !    smc(i,k) = smois(i,k)
      !    slc(i,k) = sh2o(i,k)
      !    stc(i,k) = tslb(i,k)
      !  enddo
      !enddo

      !call sleep(20)
      !stop

      return
!...................................
      end subroutine lsm_ruc_run
!-----------------------------------

!>\ingroup lsm_ruc_group
!! This subroutine contains RUC LSM initialization.
      subroutine rucinit      (restart, im, lsoil_ruc, lsoil, nlev,   & ! in
                               isot, soiltyp, vegtype, fice,          & ! in
                               land, tsurf, tsurf_ocn,                & ! in
                               tg3, smc, slc, stc,                    & ! in
                               smcref2, smcwlt2,                      & ! inout
                               lsm_ruc, lsm,                          & ! in
                               zs, sh2o, smfrkeep, tslb, smois,       & ! out
                               wetness, me, master, errmsg, errflg)

      implicit none

      logical,                                 intent(in   ) :: restart
      integer,                                 intent(in   ) :: lsm
      integer,                                 intent(in   ) :: lsm_ruc
      integer,                                 intent(in   ) :: isot
      integer,                                 intent(in   ) :: im, nlev
      integer,                                 intent(in   ) :: lsoil_ruc
      integer,                                 intent(in   ) :: lsoil
      logical,               dimension(im),    intent(in   ) :: land
      real (kind=kind_phys), dimension(im),    intent(in   ) :: tsurf, tsurf_ocn
      real (kind=kind_phys), dimension(im),    intent(inout) :: smcref2
      real (kind=kind_phys), dimension(im),    intent(inout) :: smcwlt2
      real (kind=kind_phys), dimension(im),    intent(in   ) :: tg3
      real (kind=kind_phys), dimension(im,lsoil),  intent(in   ) :: smc !  Noah
      real (kind=kind_phys), dimension(im,lsoil),  intent(in   ) :: stc !  Noah
      real (kind=kind_phys), dimension(im,lsoil),  intent(in   ) :: slc !  Noah

      integer,               dimension(im),    intent(inout) :: soiltyp
      integer,               dimension(im),    intent(inout) :: vegtype
      real (kind=kind_phys), dimension(im),    intent(inout) :: wetness
      real (kind=kind_phys), dimension(im),    intent(inout) :: fice
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: smois! ruc
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: tslb ! ruc
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: sh2o ! ruc
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: smfrkeep ! ruc

      real (kind=kind_phys), dimension(1:lsoil_ruc), intent (out)  :: zs

      integer,          intent(in ) :: me
      integer,          intent(in ) :: master
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!> local
      logical :: debug_print
      logical :: smadj ! for soil mosture adjustment
      logical :: swi_init ! for initialization in terms of SWI (soil wetness index)

      integer :: flag_soil_layers, flag_soil_levels, flag_sst
      real (kind=kind_phys),    dimension(1:lsoil_ruc) :: factorsm

      integer , dimension( 1:im , 1:1 )      :: ivgtyp
      integer , dimension( 1:im , 1:1)       :: isltyp
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: mavail
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: xice
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: sst
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: landmask
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: tsk
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: tbot
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: smtotn
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: smtotr
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: dumsm 
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: dumt
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: smfr
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: soilm
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: soiltemp
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: soilh2o

      real (kind=kind_phys) :: st_input(1:im,1:lsoil_ruc*3,1:1)
      real (kind=kind_phys) :: sm_input(1:im,1:lsoil_ruc*3,1:1)

      integer               :: ids,ide, jds,jde, kds,kde, &
                               ims,ime, jms,jme, kms,kme, &
                               its,ite, jts,jte, kts,kte, &
                               i, j, k, l, num_soil_layers, ipr

      real(kind=kind_phys), dimension(1:lsoil_ruc) :: zs2, dzs 
      integer,              dimension(1:lsoil)  :: st_levels_input ! 4 - for Noah lsm
      integer,              dimension(1:lsoil)  :: sm_levels_input ! 4 - for Noah lsm

      integer :: ii,jj
      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      debug_print = .false.

      if (lsm/=lsm_ruc) then
        write(errmsg,'(a,i0,a,i0)')                                &
              'ERROR in lsm_ruc_init: namelist variable lsm=',     &
              lsm, ' incompatible with RUC LSM, please set to ', lsm_ruc
        errflg = 1
        return
      else if (debug_print) then
        write (0,*) 'Start of RUC LSM initialization'
        write (0,*)'lsoil, lsoil_ruc =',lsoil, lsoil_ruc
      endif

      ipr = 10

      ! Set internal dimensions
      ids = 1
      ims = 1
      its = 1
      ide = im
      ime = im
      ite = im
      jds = 1
      jms = 1
      jts = 1
      jde = 1 
      jme = 1
      jte = 1
      kds = 1
      kms = 1
      kts = 1
      kde = nlev
      kme = nlev
      kte = nlev

      ! Initialize the RUC soil levels, needed for cold starts and warm starts
      CALL init_soil_depth_3 ( zs , dzs , lsoil_ruc )

      !! Check if RUC soil data (tslb, ...) is provided or not
      !if (minval(tslb)==maxval(tslb)) then
      ! For restart runs, can assume that RUC soul data is provided
      if (.not.restart) then

        flag_sst = 0

        num_soil_layers =  lsoil ! 4 - for Noah lsm

        if( lsoil_ruc == lsoil) then
          ! RUC LSM input
          smadj = .false.
          swi_init = .false.
          flag_soil_layers = 0  ! =1 for input from the Noah LSM
          flag_soil_levels = 1  ! =1 for input from RUC LSM
        else
          ! for Noah input set smadj and swi_init to .true.
          smadj = .true.
          swi_init = .true.
          flag_soil_layers = 1  ! =1 for input from the Noah LSM
          flag_soil_levels = 0  ! =1 for input from RUC LSM
        endif

        if(lsoil == 4 ) then ! for Noah input
          st_levels_input = (/ 5, 25, 70, 150/)    ! Noah centers of soil layers
          sm_levels_input = (/ 5, 25, 70, 150/)    ! Noah centers of soil layers
        elseif(lsoil /= lsoil_ruc) then
          write(errmsg,'(a,i0,a)')                                   &
                'WARNING in lsm_ruc_init: non-Noah and non-RUC input, lsoil=', lsoil
          errflg = 1
          return
        endif

      else

        ! For RUC restart data, return here
        return

      endif

      if(debug_print) then
         write (0,*)'smc(ipr,:) ==', ipr, smc(ipr,:)
         write (0,*)'stc(ipr,:) ==', ipr, stc(ipr,:)
         write (0,*)'vegtype(ipr) ==', ipr, vegtype(ipr)
         write (0,*)'soiltyp(ipr) ==', ipr, soiltyp(ipr)
         write (0,*)'its,ite,jts,jte ',its,ite,jts,jte 
      endif


        do j=jts,jte !
        do i=its,ite ! i = horizontal loop

         ! land only version
          if (land(i)) then
            tsk(i,j) = tsurf(i)
            sst(i,j) = tsurf_ocn(i)
            tbot(i,j)= tg3(i)
            ivgtyp(i,j)=vegtype(i)
            isltyp(i,j)=soiltyp(i)
            landmask(i,j)=1.
            xice(i,j)=0.
         else
            landmask(i,j)=0.
         endif ! land(i)

        enddo
        enddo

      if ( flag_soil_layers == 1 ) then
      ! Noah lsm input
        do j=jts,jte !
        do i=its,ite ! i = horizontal loop

         if (land(i)) then

          st_input(i,1,j)=tsk(i,j)
          sm_input(i,1,j)=0.

          !--- initialize smcwlt2 and smcref2 with Noah values
          smcref2 (i) = REFSMCnoah(soiltyp(i))
          smcwlt2 (i) = WLTSMCnoah(soiltyp(i))

          do k=1,lsoil
             st_input(i,k+1,j)=stc(i,k)
             ! convert volumetric soil moisture to SWI (soil wetness index)
             if(swi_init) then
               sm_input(i,k+1,j)=min(1.,max(0.,(smc(i,k) - smcwlt2(i))/  &
                                 (smcref2(i) - smcwlt2(i))))
             else
               sm_input(i,k+1,j)=smc(i,k)
             endif
          enddo
          do k=lsoil+2,lsoil_ruc * 3
             st_input(i,k,j)=0.
             sm_input(i,k,j)=0.
          enddo

         endif ! land(i)

        enddo ! i - horizontal loop
        enddo ! jme

        if(debug_print) then
          write (0,*)'st_input=',ipr, st_input(ipr,:,1)
          write (0,*)'sm_input=',ipr, sm_input(ipr,:,1)
        endif

        CALL init_soil_3_real ( tsk , tbot , dumsm , dumt ,             &
                                st_input , sm_input , landmask , sst ,  &
                                zs , dzs ,                              &
                                st_levels_input, sm_levels_input,       &
                                lsoil_ruc , num_soil_layers,            &
                                num_soil_layers,                        &
                                lsoil_ruc * 3 , lsoil_ruc * 3 ,         &
                                flag_sst,                               &
                                flag_soil_layers , flag_soil_levels ,   &
                                ids , ide , jds , jde , kds , kde ,     &
                                ims , ime , jms , jme , kms , kme ,     &
                                its , ite , jts , jte , kts , kte )

        do j=jts,jte
        do i=its,ite
         if (land(i)) then
          do k=1,lsoil_ruc
           ! convert from SWI to RUC volumetric soil moisture
           if(swi_init) then
               soilm(i,k,j)= dumsm(i,k,j) *                             &
                 (refsmc(isltyp(i,j))-drysmc(isltyp(i,j)))              &
                 + drysmc(isltyp(i,j))
           else
             soilm(i,k,j)= dumsm(i,k,j)
           endif
             soiltemp(i,k,j) = dumt(i,k,j)
          enddo
         endif ! land(i)
        enddo
        enddo

        if(debug_print) then
          write (0,*)'tsk(i,j),tbot(i,j),sst(i,j),landmask(i,j)' &
                  ,ipr,1,tsk(ipr,1),tbot(ipr,1),sst(ipr,1),landmask(ipr,1)
          write (0,*)'tsurf(ipr)=',ipr,tsurf(ipr)
          write (0,*)'stc(ipr)=',ipr,stc(ipr,:)
          write (0,*)'smc(ipr)=',ipr,smc(ipr,:)
          write (0,*)'soilt(1,:,ipr)',ipr,soiltemp(ipr,:,1)
          write (0,*)'soilm(1,:,ipr)',ipr,soilm(ipr,:,1)
        endif ! debug_print

        ! smadj should be true when the Noah LSM is used to initialize RUC
        if( smadj ) then
        ! With other LSMs as input, or when RUC soil moisture is cycled, it
        ! should be set to .false.

          do j=jts,jte
          do i=its,ite

           if (land(i)) then

            ! initialize factor
            do k=1,lsoil_ruc
               factorsm(k)=1.
            enddo
          
            ! RUC soil moisture bucket
            smtotr(i,j)=0.
            do k=1,lsoil_ruc -1
              smtotr(i,j)=smtotr(i,j) + soilm(i,k,j) *dzs(k)
            enddo
            ! Noah soil moisture bucket 
            smtotn(i,j)=smc(i,1)*0.1 + smc(i,2)*0.2 + smc(i,3)*0.7 + smc(i,4)*1.
            
            if(debug_print) then
              if(i==ipr) then
              write (0,*)'from Noah to RUC: RUC bucket and Noah bucket at',    &
                       i,j,smtotr(i,j),smtotn(i,j)
              write (0,*)'before smois=',i,j,soilm(i,:,j)
              endif
            endif
          
            ! RUC soil moisture correction to match Noah soil moisture bucket
            do k=1,lsoil_ruc-1
              soilm(i,k,j) = max(0.02,soilm(i,k,j)*smtotn(i,j)/(0.9*smtotr(i,j)))
            enddo
          
            if( soilm(i,2,j) > soilm(i,1,j) .and. soilm(i,3,j) > soilm(i,2,j)) then
            ! typical for daytime, no recent precip
              factorsm(1) = 0.75
              factorsm(2) = 0.8
              factorsm(3) = 0.85
              factorsm(4) = 0.9
              factorsm(5) = 0.95
            endif
            do k=1,lsoil_ruc
               soilm(i,k,j) = factorsm(k) * soilm(i,k,j)
            enddo
            if(debug_print) then
               if(i==ipr) write (0,*)'after smois=',i,j,soilm(i,:,j)
            endif
               smtotr(i,j) = 0.
            do k=1,lsoil_ruc - 1
               smtotr(i,j)=smtotr(i,j) + soilm(i,k,j) *dzs(k)
            enddo
            if(debug_print) then
                if(i==ipr) write (0,*)'after correction: RUC bucket and Noah bucket at',  &
                         i,j,smtotr(i,j),smtotn(i,j)
            endif

           endif ! land(i)

          enddo
          enddo

        endif ! smadj==.true.

      elseif (flag_soil_layers==0) then
      !  RUC LSM input
        if(debug_print) write (0,*)' RUC LSM input for soil variables'
        do j=jts,jte
        do i=its,ite
          do k=1,lsoil_ruc
             soilm(i,k,j)    = smc(i,k)
             soiltemp(i,k,j) = stc(i,k)
          enddo
        enddo
        enddo

      endif ! flag_soil_layers==1


      ! Initialize liquid and frozen soil moisture from total soil moisture
      ! and soil temperature, and also soil moisture availability in the top
      ! layer
      call ruclsminit( debug_print, landmask,                         &
                 lsoil_ruc, isltyp, ivgtyp, xice, mavail,             &
                 soilh2o, smfr, soiltemp, soilm,                      &
                 ims,ime, jms,jme, kms,kme,                           &
                 its,ite, jts,jte, kts,kte                            )

      do j=jts,jte
      do i=its,ite
       if (land(i)) then
        wetness(i) = mavail(i,j)
        do k = 1, lsoil_ruc
          smois(i,k) = soilm(i,k,j)
          tslb(i,k)  = soiltemp(i,k,j)
          sh2o(i,k)  = soilh2o(i,k,j)
          smfrkeep(i,k)  = smfr(i,k,j)
        enddo 
       endif ! land(i)
      enddo
      enddo

      ! For non-land points, set RUC LSM fields to input (Noah or RUC) fields
        do i=1,im
          if (.not.land(i)) then
            do k=1,min(lsoil,lsoil_ruc)
              smois(i,k) = smc(i,k)
              tslb(i,k)  = stc(i,k)
              sh2o(i,k)  = slc(i,k)
            enddo
          endif
        enddo

      if(debug_print) then
        write (0,*)'End of RUC LSM initialization'
        write (0,*)'tslb(ipr)=',ipr,tslb(ipr,:)
        write (0,*)'smois(ipr)=',ipr,smois(ipr,:)
      endif ! debug_print

      end subroutine rucinit


end module lsm_ruc
