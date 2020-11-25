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

        real(kind=kind_phys), parameter :: zero = 0.0d0, one = 1.0d0, epsln = 1.0d-10

      contains

!> This subroutine calls set_soilveg_ruc() to specify vegetation and soil parameters for 
!! a given soil and land-use classification.
!! \section arg_table_lsm_ruc_init Argument Table
!! \htmlinclude lsm_ruc_init.html
!!
      subroutine lsm_ruc_init (me, master, isot, ivegsrc, nlunit,        &
                               flag_restart, flag_init,                  &
                               im, lsoil_ruc, lsoil, kice, nlev,         & ! in
                               lsm_ruc, lsm, slmsk, stype, vtype,        & ! in 
                               tsfc_lnd, tsfc_wat,                       & ! in
                               tg3, smc, slc, stc,                       & ! in
                               zs, sh2o, smfrkeep, tslb, smois, wetness, & ! out
                               tsice, errmsg, errflg)

      implicit none
!  ---  in
      integer,              intent(in)  :: me, master, isot, ivegsrc, nlunit
      logical,              intent(in)  :: flag_restart
      logical,              intent(in)  :: flag_init
      integer,              intent(in)  :: im
      integer,              intent(in)  :: lsoil_ruc
      integer,              intent(in)  :: lsoil
      integer,              intent(in)  :: kice
      integer,              intent(in)  :: nlev
      integer,              intent(in)  :: lsm_ruc, lsm


      real (kind=kind_phys), dimension(im), intent(in) :: slmsk
      real (kind=kind_phys), dimension(im), intent(in) :: stype
      real (kind=kind_phys), dimension(im), intent(in) :: vtype
      real (kind=kind_phys), dimension(im), intent(in) :: tsfc_lnd
      real (kind=kind_phys), dimension(im), intent(in) :: tsfc_wat
      real (kind=kind_phys), dimension(im), intent(in) :: tg3

      real (kind=kind_phys), dimension(im,lsoil), intent(in) :: smc,slc,stc

!  ---  in/out:
      real (kind=kind_phys), dimension(im), intent(inout) :: wetness

!  ---  out
      real (kind=kind_phys), dimension(:),            intent(out) :: zs
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: sh2o, smfrkeep
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: tslb, smois
      real (kind=kind_phys), dimension(im,kice),      intent(out) :: tsice

      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

! --- local
      real (kind=kind_phys), dimension(lsoil_ruc) :: dzs
      integer  :: ipr, i, k
      logical  :: debug_print
      integer, dimension(im) :: soiltyp, vegtype

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0 

      ipr = 10
      debug_print = .false.

      if (ivegsrc /= 1) then
        errmsg = 'The RUC LSM expects that the ivegsrc physics namelist parameter is 1. Exiting...'
        errflg = 1
        return
      end if
      if (isot > 1) then
        errmsg = 'The RUC LSM expects that the isot physics namelist parameter is 0, or 1. Exiting...'
        errflg = 1
        return
      end if

!> - Call rucinit() to initialize soil/ice/water  variables

      if ( debug_print) then
        write (0,*) 'RUC LSM initialization'
        write (0,*) 'lsoil_ruc, lsoil',lsoil_ruc, lsoil
        write (0,*) 'me, isot, ivegsrc, nlunit ',me, isot, ivegsrc, nlunit
        write (0,*) 'noah soil temp',stc(ipr,:)
        write (0,*) 'noah mois(ipr)',ipr,smc(ipr,:)
        write (0,*) 'stype=',stype(ipr)
        write (0,*) 'vtype=',vtype(ipr)
        write (0,*) 'tsfc_lnd=',tsfc_lnd(ipr)
        write (0,*) 'tsfc_wat=',tsfc_wat(ipr)
        write (0,*) 'tg3=',tg3(ipr)
        write (0,*) 'slmsk=',slmsk(ipr)
        write (0,*) 'flag_init =',flag_init
        write (0,*) 'flag_restart =',flag_restart
      endif

      !--- initialize soil vegetation
      call set_soilveg_ruc(me, isot, ivegsrc, nlunit)

      soiltyp(:) = 0
      vegtype(:) = 0

      do i  = 1, im ! i - horizontal loop
        if (slmsk(i) == 2.) then
        !-- ice
          if (isot == 1) then
            soiltyp(i)  = 16
          else
            soiltyp(i)  = 9
          endif
          if (ivegsrc == 1) then
            vegtype(i)  = 15
          elseif(ivegsrc == 2) then
            vegtype(i)  = 13
          endif
        else
        !-- land or water
          soiltyp(i)  = int( stype(i)+0.5 )
          vegtype(i)  = int( vtype(i)+0.5 )
          if (soiltyp(i)  < 1) soiltyp(i)  = 14
          if (vegtype(i)  < 1) vegtype(i)  = 17
        endif
      enddo

      call init_soil_depth_3 ( zs , dzs , lsoil_ruc )

        call rucinit   (flag_restart, im, lsoil_ruc, lsoil, nlev,   & ! in
                        me, master, lsm_ruc, lsm, slmsk,            & ! in
                        soiltyp, vegtype,                           & ! in
                        tsfc_lnd, tsfc_wat, tg3,                    & ! in
                        zs, dzs, smc, slc, stc,                     & ! in
                        sh2o, smfrkeep, tslb, smois,                & ! out
                        wetness, errmsg, errflg)

        do i  = 1, im ! i - horizontal loop
          do k = 1, min(kice,lsoil_ruc)
          ! - at initial time set sea ice T (tsice) 
          !   equal to TSLB, initialized from the Noah STC variable
             tsice   (i,k) = tslb(i,k)
          enddo
        enddo ! i

      !-- end of initialization

      if ( debug_print) then
        write (0,*) 'ruc soil tslb',tslb(ipr,:)
        write (0,*) 'ruc soil tsice',tsice(ipr,:)
        write (0,*) 'ruc soil smois',smois(ipr,:)
        write (0,*) 'ruc wetness',wetness(ipr)
      endif

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
!> \section arg_table_lsm_ruc_run Argument Table
!! \htmlinclude lsm_ruc_run.html
!!
!>\section gen_lsmruc GSD RUC LSM General Algorithm
      subroutine lsm_ruc_run                                            & ! inputs
     &     ( iter, me, master, delt, kdt, im, nlev, lsm_ruc, lsm,       &
     &       imp_physics, imp_physics_gfdl, imp_physics_thompson,       &
     &       do_mynnsfclay, lsoil_ruc, lsoil, rdlai, zs,                &
     &       t1, q1, qc, soiltyp, vegtype, sigmaf, laixy,               &
     &       dlwflx, dswsfc, snet, tg3,                                 &
     &       land, icy,  lake,                                          &
     &       rainnc, rainc, ice, snow, graupel,                         &
     &       prsl1, zf, wind, shdmin, shdmax, alvwf, alnwf,             &
     &       srflag, snoalb, isot, ivegsrc, fice, smcwlt2, smcref2,     &
     ! --- constants
     &       con_cp, con_rd, con_rv, con_g, con_pi, con_hvap,           &
     &       con_fvirt,                                                 &
     ! for water
     &       ch_wat, tskin_wat,                                         &
     ! --- in/outs for ice and land
     &       semis_lnd, semis_ice,                                      &
     &       sncovr1_lnd, weasd_lnd, snwdph_lnd, tskin_lnd,             &
     &       sncovr1_ice, weasd_ice, snwdph_ice, tskin_ice,             &
     ! for land
     &       smois, tsice, tslb, sh2o, keepfr, smfrkeep,                & ! on RUC levels
     &       canopy, trans, tsurf_lnd, tsnow_lnd, z0rl_lnd,             &
     &       sfcqc_lnd, sfcqv_lnd,                                      &
     &       qsurf_lnd, gflux_lnd, evap_lnd, hflx_lnd,                  &
     &       runof, runoff, srunoff, drain,                             &
     &       cm_lnd, ch_lnd, evbs, evcw, stm, wetness,                  &
     &       snowfallac_lnd,                                            &
     ! for ice
     &       sfcqc_ice, sfcqv_ice,                                      &
     &       tice, tsurf_ice, tsnow_ice, z0rl_ice,                      &
     &       qsurf_ice, gflux_ice, evap_ice, ep1d_ice, hflx_ice,        &
     &       cm_ice, ch_ice, snowfallac_ice,                            &
     ! --- out
     &       rhosnf, sbsno,                                             &
     &       cmm_lnd, chh_lnd, cmm_ice, chh_ice,                        &
     !
     &       flag_iter, flag_guess, flag_init, flag_restart,            &
     &       flag_cice, frac_grid, errmsg, errflg                       &
     &     )

      implicit none

!  ---  constant parameters:
      real(kind=kind_phys), parameter :: rhoh2o  = 1000.0
      real(kind=kind_phys), parameter :: stbolt  = 5.670400e-8
      real(kind=kind_phys), parameter :: cimin   = 0.15  !--- in GFS
      !real(kind=kind_phys), parameter :: cimin   = 0.02 !--- minimum ice concentration, 0.15 in GFS
      real(kind=kind_phys), parameter :: con_tice = 271.2

!  ---  input:
      integer, intent(in) :: me, master
      integer, intent(in) :: im, nlev, iter, lsoil_ruc, lsoil, kdt, isot, ivegsrc
      integer, intent(in) :: lsm_ruc, lsm
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson

      real (kind=kind_phys), dimension(im), intent(in) ::         &
     &       t1, sigmaf, laixy, dlwflx, dswsfc, snet, tg3,        &
     &       prsl1, wind, shdmin, shdmax,                         &
     &       snoalb, alvwf, alnwf, zf, qc, q1,                    &
     ! for land
     &       cm_lnd, ch_lnd,                                      &
     ! for water
     &       ch_wat, tskin_wat,                                   &
     ! for ice
     &       cm_ice, ch_ice

      real (kind=kind_phys),  intent(in) :: delt
      real (kind=kind_phys),  intent(in) :: con_cp, con_rv, con_g,       &
                                            con_pi, con_rd,              &
                                            con_hvap, con_fvirt

      logical, dimension(im), intent(in) :: flag_iter, flag_guess
      logical, dimension(im), intent(in) :: land, icy, lake
      logical, dimension(im), intent(in) :: flag_cice
      logical,                intent(in) :: frac_grid
      logical,                intent(in) :: do_mynnsfclay

      logical,                intent(in) :: rdlai

!  ---  in/out:
      integer, dimension(im), intent(inout) :: soiltyp, vegtype
      real (kind=kind_phys), dimension(lsoil_ruc), intent(in) :: zs
      real (kind=kind_phys), dimension(im), intent(inout) :: srflag,     &
     &       canopy, trans, smcwlt2, smcref2,                            & 
     ! for land
     &       weasd_lnd, snwdph_lnd, tskin_lnd,                           &
     &       tsurf_lnd, z0rl_lnd, tsnow_lnd,                             &
     &       sfcqc_lnd, sfcqv_lnd,                                       &
     ! for ice
     &       weasd_ice, snwdph_ice, tskin_ice,                           &
     &       tsurf_ice, z0rl_ice, tsnow_ice,                             &
     &       sfcqc_ice, sfcqv_ice, fice, tice

!  ---  in
      real (kind=kind_phys), dimension(im), intent(in) ::                &
     &       rainnc, rainc, ice, snow, graupel
!  ---  in/out:
!  --- on RUC levels
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) ::   &
     &       smois, tsice, tslb, sh2o, keepfr, smfrkeep

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(inout) ::             &
     &       rhosnf, runof, drain, runoff, srunoff, evbs, evcw,          &
     &       stm, wetness, semis_lnd, semis_ice,                         &
     ! for land
     &       sncovr1_lnd, qsurf_lnd, gflux_lnd, evap_lnd,                &
     &       cmm_lnd, chh_lnd, hflx_lnd, sbsno,                          &
     &       snowfallac_lnd,                                             &
     ! for ice
     &       sncovr1_ice, qsurf_ice, gflux_ice, evap_ice, ep1d_ice,      &
     &       cmm_ice, chh_ice, hflx_ice, snowfallac_ice

      logical,          intent(in)  :: flag_init, flag_restart
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rho,                      &
     &       q0, qs1,                                                   &
     &       tprcp_old, srflag_old, sr_old, canopy_old, wetness_old,    &
     ! for land
     &       weasd_lnd_old, snwdph_lnd_old, tskin_lnd_old,              &
     &       tsnow_lnd_old, snowfallac_lnd_old,                         &
     &       sfcqv_lnd_old, sfcqc_lnd_old, z0rl_lnd_old,                &
     &       sncovr1_lnd_old,                                           &
     ! for ice
     &       weasd_ice_old, snwdph_ice_old, tskin_ice_old,              &
     &       tsnow_ice_old, snowfallac_ice_old,                         &
     &       sfcqv_ice_old, sfcqc_ice_old, z0rl_ice_old,                &
     &       sncovr1_ice_old


      real (kind=kind_phys), dimension(lsoil_ruc) :: et

      real (kind=kind_phys), dimension(im,lsoil_ruc,1) :: smsoil,       &
           slsoil, stsoil, smfrsoil, keepfrsoil, stsice

      real (kind=kind_phys), dimension(im,lsoil_ruc) :: smois_old,      &
     &       tsice_old, tslb_old, sh2o_old,                             &
     &       keepfr_old, smfrkeep_old

      real (kind=kind_phys),dimension (im,1,1)      ::                  &
     &     conflx2, sfcprs, sfctmp, q2, qcatm, rho2 
      real (kind=kind_phys),dimension (im,1)        ::                  &
     &     albbck_lnd, alb_lnd, chs_lnd, flhc_lnd, flqc_lnd,            &
     &     wet, wet_ice, smmax, cmc, drip,  ec, edir, ett,              &
     &     dew_lnd, lh_lnd, esnow_lnd, etp, qfx_lnd, acceta,            &
     &     ffrozp, lwdn, prcp, xland, xland_wat, xice, xice_lnd,        &
     &     graupelncv, snowncv, rainncv, raincv,                        &
     &     solnet_lnd, sfcexc,                                          &
     &     runoff1, runoff2, acrunoff,                                  &
     &     sfcems_lnd, hfx_lnd, shdfac, shdmin1d, shdmax1d,             &
     &     sneqv_lnd, snoalb1d_lnd, snowh_lnd, snoh_lnd, tsnav_lnd,     &
     &     snomlt_lnd, sncovr_lnd, soilw, soilm, ssoil_lnd,             &
     &     soilt_lnd, tbot,                                             &
     &     xlai, swdn, z0_lnd, znt_lnd, rhosnfr, infiltr,               &
     &     precipfr, snfallac_lnd, acsn,                                &
     &     qsfc_lnd, qsg_lnd, qvg_lnd, qcg_lnd, soilt1_lnd, chklowq
     ! ice
      real (kind=kind_phys),dimension (im,1)        ::                  &
     &     albbck_ice, alb_ice, chs_ice, flhc_ice, flqc_ice,            &
     &     dew_ice, lh_ice, esnow_ice, qfx_ice,                         &
     &     solnet_ice, sfcems_ice, hfx_ice,                             &
     &     sneqv_ice, snoalb1d_ice, snowh_ice, snoh_ice, tsnav_ice,     &
     &     snomlt_ice, sncovr_ice, ssoil_ice, soilt_ice,                &
     &     z0_ice, znt_ice, snfallac_ice,                               &
     &     qsfc_ice, qsg_ice, qvg_ice, qcg_ice, soilt1_ice


      real (kind=kind_phys) :: xice_threshold
      real (kind=kind_phys) :: fwat, qsw, evapw, hfxw

      character(len=256) :: llanduse  !< Land-use dataset.  Valid values are :
                                      !! "USGS" (USGS 24/27 category dataset) and
                                      !! "MODIFIED_IGBP_MODIS_NOAH" (MODIS 20-category dataset)

      integer :: nscat, nlcat
      real (kind=kind_phys), dimension(:,:,:), allocatable :: landusef !< fractional landuse
      real (kind=kind_phys), dimension(:,:,:), allocatable :: soilctop !< fractional soil type

      integer :: nsoil, iswater, isice
      integer, dimension (1:im,1:1) :: stype_wat, vtype_wat
      integer, dimension (1:im,1:1) :: stype_lnd, vtype_lnd
      integer, dimension (1:im,1:1) :: stype_ice, vtype_ice
      integer :: ipr

      ! local
      integer :: ims,ime, its,ite, jms,jme, jts,jte, kms,kme, kts,kte
      integer :: l, k, i, j,  fractional_seaice

      logical :: flag(im), flag_ice_uncoupled(im)
      logical :: rdlai2d, myj, frpcpn
      logical :: debug_print
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ipr = 10

      debug_print=.false.

      chklowq = 1.

      do i  = 1, im ! i - horizontal loop
        ! - Set flag for ice points for uncoupled model (islmsk(i) == 4 when coupled to CICE)
        ! - Exclude ice on the lakes if the lake model is turned on.
        flag_ice_uncoupled(i) = (icy(i) .and. .not. flag_cice(i) .and. .not.  lake(i))
        !> - Set flag for land and ice points.
        !- 10may19 - ice points are turned off.
        flag(i) = land(i) .or. flag_ice_uncoupled(i)
      enddo

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
        write (0,*)'soiltyp=',ipr,soiltyp(ipr)
        write (0,*)'vegtype=',ipr,vegtype(ipr)
        write (0,*)'kdt, iter =',kdt,iter
        write (0,*)'flag_init =',flag_init
        write (0,*)'flag_restart =',flag_restart
      endif
 
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

      !> -- number of soil categories          
      !if(isot == 1) then
      !nscat = 19 ! stasgo
      !else
      !nscat = 9  ! zobler
      !endif
      !> -- set parameters for IGBP land-use data
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
         xice_threshold = 0.02 ! HRRR value
        !xice_threshold = 0.15 ! consistent with GFS physics
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
        if (flag(i) .and. flag_guess(i)) then
          !> - Save land-related prognostic fields for guess run.
          !if(me==0 .and. i==ipr) write (0,*)'before call to RUC guess run', i
          wetness_old(i)         = wetness(i)
          canopy_old(i)          = canopy(i)
          !srflag_old(i)          = srflag(i)
          !acsnow_old(i)          = acsnow(i)
          ! for land
          weasd_lnd_old(i)       = weasd_lnd(i)
          snwdph_lnd_old(i)      = snwdph_lnd(i)
          tskin_lnd_old(i)       = tskin_lnd(i)
          tsnow_lnd_old(i)       = tsnow_lnd(i)
          snowfallac_lnd_old(i)  = snowfallac_lnd(i)
          sfcqv_lnd_old(i)       = sfcqv_lnd(i)
          sfcqc_lnd_old(i)       = sfcqc_lnd(i)
          z0rl_lnd_old(i)        = z0rl_lnd(i)
          sncovr1_lnd_old(i)     = sncovr1_lnd(i)
          ! for ice
          weasd_ice_old(i)       = weasd_ice(i)
          snwdph_ice_old(i)      = snwdph_ice(i)
          tskin_ice_old(i)       = tskin_ice(i)
          tsnow_ice_old(i)       = tsnow_ice(i)
          snowfallac_ice_old(i)  = snowfallac_ice(i)
          sfcqv_ice_old(i)       = sfcqv_ice(i)
          sfcqc_ice_old(i)       = sfcqc_ice(i)
          z0rl_ice_old(i)        = z0rl_ice(i)
          sncovr1_ice_old(i)     = sncovr1_ice(i)

          do k = 1, lsoil_ruc
            smois_old(i,k)  = smois(i,k)
            tslb_old(i,k)   = tslb(i,k)
            sh2o_old(i,k)   = sh2o(i,k)
            keepfr_old(i,k) = keepfr(i,k)
            smfrkeep_old(i,k) = smfrkeep(i,k)
            ! for ice
            tsice_old(i,k)   = tsice(i,k)
          enddo
        endif
      enddo ! im

!  --- ...  initialization block

      do j  = 1, 1
      do i  = 1, im ! i - horizontal loop
        if (flag_iter(i) .and. flag(i)) then
          evap_lnd(i)  = 0.0
          evap_ice(i)  = 0.0
          hflx_lnd (i)  = 0.0
          hflx_ice (i)  = 0.0
          gflux_lnd(i)  = 0.0
          gflux_ice(i)  = 0.0
          drain(i)  = 0.0
          canopy(i) = max(canopy(i), 0.0)

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno(i)  = 0.0

          !local i,j arrays
          dew_lnd(i,j)      = 0.0
          dew_ice(i,j)      = 0.0
          soilm(i,j)        = 0.0
          smmax(i,j)        = 0.0
          hfx_lnd(i,j)      = 0.0
          hfx_ice(i,j)      = 0.0
          qfx_lnd(i,j)      = 0.0
          qfx_ice(i,j)      = 0.0
          lh_lnd(i,j)       = 0.0
          lh_ice(i,j)       = 0.0
          acsn(i,j)         = 0.0
          sfcexc(i,j)       = 0.0
          acceta(i,j)       = 0.0
          ssoil_lnd(i,j)    = 0.0
          ssoil_ice(i,j)    = 0.0
          snomlt_lnd(i,j)   = 0.0
          snomlt_ice(i,j)   = 0.0
          infiltr(i,j)      = 0.0
          runoff1(i,j)      = 0.0
          runoff2(i,j)      = 0.0
          acrunoff(i,j)     = 0.0
          snfallac_lnd(i,j) = 0.0
          snfallac_ice(i,j) = 0.0
          rhosnfr(i,j)      = 0.0
          precipfr(i,j)     = 0.0

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
        endif ! flag_iter & flag
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
        xice(i,j)  = 0.
        if (flag_iter(i) .and. flag(i)) then

        if (frpcpn) then
          ffrozp(i,j) = srflag(i)
        else
          ffrozp(i,j) = real(nint(srflag(i)),kind_phys)
        endif

        !-- rdlai is .false. when the LAI data is not available in the
        !    - INPUT/sfc_data.nc

        rdlai2d = rdlai

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
!!\n  \a prcp    - time-step total precip (\f$kg m^{-2} \f$)
!!\n  \a raincv  - time-step convective precip (\f$kg m^{-2} \f$)
!!\n  \a rainncv - time-step non-convective precip (\f$kg m^{-2} \f$)
!!\n  \a graupelncv - time-step graupel (\f$kg m^{-2} \f$)
!!\n  \a snowncv - time-step snow (\f$kg m^{-2} \f$)
!!\n  \a precipfr - time-step precipitation in solod form (\f$kg m^{-2} \f$)
!!\n  \a shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
!!\n  \a shdmin  - minimum areal fractional coverage of green vegetation -> !shdmin1d
!!\n  \a shdmax  - maximum areal fractional coverage of green vegetation -> !shdmax1d
!!\n  \a tbot    - bottom soil temperature (local yearly-mean sfc air temp)

        lwdn(i,j)   = dlwflx(i)         !..downward lw flux at sfc in w/m2
        swdn(i,j)   = dswsfc(i)         !..downward sw flux at sfc in w/m2

        ! all precip input to RUC LSM is in [mm]
        !prcp(i,j)       = rhoh2o * tprcp(i)                   ! tprcp in [m] - convective plus explicit
        !raincv(i,j)     = rhoh2o * rainc(i)                   ! total time-step convective precip
        !rainncv(i,j)    = rhoh2o * max(rain(i)-rainc(i),0.0)  ! total time-step explicit precip 
        !graupelncv(i,j) = rhoh2o * graupel(i)
        !snowncv(i,j)    = rhoh2o * snow(i)
        prcp(i,j)       = rhoh2o * (rainc(i)+rainnc(i))        ! tprcp in [m] - convective plus explicit
        raincv(i,j)     = rhoh2o * rainc(i)                    ! total time-step convective precip
        rainncv(i,j)    = rhoh2o * rainnc(i)                   ! total time-step explicit precip 
        graupelncv(i,j) = rhoh2o * graupel(i)
        snowncv(i,j)    = rhoh2o * snow(i)
        ! ice precipitation is not used
        ! precipfr(i,j)   = rainncv(i,j) * ffrozp(i,j)

        ! ice not used
        ! precipfr(i,j)   = rainncv(i,j) * ffrozp(i,j)
        !acsn(i,j)         = acsnow(i)
        acsn(i,j) = 0.0

        ! --- units %
        shdfac(i,j) = sigmaf(i)*100.
        shdmin1d(i,j) = shdmin(i)*100.
        shdmax1d(i,j) = shdmax(i)*100.

        tbot(i,j) = tg3(i)

!>  -   3. canopy/soil characteristics (s):
!!\n \a vegtyp  - vegetation type (integer index)                   -> vtype
!!\n \a soiltyp - soil type (integer index)                         -> stype
!!\n \a sfcems  -  surface emmisivity                                   -> sfcemis
!!\n \a 0.5*(alvwf + alnwf) - backround snow-free surface albedo (fraction)         -> albbck
!!\n \a snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d

        if(ivegsrc == 1) then   ! IGBP - MODIS
            vtype_wat(i,j) = 17 ! 17 - water (oceans and lakes) in MODIS
            stype_wat(i,j) = 14
            xland_wat(i,j) = 2. ! xland = 2 for water
            vtype_lnd(i,j) = vegtype(i)
            stype_lnd(i,j) = soiltyp(i)
            vtype_ice(i,j) = 15 ! MODIS
            if(isot == 0) then
              stype_ice(i,j) = 9  ! ZOBLER
            else
              stype_ice(i,j) = 16 ! STASGO
            endif
        !> - Prepare land/ice/water masks for RUC LSM
        !SLMSK0   - SEA(0),LAND(1),ICE(2) MASK
          !if(islmsk(i) == 0.) then
          !elseif(islmsk(i) == 1.) then ! land

          if(land(i)) then ! some land
            xland(i,j) = 1.
            xice_lnd(i,j) = 0.
          elseif(flag_ice_uncoupled(i)) then  ! some ice
            xland(i,j) = 1.
            xice(i,j)  = fice(i)  ! fraction of sea-ice
          endif
        else
          write (0,*)'MODIS landuse is not available'
        endif

        if(rdlai2d) then
          xlai(i,j) = laixy(i)
        else
          xlai(i,j) = 0.
        endif

   if (land(i)) then ! at least some land in the grid cell

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
!!\n \a snowh_lnd  - actual snow depth (\f$m\f$)
!!\n \a sneqv_lnd  - liquid water-equivalent snow depth (\f$m\f$)
!!\n \a sncovr_lnd - fraction of snow in the grid cell
!!\n \a chh_lnd    - surface exchange coefficient for heat (\f$m s^{-1}\f$)      -> chs
!!\n \a z0_lnd     - surface roughness (\f$m\f$)     -> zorl(\f$cm\f$)
!!\n \a qsfc_lnd   - specific humidity at surface (\f$kg kg^{-1}\f$)
!!\n \a qvg_lnd    - water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n \a qsg_lnd    - saturated water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n \a qcg_lnd    - cloud water mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n \a solnet_lnd - net sw radiation flux (dn-up) (\f$W m^{-2}\f$)

        qvg_lnd(i,j)    = sfcqv_lnd(i)
        qsfc_lnd(i,j)   = sfcqv_lnd(i)/(1.+sfcqv_lnd(i))
        qsg_lnd(i,j)    = rslf(prsl1(i),tsurf_lnd(i))
        qcg_lnd(i,j)    = sfcqc_lnd(i) 
        sfcems_lnd(i,j) = semis_lnd(i)
        sncovr_lnd(i,j) = sncovr1_lnd(i)
        snoalb1d_lnd(i,j) = snoalb(i)
        albbck_lnd(i,j)   = max(0.01, 0.5 * (alvwf(i) + alnwf(i))) 
        ! alb_lnd takes into account snow on the ground
        if (sncovr_lnd(i,j) > 0.) then
        !- averaged of snow-free and snow-covered
          alb_lnd(i,j) = albbck_lnd(i,j) * (1.-sncovr_lnd(i,j)) + snoalb(i) * sncovr_lnd(i,j)
        else
          alb_lnd(i,j) = albbck_lnd(i,j)
        endif
        solnet_lnd(i,j) = dswsfc(i)*(1.-alb_lnd(i,j)) !snet(i) !..net sw rad flx (dn-up) at sfc in w/m2

        cmc(i,j) = canopy(i)            !  [mm] 
        soilt_lnd(i,j) = tsurf_lnd(i)            ! clu_q2m_iter
        ! sanity check for snow temperature tsnow
        if (tsnow_lnd(i) > 0. .and. tsnow_lnd(i) < 273.15) then
          soilt1_lnd(i,j) = tsnow_lnd(i)
        else
          soilt1_lnd(i,j) = tsurf_lnd(i)
        endif
        tsnav_lnd(i,j) = 0.5*(soilt_lnd(i,j) + soilt1_lnd(i,j)) - 273.15
        do k = 1, lsoil_ruc
          smsoil  (i,k,j) = smois(i,k)
          slsoil  (i,k,j) = sh2o(i,k)
          stsoil  (i,k,j) = tslb(i,k)
          smfrsoil(i,k,j) = smfrkeep(i,k)
          keepfrsoil(i,k,j) = keepfr(i,k)
        enddo
        ! land
        if (wetness(i) > 0.) then
          wet(i,j) = wetness(i)
        else
          wet(i,j) = max(0.0001,smsoil(i,1,j)/0.3)
        endif

        chs_lnd (i,j)   = ch_lnd(i) * wind(i) ! compute conductance 
        flhc_lnd(i,j)   = chs_lnd(i,j) * rho(i) * con_cp ! * (1. + 0.84*q2(i,1,j))
        flqc_lnd(i,j)   = chs_lnd(i,j) * rho(i) * wet(i,j)
        ! for output
        cmm_lnd(i) = cm_lnd(i) * wind(i)
        chh_lnd(i) = chs_lnd(i,j) * rho(i)
        !
        snowh_lnd(i,j) = snwdph_lnd(i) * 0.001         ! convert from mm to m
        sneqv_lnd(i,j) = weasd_lnd(i)                  ! [mm]
        snfallac_lnd(i,j) = snowfallac_lnd(i)
        !> -- sanity checks on sneqv and snowh
        if (sneqv_lnd(i,j) /= 0.0 .and. snowh_lnd(i,j) == 0.0) then
          snowh_lnd(i,j) = 0.003 * sneqv_lnd(i,j) ! snow density ~300 kg m-3 
        endif

        if (snowh_lnd(i,j) /= 0.0 .and. sneqv_lnd(i,j) == 0.0) then
          sneqv_lnd(i,j) = 300. * snowh_lnd(i,j) ! snow density ~300 kg m-3 
        endif

        if (sneqv_lnd(i,j) > 0. .and. snowh_lnd(i,j) > 0.) then
          if(sneqv_lnd(i,j)/snowh_lnd(i,j) > 950.) then
            sneqv_lnd(i,j) = 300. * snowh_lnd(i,j)
          endif
        endif
        !  ---- ... outside sflx, roughness uses cm as unit
        z0_lnd(i,j)  = z0rl_lnd(i)/100.
        znt_lnd(i,j) = z0rl_lnd(i)/100.

        if(debug_print) then
          if(me==0 ) then
            write (0,*)'before LSMRUC for land'
            write (0,*)'sfcems(i,j) =',i,j,sfcems_lnd(i,j)
            write (0,*)'chklowq(i,j) =',i,j,chklowq(i,j)
            write (0,*)'chs(i,j) =',i,j,chs_lnd(i,j)
            write (0,*)'flqc(i,j) =',i,j,flqc_lnd(i,j)
            write (0,*)'flhc(i,j) =',i,j,flhc_lnd(i,j)
            write (0,*)'wet(i,j) =',i,j,wet(i,j)
            write (0,*)'cmc(i,j) =',i,j,cmc(i,j)
            write (0,*)'shdfac(i,j) =',i,j,shdfac(i,j)
            write (0,*)'alb(i,j) =',i,j,alb_lnd(i,j)
            write (0,*)'znt(i,j) =',i,j,znt_lnd(i,j)
            write (0,*)'z0(i,j) =',i,j,z0_lnd(i,j)
            write (0,*)'snoalb1d(i,j) =',i,j,snoalb1d_lnd(i,j)
            write (0,*)'landusef(i,:,j) =',i,j,landusef(i,:,j)
            write (0,*)'soilctop(i,:,j) =',i,j,soilctop(i,:,j)
            write (0,*)'nlcat=',nlcat
            write (0,*)'nscat=',nscat
            write (0,*)'qsfc(i,j) =',i,j,qsfc_lnd(i,j)
            write (0,*)'qvg(i,j) =',i,j,qvg_lnd(i,j)
            write (0,*)'qsg(i,j) =',i,j,qsg_lnd(i,j)
            write (0,*)'qcg(i,j) =',i,j,qcg_lnd(i,j)
            write (0,*)'dew(i,j) =',i,j,dew_lnd(i,j)
            write (0,*)'soilt(i,j) =',i,j,soilt_lnd(i,j)
            write (0,*)'tskin(i) =',i,j,tskin_lnd(i)
            write (0,*)'soilt1(i,j) =',i,j,soilt1_lnd(i,j)
            write (0,*)'tsnav(i,j) =',i,j,tsnav_lnd(i,j)
            write (0,*)'tbot(i,j) =',i,j,tbot(i,j)
            write (0,*)'vtype(i,j) =',i,j,vtype_lnd(i,j)
            write (0,*)'stype(i,j) =',i,j,stype_lnd(i,j)
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
            write (0,*)'acrunoff(i,j) =',i,j,acrunoff(i,j)
            write (0,*)'acsn(i,j) =',i,j,acsn(i,j)
            write (0,*)'shdmin1d(i,j) =',i,j,shdmin1d(i,j)
            write (0,*)'shdmax1d(i,j) =',i,j,shdmax1d(i,j)
            write (0,*)'rdlai2d =',rdlai2d
          endif
        endif

!> - Call RUC LSM lsmruc() for land. 
      call lsmruc(                                                           &
     &          delt, flag_init, flag_restart, kdt, iter, nsoil,             &
     &          graupelncv(i,j), snowncv(i,j), rainncv(i,j), raincv(i,j),    &
     &          zs, prcp(i,j), sneqv_lnd(i,j), snowh_lnd(i,j),               &
     &          sncovr_lnd(i,j),                                             &
     &          ffrozp(i,j), frpcpn,                                         &
     &          rhosnfr(i,j), precipfr(i,j),                                 &
!  ---  inputs:
     &          conflx2(i,1,j), sfcprs(i,1,j), sfctmp(i,1,j), q2(i,1,j),     &
     &          qcatm(i,1,j), rho2(i,1,j),                                   &
     &          lwdn(i,j), solnet_lnd(i,j), sfcems_lnd(i,j), chklowq(i,j),   &
     &          chs_lnd(i,j), flqc_lnd(i,j), flhc_lnd(i,j),                  &
!  ---  input/outputs:
     &          wet(i,j), cmc(i,j), shdfac(i,j), alb_lnd(i,j), znt_lnd(i,j), &
     &          z0_lnd(i,j), snoalb1d_lnd(i,j), albbck_lnd(i,j),             &
     &          xlai(i,j), landusef(i,:,j), nlcat,                           &
!  --- mosaic_lu and mosaic_soil are moved to the namelist
!     &          mosaic_lu, mosaic_soil,                                      &
     &          soilctop(i,:,j), nscat,                                      &
     &          qsfc_lnd(i,j), qsg_lnd(i,j), qvg_lnd(i,j), qcg_lnd(i,j),     &
     &          dew_lnd(i,j), soilt1_lnd(i,j),                               &
     &          tsnav_lnd(i,j), tbot(i,j), vtype_lnd(i,j), stype_lnd(i,j),   &
     &          xland(i,j), iswater, isice, xice_lnd(i,j), xice_threshold,   & !  xice=0. for the land portion of grid area
!  ---  constants
     &          con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, stbolt,     &
!  ---  input/outputs:
     &          smsoil(i,:,j), slsoil(i,:,j), soilm(i,j), smmax(i,j),        &
     &          stsoil(i,:,j), soilt_lnd(i,j),                               &
     &          hfx_lnd(i,j), qfx_lnd(i,j), lh_lnd(i,j),                     &
     &          infiltr(i,j), runoff1(i,j), runoff2(i,j), acrunoff(i,j),     &
     &          sfcexc(i,j), acceta(i,j), ssoil_lnd(i,j),                    &
     &          snfallac_lnd(i,j), acsn(i,j), snomlt_lnd(i,j),               &
     &          smfrsoil(i,:,j),keepfrsoil(i,:,j), .false.,                  &
     &          shdmin1d(i,j), shdmax1d(i,j), rdlai2d,                       &
     &          ims,ime, jms,jme, kms,kme,                                   &
     &          its,ite, jts,jte, kts,kte                                    )
        if(debug_print) then
          write (0,*)'after LSMRUC for land'
          write (0,*)'after sneqv(i,j) =',i,j,sneqv_lnd(i,j)
          write (0,*)'after snowh(i,j) =',i,j,snowh_lnd(i,j)
          write (0,*)'after sncovr(i,j) =',i,j,sncovr_lnd(i,j)
          write (0,*)'after vtype(i,j) =',i,j,vtype_lnd(i,j)
          write (0,*)'after stype(i,j) =',i,j,stype_lnd(i,j)
          write (0,*)'after wet(i,j) =',i,j,wet(i,j)
          write (0,*)'after cmc(i,j) =',i,j,cmc(i,j)
          write (0,*)'after qsfc(i,j) =',i,j,qsfc_lnd(i,j)
          write (0,*)'after qvg(i,j) =',i,j,qvg_lnd(i,j)
          write (0,*)'after qsg(i,j) =',i,j,qsg_lnd(i,j)
          write (0,*)'after qcg(i,j) =',i,j,qcg_lnd(i,j)
          write (0,*)'after dew(i,j) =',i,j,dew_lnd(i,j)
          write (0,*)'after soilt(i,j) =',i,j,soilt_lnd(i,j)
          write (0,*)'after tskin(i) =',i,j,tskin_lnd(i)
          write (0,*)'after soilt1(i,j) =',i,j,soilt1_lnd(i,j)
          write (0,*)'after tsnav(i,j) =',i,j,tsnav_lnd(i,j)
          write (0,*)'after smsoil(i,:,j)=',i,j,smsoil(i,:,j)
          write (0,*)'after slsoil(i,:,j)=',i,j,slsoil(i,:,j)
          write (0,*)'after stsoil(i,:,j)=',i,j,stsoil(i,:,j)
          write (0,*)'after smfrsoil(i,:,j)=',i,j,smfrsoil(i,:,j)
          write (0,*)'after keepfrsoil(i,:,j)=',i,j,keepfrsoil(i,:,j)
          write (0,*)'after soilm(i,j) =',i,j,soilm(i,j)
          write (0,*)'after smmax(i,j) =',i,j,smmax(i,j)
          write (0,*)'after hfx(i,j) =',i,j,hfx_lnd(i,j)
          write (0,*)'after qfx(i,j) =',i,j,qfx_lnd(i,j)
          write (0,*)'after lh(i,j) =',i,j,lh_lnd(i,j)
          write (0,*)'after infiltr(i,j) =',i,j,infiltr(i,j)
          write (0,*)'after runoff1(i,j) =',i,j,runoff1(i,j)
          write (0,*)'after runoff2(i,j) =',i,j,runoff2(i,j)
          write (0,*)'after ssoil(i,j) =',i,j,ssoil_lnd(i,j)
          write (0,*)'after snfallac(i,j) =',i,j,snfallac_lnd(i,j)
          write (0,*)'after acsn(i,j) =',i,j,acsn(i,j)
          write (0,*)'after snomlt(i,j) =',i,j,snomlt_lnd(i,j)
        endif


!> - RUC LSM: prepare variables for return to parent model and unit conversion.
!>  -   6. output (o):
!!\n \a lh     - actual latent heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a hfx    - sensible heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a ssoil   - soil heat flux (\f$W m^{-2}\f$: negative if downward from surface)
!!\n \a runoff1 - surface runoff (\f$m s^{-1}\f$), not infiltrating the surface
!!\n \a runoff2 - subsurface runoff (\f$m s^{-1}\f$), drainage out bottom
!!\n \a snoh    - phase-change heat flux from snowmelt (w m-2)
!!\n \a lh     - actual latent heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a hfx    - sensible heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a ssoil   - soil heat flux (\f$W m^{-2}\f$: negative if downward from surface)
!!\n \a runoff1 - surface runoff (\f$m s^{-1}\f$), not infiltrating the surface
!!\n \a runoff2 - subsurface runoff (\f$m s^{-1}\f$), drainage out bottom
!!\n \a snoh    - phase-change heat flux from snowmelt (w m-2)
!
!  --- ...  do not return the following output fields to parent model
!    ec      - canopy water evaporation (m s-1)
!    edir    - direct soil evaporation (m s-1)
!    et(nsoil)-plant transpiration from a particular root layer (m s-1)
!    ett     - total plant transpiration (m s-1)
!    esnow   - sublimation from (or deposition to if <0) snowpack (m s-1)
!    drip    - through-fall of precip and/or dew in excess of canopy
!              water-holding capacity (m)
!    snomlt  - snow melt (m) (water equivalent)
!    xlai    - leaf area index (dimensionless)
!    soilw   - available soil moisture in root zone (unitless fraction
!              between smcwlt and smcmax)
!    soilm   - total soil column moisture content (frozen+unfrozen) (m)
!    nroot   - number of root layers, a function of veg type, determined
!              in subroutine redprm.


        !evbs(i)  = edir(i,j)
        !evcw(i)  = ec(i,j)
        !trans(i) = ett(i,j)
        !sbsno(i) = esnow(i,j)
        !snohf(i) = snoh(i,j)

        ! Interstitial
        evap_lnd(i)   = qfx_lnd(i,j) / rho(i)           ! kinematic
        hflx_lnd(i)   = hfx_lnd(i,j) / (con_cp*rho(i))  ! kinematic
        gflux_lnd(i)  = ssoil_lnd(i,j)
        qsurf_lnd(i)   = qsfc_lnd(i,j)
        tsurf_lnd(i)   = soilt_lnd(i,j)
        stm(i)         = soilm(i,j) * 1.e-3 ! convert to [m]

        runof (i)  = runoff1(i,j)
        drain (i)  = runoff2(i,j)

        wetness(i) = wet(i,j)

        ! tsnow(i)   = soilt1(i,j)
        sfcqv_lnd(i)  = qvg_lnd(i,j)
        sfcqc_lnd(i)  = qcg_lnd(i,j)
        !  --- ...  units [m/s] = [g m-2 s-1] 
        rhosnf(i) = rhosnfr(i,j)
        !acsnow(i) = acsn(i,j)     ! kg m-2

        ! --- ... accumulated total runoff and surface runoff
        runoff(i)  = runoff(i)  + (drain(i)+runof(i)) * delt * 0.001 ! kg m-2
        srunoff(i) = srunoff(i) + runof(i) * delt * 0.001            ! kg m-2

        ! --- ... accumulated frozen precipitation (accumulation in lsmruc)
        snowfallac_lnd(i) = snfallac_lnd(i,j) ! kg m-2
        !  --- ...  unit conversion (from m to mm)
        snwdph_lnd(i)  = snowh_lnd(i,j) * 1000.0

        canopy(i)      = cmc(i,j)   ! mm
        weasd_lnd(i)   = sneqv_lnd(i,j) ! mm
        sncovr1_lnd(i) = sncovr_lnd(i,j)
        !  ---- ... outside RUC LSM, roughness uses cm as unit 
        !  (update after snow's effect)
        z0rl_lnd(i) = znt_lnd(i,j)*100.

        do k = 1, lsoil_ruc
          smois(i,k)  = smsoil(i,k,j)
          sh2o(i,k)   = slsoil(i,k,j)
          tslb(i,k)   = stsoil(i,k,j)
          keepfr(i,k)   = keepfrsoil(i,k,j)
          smfrkeep(i,k) = smfrsoil(i,k,j)
        enddo
     if(debug_print) then
       write (0,*)'LAND -i,j,stype_lnd,vtype_lnd',i,j,stype_lnd(i,j),vtype_lnd(i,j)
       write (0,*)'i,j,tsurf_lnd(i)',i,j,tsurf_lnd(i)
       write (0,*)'kdt,iter,stsoil(i,:,j)',kdt,iter,stsoil(i,:,j)
     endif
   endif ! end of land

   if (flag_ice_uncoupled(i)) then ! at least some ice in the grid cell
   !-- ice point

        sncovr_ice(i,j)   = sncovr1_ice(i)
        snoalb1d_ice(i,j) = 0.75 ! RAP value for max snow alb on ice
        albbck_ice(i,j)   = 0.55 ! RAP value for ice alb
        if (sncovr_ice(i,j) > 0.) then
        !- averaged of snow-free and snow-covered ice
          alb_ice(i,j) = albbck_ice(i,j) * (1.-sncovr_ice(i,j)) + snoalb1d_ice(i,j) * sncovr_ice(i,j)
        else
        ! snow-free ice
          alb_ice(i,j) = albbck_ice(i,j)
        endif

        solnet_ice(i,j) = dswsfc(i)*(1.-alb_ice(i,j))
        qvg_ice(i,j)    = sfcqv_ice(i)
        qsfc_ice(i,j)   = sfcqv_ice(i)/(1.+sfcqv_ice(i))
        qsg_ice(i,j)    = rslf(prsl1(i),tsurf_ice(i))
        qcg_ice(i,j)    = sfcqc_ice(i)
        sfcems_ice(i,j) = semis_ice(i)

        cmc(i,j) = canopy(i)                     ! [mm]
        soilt_ice(i,j) = tsurf_ice(i)            ! clu_q2m_iter
        if (tsnow_ice(i) > 0. .and. tsnow_ice(i) < 273.15) then
          soilt1_ice(i,j) = tsnow_ice(i)
        else
          soilt1_ice(i,j) = tsurf_ice(i)
        endif
        tsnav_ice(i,j) = 0.5*(soilt_ice(i,j) + soilt1_ice(i,j)) - 273.15
        do k = 1, lsoil_ruc
          stsice  (i,k,j) = tsice(i,k)
          smsoil  (i,k,j) = 1.
          slsoil  (i,k,j) = 0.
          smfrsoil(i,k,j) = 1.
          keepfrsoil(i,k,j) = 1.
        enddo

        wet_ice(i,j) = 1.

        chs_ice (i,j)   = ch_ice(i) * wind(i) ! compute conductance 
        flhc_ice(i,j)   = chs_ice(i,j) * rho(i) * con_cp ! * (1. + 0.84*q2(i,1,j))
        flqc_ice(i,j)   = chs_ice(i,j) * rho(i) * wet_ice(i,j)
        ! for output
        cmm_ice(i) = cm_ice (i) * wind(i)
        chh_ice(i) = chs_ice(i,j) * rho(i)


        snowh_ice(i,j) = snwdph_ice(i) * 0.001         ! convert from mm to m
        sneqv_ice(i,j) = weasd_ice(i)                  ! [mm]
        snfallac_ice(i,j) = snowfallac_ice(i)

        !> -- sanity checks on sneqv and snowh
        if (sneqv_ice(i,j) /= 0.0 .and. snowh_ice(i,j) == 0.0) then
          snowh_ice(i,j) = 0.003 * sneqv_ice(i,j) ! snow density ~300 kg m-3 
        endif

        if (snowh_ice(i,j) /= 0.0 .and. sneqv_ice(i,j) == 0.0) then
          sneqv_ice(i,j) = 300. * snowh_ice(i,j) ! snow density ~300 kg m-3 
        endif

        if (sneqv_ice(i,j) > 0. .and. snowh_ice(i,j) > 0.) then
          if(sneqv_ice(i,j)/snowh_ice(i,j) > 950.) then
            sneqv_ice(i,j) = 300. * snowh_ice(i,j)
          endif
        endif

        z0_ice(i,j)  = z0rl_ice(i)/100.
        znt_ice(i,j) = z0rl_ice(i)/100.

!> - Call RUC LSM lsmruc() for ice. 
      call lsmruc(                                                           &
     &          delt, flag_init, flag_restart, kdt, iter, nsoil,             &
     &          graupelncv(i,j), snowncv(i,j), rainncv(i,j), raincv(i,j),    &
     &          zs, prcp(i,j), sneqv_ice(i,j), snowh_ice(i,j),               &
     &          sncovr_ice(i,j),                                             &
     &          ffrozp(i,j), frpcpn,                                         &
     &          rhosnfr(i,j), precipfr(i,j),                                 &
!  ---  inputs:
     &          conflx2(i,1,j), sfcprs(i,1,j), sfctmp(i,1,j), q2(i,1,j),     &
     &          qcatm(i,1,j), rho2(i,1,j),                                   &
     &          lwdn(i,j), solnet_ice(i,j), sfcems_ice(i,j), chklowq(i,j),   &
     &          chs_ice(i,j), flqc_ice(i,j), flhc_ice(i,j),                  &
!  ---  input/outputs:
     &          wet_ice(i,j), cmc(i,j), shdfac(i,j), alb_ice(i,j),           &
     &          znt_ice(i,j), z0_ice(i,j), snoalb1d_ice(i,j),                &
     &          albbck_ice(i,j), xlai(i,j),landusef(i,:,j), nlcat,           &
!  --- mosaic_lu and mosaic_soil are moved to the namelist
!     &          mosaic_lu, mosaic_soil,                                      &
     &          soilctop(i,:,j), nscat,                                      &
     &          qsfc_ice(i,j), qsg_ice(i,j), qvg_ice(i,j), qcg_ice(i,j),     &
     &          dew_ice(i,j), soilt1_ice(i,j),                               &
     &          tsnav_ice(i,j), tbot(i,j), vtype_ice(i,j), stype_ice(i,j),   &
     &          xland(i,j), iswater, isice, xice(i,j), xice_threshold,       &
!  ---  constants
     &          con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, stbolt,     &
!  ---  input/outputs:
     &          smsoil(i,:,j), slsoil(i,:,j), soilm(i,j), smmax(i,j),        &
     &          stsice(i,:,j), soilt_ice(i,j),                               &
     &          hfx_ice(i,j), qfx_ice(i,j), lh_ice(i,j),                     &
     &          infiltr(i,j), runoff1(i,j), runoff2(i,j), acrunoff(i,j),     &
     &          sfcexc(i,j), acceta(i,j), ssoil_ice(i,j),                    &
     &          snfallac_ice(i,j), acsn(i,j), snomlt_ice(i,j),                 &
     &          smfrsoil(i,:,j),keepfrsoil(i,:,j), .false.,                  &
     &          shdmin1d(i,j), shdmax1d(i,j), rdlai2d,                       &
     &          ims,ime, jms,jme, kms,kme,                                   &
     &          its,ite, jts,jte, kts,kte                                    )

        ! Interstitial
        evap_ice(i)   = qfx_ice(i,j) / rho(i)           ! kinematic
        ep1d_ice(i)   = qfx_ice(i,j) * con_hvap
        hflx_ice(i)   = hfx_ice(i,j) / (con_cp*rho(i))  ! kinematic
        gflux_ice(i)  = ssoil_ice(i,j)

        qsurf_ice(i)   = qsfc_ice(i,j)
        tsurf_ice(i)   = soilt_ice(i,j)

        sfcqv_ice(i)  = qvg_ice(i,j)
        sfcqc_ice(i)  = qcg_ice(i,j)

        snowfallac_ice(i) = snfallac_ice(i,j) ! kg m-2
        !  --- ...  unit conversion (from m to mm)
        snwdph_ice(i)  = snowh_ice(i,j) * 1000.0
        weasd_ice(i)   = sneqv_ice(i,j) ! mm
        sncovr1_ice(i) = sncovr_ice(i,j)
        z0rl_ice(i) = znt_ice(i,j)*100.

        do k = 1, lsoil_ruc
          tsice(i,k)  = stsice(i,k,j)
          if(.not. frac_grid) then
            smois(i,k)  = 1.
            sh2o(i,k)   = 0.
            tslb(i,k)   = stsice(i,k,j)
            keepfr(i,k)   = 1.
            smfrkeep(i,k) = 1.
          endif
        enddo
     if(debug_print) then
       write (0,*)'ICE - i,j,stype_ice,vtype_ice)',i,j,stype_ice(i,j),vtype_ice(i,j)
       write (0,*)'i,j,tsurf_ice(i)',i,j,tsurf_ice(i)
       write (0,*)'kdt,iter,stsice(i,:,j)',kdt,iter,stsice(i,:,j)
     endif

   endif ! ice


        endif   ! end if_flag_iter_and_flag
      enddo   ! j
      enddo   ! i

      !-- Take care of fractional sea ice for uncoupled run with frac_grid=.false.
      !-- When frac_grid=.true. GFS_surface_composite will take care of this.
      do i  = 1, im   ! i - horizontal loop
        if ( flag_iter(i) .and. flag(i) ) then
        ! Do this only when the fractional grid is not turned on!
        ! Compute composite for a fractional sea ice: fice(i) < 1.
        ! This is needed for the 2-way coupling 
        ! in the upcoupled case (when sfc_cice is not used).
          if(.not. frac_grid) then
            if( flag_ice_uncoupled(i) .and. fice(i) < 1.) then
              !write (0,*)'Fractional sea ice at i', i, fice(i)
              fwat =  1.0 - fice(i)
             ! Check if ice fraction is below the minimum value: 15% in GFS
             ! physics.
              if (fice(i) < cimin) then ! cimin - minimal ice fraction
                        write (0,*)'warning: ice fraction is low:', fice(i)
                fice(i) = cimin
                fwat    = 1.0 - cimin
                write (0,*)'fix ice fraction: reset it to:', fice(i), tskin_wat(i)
              endif

            ! Compute the composite of ice and open water for 2-way coupling in the
            ! uncoupled sea-ice model. Use ice variables for the composite.
              tsurf_ice(i) = tsurf_ice(i) * fice(i) + min(con_tice,tskin_wat(i)) * fwat
              chh_ice(i)   = chh_ice(i) * fice(i) + ch_wat(i) * wind(i) * rho(i) * fwat
              hfxw         = ch_wat(i) * wind(i) * (min(con_tice,tskin_wat(i)) - t1(i))
              hflx_ice(i)  = hflx_ice(i) * fice(i) + hfxw * fwat
              qsw          = rslf(prsl1(i),min(con_tice,tskin_wat(i)))
              evapw        = ch_wat(i) * wind(i) * (qsw - q0(i))
              evap_ice(i)  = evap_ice(i) * fice(i) + evapw * fwat
              qsurf_ice(i) = q1(i) + evap_ice(i) * rho(i) / chh_ice(i) 
            endif ! flag_ice_uncoupled(i) .and. fice(i) < 1.
          endif ! flag_iter, icy, not frac_grid
        endif
      enddo ! i

!> - Restore land-related prognostic fields for guess run.
      do j  = 1, 1
      do i  = 1, im
        if (flag(i)) then
          if(debug_print) write (0,*)'end ',i,flag_guess(i),flag_iter(i)
          if (flag_guess(i)) then
            if(debug_print) write (0,*)'guess run'

            weasd_lnd(i)       = weasd_lnd_old(i)
            snwdph_lnd(i)      = snwdph_lnd_old(i)
            tskin_lnd(i)       = tskin_lnd_old(i)
            canopy(i)          = canopy_old(i)
            !srflag(i)          = srflag_old(i)
            tsnow_lnd(i)       = tsnow_lnd_old(i)
            snowfallac_lnd(i)  = snowfallac_lnd_old(i)
            !acsnow(i)          = acsnow_old(i)
            sfcqv_lnd(i)       = sfcqv_lnd_old(i)
            sfcqc_lnd(i)       = sfcqc_lnd_old(i)
            wetness(i)         = wetness_old(i)
            z0rl_lnd(i)        = z0rl_lnd_old(i)
            sncovr1_lnd(i)     = sncovr1_lnd_old(i)
            !ice
            weasd_ice(i)       = weasd_ice_old(i)
            snwdph_ice(i)      = snwdph_ice_old(i)
            tskin_ice(i)       = tskin_ice_old(i)
            tsnow_ice(i)       = tsnow_ice_old(i)
            snowfallac_ice(i)  = snowfallac_ice_old(i)
            sfcqv_ice(i)       = sfcqv_ice_old(i)
            sfcqc_ice(i)       = sfcqc_ice_old(i)
            z0rl_ice(i)        = z0rl_ice_old(i)
            sncovr1_ice(i)     = sncovr1_ice_old(i)

            do k = 1, lsoil_ruc
              smois(i,k)    = smois_old(i,k)
              tslb(i,k)     = tslb_old(i,k)
              tsice(i,k)    = tsice_old(i,k)
              sh2o(i,k)     = sh2o_old(i,k)
              keepfr(i,k)   = keepfr_old(i,k)
              smfrkeep(i,k) = smfrkeep_old(i,k)
            enddo
          else ! flag_guess
            if(debug_print) write (0,*)'iter run', i,j, tskin_ice(i),tsurf_ice(i)
            tskin_lnd(i) = tsurf_lnd(i)
            tskin_ice(i) = tsurf_ice(i)
            tice(i)      = tsurf_ice(i)
          endif ! flag_guess
        endif ! flag
      enddo  ! i
      enddo  ! j
!
      deallocate(soilctop)
      deallocate(landusef)
!
      return
!...................................
      end subroutine lsm_ruc_run
!-----------------------------------

!>\ingroup lsm_ruc_group
!! This subroutine contains RUC LSM initialization.
      subroutine rucinit        (restart, im, lsoil_ruc, lsoil, nlev,   & ! in
                                 me, master, lsm_ruc, lsm, slmsk,       & ! in
                                 soiltyp, vegtype,                      & ! in
                                 tskin_lnd, tskin_wat, tg3,             & ! !in
                                 zs, dzs, smc, slc, stc,                & ! in
                                 sh2o, smfrkeep, tslb, smois,           & ! out
                                 wetness, errmsg, errflg)

      implicit none

      logical,                                        intent(in   ) :: restart
      integer,                                        intent(in   ) :: lsm
      integer,                                        intent(in   ) :: lsm_ruc
      integer,                                        intent(in   ) :: im, nlev
      integer,                                        intent(in   ) :: lsoil_ruc
      integer,                                        intent(in   ) :: lsoil
      real (kind=kind_phys), dimension(im),           intent(in   ) :: slmsk
      real (kind=kind_phys), dimension(im),           intent(in   ) :: tskin_lnd, tskin_wat
      real (kind=kind_phys), dimension(im),           intent(in   ) :: tg3
      real (kind=kind_phys), dimension(1:lsoil_ruc),  intent(in   ) :: zs
      real (kind=kind_phys), dimension(1:lsoil_ruc),  intent(in   ) :: dzs
      real (kind=kind_phys), dimension(im,lsoil),     intent(in   ) :: smc !  Noah
      real (kind=kind_phys), dimension(im,lsoil),     intent(in   ) :: stc !  Noah
      real (kind=kind_phys), dimension(im,lsoil),     intent(in   ) :: slc !  Noah

      integer,               dimension(im),    intent(inout) :: soiltyp
      integer,               dimension(im),    intent(inout) :: vegtype
      real (kind=kind_phys), dimension(im),    intent(inout) :: wetness
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: smois! ruc
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: tslb ! ruc
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: sh2o ! ruc
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: smfrkeep ! ruc

      integer,          intent(in ) :: me
      integer,          intent(in ) :: master
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!> local
      logical :: debug_print
      logical :: smadj ! for soil mosture adjustment
      logical :: swi_init ! for initialization in terms of SWI (soil wetness index)

      integer :: flag_soil_layers, flag_soil_levels, flag_sst
      real (kind=kind_phys),    dimension(1:lsoil_ruc)  :: factorsm
      real (kind=kind_phys),    dimension(im)           :: smcref2
      real (kind=kind_phys),    dimension(im)           :: smcwlt2

      integer , dimension( 1:im , 1:1 )      :: ivgtyp
      integer , dimension( 1:im , 1:1)       :: isltyp
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: mavail
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
        write (0,*)'restart = ',restart
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
         write (0,*)'tskin_lnd(:)=',tskin_lnd(:)
         write (0,*)'tskin_wat(:)=',tskin_wat(:)
         write (0,*)'vegtype(ipr) ==', ipr, vegtype(ipr)
         write (0,*)'soiltyp(ipr) ==', ipr, soiltyp(ipr)
         write (0,*)'its,ite,jts,jte ',its,ite,jts,jte 
      endif


        do j=jts,jte !
        do i=its,ite ! i = horizontal loop

            sst(i,j) = tskin_wat(i)
            tbot(i,j) = tg3(i)
            ivgtyp(i,j) = vegtype(i)
            isltyp(i,j) = soiltyp(i)
          if (slmsk(i) == 0.) then
          !-- water
            tsk(i,j) = tskin_wat(i)
            landmask(i,j)=0.
          else
          !-- land or ice
            tsk(i,j) = tskin_lnd(i)
            landmask(i,j)=1.
          endif ! land(i)

        enddo
        enddo

      if ( flag_soil_layers == 1 ) then
      ! Noah lsm input
        do j=jts,jte !
        do i=its,ite ! i = horizontal loop

          st_input(i,1,j)=tsk(i,j)
          sm_input(i,1,j)=0.

          !--- initialize smcwlt2 and smcref2 with Noah values
          if(slmsk(i) == 1.) then
            smcref2 (i) = REFSMCnoah(soiltyp(i))
            smcwlt2 (i) = WLTSMCnoah(soiltyp(i))
          else
            smcref2 (i) = 1.
            smcwlt2 (i) = 0.
          endif

          do k=1,lsoil
             st_input(i,k+1,j)=stc(i,k)
             ! convert volumetric soil moisture to SWI (soil wetness index)
             if(slmsk(i) == 1. .and. swi_init) then
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
         if (slmsk(i) == 1.) then
         !-- land
           do k=1,lsoil_ruc
           ! convert from SWI to RUC volumetric soil moisture
             if(swi_init) then
               soilm(i,k,j) = dumsm(i,k,j) *                            &
                 (refsmc(isltyp(i,j))-drysmc(isltyp(i,j)))              &
                 + drysmc(isltyp(i,j))
             else
               soilm(i,k,j) = dumsm(i,k,j)
             endif
             soiltemp(i,k,j) = dumt(i,k,j)
           enddo ! k
         else
         !-- ice or water
           do k=1,lsoil_ruc
             soilm(i,k,j) = 1.
             soiltemp(i,k,j) = dumt(i,k,j)
           enddo ! k
         endif ! land
        enddo
        enddo

        if(debug_print) then
          write (0,*)'tsk(i,j),tbot(i,j),sst(i,j),landmask(i,j)' &
                  ,ipr,1,tsk(ipr,1),tbot(ipr,1),sst(ipr,1),landmask(ipr,1)
          write (0,*)'tskin_lnd(ipr)=',ipr,tskin_lnd(ipr)
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

           if (slmsk(i) == 1.) then

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

      call ruclsminit( debug_print, slmsk,                            &
                 lsoil_ruc, isltyp, ivgtyp, mavail,                   &
                 soilh2o, smfr, soiltemp, soilm,                      &
                 ims,ime, jms,jme, kms,kme,                           &
                 its,ite, jts,jte, kts,kte                            )

      do j=jts,jte
      do i=its,ite
        wetness(i) = mavail(i,j)
        do k = 1, lsoil_ruc
          smois(i,k) = soilm(i,k,j)
          tslb(i,k)  = soiltemp(i,k,j)
          sh2o(i,k)  = soilh2o(i,k,j)
          smfrkeep(i,k)  = smfr(i,k,j)
        enddo 
      enddo
      enddo

        !do i=1,im
        !    wetness (i) = 1.
        !    do k=1,min(lsoil,lsoil_ruc)
        !      smois(i,k) = smc(i,k)
        !      tslb(i,k)  = stc(i,k)
        !      sh2o(i,k)  = slc(i,k)
        !    enddo
        !enddo

      if(debug_print) then
        do i=1,im
        write (0,*)'End of RUC LSM initialization'
        write (0,*)'tslb(i)=',i,tslb(i,:)
        write (0,*)'smois(i)=',i,smois(i,:)
        write (0,*)'wetness(i)=',i,wetness(i)
        enddo
      endif ! debug_print

      end subroutine rucinit


end module lsm_ruc
