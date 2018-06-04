!>  \file sfc_drv_ruc.f
!!  This file contains the RUC land surface scheme driver.

      module lsm_ruc_pre
        use module_soil_pre
        use module_sf_ruclsm,  only: ruclsminit
        use machine,           only: kind_phys
      contains


#if 0
!! \section arg_table_lsm_ruc_pre_init Argument Table
!! | local_name     | standard_name                                                    | long_name                                               | units         | rank | type      |    kind   | intent | optional |
!|------------------|------------------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                           | horizontal loop extent                                  | count         |    0 | integer   |           | in     | F        |
!! | nlev           | vertical_dimension                                               | number of vertical levels                               | count         |    0 | integer   |           | in     | F        |
!! | km             | soil_vertical_dimension                                          | soil vertical layer dimension                           | count         |    0 | integer   |           | in     | F        |
!! | con_hvap       | latent_heat_of_vaporization_of_water_at_0C                       | latent heat of vaporization/sublimation (hvap)          | J kg-1        |    0 | real      | kind_phys | in     | F        |
!! | vegtype        | cell_vegetation_type                                             | vegetation type at each grid cell                       | index         |    1 | integer   |           | in     | F        |
!! | soiltyp        | cell_soil_type                                                   | soil type at each grid cell                             | index         |    1 | integer   |           | in     | F        |
!! | fice           | sea_ice_concentration                                            | ice fraction over open water                            | frac          |    1 | real      | kind_phys | in     | F        |
!! | islimsk        | sea_land_ice_mask                                                | landmask: sea/land/ice=0/1/2                            | flag          |    1 | integer   |           | in     | F        |
!! | tsk            | surface_skin_temperature                                         | surface skin temperature                                | K             |    1 | real      | kind_phys | in     | F        |
!! | tg3            | deep_soil_temperature                                            | bottom soil temperature                                 | K             |    1 | real      | kind_phys | in     | F        |
!! | smc            | volume_fraction_of_soil_moisture                                 | volumetric fraction of soil moisture                    | frac          |    2 | real      | kind_phys | in     | F        |
!! | slc            | volume_fraction_of_unfrozen_soil_moisture                        | volume fraction of unfrozen soil moisture               | frac          |    2 | real      | kind_phys | in     | F        |
!! | stc            | soil_temperature                                                 | soil temperature                                        | K             |    2 | real      | kind_phys | in     | F        |
!! | zs             | depth_of_soil_levels_for_land_surface_model                      | depth of soil levels for lsm                            | count         |    0 | real      | kind_phys | inout  | F        |
!! | smois          | volume_fraction_of_soil_moisture_for_land_surface_model          | volumetric fraction of soil moisture for lsm            | frac          |    2 | real      | kind_phys | inout  | F        |
!! | sh2o           | volume_fraction_of_unfrozen_soil_moisture_for_land_surface_model | volume fraction of unfrozen soil moisture for lsm       | frac          |    2 | real      | kind_phys | inout  | F        |
!! | tslb           | soil_temperature_for_land_surface_model                          | soil temperature for land surface model                 | K             |    2 | real      | kind_phys | inout  | F        |
!! | wet1           | normalized_soil_wetness                                          | normalized soil wetness                                 | frac          |    1 | real      | kind_phys | inout  | F        |
!! | lsm_ruc        | flag_for_ruc_land_surface_scheme                                 | flag for RUC land surface model                         | flag          |    0 | integer   |           | none   | F        |
!! | lsm            | flag_for_land_surface_scheme                                     | flag for land surface model                             | flag          |    0 | integer   |           | none   | F        |
!! | errmsg         | error_message                                                    | error message for error handling in CCPP                | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                       | error flag for error handling in CCPP                   | flag          |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine lsm_ruc_pre_init (im, km, nlev, soiltyp, vegtype, fice,  & ! in
                               islimsk, tsurf, tg3,                   & ! in
                               smc, slc, stc,                         & ! in
                               zs, sh2o, tslb, smois, wet1,           & ! out
                               errmsg, errflg)

      implicit none


      integer,                                 intent(in   ) :: im, nlev   
      integer,                                 intent(in   ) :: km   
      integer,               dimension(im),    intent(in   ) :: islimsk
      real (kind=kind_phys), dimension(im),    intent(in   ) :: tsurf
      real (kind=kind_phys), dimension(im),    intent(in   ) :: tg3
      real (kind=kind_phys), dimension(im,4),  intent(in   ) :: smc !  Noah
      real (kind=kind_phys), dimension(im,4),  intent(in   ) :: stc !  Noah
      real (kind=kind_phys), dimension(im,4),  intent(in   ) :: slc !  Noah

      integer,               dimension(im),    intent(inout) :: soiltyp
      integer,               dimension(im),    intent(inout) :: vegtype
      real (kind=kind_phys), dimension(im),    intent(inout) :: wet1
      real (kind=kind_phys), dimension(im),    intent(inout) :: fice
      real (kind=kind_phys), dimension(1,km,im), intent(inout) :: smois! ruc
      real (kind=kind_phys), dimension(1,km,im), intent(inout) :: tslb ! ruc
      real (kind=kind_phys), dimension(1,km,im), intent(inout) :: sh2o ! ruc

      real (kind=kind_phys), dimension(1:km), intent (out)  :: zs

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!> local
      logical :: debug_print

      integer :: flag_soil_layers, flag_soil_levels, flag_sst
      real (kind=kind_phys),    dimension(1:km) :: factorsm

      real (kind=kind_phys),    dimension( 1:1 , 1:im )       :: sst, landmask, tsk
      real (kind=kind_phys),    dimension( 1:1 , 1:im )       :: smtotn,smtotr
      real (kind=kind_phys),    dimension( 1:1 , 1:km, 1:im ) :: dumsm , dumt
      real (kind=kind_phys),    dimension( 1:1 , 1:km, 1:im ) :: smfr

      real (kind=kind_phys) :: st_input(1:1,1:km*3,1:im)        
      real (kind=kind_phys) :: sm_input(1:1,1:km*3,1:im)       

      integer               :: ids,ide, jds,jde, kds,kde, &
                               ims,ime, jms,jme, kms,kme, &
                               its,ite, jts,jte, kts,kte, &
                               i, j, k, l, num_soil_layers

      real(kind=kind_phys), dimension(1:km) :: zs2, dzs 
      integer,              dimension(1:4)  :: st_levels_input, sm_levels_input ! 4 - for Noah lsm

       ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

      if (lsm/=lsm_ruc) then
        write(errmsg,’(a,i0,a,i0)’) ‘ERROR in lsm_ruc_pre_init: namelist variable lsm=‘, lsm, ’ incompatible with RUC LSM, please set to ’, lsm_ruc
        errflg = 1
        return
      endif

       debug_print = .false.

         st_levels_input = (/ 5, 25, 70, 150/)    ! Noah soil levels
         sm_levels_input = (/ 5, 25, 70, 150/)    ! Noah soil levels

! Set internal dimensions
         ids = 1
         ims = 1
         its = 1
         ide = 1
         ime = 1
         ite = 1
         jds = 1
         jms = 1
         jts = 1
         jde = im 
         jme = im 
         jte = im 
         kds = 1
         kms = 1
         kts = 1
         kde = nlev
         kme = nlev
         kte = nlev

         num_soil_layers =  4 ! for Noah lsm

         CALL init_soil_depth_3 ( zs , dzs , km )

         flag_soil_layers = 1  ! 1 - for input from the Noah LSM
         flag_soil_levels = 0  ! 1 - for input from RUC LSM
         flag_sst=1

      IF ( islimsk(i) == 1 ) then
          landmask = 1 ! land
      else
          landmask = 0 ! water
          sst = tsk
      endif

       do i=its,ite
       do j=jts,jte
            st_input(i,1,j)=tsurf(j)
            sm_input(i,1,j)=0.
         do k=1,4
            st_input(i,k+1,j)=stc(j,k)
            sm_input(i,k+1,j)=smc(j,k)
         enddo
         do k=6,km * 3
            st_input(i,k,j)=0.
            sm_input(i,k,j)=0.
         enddo
! Noah soil moisture bucket 
         smtotn(i,j)=sm_input(i,2,j)*0.1 + sm_input(i,3,j)*0.2          &
                + sm_input(i,4,j)*0.7 + sm_input(i,5,j)*1.
       enddo
       enddo
         CALL init_soil_3_real ( tsk , tg3 , dumsm , dumt ,             &
                                 st_input , sm_input , landmask , sst , &
                                 zs , dzs ,                             &
                                 st_levels_input, sm_levels_input,      &
                                 km , num_soil_layers, num_soil_layers, &
                                 km * 3 , km * 3 , flag_sst,            &
                                 flag_soil_layers , flag_soil_levels ,  &
                                 ids , ide , jds , jde , kds , kde ,    &
                                 ims , ime , jms , jme , kms , kme ,    &
                                 its , ite , jts , jte , kts , kte )

      do i=its,ite
      do j=jts,jte

         do k=1,km
          smois(i,k,j)= dumsm(i,k,j)
          tslb(i,k,j) = dumt(i,k,j)
         enddo
      if(debug_print) then
         do k=1,km
           print *,'tslb(1,k,j)',j,k,tslb(1,k,j)
           print *,'smois(1,k,j)',j,k,smois(1,k,j)
         enddo
         do k=1,4
           print *,'stc(i,k,1)',j,k,stc(j,k)
           print *,'smc(i,k,1)',j,k,smc(j,k)
         enddo
      endif

      IF ( islimsk(i) == 1 ) then  ! Land
! initialize factor
        do k=1,km
           factorsm(k)=1.
        enddo

! RUC soil moisture bucket
           smtotr(i,j)=0.
        do k=1,km -1
          smtotr(i,j)=smtotr(i,j) + smois(i,k,j) *dzs(k)
        enddo

        if(debug_print) then
          print *,'from Noah to RUC: RUC bucket and Noah bucket at',i,j,smtotr(i,j),smtotn(i,j)
          print *,'before smois=',i,j,smois(i,:,j)
        endif

! RUC soil moisture correction to match Noah soil moisture bucket
        do k=1,km-1
           smois(i,k,j) = max(0.02,smois(i,k,j)*smtotn(i,j)/(0.9*smtotr(i,j)))
        enddo

        if( smois(i,2,j) > smois(i,1,j) .and. smois(i,3,j) > smois(i,2,j)) then
! typical for daytime, no recent precip
          factorsm(1) = 0.75
          factorsm(2) = 0.8
          factorsm(3) = 0.85
          factorsm(4) = 0.9
          factorsm(5) = 0.95
        endif
        do k=1,km
           smois(i,k,j) = factorsm(k) * smois(i,k,j)
        enddo
        if(debug_print) then
           print *,'after smois=',i,j,smois(i,:,j)
        endif
           smtotr(i,j) = 0.
        do k=1,km - 1
           smtotr(i,j)=smtotr(i,j) + smois(i,k,j) *dzs(k)
        enddo
        if(debug_print) then
            print *,'after correction: RUC bucket and Noah bucket at',i,j,smtotr(i,j),smtotn(i,j)
        endif
      ENDIF ! land
       enddo
       enddo

!> Initialize soil and vegetation parameters
        call ruclsminit( debug_print,                                   &
                   km, soiltyp, vegtype, fice, wet1,                    &
                   sh2o, smfr, tslb, smois,                             &
                   ims,ime, jms,jme, kms,kme,                           &
                   its,ite, jts,jte, kts,kte                            )

          if (errflg /= 0) return

      end subroutine lsm_ruc_pre_init

!! \brief Brief description of the subroutine
!!
#if 0
!! \section arg_table_lsm_ruc_pre_run Argument Table
!! | local_name     | standard_name                                               | long_name                                  | units      | rank | type      |    kind   | intent | optional |
!! |----------------|-------------------------------------------------------------|--------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                      | horizontal loop extent =1                  | count      |    0 | integer   |           | in     | F        |
!! | km             | soil_vertical_dimension                                     | soil vertical layer dimension              | count      |    0 | integer   |           | in     | F        |
!! | drain          | subsurface_runoff_flux                                      | subsurface runoff flux                     | g m-2 s-1  |    1 | real      | kind_phys | out    | F        |
!! | runof          | surface_runoff_flux                                         | surface runoff flux                        | g m-2 s-1  |    1 | real      | kind_phys | out    | F        |
!! | evbs           | soil_upward_latent_heat_flux                                | soil upward latent heat flux               | W m-2      |    1 | real      | kind_phys | out    | F        |
!! | evcw           | canopy_upward_latent_heat_flux                              | canopy upward latent heat flux             | W m-2      |    1 | real      | kind_phys | out    | F        |
!! | trans          | transpiration_flux                                          | total plant transpiration rate             | kg m-2 s-1 |    1 | real      | kind_phys | out    | F        |
!! | sbsno          | snow_deposition_sublimation_upward_latent_heat_flux         | latent heat flux from snow depo/subl       | W m-2      |    1 | real      | kind_phys | out    | F        |
!! | snowc          | surface_snow_area_fraction                                  | surface snow area fraction                 | frac       |    1 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                               | error message for error handling in CCPP   | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                  | error flag for error handling in CCPP      | flag       |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine lsm_ruc_pre_run                                       &
     &  (im,km,drain,runof,evbs,evcw,trans,sbsno,snowc,snohf,          &
     &   errmsg,errflg                                                 &
     &  )

      implicit none

!  ---  interface variables
      integer, intent(in) :: im, km

      real(kind=kind_phys), dimension(im), intent(inout)  ::            &
     &    drain,runof,evbs,evcw,trans,sbsno,snowc,snohf

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      drain(:)   = 0.0
      runof(:)   = 0.0
      evbs(:)    = 0.0
      evcw(:)    = 0.0
      trans(:)   = 0.0
      sbsno(:)   = 0.0
      snowc(:)   = 0.0

      end subroutine lsm_ruc_pre_run

      subroutine lsm_ruc_pre_finalize
      end subroutine lsm_ruc_pre_finalize

      end module lsm_ruc_pre


      module lsm_ruc_post
        use machine,           only: kind_phys
      contains

      subroutine lsm_ruc_post_init
      end subroutine lsm_ruc_post_init

      subroutine lsm_ruc_post_finalize
      end subroutine lsm_ruc_post_finalize

!> \brief Brief description of the subroutine
!!
#if 0
!! \section arg_table_lsm_ruc_post_run Argument Table
!! | local_name     | standard_name                                               | long_name                                  | units      | rank | type      |    kind   | intent | optional |
!! |----------------|-------------------------------------------------------------|--------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                      | horizontal loop extent                     | count      |    0 | integer   |           | in     | F        |
!! | km             | soil_vertical_dimension                                     | soil vertical layer dimension              | count      |    0 | integer   |           | in     | F        |
!! | flag_lssav     | flag_diagnostics                                            | flag for calculating diagnostic fields     | flag       |    0 | logical   |           | in     | F        |
!! | dtf            | time_step_for_dynamics                                      | dynamics time step                         | s          |    0 | real      | kind_phys | in     | F        |
!! | drain          | subsurface_runoff_flux                                      | subsurface runoff flux                     | g m-2 s-1  |    1 | real      | kind_phys | in     | F        |
!! | runof          | surface_runoff_flux                                         | surface runoff flux                        | g m-2 s-1  |    1 | real      | kind_phys | in     | F        |
!! | runoff         | total_runoff                                                | total runoff                               | kg m-2     |    1 | real      | kind_phys | inout  | F        |
!! | srunoff        | surface_runoff                                              | surface runoff                             | kg m-2     |    1 | real      | kind_phys | inout  | F        |
!! | errmsg         | error_message                                               | error message for error handling in CCPP   | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                  | error flag for error handling in CCPP      | flag       |    0 | integer   |           | out    | F        |
!!
#endif
!!  \section lsm_post_general General Algorithm
!!  \section lsm_post_detailed Detailed Algorithm
!!  @{

      subroutine lsm_ruc_post_run                                      &
     &  (im,km, flag_lssav,dtf,drain,runof,runoff,srunoff,errmsg,errflg &
     &  )
      implicit none

!  ---  interface variables
      integer, intent(in) :: im, km

      logical, intent(in) :: flag_lssav
      real(kind=kind_phys), intent (in)   :: dtf

      real(kind=kind_phys), dimension(im), intent(in   )  ::            &
     &    drain, runof

      real(kind=kind_phys), dimension(im), intent(inout)  ::            &
     &    runoff, srunoff

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if(flag_lssav) then
        runoff(:)  = runoff(:)  + (drain(:)+runof(:)) * dtf * 1000. ! kg m-2
        srunoff(:) = srunoff(:) + runof(:) * dtf * 1000.            ! kg m-2
      end if

      end subroutine lsm_ruc_post_run

!! @}
      end module lsm_ruc_post
!! @}

      module lsm_ruc

      use machine,           only: kind_phys
      use module_sf_ruclsm
      use funcphys, only : fpvs

      contains

      subroutine lsm_ruc_init
      end subroutine lsm_ruc_init

      subroutine lsm_ruc_finalize
      end subroutine lsm_ruc_finalize

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!      call sfc_drv_ruc                                                 !
! --- inputs                                                            !
!         ( kdt, im, km, u1, v1, t1, q1, qc, soiltyp, vegtype, sigmaf,  !
!           sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,           !
!           prsl1, prslki, zf, islimsk, shdmin, shdmax,                 ! 
!           snoalb, sfalb, flag_iter, flag_guess, isot, ivegsrc, fice,  !
! --- constants                                                         !
!           con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, con_fvirt, !
! --- in/outs                                                           !
!           weasd, snwdph, tskin, tprcp, rain, rainc, snow,             !
!           graupel, srflag, sr, smc, stc, slc, keepfr,                 !
!           canopy, trans, tsurf, tsnow, zorl,                          !
! --- ourputs                                                           !
!           sncovr1, qsurf, gflux, drain, dtsfc1, evap, dqsfc1, ep,     !
!           runoff, evbs, evcw, sbsno, snowc, stm, wet1,                !
!           errmsg, errflg                                              !
!         )                                                             !
!                                                                       !
!  program history log:                                                 !
!    may  2018  -- tanya smirnova                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horiz dimention and num of used pts      1    !
!     km       - integer, vertical soil layer dimension            9    !
!     ps       - real, surface pressure (pa)                       im   !
!     u1, v1   - real, u/v component of surface layer wind         im   !
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
!     islimsk  - integer, sea/land/ice mask (=0/1/2)               im   !
!     slopetyp - integer, class of sfc slope (integer index)       im   !
!     shdmin   - real, min fractional coverage of green veg        im   !
!     shdmax   - real, max fractnl cover of green veg (not used)   im   !
!     snoalb   - real, upper bound on max albedo over deep snow    im   !
!     sfalb    - real, mean sfc diffused sw albedo (fractional)    im   !
!     flag_iter- logical,                                          im   !
!     flag_guess-logical,                                          im   !
!     isot     - integer, sfc soil type data source zobler or statsgo   !
!     ivegsrc  - integer, sfc veg type data source umd or igbp          !
!     smc      - real, total soil moisture content (fractional)   im,km !
!                                                                       !
!  input/outputs:                                                       !
!     weasd    - real, water equivalent accumulated snow depth (mm) im  !
!     snwdph   - real, snow depth (water equiv) over land          im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     sr       - real, mixed-phase precipitation fraction          im   !
!     stc      - real, soil temp (k)                              im,km !
!     slc      - real, liquid soil moisture                       im,km !
!     canopy   - real, canopy moisture content (m)                 im   !
!     trans    - real, total plant transpiration (m/s)             im   !
!     tsurf    - real, surface skin temperature (after iteration)  im   !
!                                                                       !
!  outputs:                                                             !
!     sncovr1  - real, snow cover over land (fractional)           im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     drain    - real, subsurface runoff (m/s)                     im   !
!     evap     - real, latent heat flux in kg kg-1 m s-1           im   !
!     dtsfc1   - real, sensible heat flux in W m-2                 im   !
!     dqsfc1   - real, sensible heat flux in W m-2                 im   !
!     runof    - real, surface runoff (m/s)                        im   !
!     evbs     - real, direct soil evaporation (m/s)               im   !
!     evcw     - real, canopy water evaporation (m/s)              im   !
!     sbsno    - real, sublimation/deposit from snopack (m/s)      im   !
!     snowc    - real, fractional snow cover                       im   !
!     stm      - real, total soil column moisture content (m)      im   !
!     zorl     - real, surface roughness                           im   !
!     wet1     - real, normalized soil wetness                     im   !
!                                                                       !
!  ====================    end of description    =====================  !

!-----------------------------------
!      subroutine sfc_drv_ruc                                                &
! \defgroup RUC Surface Model - wrf4.0
!> \defgroup RUC_drv  RUC LSM Driver
!! \brief This is the RUC LSM driver module, with the functionality of 
!! preparing variables to run RUC LSM lsmruc(), calling lsmruc() and post-processing
!! variables for return to the parent model suite including unit conversion, as well 
!! as diagnotics calculation. 
#if 0
!! \section arg_table_lsm_ruc_run Argument Table
!! | local_name      | standard_name                                                                | long_name                                                       | units         | rank | type      |    kind   | intent | optional |
!|-------------------|------------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | delt            | time_step_for_dynamics                                                       | physics time step                                               | s             |    0 | real      | kind_phys | in     | F        |
!! | kdt             | index_of_time_step                                                           | current number of time steps                                    | index         |    0 | integer   |           | in     | F        |
!! | im              | horizontal_loop_extent                                                       | horizontal loop extent                                          | count         |    0 | integer   |           | in     | F        |
!! | km              | soil_vertical_dimension                                                      | soil vertical layer dimension                                   | count         |    0 | integer   |           | in     | F        |
!! | zs              | depth_of_soil_levels_for_land_surface_model                                  | depth of soil levels for lsm                                    | count         |    0 | real      | kind_phys | inout  | F        |
!! | con_cp          | specific_heat_of_dry_air_at_constant_pressure                                | specific heat !of dry air at constant pressure                  | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | con_g           | gravitational_acceleration                                                   | gravitational acceleration                                      | m s-2         |    0 | real      | kind_phys | in     | F        |
!! | con_pi          | pi                                                                           | ratio of a circle's circumference to its diameter               | radians       |    0 | real      | kind_phys | in     | F        |
!! | con_rd          | gas_constant_dry_air                                                         | ideal gas constant for dry air                                  | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | con_rv          | gas_constant_water_vapor                                                     | ideal gas constant for water vapor                              | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | con_hvap        | latent_heat_of_vaporization_of_water_at_0C    | latent heat of vaporization/sublimation (hvap)                  | J kg-1        |    0 | real      | kind_phys | in     | F        |      
!! | con_fvirt       | ratio_of_vapor_to_dry_air_gas_constants_minus_one                            | rv/rd - 1 (rv = ideal gas constant for water vapor)             | none          |    0 | real      | kind_phys | in     | F        |
!! | tprcp           | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep       | nonnegative precipitation amount in one dynamics time step      | m             |    1 | real      | kind_phys | inout  | F        | 
!! | rain            | lwe_thickness_of_precipitation_amount_on_dynamics_timestep                   | total rain at this time step                                    | m             |    1 | real      | kind_phys | in     | F        |
!! | rainc           | lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep        | convective rain at this time step                               | m             |    1 | real      | kind_phys | in     | F        |
!! | graupel         | lwe_thickness_of_graupel_amount_on_dynamics_timestep                         | graupel fall at this time step                                  | m             |    1 | real      | kind_phys | in     | F        |
!! | snow            | lwe_thickness_of_graupel_amount_on_dynamics_timestep                         | graupel fall at this time step                                  | m             |    1 | real      | kind_phys | in     | F        |
!! | sncovr1         | surface_snow_area_fraction_for_diagnostics                                   | surface snow area fraction                                      | frac          |    1 | real      | kind_phys | inout  | F        | 
!! | weasd           | water_equivalent_accumulated_snow_depth                                      | water equivalent accumulated snow depth                         | mm            |    1 | real      | kind_phys | inout  | F        |
!! | snwdph          | surface_snow_thickness_water_equivalent                                      | water equivalent snow depth over land                           | mm            |    1 | real      | kind_phys | inout  | F        | 
!! | sr              | ratio_of_snowfall_to_rainfall                                                | snow ratio: ratio of snow to total precipitation                | frac          |    1 | real      | kind_phys | in     | F        | 
!! | rhosnf          | density_of_frozen_precipitation                                              | density of frozen precipitation                                 | kg m-3        |    1 | real      | kind_phys | out    | F        |
!! | zf              | height_above_ground_at_lowest_model_layer                                    | layer 1 height above ground (not MSL)                           | m             |    1 | real      | kind_phys | in     | F        |
!! | u1              | x_wind_at_lowest_model_layer                                                 | x component of 1st model layer wind                             | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v1              | y_wind_at_lowest_model_layer                                                 | y component of 1st model layer wind                             | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | prsl1           | air_pressure_at_lowest_model_layer                                           | Model layer 1 mean pressure                                     | Pa            |    1 | real      | kind_phys | in     | F        | 
!! | t1              | air_temperature_at_lowest_model_layer                                        | 1st model layer air temperature                                 | K             |    1 | real      | kind_phys | in     | F        | 
!! | q1              | specific_humidity_at_lowest_model_layer                                      | 1st model layer specific humidity                               | kg kg-1       |    1 | real      | kind_phys | in     | F        | 
!! | qc              | cloud_condensed_water_concentration_at_lowest_model_layer                    | concentration of cloud water at lowest model layer              | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | dlwflx          | surface_downwelling_longwave_flux_on_radiation_time_step                     | total sky sfc downward lw flux                                  | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | dswsfc          | surface_downwelling_shortwave_flux                                           | total sky surface downward shortwave flux                       | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | snet            | surface_net_downwelling_shortwave_flux                                       | total sky surface net shortwave flux                            | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | sfcemis         | surface_longwave_emissivity                                                  | surface longwave emissivity                                     | frac          |    1 | real      | kind_phys | inout  | F        |
!! | cm              | surface_drag_coefficient_for_momentum_in_air                                 | surface exchange coeff for momentum                             | none          |    1 | real      | kind_phys | in     | F        |
!! | ch              | surface_drag_coefficient_for_heat_and_moisture_in_air                        | surface exchange coeff heat & moisture                          | none          |    1 | real      | kind_phys | in     | F        | 
!! | wet1            | normalized_soil_wetness                                                      | normalized soil wetness                                         | frac          |    1 | real      | kind_phys | inout  | F        | 
!! | canopy          | canopy_water_amount                                                          | canopy moisture content                                         | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | sigmaf          | vegetation_area_fraction                                                     | areal fractional cover of green vegetation                      | frac          |    1 | real      | kind_phys | in     | F        | 
!! | sfalb           | surface_diffused_shortwave_albedo                                            | mean surface diffused shortwave albedo                          | frac          |    1 | real      | kind_phys | inout  | F        |
!! | zorl            | surface_roughness_length                                                     | surface roughness length                                        | cm            |    1 | real      | kind_phys | inout  | F        | 
!! | snoalb          | upper_bound_on_max_albedo_over_deep_snow                                     | maximum snow albedo                                             | frac          |    1 | real      | kind_phys | in     | F        |
!! | qsurf           | surface_specific_humidity                                                    | surface specific humidity                                       | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcqs           | saturation_water_vapor_concentration_at_surface                              | saturation water vapor concentration at surface                 | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcqc           | cloud_condensed_water_concentration_at_surface                               | cloud water concentration at surface                            | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcdew          | surface_condensation_mass                                                    | mass of condensed water at surface                              | kg m-2        |    1 | real      | kind_phys | out    | F        |
!! | tg3             | deep_soil_temperature                                                        | bottom soil temperature                                         | K             |    1 | real      | kind_phys | in     | F        |
!! | vegtype         | cell_vegetation_type                                                         | vegetation type at each grid cell                               | index         |    1 | integer   |           | in     | F        |
!! | soiltyp         | cell_soil_type                                                               | soil type at each grid cell                                     | index         |    1 | integer   |           | in     | F        |
!! | fice            | sea_ice_concentration                                                        | ice fraction over open water                                    | frac          |    1 | real      | kind_phys | in     | F        |
!! | keepfr          | flag_for_frozen_soil_physics                                                 | flag for processes in frozen soil: 0, 1-limit on ice increase   | flag          |    2 | real      | kind_phys | inout  | F        |
!! | smc             | volume_fraction_of_soil_moisture_for_land_surface_model                      | volumetric fraction of soil moisture for lsm                    | frac          |    2 | real      | kind_phys | inout  | F        |
!! | slc             | volume_fraction_of_unfrozen_soil_moisture_for_land_surface_model             | volume fraction of unfrozen soil moisture for lsm               | frac          |    2 | real      | kind_phys | inout  | F        |
!! | stc             | soil_temperature_for_land_surface_model                                      | soil temperature for land surface model                         | K             |    2 | real      | kind_phys | inout  | F        |
!! | tsurf           | surface_skin_temperature                                                     | surface skin temperature                                        | K             |    1 | real      | kind_phys | inout  | F        |
!! | tsnow           | snow_temperature_bottom_first_layer                                          | snow temperature at the bottom of first snow layer              | K             |    1 | real      | kind_phys | inout  | F        |
!! | dqsfc1          | instantaneous_surface_upward_latent_heat_flux                                | surface upward latent heat flux                                 | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | dtsfc1          | instantaneous_surface_upward_sensible_heat_flux                              | surface upward sensible heat flux                               | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | evap            | kinematic_surface_upward_latent_heat_flux                                    | surface upward evaporation flux                                 | kg kg-1 m s-1 |    1 | real      | kind_phys | out    | F        |
!! | runof           | surface_runoff_flux                                                          | surface runoff flux                                             | g m-2 s-1     |    1 | real      | kind_phys | out    | F        |
!! | drain           | subsurface_runoff_flux                                                       | subsurface runoff flux                                          | g m-2 s-1     |    1 | real      | kind_phys | out    | F        |
!! | runoff          | total_runoff                                                                 | total water runoff                                              | kg m-2        |    1 | real      | kind_phys | none   | F        |
!! | gflux           | upward_heat_flux_in_soil                                                     | upward soil heat flux                                           | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | shdmin          | minimum_vegetation_area_fraction                                             | min fractional coverage of green veg                            | frac          |    1 | real      | kind_phys | in     | F        |
!! | shdmax          | maximum_vegetation_area_fraction                                             | max fractional coverage of green vegetation                     | frac          |    1 | real      | kind_phys | in     | F        |
!! | errmsg          | error_message                                                                | error message for error handling in CCPP                        | none          |    0 | character | len=*     | out    | F        |
!! | errflg          | error_flag                                                                   | error flag for error handling in CCPP                           | flag          |    0 | integer   |           | out    | F        |
!!
#endif
!!  \section general_ruc_drv RUC Driver General Algorithm
!!  @{
!  \section detailed_ruc RUC Driver Detailed Algorithm
!  @{


      subroutine lsm_ruc_run                                            &
! --- inputs
     &     ( kdt, im, km, zs,                                           &
     &       u1, v1, t1, q1, qc, soiltyp, vegtype, sigmaf,              &
     &       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
     &       prsl1, prslki, zf, islimsk, shdmin, shdmax,                &
     &       snoalb, sfalb, flag_iter, flag_guess, isot, ivegsrc, fice, &
! --- constants
     &       con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, con_fvirt,&
! --- in/outs
     &       weasd, snwdph, tskin, tprcp, rain, rainc, snow,            &
     &       graupel, srflag, sr, smc, stc, slc, keepfr,                &
     &       canopy, trans, tsurf, tsnow, zorl,                         &
     &       sfcqs, sfcqc, sfcdew,                                      &
! --- outputs
     &       sncovr1, qsurf, gflux, drain, dtsfc1, evap, dqsfc1, ep,    &
     &       runoff, evbs, evcw, sbsno, snowc, stm, wet1,               &
     &       errmsg, errflg                                             &
     &     )

      implicit none

!  ---  constant parameters:
      integer,              parameter :: ime=1
      integer,              parameter :: jme=1
      real(kind=kind_phys), parameter :: stbolt  = 5.670400e-8

!      real(kind=kind_phys), save         :: zsoil_ruc(km)
!      data zsoil_ruc / 0., 0.05, 0.1, 0.2, 0.4, 0.6, 1.0, 1.6, 3.0 /   ! in [m]

!  ---  input:
      integer, intent(in) :: im, km, kdt, isot, ivegsrc

      integer, dimension(im), intent(in) :: soiltyp, vegtype
      real (kind=kind_phys), dimension(km), intent(in   ) :: zs

      real (kind=kind_phys), dimension(im), intent(in) :: u1, v1,&
     &       t1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,    &
     &       ch, prsl1, prslki, shdmin, shdmax,                  &
     &       snoalb, sfalb, zf, fice, qc, q1 

      integer, dimension(im), intent(in) :: islimsk
      real (kind=kind_phys),  intent(in) :: delt
      real (kind=kind_phys),  intent(in) :: con_cp, con_rv, con_g,      &
                                            con_pi, con_rd,             &
                                            con_hvap, con_fvirt

      logical, dimension(im), intent(in) :: flag_iter, flag_guess

!  ---  in/out:
      real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
     &       snwdph, tskin, tprcp, rain, rainc, graupel, snow,          &
             srflag, sr, canopy, trans, tsurf, zorl, tsnow,             &
             sfcqs, sfcqc, sfcdew

      real (kind=kind_phys), dimension(im,km), intent(inout) ::         &
     &       smc, stc, slc, keepfr

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(inout) :: sncovr1,   &
     &       qsurf, gflux, dtsfc1, evap, dqsfc1, ep, runoff, drain,     &
     &       evbs, evcw, sbsno, snowc, stm, wet1

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rch, rho,                 &
     &       q0, qs1, theta1, weasd_old, snwdph_old,              &
     &       tprcp_old, srflag_old, sr_old, tskin_old, canopy_old

      real (kind=kind_phys), dimension(km) :: et, sldpth

      real (kind=kind_phys), dimension(1:ime,km,1:jme) :: stsoil, & 
     &       smsoil, slsoil, smfrsoil, keepfrsoil

      real (kind=kind_phys), dimension(im,km) :: zsoil, smc_old,        &
     &       stc_old, slc_old, keepfr_old

      real (kind=kind_phys),dimension (1:ime,1:jme) ::                  &
     &     alb, albedo, chs, flhc, flqc, wet, smmax, cmc,               &
     &     dew, drip,  ec, edir, ett, eta, esnow, etp, qfx,             &
     &     acceta, ffrozp, lwdn, prcp, q2, xland, xice,                 &
     &     graupelncv, snowncv, rainncv, raincv,                        &
     &     qcatm, solnet, sfcexc, conflx2,                              &
     &     runoff1, runoff2, acrunoff,sfcprs, sfctmp,                   &
     &     sfcems, sheat, shdfac, shdmin1d, shdmax1d,                   &
     &     sneqv, snoalb1d, snowh, snoh, tsnav,                         &
     &     snomlt, sncovr, soilw, soilm, ssoil, soilt, th2, tbot,       &
     &     xlai, swdn, z0, znt, rho2, rhosnf, infiltr,                  &
     &     precipfr, snowfallac, acsnow,                                &
     &     qsfc, qsg, qvg, qcg, soilt1, chklowq

      real (kind=kind_phys) :: xice_threshold, eps, epsm1, rhoh2o

      character(len=256) :: llanduse  ! Land-use dataset.  Valid values are :
                                      ! "USGS" (USGS 24/27 category dataset) and
                                      ! "MODIFIED_IGBP_MODIS_NOAH" (MODIS 20-category dataset)

      real (kind=kind_phys), dimension(21) :: landusef ! fractional landuse
      real (kind=kind_phys), dimension(19) ::  soilctop ! fractional soil type

      integer :: nsoil, iswater, isice
      integer, dimension (1:ime,1:jme) :: stype, vtype

      integer :: l, i, k, i1, j, nlcat, nscat, fractional_seaice

      logical :: flag(im)
      logical :: rdlai2d, myj, frpcpn
!
!===> ...  begin here
!
      real(kind=kind_phys) :: st_levels_input(4), sm_levels_input(4)
            st_levels_input = (/ 5, 25, 70, 150/)    ! Noah soil levels
            sm_levels_input = (/ 5, 25, 70, 150/)    ! Noah soil levels


      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!> - Set flag for land and ice points.
      do i = 1, im
        flag(i) = (islimsk(i) == 1 .or. islimsk(i) == 2)
      enddo

           chklowq = 1.
           eps   = con_rd/con_rv
           epsm1 = con_rd/con_rv-1.
!> -- number of soil categories          
          if(isot == 1) then
            nscat = 19 ! stasgo
          else
            nscat = 9  ! zobler
          endif
!> -- set parameters for IGBP land-use data
          if(ivegsrc == 2) then
            llanduse = 'MODI-RUC'  ! IGBP
            nlcat = 20  ! IGBP - "MODI-RUC"
            iswater = 17
            isice = 15
          endif

          fractional_seaice = 1
        if ( fractional_seaice == 0 ) then
          xice_threshold = 0.5
        else if ( fractional_seaice == 1 ) then
          xice_threshold = 0.02
        endif

!> - Save land-related prognostic fields for guess run.

      do i = 1, im
        if (flag(i) .and. flag_guess(i)) then
          weasd_old(i)  = weasd(i)
          snwdph_old(i) = snwdph(i)
          tskin_old(i)  = tskin(i)
          canopy_old(i) = canopy(i)
          tprcp_old(i)  = tprcp(i)
          srflag_old(i) = srflag(i)
          sr_old(i)     = sr(i)

          do k = 1, km
            smc_old(i,k) = smc(i,k)
            stc_old(i,k) = stc(i,k)
            slc_old(i,k) = slc(i,k)
            keepfr_old(i,k) = keepfr(i,k)
          enddo
        endif
      enddo

!  --- ...  initialization block

      do i = 1, im
        if (flag_iter(i)) then
          ep(i)     = 0.0
          evap (i)  = 0.0
          dqsfc1 (i)= 0.0
          gflux(i)  = 0.0
          drain(i)  = 0.0
          canopy(i) = max(canopy(i), 0.0)

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno(i)  = 0.0
          snowc(i)  = 0.0
        endif
      enddo

!  --- ...  initialize variables

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          q0(i)   = max(q1(i)/(1.-q1(i)), 1.e-8)   !* q1=specific humidity at level 1 (kg/kg)
          theta1(i) = t1(i) * prslki(i) !* adiabatic temp at level 1 (k)

          rho(i) = prsl1(i) / (con_rd*t1(i)*(1.0+con_fvirt*q0(i)))
          qs1(i) = fpvs( t1(i) )        !* qs1=sat. humidity at level 1 (kg/kg)
          qs1(i) = max(eps*qs1(i) / (prsl1(i)+epsm1*qs1(i)), 1.e-8)
          q0 (i) = min(qs1(i), q0(i))
        endif
      enddo

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          do k = 1, km
            zsoil(i,k) = zs(k)   ! [m]
          enddo
        endif
      enddo

      do i1 = 1,ime
      do j  = 1,jme

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then

!> - Prepare variables to run RUC LSM: 
!!  -   1. configuration information (c):
!!\n  ----------------------------------------
!!\n  \a ffrozp  - fraction of frozen precipitation
!!\n  \a frpcpn  - .true. if mixed phase precipitation available
!!\n  \a ids,ide  - horizontal_loop_extent
!!\n  \a jds,jde  - horizontal_loop_extent
!!\n  \a fice    - fraction of sea-ice in the grid cell
!!\n  \a delt    - timestep (sec) (dt should not exceed 3600 secs) 
!!\n  \a zlvl    - height (\f$m\f$) above ground of atmospheric forcing variables
!!\n  \a nsoil   - number of soil layers (= 6 or 9)
!!\n  \a sldpth  - the depth of each soil level (\f$m\f$)

          frpcpn = .true.                 ! .true. if mixed phase precipitation available (Thompson)

        if(.not.frpcpn) then ! no mixed-phase precipitation available
          if     (srflag(i) == 1.0) then  ! snow phase
            ffrozp(i1,j) = 1.0
          elseif (srflag(i) == 0.0) then  ! rain phase
            ffrozp(i1,j) = 0.0
          endif
        else ! mixed-phase precipitation is available
            ffrozp(i1,j) = sr(i)
        endif ! frpcpn

! zf(i)  - first atm. level above ground surface 
          conflx2 = zf(i) * 2. ! multiplied by 2., because inside RUC LSM surface layer CONFLX=0.5*zlvl

          nsoil = km
          sldpth(1) =  zsoil(i,1)
          do k = 2, km
            sldpth(k) = zsoil(i,k) - zsoil(i,k-1) ! [m]
          enddo

!>  -   2. forcing data (f):
!!\n  ---------------------------------------
!!\n  \a lwdn    - lw dw radiation flux at surface (\f$W m^{-2}\f$)
!!\n  \a swdn    - sw dw radiation flux at surface (\f$W m^{-2}\f$)
!!\n  \a solnet  - net sw radiation flux (dn-up) (\f$W m^{-2}\f$)
!!\n  \a sfcprs  - pressure at height zlvl above ground (pascals)
!!\n  \a prcp    - time-step total precip (\f$kg m^{-2} \f$)
!!\n  \a raincv  - time-step convective precip (\f$kg m^{-2} \f$)
!!\n  \a rainncv - time-step non-convective precip (\f$kg m^{-2} \f$)
!!\n  \a graupelncv - time-step graupel (\f$kg m^{-2} \f$)
!!\n  \a snowncv - time-step snow (\f$kg m^{-2} \f$)
!!\n  \a sfctmp  - air temperature (\f$K\f$) at height zlvl above ground
!!\n  \a th2     - air potential temperature (\f$K\f$) at height zlvl above ground
!!\n  \a q2      - mixing ratio at height zlvl above ground (\f$kg kg^{-1}\f$)
!!\n  \a qcatm   - cloud water at the first atm. level

          lwdn   = dlwflx(i)         !..downward lw flux at sfc in w/m2
          swdn   = dswsfc(i)         !..downward sw flux at sfc in w/m2
          solnet = snet(i)           !..net sw rad flx (dn-up) at sfc in w/m2
          sfcems = sfcemis(i)

          sfcprs = prsl1(i)

          prcp       = rhoh2o * tprcp(i) 
          rainncv    = rhoh2o * rain(i) 
          raincv     = rhoh2o * rainc(i)
          graupelncv = rhoh2o * graupel(i)
          snowncv    = rhoh2o * snow(i)
          precipfr   = rainncv * ffrozp
          sfctmp = t1(i)
          th2    = theta1(i)
          q2     = q0(i)
          qcatm  = qc(i)
          rho2   = rho(i)
          qsfc   = qsurf(i)
          qvg    = qsurf(i)/(1.-qsurf(i))

!>  -   3. canopy/soil characteristics (s):
!!\n      ------------------------------------
!!\n \a vegtyp  - vegetation type (integer index)                   -> vtype
!!\n \a soiltyp - soil type (integer index)                         -> stype
!!\n \a shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
!!\n \a shdmin  - minimum areal fractional coverage of green vegetation -> shdmin1d
!!\n \a shdmax  - maximum areal fractional coverage of green vegetation -> shdmax1d
!!\n \a alb     - backround snow-free surface albedo (fraction)
!!\n \a snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d
!!\n \a tbot    - bottom soil temperature (local yearly-mean sfc air temp)

      if(ivegsrc == 1) then   ! IGBP - MODIS
!> - Prepare land/ice/water masks for RUC LSM
!SLMSK   - SEA(0),LAND(1),ICE(2) MASK
       if(islimsk(i) == 0.) then
           vtype = 17 ! 17 - water (oceans and lakes) in MODIS
           stype = 14
           xland = 2. ! xland = 2 for water
           xice  = 0.
!           landmask = 0.
       elseif(islimsk(i) == 1.) then ! land
           vtype = vegtype(i)
           stype = soiltyp(i)
           xland  = 1.
           xice   = 0.
!           landmask = 1.
       elseif(islimsk(i) == 2) then  ! ice
           vtype = 15 ! MODIS
          if(isot == 0) then
              stype = 9  ! ZOBLER
          else
              stype = 16 ! STASGO
          endif
           xland    = 1.
!           landmask = 1.
!tgs - check the name of sea ice fraction
           xice = fice(i)  ! fraction of sea-ice
       endif
      else
        print *,'MODIS landuse is not available'
      endif

! --- units %
          shdfac = sigmaf(i)*100.
          shdmin1d = shdmin(i)*100.
          shdmax1d = shdmax(i)*100.

          snoalb1d = snoalb(i)
          alb  = sfalb(i)
          tbot = tg3(i)

!tgs - for now set rdlai2d to .false., WRF has LAI maps, and RUC LSM
!      uses rdlai2d = .true.
!          rdlai2d = .false.
!       if( .not. rdlai2d) xlai = lai_data(vtype)
!
!>  -   4. history (state) variables (h):
!!\n      ------------------------------
!!\n \a cmc        - canopy moisture content (\f$m\f$)
!!\n \a t1         - ground/canopy/snowpack effective skin temperature (\f$K\f$)   -> tskin
!!\n \a stc(nsoil) - soil temp (\f$K\f$)                                         -> stsoil
!!\n \a smc(nsoil) - total soil moisture content (volumetric fraction)     -> smsoil
!!\n \a sh2o(nsoil)- unfrozen soil moisture content (volumetric fraction)  -> slsoil
!!\n \a smfrsoil(nsoil)- frozen soil moisture content (volumetric fraction)  -> smfrsoil
!!\n \a keepfrflag(nsoil) - flag for frozen soil physics: 0. or 1.
!!\n \a snowh      - actual snow depth (\f$m\f$)
!!\n \a sneqv      - liquid water-equivalent snow depth (\f$m\f$)
!!\n \a albedo     - surface albedo including snow effect (unitless fraction)
!!\n \a ch         - surface exchange coefficient for heat and moisture (\f$m s^{-1}\f$) -> chx
!!\n \a cm         - surface exchange coefficient for momentum (\f$m s^{-1}\f$)          -> cmx
!!\n \a z0         - surface roughness (\f$m\f$)     -> zorl(\f$cm\f$)
!!\n \a rnet       - Residual of the surface energy balance equation (\f$ w m^{-2} \f$)

          cmc = canopy(i)            !  [mm] 
          soilt = tsurf(i)            ! clu_q2m_iter
          soilt1 = tsnow(i)

          do k = 1, km
            stsoil(:,k,:) = stc(i,k)
            smsoil(:,k,:) = smc(i,k)
            slsoil(:,k,:) = slc(i,k)
            smfrsoil(:,k,:) = (smsoil(i1,k,j)-slsoil(i1,k,j))/0.9
            keepfrsoil(:,k,:) = keepfr(i,k)
          enddo

          wet = wet1(i)

          snowh = snwdph(i) * 0.001         ! convert from mm to m
          sneqv = weasd(i)                  ! [mm]

!> -- sanity checks on sneqv and snowh
          if (sneqv(i1,j) /= 0.0 .and. snowh(i1,j) == 0.0) then
             snowh(i1,j) = 0.003 * sneqv(i1,j) ! snow density ~300 kg m-3 
          endif

          if (snowh(i1,j) /= 0.0 .and. sneqv(i1,j) == 0.0) then
             sneqv(i1,j) = 300. * snowh(i1,j) ! snow density ~300 kg m-3 
          endif

          if (sneqv(i1,j) > 0. .and. snowh(i1,j) > 0.) then
            if(sneqv(i1,j)/snowh(i1,j) > 950.) then
              sneqv(i1,j) = 300. * snowh(i1,j)
            endif
          endif

!> -- for now set fractions of differnet landuse and soil types in the grid
!     cell to zero
          landusef(:) = 0.0
          soilctop(:) = 0.0

          sncovr = sncovr1(i)

          chs    = ch(i)  
          flhc   = ch(i) * rho(i) * con_cp * (1. + 0.84*q2(i1,j))
          flqc   = ch(i) * rho(i) * wet(i1,j)

!  ---- ... outside sflx, roughness uses cm as unit
          z0  = zorl(i)/100.
          znt = zorl(i)/100.

!> - Call RUC LSM lsmruc(). 

      call lsmruc( delt, kdt, nsoil,                             &
     &          graupelncv, snowncv, rainncv, raincv,                   &
     &          sldpth, prcp, sneqv, snowh, sncovr, ffrozp, frpcpn,             &
     &          rhosnf, precipfr,                                       &
!  ---  inputs:
     &          conflx2, sfcprs, sfctmp, q2, qcatm, rho2,               &
     &          lwdn, solnet, sfcems, chklowq, chs, flqc, flhc,                  &
!  ---  input/outputs:
     &          wet, cmc, shdfac, albedo, znt,                      &
     &          z0, snoalb1d, alb,                                      &
!     &          z0, snoalb1d, alb, xlai,                                &
     &          landusef, nlcat,                              &
!  --- mosaic_lu and mosaic_soil are moved to the namelist
!     &          mosaic_lu, mosaic_soil,                                 &
     &          soilctop, nscat,                                        &
     &          qsfc, qsg, qvg, qcg, dew, soilt1,                       &
     &          tsnav, tbot, vtype, stype, xland,                              &
     &          iswater, isice, xice, xice_threshold,                   &
!  ---  constants
     &          con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, stbolt,&
!  ---  input/outputs:
     &          smsoil, slsoil, soilm, smmax,                           &
     &          stsoil, soilt, sheat, qfx, eta, infiltr,                        &
     &          runoff1, runoff2, acrunoff, sfcexc,                     &
     &          acceta, ssoil, snowfallac, acsnow, snomlt,        &
     &          smfrsoil,keepfrsoil, .false.,                           &
     &          shdmin1d, shdmax1d, rdlai2d,                            &
     &          1,1, 1,1, 1,1, 1,1, 1,1, 1,1 )


!> - RUC LSM: prepare variables for return to parent model and unit conversion.
!>  -   6. output (o):
!!\n  ------------------------------
!!\n \a eta     - actual latent heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a sheat   - sensible heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a ssoil   - soil heat flux (\f$W m^{-2}\f$: negative if downward from surface)
!!\n \a runoff1 - surface runoff (\f$m s^{-1}\f$), not infiltrating the surface
!!\n \a runoff2 - subsurface runoff (\f$m s^{-1}\f$), drainage out bottom
!!\n \a snoh    - phase-change heat flux from snowmelt (w m-2)
!

          evap(i)   = qfx(i1,j)
          dqsfc1(i) = eta(i1,j)
          dtsfc1(i) = sheat(i1,j)
          gflux(i)  = ssoil(i1,j)

!          evbs(i)  = edir(i1,j)
!          evcw(i)  = ec(i1,j)
!          trans(i) = ett(i1,j)
!          sbsno(i) = esnow(i1,j)
!          snohf(i) = snoh(i1,j)

          sfcqs(i)  = qsg(i1,j)
          sfcqc(i)  = qcg(i1,j)
          sfcdew(i) = dew(i1,j)

          snowc(i) = sncovr(i1,j)
          stm(i)   = soilm(i1,j)

          tsurf(i)   = soilt(i1,j)

          do k = 1, km
            stc(i,k) = stsoil(i1,k,j)
            smc(i,k) = smsoil(i1,k,j)
            slc(i,k) = slsoil(i1,k,j)
            keepfr(i,k) = keepfrsoil(i1,k,j)
          enddo

          wet1(i) = wet(i1,j) 

!  --- ...  units [m/s] = [g m-2 s-1] 
          runoff (i)  = runoff1(i1,j) 
          drain (i)  = runoff2(i1,j)

!  --- ...  unit conversion (from m to mm)
          snwdph(i)  = snowh(i1,j) * 1000.0

          canopy(i)  = cmc(i1,j)   ! mm
          weasd(i)   = sneqv(i1,j) ! mm
          sncovr1(i) = sncovr(i1,j)
!  ---- ... outside sflx, roughness uses cm as unit (update after snow's
!  effect)
          zorl(i) = znt(i1,j)*100.

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
      enddo   ! end do_i_loop

!> - Compute specific humidity at surface (\a qsurf).

      do i = 1, im
        if (flag_iter(i)  .and. flag(i)) then
!          rch(i)   = rho(i) * cp * ch(i) * wind(i)
!          qsurf(i) = q1(i)  + evap(i) / (elocp * rch(i))
          qsurf(i) = qsfc(i1,j)
        endif
      enddo

!> - Compute surface upward sensible heat flux (\a hflx) and evaporation
!! flux (\a evap). 
!      do i = 1, im
!        if (flag_iter(i) .and. flag(i)) then
!          tem     = 1.0 / rho(i)
!          hflx(i) = hflx(i) * tem * cpinv
!          evap(i) = evap(i) * tem * hvapi
!        endif
!      enddo

      enddo  ! i1
      enddo  ! j

!> - Restore land-related prognostic fields for guess run.

      do i = 1, im
        if (flag(i)) then
          if (flag_guess(i)) then
            weasd(i)  = weasd_old(i)
            snwdph(i) = snwdph_old(i)
            tskin(i)  = tskin_old(i)
            canopy(i) = canopy_old(i)
            tprcp(i)  = tprcp_old(i)
            srflag(i) = srflag_old(i)

            do k = 1, km
              smc(i,k) = smc_old(i,k)
              stc(i,k) = stc_old(i,k)
              slc(i,k) = slc_old(i,k)
            enddo
          else
            tskin(i) = tsurf(i)
          endif
        endif
      enddo
!
      return
!...................................
      end subroutine lsm_ruc_run
!-----------------------------------
!! @}

      end module lsm_ruc
