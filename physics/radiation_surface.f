!>  \file radiation_surface.f
!!  This file contains routines that set up surface albedo for SW
!!  radiation and surface emissivity for LW radiation.

!  ==========================================================  !!!!!
!            'module_radiation_surface' description            !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!    this module sets up surface albedo for sw radiation and surface   !
!    emissivity for lw radiation.                                      !
!                                                                      !
!                                                                      !
!    in the module, the externally callabe subroutines are :           !
!                                                                      !
!      'sfc_init'   -- initialization radiation surface data           !
!         inputs:                                                      !
!           ( me )                                                     !
!         outputs:                                                     !
!           (none)                                                     !
!                                                                      !
!      'setalb'     -- set up four-component surface albedoes          !
!         inputs:                                                      !
!           (slmsk,snowf,sncovr,snoalb,zorlf,coszf,tsknf,tairf,hprif,  !
!            alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc            !
!            IMAX)                                                     !
!         outputs:                                                     !
!           (sfcalb)                                                   !
!                                                                      !
!      'setemis'    -- set up surface emissivity for lw radiation      !
!         inputs:                                                      !
!           (xlon,xlat,slmsk,snowf,sncovr,zorlf,tsknf,tairf,hprif,     !
!            IMAX)                                                     !
!         outputs:                                                     !
!           (sfcemis)                                                  !
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module machine'             in 'machine.f'                    !
!       'module physcons'            in 'physcons.f'                   !
!       'module module_iounitdef'    in 'iounitdef.f'                  !
!                                                                      !
!                                                                      !
!    program history log:                                              !
!           1995   y.t. hou     - created albaer.f (include albedo     !
!                                 and aerosols calculations)           !
!      nov  1997   y.t. hou     - modified snow albedo                 !
!      jan  1998   y.t. hou     - included grumbine's sea-ice scheme   !
!      feb  1998   h.l. pan     - seasonal interpolation in cycle      !
!      mar  2000   y.t. hou     - modified to use opac aerosol data    !
!      apr  2003   y.t. hou     - seperate albedo and aerosols into    !
!                    two subroutines, rewritten in f90 modulized form  !
!      jan  2005   s. moorthi   - xingren's sea-ice fraction added     !
!      apr  2005   y.t. hou     - revised module structure             !
!      feb  2006   y.t. hou     - add varying surface emissivity,      !
!                    modified sfc albedo structure for modis shceme    !
!      Mar  2006   s. moorthi   - added surface temp over ice fraction !
!      mar  2007   c. marshall & h. wei                                !
!                               - added modis based sfc albedo scheme  !
!      may  2007   y. hou & s. moorthi                                 !
!                               - fix bug in modis albedo over ocean   !
!      aug  2007   h. wei & s. moorthi                                 !
!                               - fix bug in modis albedo over sea-ice !
!      aug  2007   y. hou       - fix bug in emissivity over ocean in  !
!                                 the modis scheme option              !
!      dec  2008   f. yang      - modified zenith angle dependence on  !
!                                 surface albedo over land. (2008 jamc)!
!      aug  2012   y. hou       - minor modification in initialization !
!                                 subr 'sfc_init'.                     !
!      nov  2012   y. hou       - modified control parameters through  !
!                    module 'physparam'.                               !
!      jun  2018   h-m lin/y-t hou - correct error in clim-scheme of   !
!                    weak/strong factor and restore to the orig form   !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!


!>\ingroup RRTMG
!> \defgroup module_radiation_surface RRTMG Surface Module
!> This module sets up surface albedo for SW radiation and surface
!! emissivity for LW radiation.
!!
!! In the module, the externally callable subroutines are :
!! - sfc_init(): initialization radiation surface data
!! - setalb(): set up four-component surface albedoes
!! - setemis(): set up surface emissivity for lw radiation
!!
!! SW surface albedo (namelist control parameter - \b IALB=1)
!!\n IALB=1: MODIS retrievals based monthly mean climatology
!!\n IALB=2: use surface albedo from land model
!!
!! LW surface emissivity (namelist control parameter - \b IEMS=1)
!!\n IEMS=1: surface type based climatology in \f$1^o\f$ horizontal resolution
!!\n IEMS=2: use surface emissivity from land model
!!
!!\version NCEP-Radiation_surface   v5.1  Nov 2012

!> This module sets up surface albedo for SW radiation and surface  
!! emissivity for LW radiation.  
      module module_radiation_surface   
!
!! \section arg_table_module_radiation_surface
!! \htmlinclude module_radiation_surface.html
!!
      use physparam,         only : ialbflg, iemsflg, semis_file,       &
     &                              kind_phys
      use physcons,          only : con_t0c, con_ttp, con_pi, con_tice
      use module_iounitdef,  only : NIRADSF
      use surface_perturbation, only : ppfbet
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGSFC='NCEP-Radiation_surface   v5.1  Nov 2012 '
!    &   VTAGSFC='NCEP-Radiation_surface   v5.0  Aug 2012 '

!  ---  constant parameters
      integer, parameter, public :: NF_ALBD = 4     !< number of surface albedo components
      integer, parameter, public :: IMXEMS = 360    !< number of longtitude points in global emis-type map
      integer, parameter, public :: JMXEMS = 180    !< number of latitude points in global emis-type map
      real (kind=kind_phys), parameter :: f_zero = 0.0
      real (kind=kind_phys), parameter :: f_one  = 1.0
      real (kind=kind_phys), parameter :: epsln  = 1.0e-6
      real (kind=kind_phys), parameter :: rad2dg= 180.0 / con_pi
      integer, allocatable  ::  idxems(:,:)         !< global surface emissivity index array
      integer :: iemslw = 1                         !< global surface emissivity control flag set up in 'sfc_init'
!
      public  sfc_init, setalb, setemis
      public  f_zero, f_one, epsln

! =================
      contains
! =================

!> \ingroup module_radiation_surface
!> This subroutine is the initialization program for surface radiation
!! related quantities (albedo, emissivity, etc.)
!!\param me       print control flag
!>\section gen_sfc_init sfc_init General Algorithm
!! @{
!-----------------------------------
      subroutine sfc_init                                               &
     &     ( me, errmsg, errflg )!  ---  inputs/outputs:
!
!  ===================================================================  !
!                                                                       !
!  this program is the initialization program for surface radiation     !
!  related quantities (albedo, emissivity, etc.)                        !
!                                                                       !
! usage:         call sfc_init                                          !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!      me           - print control flag                                !
!                                                                       !
!  outputs: (none) to module variables only                             !
!                                                                       !
!  external module variables:                                           !
!     ialbflg       - control flag for surface albedo schemes           !
!                     =1: use modis based surface albedo                !
!                     =2: use surface albedo from land model            !
!     iemsflg       - control flag for sfc emissivity schemes (ab:2-dig)!
!                     a:=0 set sfc air/ground t same for lw radiation   !
!                       =1 set sfc air/ground t diff for lw radiation   !
!                     b:=1 use varying climtology sfc emiss (veg based) !
!                       =2 use surface emissivity from land model       !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: me

!  ---  outputs: ( none )
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      integer    :: i, k
!     integer    :: ia, ja
      logical    :: file_exist
      character  :: cline*80
!
!===> ...  begin here
!
      errmsg = ''
      errflg = 0
!
      if ( me == 0 ) print *, VTAGSFC   ! print out version tag

!> - Initialization of surface albedo section
!! \n physparam::ialbflg
!!  - = 1: using MODIS based land surface albedo for SW
!!  - = 2: using albedo from land model

      if ( ialbflg == 1 ) then

        if ( me == 0 ) then
          print *,' - Using MODIS based land surface albedo for sw'
        endif

      elseif ( ialbflg == 2 ) then      ! use albedo from land model

        if ( me == 0 ) then
          print *,' - Using Albedo From Land Model'
        endif

      else

        errmsg = 'module_radiation_surface: invalid ialbflg option'
        errflg = 1
        return

      endif    ! end if_ialbflg_block

!> - Initialization of surface emissivity section
!! \n physparam::iemsflg
!!  - = 1: input SFC emissivity type map from "semis_file"
!!  - = 2: input SFC emissivity from land model

      iemslw = mod(iemsflg, 10)          ! emissivity control

      if ( iemslw == 1 ) then        ! input sfc emiss type map

!  ---  allocate data space
        if ( .not. allocated(idxems) ) then
          allocate ( idxems(IMXEMS,JMXEMS) )
        endif

!  ---  check to see if requested emissivity data file existed

        inquire (file=semis_file, exist=file_exist)

        if ( .not. file_exist ) then
          if ( me == 0 ) then
            print *,' - Using Varying Surface Emissivity for lw'
            print *,'   Requested data file "',semis_file,'" not found!'
            print *,'   Change to fixed surface emissivity = 1.0 !'
          endif

          iemslw = 0
        else
          close(NIRADSF)
          open (NIRADSF,file=semis_file,form='formatted',status='old')
          rewind NIRADSF

          read (NIRADSF,12) cline
  12      format(a80)

          read (NIRADSF,14) idxems
  14      format(80i1)

          if ( me == 0 ) then
            print *,' - Using Varying Surface Emissivity for lw'
            print *,'   Opened data file: ',semis_file
            print *, cline
!check      print *,' CHECK: Sample emissivity index data'
!           ia = IMXEMS / 5
!           ja = JMXEMS / 5
!           print *, idxems(1:IMXEMS:ia,1:JMXEMS:ja)
          endif

          close(NIRADSF)
        endif    ! end if_file_exist_block

      elseif ( iemslw == 2 ) then        ! use emiss from land model

        if ( me == 0 ) then
          print *,' - Using Surface Emissivity From Land Model'
        endif

      else

         errmsg = 'module_radiation_surface: invalid iemslw option'
         errflg = 1
         return

      endif   ! end if_iemslw_block

!
      return
!...................................
      end subroutine sfc_init
!-----------------------------------
!! @}

!> \ingroup module_radiation_surface
!> This subroutine computes four components of surface albedos (i.e.,
!! vis-nir, direct-diffused) according to control flag ialbflg.
!! \n 1) climatological surface albedo scheme (\cite briegleb_1992)
!! \n 2) MODIS retrieval based scheme from Boston univ.
!!\param slmsk      (IMAX), sea(0),land(1),ice(2) mask on fcst model grid
!!\param snowf      (IMAX), snow depth water equivalent in mm
!!\param sncovr     (IMAX), snow cover over land
!!\param snoalb     (IMAX), maximum snow albedo over land (for deep snow)
!!\param zorlf      (IMAX), surface roughness in cm
!!\param coszf      (IMAX), cosin of solar zenith angle
!!\param tsknf      (IMAX), ground surface temperature in K
!!\param tairf      (IMAX), lowest model layer air temperature in K
!!\param hprif      (IMAX), topographic sdv in m
!!\n ---  for ialbflg=0 climtological albedo scheme  ---
!!\param alvsf      (IMAX), 60 degree vis albedo with strong cosz dependency
!!\param alnsf      (IMAX), 60 degree nir albedo with strong cosz dependency
!!\param alvwf      (IMAX), 60 degree vis albedo with weak cosz dependency
!!\param alnwf      (IMAX), 60 degree nir albedo with weak cosz dependency
!!\n ---  for ialbflg=1 MODIS based land albedo scheme ---
!!\param alvsf      (IMAX), visible black sky albedo at zenith 60 degree
!!\param alnsf      (IMAX), near-ir black sky albedo at zenith 60 degree
!!\param alvwf      (IMAX), visible white sky albedo
!!\param alnwf      (IMAX), near-ir white sky albedo
!!\param facsf      (IMAX), fractional coverage with strong cosz dependency
!!\param facwf      (IMAX), fractional coverage with weak cosz dependency
!!\param fice       (IMAX), sea-ice fraction
!!\param tisfc      (IMAX), sea-ice surface temperature
!!\param IMAX       array horizontal dimension
!!\param albppert   (IMAX), a probability value in the interval [0,1]
!!\param pertalb    (5), magnitude of perturbation of surface albedo
!!\param sfcalb     (IMAX,NF_ALBD), mean sfc albedo
!!\n                    ( :, 1) -     near ir direct beam albedo
!!\n                    ( :, 2) -     near ir diffused albedo
!!\n                    ( :, 3) -     uv+vis direct beam albedo
!!\n                    ( :, 4) -     uv+vis diffused albedo
!>\section general_setalb setalb General Algorithm
!! @{
!-----------------------------------
      subroutine setalb                                                 &
     &     ( slmsk,lsm,lsm_noahmp,lsm_ruc,snowf,                        & !  ---  inputs:
     &       sncovr,sncovr_ice,snoalb,zorlf,coszf,                      &
     &       tsknf,tairf,hprif,landfrac,frac_grid,min_seaice,           & 
     &       alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,            &
     &       lsmalbdvis, lsmalbdnir, lsmalbivis, lsmalbinir,            &
     &       icealbdvis, icealbdnir, icealbivis, icealbinir,            &
     &       IMAX, albPpert, pertalb, fracl, fraco, fraci, icy,         &
     &       sfcalb                                                     & !  ---  outputs:
     &     )

!  ===================================================================  !
!                                                                       !
!  this program computes four components of surface albedos (i.e.       !
!  vis-nir, direct-diffused) according to controflag ialbflg.           !
!   1) climatological surface albedo scheme (briegleb 1992)             !
!   2) modis retrieval based scheme from boston univ.                   !
!                                                                       !
!                                                                       !
! usage:         call setalb                                            !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     slmsk (IMAX)  - sea(0),land(1),ice(2) mask on fcst model grid     !
!     snowf (IMAX)  - snow depth water equivalent in mm                 !
!     sncovr(IMAX)  - ialgflg=0: not used                               !
!                     ialgflg=1: snow cover over land in fraction       !
!     sncovr_ice(IMAX)  - ialgflg=0: not used                           !
!                     ialgflg=1: snow cover over ice in fraction        !
!     snoalb(IMAX)  - ialbflg=0: not used                               !
!                     ialgflg=1: max snow albedo over land in fraction  !
!     zorlf (IMAX)  - surface roughness in cm                           !
!     coszf (IMAX)  - cosin of solar zenith angle                       !
!     tsknf (IMAX)  - ground surface temperature in k                   !
!     tairf (IMAX)  - lowest model layer air temperature in k           !
!     hprif (IMAX)  - topographic sdv in m                              !
!           ---  for ialbflg=0 climtological albedo scheme  ---         !
!     alvsf (IMAX)  - 60 degree vis albedo with strong cosz dependency  !
!     alnsf (IMAX)  - 60 degree nir albedo with strong cosz dependency  !
!     alvwf (IMAX)  - 60 degree vis albedo with weak cosz dependency    !
!     alnwf (IMAX)  - 60 degree nir albedo with weak cosz dependency    !
!           ---  for ialbflg=1 modis based land albedo scheme ---       !
!     alvsf (IMAX)  - visible black sky albedo at zenith 60 degree      !
!     alnsf (IMAX)  - near-ir black sky albedo at zenith 60 degree      !
!     alvwf (IMAX)  - visible white sky albedo                          !
!     alnwf (IMAX)  - near-ir white sky albedo                          !
!                                                                       !
!     facsf (IMAX)  - fractional coverage with strong cosz dependency   !
!     facwf (IMAX)  - fractional coverage with weak cosz dependency     !
!     fice  (IMAX)  - sea-ice fraction                                  !
!     tisfc (IMAX)  - sea-ice surface temperature                       !
!     IMAX          - array horizontal dimension                        !
!                                                                       !
!  outputs:                                                             !
!     sfcalb(IMAX,NF_ALBD)                                              !
!           ( :, 1) -     near ir direct beam albedo                    !
!           ( :, 2) -     near ir diffused albedo                       !
!           ( :, 3) -     uv+vis direct beam albedo                     !
!           ( :, 4) -     uv+vis diffused albedo                        !
!                                                                       !
!  module internal control variables:                                   !
!     ialbflg       - =0 use the default climatology surface albedo     !
!                     =1 use modis retrieved albedo and input snow cover!
!                        for land areas                                 !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer, intent(in) :: IMAX
      integer, intent(in) :: lsm, lsm_noahmp, lsm_ruc
      logical, intent(in) :: frac_grid

      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       slmsk, snowf, zorlf, coszf, tsknf, tairf, hprif, landfrac, &
     &       alvsf, alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc,     &
     &       lsmalbdvis, lsmalbdnir, lsmalbivis, lsmalbinir,            &
     &       icealbdvis, icealbdnir, icealbivis, icealbinir,            &
     &       sncovr, sncovr_ice, snoalb, albPpert           ! sfc-perts, mgehne
      real (kind=kind_phys),  intent(in) :: pertalb         ! sfc-perts, mgehne
      real (kind=kind_phys),  intent(in) :: min_seaice
      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       fracl, fraco, fraci
      logical, dimension(:), intent(in) ::                              &
     &       icy

!  ---  outputs
      real (kind=kind_phys), dimension(IMAX,NF_ALBD), intent(out) ::    &
     &       sfcalb

!  ---  locals:
      real (kind=kind_phys) :: asnvb, asnnb, asnvd, asnnd, asevb        &
     &,     asenb, asevd, asend, fsno,  fsea,  rfcs,  rfcw,  flnd       &
     &,     asnow, argh,  hrgh,  fsno0, fsno1, flnd0, fsea0, csnow      &
     &,     a1, a2, b1, b2, b3, ab1bm, ab2bm, m, s, alpha, beta, albtmp

      real (kind=kind_phys) :: asevb_wat,asenb_wat,asevd_wat,asend_wat, &
     &                         asevb_ice,asenb_ice,asevd_ice,asend_ice

      real (kind=kind_phys) :: alndnb, alndnd, alndvb, alndvd

      real (kind=kind_phys) ffw, dtgd

      integer :: i, k, kk, iflag

!
!===> ...  begin here
!
!> - Use modis based albedo for land area:
      if ( ialbflg == 1 ) then

        do i = 1, IMAX

          !-- water albedo
          asevd_wat = 0.06
          asend_wat = 0.06
          asevb_wat = asevd_wat
          asenb_wat = asevd_wat

          ! direct albedo CZA dependence over water
          if (fraco(i) > f_zero .and. coszf(i) > 0.0001) then
            asevb_wat = max (asevd_wat, 0.026/(coszf(i)**1.7 + 0.065)   &
     &                  + 0.15 * (coszf(i)-0.1) * (coszf(i)-0.5)        &
     &                  * (coszf(i)-f_one))
            asenb_wat = asevb_wat
          endif

          if (icy(i)) then
          !-- Computation of ice albedo
            asnow = 0.02*snowf(i)
            argh  = min(0.50, max(.025, 0.01*zorlf(i)))
            hrgh  = min(f_one,max(0.20,1.0577-1.1538e-3*hprif(i)))
            fsno0 = asnow / (argh + asnow) * hrgh ! snow fraction on ice
            ! diffused
            if (tsknf(i) > 271.1 .and. tsknf(i) < 271.5) then
            !tgs: looks like albedo reduction from puddles on ice
              a1 = (tsknf(i) - 271.1)**2
              asevd_ice = 0.7 - 4.0*a1
              asend_ice = 0.65 - 3.6875*a1
            else
              asevd_ice = 0.70
              asend_ice = 0.65
            endif
            ! direct
            asevb_ice = asevd_ice
            asenb_ice = asend_ice

            if (fsno0 > f_zero) then
            ! Snow on ice
              dtgd = max(f_zero, min(5.0, (con_ttp-tisfc(i)) ))
              b1   = 0.03 * dtgd
              asnvd = (asevd_ice + b1) ! diffused snow albedo
              asnnd = (asend_ice + b1)
              if (coszf(i) > 0.0001 .and. coszf(i) < 0.5) then ! direct snow albedo
                csnow = 0.5 * (3.0 / (f_one+4.0*coszf(i)) - f_one)
                asnvb = min( 0.98, asnvd+(f_one-asnvd)*csnow )
                asnnb = min( 0.98, asnnd+(f_one-asnnd)*csnow )
              else
                asnvb = asnvd
                asnnb = asnnd
              endif

              ! composite ice and snow albedos
              asevd_ice = asevd_ice * (1. - fsno0) + asnvd * fsno0
              asend_ice = asend_ice * (1. - fsno0) + asnnd * fsno0
              asevb_ice = asevb_ice * (1. - fsno0) + asnvb * fsno0
              asenb_ice = asenb_ice * (1. - fsno0) + asnnb * fsno0
            endif ! snow
          else
          ! icy = false, fill in values
            asevd_ice = 0.70
            asend_ice = 0.65
            asevb_ice = 0.70
            asenb_ice = 0.65
          endif ! end icy

          if (fracl(i) > f_zero) then
!>  - Use snow cover input directly for land model, no
!!      conversion needed.

            fsno0 = sncovr(i) ! snow fraction on land

            fsno1 = f_one - fsno0 
            flnd0 = min(f_one, facsf(i)+facwf(i))
            flnd  = flnd0 * fsno1 ! snow-free fraction
            fsno  = f_one - flnd  ! snow-covered fraction

            !>  - use Fanglin's zenith angle treatment.
            if (coszf(i) > 0.0001) then
              rfcs = 1.775/(1.0+1.55*coszf(i))
            else
            !- no sun
              rfcs  = f_one
            endif
            !- zenith dependence is applied only to direct beam albedo
            ab1bm = min(0.99, alnsf(i)*rfcs)
            ab2bm = min(0.99, alvsf(i)*rfcs)

            alndnb = ab1bm   *flnd + snoalb(i) * fsno
            alndnd = alnwf(i)*flnd + snoalb(i) * fsno
            alndvb = ab2bm   *flnd + snoalb(i) * fsno
            alndvd = alvwf(i)*flnd + snoalb(i) * fsno
          else
          !-- fill in values for land albedo
            alndnb = 0.
            alndnd = 0.
            alndvb = 0.
            alndvd = 0.
          endif ! end land

          !-- Composite mean surface albedo from land, open water and
          !-- ice fractions
          sfcalb(i,1) = min(0.99,max(0.01,alndnb))*fracl(i)             & ! direct beam NIR
     &                  + asenb_wat*fraco(i) + asenb_ice*fraci(i)
          sfcalb(i,2) = min(0.99,max(0.01,alndnd))*fracl(i)             & ! diffuse NIR
     &                  + asend_wat*fraco(i) + asend_ice*fraci(i)
          sfcalb(i,3) = min(0.99,max(0.01,alndvb))*fracl(i)             & ! direct beam visible
     &                  + asevb_wat*fraco(i) + asevb_ice*fraci(i)
          sfcalb(i,4) = min(0.99,max(0.01,alndvd))*fracl(i)             & ! diffuse visible
     &                  + asevd_wat*fraco(i) + asevd_ice*fraci(i)

        enddo    ! end_do_i_loop

!> -# use land model output for land area: Noah MP, RUC (land and ice).
      elseif ( ialbflg == 2 ) then
        do i = 1, IMAX

          !-- water albedo
          asevd_wat = 0.06
          asend_wat = 0.06
          asevb_wat = asevd_wat
          asenb_wat = asevd_wat

          ! direct albedo CZA dependence over water
          if (fraco(i) > f_zero .and. coszf(i) > 0.0001) then
            asevb_wat = max (asevd_wat, 0.026/(coszf(i)**1.7 + 0.065)   &
     &                  + 0.15 * (coszf(i)-0.1) * (coszf(i)-0.5)        &
     &                  * (coszf(i)-f_one))
            asenb_wat = asevb_wat
          endif

          !-- ice albedo
          !tgs: this part of the code needs the input from the ice
          !     model. Otherwise it uses the backup albedo computation 
          !     from ialbflg = 1.
          if (icy(i)) then
            if(lsm == lsm_ruc ) then
            !-- use ice albedo from the RUC ice model
              asevd_ice = icealbivis(i)
              asend_ice = icealbinir(i)
              asevb_ice = icealbdvis(i)
              asenb_ice = icealbdnir(i)
            else
            !-- Computation of ice albedo
              asnow = 0.02*snowf(i)
              argh  = min(0.50, max(.025, 0.01*zorlf(i)))
              hrgh  = min(f_one,max(0.20,1.0577-1.1538e-3*hprif(i)))
              fsno0 = asnow / (argh + asnow) * hrgh
              ! diffused
              if (tsknf(i) > 271.1 .and. tsknf(i) < 271.5) then
              !tgs: looks like albedo reduction from puddles on ice
                a1 = (tsknf(i) - 271.1)**2
                asevd_ice = 0.7 - 4.0*a1
                asend_ice = 0.65 - 3.6875*a1
              else
                asevd_ice = 0.70
                asend_ice = 0.65
              endif
              ! direct
              asevb_ice = asevd_ice
              asenb_ice = asend_ice
  
              if (fsno0 > f_zero) then 
              ! Snow on ice
                dtgd = max(f_zero, min(5.0, (con_ttp-tisfc(i)) ))
                b1   = 0.03 * dtgd
                asnvd = (asevd_ice + b1)                                ! diffused snow albedo
                asnnd = (asend_ice + b1)

                if (coszf(i) > 0.0001 .and. coszf(i) < 0.5) then        ! direct snow albedo
                  csnow = 0.5 * (3.0 / (f_one+4.0*coszf(i)) - f_one)
                  asnvb = min( 0.98, asnvd+(f_one-asnvd)*csnow )
                  asnnb = min( 0.98, asnnd+(f_one-asnnd)*csnow )
                else
                  asnvb = asnvd
                  asnnb = asnnd
                endif
  
                ! composite ice and snow albedos
                asevd_ice = asevd_ice * (1. - fsno0) + asnvd * fsno0
                asend_ice = asend_ice * (1. - fsno0) + asnnd * fsno0
                asevb_ice = asevb_ice * (1. - fsno0) + asnvb * fsno0
                asenb_ice = asenb_ice * (1. - fsno0) + asnnb * fsno0
              endif ! snow
            endif ! ice option from LSM or otherwise
          else
          ! icy = false, fill in values
            asevd_ice = 0.70
            asend_ice = 0.65
            asevb_ice = 0.70
            asenb_ice = 0.65
          endif ! end icy
 
          !-- Composite mean surface albedo from land, open water and
          !-- ice fractions
          sfcalb(i,1) = min(0.99,max(0.01,lsmalbdnir(i)))*fracl(i)      & ! direct beam NIR
     &                  + asenb_wat*fraco(i) + asenb_ice*fraci(i)
          sfcalb(i,2) = min(0.99,max(0.01,lsmalbinir(i)))*fracl(i)      & ! diffuse NIR
     &                  + asend_wat*fraco(i) + asend_ice*fraci(i)
          sfcalb(i,3) = min(0.99,max(0.01,lsmalbdvis(i)))*fracl(i)      & ! direct beam visible
     &                  + asevb_wat*fraco(i) + asevb_ice*fraci(i)
          sfcalb(i,4) = min(0.99,max(0.01,lsmalbivis(i)))*fracl(i)      & ! diffuse visible
     &                  + asevd_wat*fraco(i) + asevd_ice*fraci(i)

        enddo    ! end_do_i_loop

      endif   ! end if_ialbflg
!

! sfc-perts, mgehne ***
!> - Call ppebet () to perturb all 4 elements of surface albedo sfcalb(:,1:4).
      if (pertalb>0.0) then
        do i = 1, imax
          do kk=1, 4
            ! compute beta distribution parameters for all 4 albedos
            m = sfcalb(i,kk)
            s = pertalb*m*(1.-m)
            alpha = m*m*(1.-m)/(s*s)-m
            beta  = alpha*(1.-m)/m
            ! compute beta distribution value corresponding
            ! to the given percentile albPpert to use as new albedo
            call ppfbet(albPpert(i),alpha,beta,iflag,albtmp)
            sfcalb(i,kk) = albtmp
          enddo
        enddo     ! end_do_i_loop
      endif

! *** sfc-perts, mgehne


      return
!...................................
      end subroutine setalb
!-----------------------------------
!! @}

!> \ingroup module_radiation_surface
!> This subroutine computes surface emissivity for LW radiation.
!!\param xlon      (IMAX), longitude in radiance, ok for both 0->2pi
!!                  or -pi -> +pi ranges
!!\param xlat      (IMAX), latitude  in radiance, default to pi/2 ->
!!                  -pi/2 range, otherwise see in-line comment
!!\param landfrac  (IMAX), fraction of grid that is land
!!\param snowf     (IMAX), snow depth water equivalent in mm
!!\param sncovr    (IMAX), snow cover over land
!!\param zorlf     (IMAX), surface roughness in cm
!!\param tsknf     (IMAX), ground surface temperature in K
!!\param tairf     (IMAX), lowest model layer air temperature in K
!!\param hprif     (IMAX), topographic standard deviation in m
!!\param IMAX       array horizontal dimension
!!\param sfcemis  (IMAX), surface emissivity
!>\section general_setemis setemis General Algorithm
!! @{
!-----------------------------------
      subroutine setemis                                                &
     &     ( lsm,lsm_noahmp,lsm_ruc,vtype,landfrac,frac_grid,           &  !  ---  inputs:
     &       min_seaice,xlon,xlat,slmsk,snowf,sncovr,sncovr_ice,        &
     &       zorlf,tsknf,tairf,hprif,                                   &
     &       semis_lnd,semis_ice,IMAX,fracl,fraco,fraci,icy,            &
     &       semisbase, sfcemis                                         &  !  ---  outputs:
     &     )

!  ===================================================================  !
!                                                                       !
!  this program computes surface emissivity for lw radiation.           !
!                                                                       !
!  usage:         call setemis                                          !
!                                                                       !
!  subprograms called:  none                                            !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     xlon  (IMAX)  - longitude in radiance, ok for both 0->2pi or      !
!                     -pi -> +pi ranges                                 !
!     xlat  (IMAX)  - latitude  in radiance, default to pi/2 -> -pi/2   !
!                     range, otherwise see in-line comment              !
!     slmsk (IMAX)  - sea(0),land(1),ice(2) mask on fcst model grid     !
!     landfrac (IMAX) - fraction of land on on fcst model grid          !
!     snowf (IMAX)  - snow depth water equivalent in mm                 !
!     sncovr(IMAX)  - ialbflg=1: snow cover over land in fraction       !
!     sncovr_ice(IMAX) - snow cover over ice in fraction                !
!     zorlf (IMAX)  - surface roughness in cm                           !
!     tsknf (IMAX)  - ground surface temperature in k                   !
!     tairf (IMAX)  - lowest model layer air temperature in k           !
!     hprif (IMAX)  - topographic sdv in m                              !
!     semis_lnd (IMAX) - emissivity from lsm                            !
!     IMAX          - array horizontal dimension                        !
!                                                                       !
!  outputs:                                                             !
!     sfcemis(IMAX)   - surface emissivity                              !
!                                                                       !
!  -------------------------------------------------------------------  !
!                                                                       !
!  surface type definations:                                            !
!     1. open water                   2. grass/wood/shrub land          !
!     3. tundra/bare soil             4. sandy desert                   !
!     5. rocky desert                 6. forest                         !
!     7. ice                          8. snow                           !
!                                                                       !
!  input index data lon from 0 towards east, lat from n to s            !
!                                                                       !
!  ====================    end of description    =====================  !
!
      use set_soilveg_ruc_mod,  only: set_soilveg_ruc
      use namelist_soilveg_ruc

      implicit none

!  ---  inputs
      integer, intent(in) :: IMAX
      integer, intent(in) :: lsm, lsm_noahmp, lsm_ruc
      logical, intent(in) :: frac_grid
      real (kind=kind_phys), dimension(:), intent(in) :: vtype
      real (kind=kind_phys), dimension(:), intent(in) :: landfrac
      real (kind=kind_phys), intent(in) :: min_seaice

      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       xlon,xlat, slmsk, snowf,sncovr, sncovr_ice,                &
     &       zorlf, tsknf, tairf, hprif, semis_lnd, semis_ice
      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       fracl, fraco, fraci
      logical, dimension(:), intent(in) ::                              &
     &       icy

!  ---  outputs
      real (kind=kind_phys), dimension(:), intent(out) :: semisbase
      real (kind=kind_phys), dimension(:), intent(out) :: sfcemis

!  ---  locals:
      integer :: i, i1, i2, j1, j2, idx
      integer :: ivgtyp

      real (kind=kind_phys) :: dltg, hdlt, tmp1, tmp2,                  &
     &      asnow, argh, hrgh, fsno
      real (kind=kind_phys) :: sfcemis_land, sfcemis_ice

!  ---  reference emiss value for diff surface emiss index
!       1-open water, 2-grass/shrub land, 3-bare soil, tundra,
!       4-sandy desert, 5-rocky desert, 6-forest, 7-ice, 8-snow

      real (kind=kind_phys) ::  emsref(8)
      data  emsref / 0.97, 0.95, 0.94, 0.90, 0.93, 0.96, 0.96, 0.99 /

!
!===> ...  begin here
!
!> -# Set emissivity by surface type and conditions
      if ( iemslw == 1 ) then

        dltg = 360.0 / float(IMXEMS)
        hdlt = 0.5 * dltg

!  --- ...  mapping input data onto model grid
!           note: this is a simple mapping method, an upgrade is needed if
!           the model grid is much coarser than the 1-deg data resolution

        lab_do_IMAX : do i = 1, IMAX

          if (fracl(i) < epsln) then                    ! no land
            if ( abs(fraco(i)-f_one) < epsln ) then     ! open water point
              sfcemis(i) = emsref(1)
            elseif ( abs(fraci(i)-f_one) > epsln ) then ! complete sea/lake ice
              sfcemis(i) = emsref(7)
            else
            !-- fractional sea ice
              sfcemis(i) = fraco(i)*emsref(1) + fraci(i)*emsref(7)
            endif

          else                                     ! land or fractional grid

!  ---  map grid in longitude direction
            i2 = 1
            j2 = 1

            tmp1 = xlon(i) * rad2dg
            if (tmp1 < f_zero) tmp1 = tmp1 + 360.0

            lab_do_IMXEMS : do i1 = 1, IMXEMS
              tmp2 = dltg * (i1 - 1) + hdlt

              if (abs(tmp1-tmp2) <= hdlt) then
               i2 = i1
                exit lab_do_IMXEMS
              endif
            enddo  lab_do_IMXEMS

!  ---  map grid in latitude direction
            tmp1 = xlat(i) * rad2dg           ! if xlat in pi/2 -> -pi/2 range
!           tmp1 = 90.0 - xlat(i)*rad2dg      ! if xlat in 0 -> pi range

            lab_do_JMXEMS : do j1 = 1, JMXEMS
              tmp2 = 90.0 - dltg * (j1 - 1)

              if (abs(tmp1-tmp2) <= hdlt) then
                j2 = j1
                exit lab_do_JMXEMS
              endif
            enddo  lab_do_JMXEMS

            idx = max( 2, idxems(i2,j2) )
            if ( idx >= 7 ) idx = 2
            if (abs(fracl(i)-f_one) < epsln) then
              sfcemis(i) = emsref(idx)
            else
              sfcemis(i) = fracl(i)*emsref(idx) + fraco(i)*emsref(1)          &
     &                                       + fraci(i)*emsref(7)
            endif
            semisbase(i) = sfcemis(i)

          endif   ! end if_slmsk_block

!> - Check for snow covered area.
          if ( sncovr(i) > f_zero ) then ! input land/ice area snow cover

            fsno = sncovr(i)
            sfcemis(i) = sfcemis(i)*(f_one - fsno) + emsref(8)*fsno

          else                                           ! compute snow cover from snow depth
            if ( snowf(i) > f_zero ) then
              asnow = 0.02*snowf(i)
              argh  = min(0.50, max(.025, 0.01*zorlf(i)))
              hrgh  = min(f_one, max(0.20, 1.0577-1.1538e-3*hprif(i) ) )
              fsno = asnow / (argh + asnow) * hrgh

              if (abs(fraco(i)-f_one) < epsln) fsno = f_zero         ! no snow over open water
              sfcemis(i) = sfcemis(i)*(f_one - fsno) + emsref(8)*fsno
            endif

          endif                                          ! end if_ialbflg

        enddo  lab_do_IMAX

      elseif ( iemslw == 2 ) then   ! sfc emiss updated in land model: Noah MP or RUC

        do i = 1, IMAX

          !-- ice emissivity
          sfcemis_ice = emsref(7)

          if ( icy(i) ) then 
          !-- complete or fractional ice
            if (lsm == lsm_noahmp) then
              if ( snowf(i) > f_zero ) then
                asnow = 0.02*snowf(i)
                argh  = min(0.50, max(.025,0.01*zorlf(i)))
                hrgh  = min(f_one,max(0.20,1.0577-1.1538e-3*hprif(i)))
                fsno = asnow / (argh + asnow) * hrgh
                sfcemis_ice = sfcemis_ice*(f_one-fsno)+emsref(8)*fsno
              endif
            elseif (lsm == lsm_ruc) then
              sfcemis_ice = semis_ice(i) ! output from lsm (with snow effect)
            endif ! lsm check
          endif ! icy

          !-- land emissivity 
          !-- from Noah MP or RUC lsms
          sfcemis_land = semis_lnd(i) ! albedo with snow effect from LSM

          !-- Composite emissivity from land, water and ice fractions.
          sfcemis(i) = fracl(i)*sfcemis_land + fraco(i)*emsref(1)             &
     &                                    + fraci(i)*sfcemis_ice

         enddo  ! i

      endif   ! end if_iemslw_block

!chk  print *,' In setemis, iemsflg, sfcemis =',iemsflg,sfcemis

!
      return
!...................................
      end subroutine setemis
!! @}
!-----------------------------------

!.........................................!
      end module module_radiation_surface !
!=========================================!
