!> \file GFS_RRTMG_pre.f90
!! This file contains
      module GFS_RRTMG_pre

      public GFS_RRTMG_pre_run

      contains

!> \defgroup GFS_RRTMG_pre GFS RRTMG Scheme Pre
!! @{
!!\section arg_table_GFS_RRTMG_pre_init Argument Table
!!
      subroutine GFS_RRTMG_pre_init
      end subroutine GFS_RRTMG_pre_init      

!!\section arg_table_GFS_RRTMG_pre_run Argument Table
!!| local var name | longname                                      | description                                                          | units       | rank |  type                         |   kind    | intent | optional |
!!|----------------|----------------------------------- -----------|----------------------------------------------------------------------|-------------|------|-------------------------------|-----------|--------|----------|
!!|   Model        | FV3-GFS_Control_type                          | Fortran DDT containing FV3-GFS model control parameters              |  DDT        |  0   | GFS_typedefs%GFS_control_type |           | in     | F        |
!!|   Grid         | FV3-GFS_Grid_type                             | Fortran DDT containing FV3-GFS grid and interpolation related data   |  DDT        |  0   | GFS_typedefs%GFS_grid_type    |           | in     | F        |
!!|   lm
!!|   me
!!|   im
!!|   ntrac
!!|  lmk
      subroutine GFS_RRTMG_pre_run (Model, Grid, lm, me, im,  ntrac, &
          lmk, lmp, kd, kt, kb, lla, llb, lya, lyb, lp1, raddt,  &
           tskn, tsfg, Sfcprop,  Statein, plvl, plyr,            &
          tlyr, prslk1, rhly, qstl, tracer1, olyr, Radtend,      &
          gasvmr, tlvl, tsfa, tvly, qlyr, nday, idxday, faersw,  &
          faerlw, aerodp, Tbd, Cldprop, deltaq, clouds, cldsa,   &
          mtopa, mbota, sfcalb)



!zhang        implicit none
      use machine,                   only: kind_phys
      use GFS_typedefs,              only: GFS_statein_type,   &
                                           GFS_stateout_type,  &
                                           GFS_sfcprop_type,   &
                                           GFS_coupling_type,  &
                                           GFS_control_type,   &
                                           GFS_grid_type,      &
                                           GFS_tbd_type,       &
                                           GFS_cldprop_type,   &
                                           GFS_radtend_type,   &
                                           GFS_diag_type
      use physparam
      use physcons,                  only: eps   => con_eps,         &
     &                                     epsm1 => con_epsm1,       &
     &                                     fvirt => con_fvirt        &
     &,                                    rocp  => con_rocp
      use funcphys,                  only: fpvs

      use module_radiation_astronomy,only: sol_init, sol_update, coszmn
      use module_radiation_gases,    only: NF_VGAS, getgases, getozn,  &
     &                                     gas_init, gas_update
      use module_radiation_aerosols, only: NF_AESW, NF_AELW, setaer,   &
     &                                     aer_init, aer_update,       &
     &                                     NSPC1
      use module_radiation_surface,  only: NF_ALBD, sfc_init, setalb,  &
     &                                     setemis
      use module_radiation_clouds,   only: NF_CLDS, cld_init,          &
     &                                     progcld1, progcld2,progcld3,&
     &                                     progclduni, diagcld1

      use module_radsw_parameters,   only: topfsw_type, sfcfsw_type,   &
     &                                     profsw_type,cmpfsw_type,NBDSW
!zhang      use module_radsw_main,         only: rswinit,  swrad

      use module_radlw_parameters,   only: topflw_type, sfcflw_type,    &
     &                                     proflw_type, NBDLW
!zhang      use module_radlw_main,         only: rlwinit,  lwrad


      implicit none
        !integer, intent(inout) :: me, lm, im, lp1, ntrac
        !integer, intent(inout) :: lmk, lmp, kd, kt, kb, lla, llb, lya, lyb
        type(GFS_control_type),   intent(in) :: Model
        type(GFS_grid_type),      intent(in) :: Grid
        type(GFS_sfcprop_type),         intent(in)    :: Sfcprop
        type(GFS_statein_type), intent(in) :: Statein
        type(GFS_radtend_type), intent(in) :: Radtend
        type(GFS_tbd_type),     intent(in) :: Tbd
        type(GFS_cldprop_type), intent(in) :: Cldprop

        !integer, intent(out) :: nday
        !integer, dimension(Size (Grid%xlon, 1)), intent(inout) :: idxday
        !real(kind=kind_phys), intent(out)    :: raddt
        !real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(inout) :: tsfg, tskn
        !real(kind=kind_phys), dimension(Size(Grid%xlon, 1), Model%levr+1+LTP), intent(inout) :: plvl
        !real(kind=kind_phys), dimension(size(Grid%xlon, 1), Model%levr+LTP), intent(inout) :: plyr, tlyr, prslk1, rhly, qstl
        !real(kind=kind_phys), dimension(Size (Grid%xlon, 1), Model%levr + &
        !    LTP, 2:Model%ntrac), intent(inout) :: tracer1
        !real(kind=kind_phys), dimension(Size (Grid%xlon, 1), Model%levr + &
        !     LTP), intent(inout) :: olyr
        !real(kind = kind_phys), dimension(Size (Grid%xlon, 1), Model%levr + &
        !     LTP, NF_VGAS), intent(inout) :: gasvmr
        !real(kind = kind_phys), dimension(Size (Grid%xlon, 1), Model%levr + &
        !    1 + LTP), intent(inout) :: tlvl
        !real(kind = kind_phys), dimension(Size (Grid%xlon, 1)) :: tsfa,tem1d
        !real(kind = kind_phys), dimension(Size (Grid%xlon, 1), Model%levr+LTP), intent(inout) :: qlyr, tvly
        !real(kind = kind_phys), dimension(Size (Grid%xlon, 1), Model%levr + &
        !    LTP, NBDSW, NF_AESW), intent(inout) :: faersw
        !real(kind = kind_phys), dimension(Size (Grid%xlon, 1), Model%levr + &
        !    LTP, NBDLW, NF_AELW), intent(inout) :: faerlw
        !real(kind = kind_phys), dimension(Size (Grid%xlon, 1), NSPC1), intent(inout) :: aerodp
        !real(kind = kind_phys), dimension(size(Grid%xlon, 1), Model%levr + &
        !    LTP), intent(out) :: deltaq
        !real(kind = kind_phys), dimension(Size (Grid%xlon, 1), Model%levr + &
        !    LTP, NF_CLDS), intent(inout) :: clouds
        !real(kind = kind_phys), dimension(Size (Grid%xlon, 1), 5), intent(out) :: cldsa
        !integer, dimension(size(Grid%xlon, 1), 3), intent(out) :: mbota, mtopa
        !real (kind = kind_phys), dimension(im, NF_ALBD), intent(out) :: sfcalb
        !real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP) :: &
        !   htswc, htlwc, gcice, grain, grime, htsw0, htlw0, plyr, tlyr,     &
        !   qlyr, olyr, rhly, tvly,qstl, vvel, clw, ciw, prslk1, tem2da,     &
        !   tem2db, cldcov, deltaq, cnvc, cnvw


!  ---  version tag and last revision date
      character(40), parameter ::                                    &
     &   VTAGRAD='NCEP-Radiation_driver    v5.2  Jan 2013 '
!    &   VTAGRAD='NCEP-Radiation_driver    v5.1  Nov 2012 '
!    &   VTAGRAD='NCEP-Radiation_driver    v5.0  Aug 2012 '

!>\name Constant values

!> lower limit of saturation vapor pressure (=1.0e-10)
      real (kind=kind_phys) :: QMIN
!> lower limit of specific humidity (=1.0e-7)
      real (kind=kind_phys) :: QME5
!> lower limit of specific humidity (=1.0e-7)
      real (kind=kind_phys) :: QME6
!> EPSQ=1.0e-12
      real (kind=kind_phys) :: EPSQ
!     parameter (QMIN=1.0e-10, QME5=1.0e-5,  QME6=1.0e-6,  EPSQ=1.0e-12)
      parameter (QMIN=1.0e-10, QME5=1.0e-7,  QME6=1.0e-7,  EPSQ=1.0e-12)
!     parameter (QMIN=1.0e-10, QME5=1.0e-20, QME6=1.0e-20, EPSQ=1.0e-12)

!> lower limit of toa pressure value in mb
      real, parameter :: prsmin = 1.0e-6

!> control flag for LW surface temperature at air/ground interface
!! (default=0, the value will be set in subroutine radinit)
      integer :: itsfc  =0

!> new data input control variables (set/reset in subroutines
!radinit/radupdate):
      integer :: month0=0,   iyear0=0,   monthd=0

!> control flag for the first time of reading climatological ozone data
!! (set/reset in subroutines radinit/radupdate, it is used only if the
!! control parameter ioznflg=0)
      logical :: loz1st =.true.

!> optional extra top layer on top of low ceiling models
!!\n LTP=0: no extra top layer
      integer, parameter :: LTP = 0   ! no extra top layer
!     integer, parameter :: LTP = 1   ! add an extra top layer

!> control flag for extra top layer
      logical, parameter :: lextop = (LTP > 0)

!
!  ---  local variables: (horizontal dimensioned by IM)
      !--- INTEGER VARIABLES
      integer :: me, im, lm, nfxr, ntrac
      integer :: i, j, k, k1, lv, itop, ibtc, nday, LP1, LMK, LMP, kd, &
                 lla, llb, lya, lyb, kt, kb
      integer, dimension(size(Grid%xlon,1)) :: idxday
      integer, dimension(size(Grid%xlon,1),3) :: mbota, mtopa

      !--- REAL VARIABLES
      real(kind=kind_phys) :: raddt, es, qs, delt, tem0d

      real(kind=kind_phys), dimension(size(Grid%xlon,1)) :: &
           tsfa, cvt1, cvb1, tem1d, tsfg, tskn

      real(kind=kind_phys), dimension(size(Grid%xlon,1),5)       :: cldsa
      real(kind=kind_phys), dimension(size(Grid%xlon,1),NSPC1)   :: aerodp
      real(kind=kind_phys), dimension(size(Grid%xlon,1),NF_ALBD) :: sfcalb

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP) :: &
           htswc, htlwc, gcice, grain, grime, htsw0, htlw0, plyr, tlyr, &
           qlyr, olyr, rhly, tvly,qstl, vvel, clw, ciw, prslk1, tem2da, &
           tem2db, cldcov, deltaq, cnvc, cnvw

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+1+LTP) :: plvl, tlvl

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,2:Model%ntrac) :: tracer1
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NF_CLDS) :: clouds
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NF_VGAS) :: gasvmr

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDSW,NF_AESW)::faersw
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDLW,NF_AELW)::faerlw

      !--- TYPED VARIABLES
      type (cmpfsw_type),    dimension(size(Grid%xlon,1)) :: scmpsw



          !pedro Set commonly used integers
          !pedro call Set_common_int (Model, Grid, lm, me, im, lp1, ntrac)

!
!===> ...  begin here
!
      !--- set commonly used integers
      me = Model%me
      LM = Model%levr
      IM = size(Grid%xlon,1)
      NFXR = Model%nfxr
      NTRAC = Model%ntrac        ! tracers in grrad strip off sphum - start tracer1(2:NTRAC)

      LP1 = LM + 1               ! num of in/out levels


          !pedro Set local /level/layer indexes corresponding
          !pedro  to in/out variables
        !pedro call Set_local_int (lmk, lm, lmp, kd, kt, &
        !pedro    kb, lla, llb, lya, lyb, lp1, raddt, Model)

!  --- ...  set local /level/layer indexes corresponding to in/out
!  variables

      LMK = LM + LTP             ! num of local layers
      LMP = LMK + 1              ! num of local levels

      if ( lextop ) then
        if ( ivflip == 1 ) then    ! vertical from sfc upward
          kd = 0                   ! index diff between in/out and local
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
          lla = LMK                ! local index at the 2nd level from top
          llb = LMP                ! local index at toa level
          lya = LM                 ! local index for the 2nd layer from top
          lyb = LP1                ! local index for the top layer
        else                       ! vertical from toa downward
          kd = 1                   ! index diff between in/out and local
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
          lla = 2                  ! local index at the 2nd level from top
          llb = 1                  ! local index at toa level
          lya = 2                  ! local index for the 2nd layer from top
          lyb = 1                  ! local index for the top layer
        endif                    ! end if_ivflip_block
      else
        kd = 0
        if ( ivflip == 1 ) then  ! vertical from sfc upward
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
        else                     ! vertical from toa downward
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
        endif                    ! end if_ivflip_block
      endif   ! end if_lextop_block

      raddt = min(Model%fhswr, Model%fhlwr)
!     print *,' in grrad : raddt=',raddt

          !pedro Setup surface ground temperature and
          !pedro ground/air skin temperature if required.
          !pedro call Set_sfc_vars (im, tskn, tsfg, Sfcprop, Grid)

!> -# Setup surface ground temperature and ground/air skin temperature
!! if required.

      if ( itsfc == 0 ) then            ! use same sfc skin-air/ground temp
        do i = 1, IM
          tskn(i) = Sfcprop%tsfc(i)
          tsfg(i) = Sfcprop%tsfc(i)
        enddo
      else                              ! use diff sfc skin-air/ground temp
        do i = 1, IM
          tskn(i) = Sfcprop%tsfc(i)
          tsfg(i) = Sfcprop%tsfc(i)
        enddo
      endif



          !pedro Prepare atmospheric profiles.
          !pedro Convert pressure unit from pa to mb
          !pedro call Prep_profiles (lm, kd, im, Statein, plvl, plyr, tlyr, &
          !  prslk1, rhly, qstl, Model, Grid)

!> -# Prepare atmospheric profiles for radiation input.
!
!           convert pressure unit from pa to mb
      do k = 1, LM
        k1 = k + kd
        do i = 1, IM
          plvl(i,k1)   = 0.01 * Statein%prsi(i,k)   ! pa to mb (hpa)
          plyr(i,k1)   = 0.01 * Statein%prsl(i,k)   ! pa to mb (hpa)
          tlyr(i,k1)   = Statein%tgrs(i,k)
          prslk1(i,k1) = Statein%prslk(i,k)

!>  - Compute relative humidity.
!         es  = min( Statein%prsl(i,k), 0.001 * fpvs( Statein%tgrs(i,k)   ) )   ! fpvs in pa
          es  = min( Statein%prsl(i,k),  fpvs( Statein%tgrs(i,k) ) )  ! fpvs and prsl in pa
          qs  = max( QMIN, eps * es / (Statein%prsl(i,k) + epsm1*es) )
          rhly(i,k1) = max( 0.0, min( 1.0, max(QMIN, Statein%qgrs(i,k,1))/qs ) )
          qstl(i,k1) = qs
        enddo
      enddo



          !pedro Recast remaining all tracers (except sphum)
          !pedro  forcing them all to be positive
        !pedro call Recast_tracers (tracer1, plvl, plyr, tlyr, prslk1, rhly, &
        !    qstl, Statein, Grid, Model, ntrac, lm, im, kd, lp1, llb, &
        !    lla, lya, lyb)

      !--- recast remaining all tracers (except sphum) forcing them all
      !to be positive
      do j = 2, NTRAC
        do k = 1, LM
          k1 = k + kd
          tracer1(:,k1,j) = max(0.0,Statein%qgrs(:,k,j))
        enddo
      enddo

      do i = 1, IM
        plvl(i,LP1+kd) = 0.01 * Statein%prsi(i,LP1)  ! pa to mb (hpa)
      enddo

      if ( lextop ) then                 ! values for extra top layer
        do i = 1, IM
          plvl(i,llb) = prsmin
          if ( plvl(i,lla) <= prsmin ) plvl(i,lla) = 2.0*prsmin
          plyr(i,lyb)   = 0.5 * plvl(i,lla)
          tlyr(i,lyb)   = tlyr(i,lya)
          prslk1(i,lyb) = (plyr(i,lyb)*0.00001) ** rocp ! plyr in Pa
          rhly(i,lyb)   = rhly(i,lya)
          qstl(i,lyb)   = qstl(i,lya)
        enddo

!  ---  note: may need to take care the top layer amount
       tracer1(:,lyb,:) = tracer1(:,lya,:)
      endif


          !pedro Get layer ozone mass mixing ratio
        !pedro call Prep_ozone  (Model, Grid, im, lmk, tracer1, olyr, prslk1)

!>  - Get layer ozone mass mixing ratio (if use ozone climatology data,
!!    call getozn()).

      if (Model%ntoz > 0) then            ! interactive ozone generation
        olyr(:,:) = max( QMIN, tracer1(:,1:LMK,Model%ntoz) )
      else                                ! climatological ozone
        call getozn (prslk1, Grid%xlat, IM, LMK,    &     !  ---  inputs
                     olyr)                                !  ---  outputs
      endif                               ! end_if_ntoz



         !pedro Compute cosine of zenith angle.
        !pedro call coszmn (Grid%xlon,Grid%sinlat, Grid%coslat, Model%solhr, &
        !pedor    im, me, Radtend%coszen, Radtend%coszdg)
!>  - Call coszmn(), to compute cosine of zenith angle.
      call coszmn (Grid%xlon,Grid%sinlat,           &     !  ---  inputs
                   Grid%coslat,Model%solhr, IM, me, &
                   Radtend%coszen, Radtend%coszdg)        !  --- outputs


          !pedro Set up non-prognostic gas volume mixing ratioes
        !pedro call getgases (plvl, Grid%xlon, Grid%xlat, im, lmk, gasvmr)

!>  - Call getgases(), to set up non-prognostic gas volume mixing
!!    ratioes (gasvmr).
!  - gasvmr(:,:,1)  -  co2 volume mixing ratio
!  - gasvmr(:,:,2)  -  n2o volume mixing ratio
!  - gasvmr(:,:,3)  -  ch4 volume mixing ratio
!  - gasvmr(:,:,4)  -  o2  volume mixing ratio
!  - gasvmr(:,:,5)  -  co  volume mixing ratio
!  - gasvmr(:,:,6)  -  cf11 volume mixing ratio
!  - gasvmr(:,:,7)  -  cf12 volume mixing ratio
!  - gasvmr(:,:,8)  -  cf22 volume mixing ratio
!  - gasvmr(:,:,9)  -  ccl4 volume mixing ratio

!  --- ...  set up non-prognostic gas volume mixing ratioes

      call getgases (plvl, Grid%xlon, Grid%xlat, IM, LMK,  & !  --- inputs
                     gasvmr)                                 !  --- outputs



       !pedro Get temperature at layer interface, and layer moisture.
       !pedro call Prep_t_and_moist (Grid, Model, Statein, lmp, kd, lmk, lm, &
       !pedro     im, lya, lyb, plyr, tlyr, tlvl, plvl, tsfa, tskn, tvly, qlyr)

!>  - Get temperature at layer interface, and layer moisture.
      do k = 2, LMK
        do i = 1, IM
          tem2da(i,k) = log( plyr(i,k) )
          tem2db(i,k) = log( plvl(i,k) )
        enddo
      enddo

      if (ivflip == 0) then              ! input data from toa to sfc
        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = 1.0
          tsfa  (i)   = tlyr(i,LMK)                  ! sfc layer air temp
          tlvl(i,1)   = tlyr(i,1)
          tlvl(i,LMP) = tskn(i)
        enddo

        do k = 1, LM
          k1 = k + kd
          do i = 1, IM
            qlyr(i,k1) = max( tem1d(i), Statein%qgrs(i,k,1) )
            tem1d(i)   = min( QME5, qlyr(i,k1) )
            tvly(i,k1) = Statein%tgrs(i,k) * (1.0 + fvirt*qlyr(i,k1)) ! virtual T (K)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        do k = 2, LMK
          do i = 1, IM
            tlvl(i,k) = tlyr(i,k) + (tlyr(i,k-1) - tlyr(i,k))           &
     &                * (tem2db(i,k)   - tem2da(i,k))                   &
     &                / (tem2da(i,k-1) - tem2da(i,k))
          enddo
        enddo

      else                               ! input data from sfc to toa

        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = log( plvl(i,1) )
          tsfa  (i)   = tlyr(i,1)                    ! sfc layer air temp
          tlvl(i,1)   = tskn(i)
          tlvl(i,LMP) = tlyr(i,LMK)
        enddo

        do k = LM, 1, -1
          do i = 1, IM
            qlyr(i,k) = max( tem1d(i), Statein%qgrs(i,k,1) )
            tem1d(i)  = min( QME5, qlyr(i,k) )
            tvly(i,k) = Statein%tgrs(i,k) * (1.0 + fvirt*qlyr(i,k)) ! virtual T (K)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        do k = 1, LMK-1
          do i = 1, IM
            tlvl(i,k+1) = tlyr(i,k) + (tlyr(i,k+1) - tlyr(i,k))         &
     &                  * (tem2db(i,k+1) - tem2da(i,k))                 &
     &                  / (tem2da(i,k+1) - tem2da(i,k))
          enddo
        enddo

      endif                              ! end_if_ivflip

       !pedro Check for daytime points for SW radiation.
       !pedro call Find_daytime (im, Radtend, Grid, nday, idxday)

!>  - Check for daytime points for SW radiation.

      nday = 0
      do i = 1, IM
        if (Radtend%coszen(i) >= 0.0001) then
          nday = nday + 1
          idxday(nday) = i
        endif
      enddo

      !pedro Setup aerosols
      !pedro  call setaer (plvl, plyr, prslk1, tvly, rhly, Sfcprop%slmsk,   &
      !pedro      tracer1, Grid%xlon, Grid%xlat, im, lmk, lmp, Model%lsswr, &
      !pedro      Model%lslwr, faersw,faerlw,aerodp)

!>  - Call module_radiation_aerosols::setaer(),to setup aerosols
!! property profile for radiation.

!check  print *,' in grrad : calling setaer '

      call setaer (plvl, plyr, prslk1, tvly, rhly, Sfcprop%slmsk, &  !  ---  inputs
                   tracer1, Grid%xlon, Grid%xlat, IM, LMK, LMP,    &
                   Model%lsswr,Model%lslwr,                        &
                   faersw,faerlw,aerodp)                              !  ---  outputs



          !pedro  Obtain cloud information
        !pedro call Get_cloud_info (Model, Grid, Tbd, Sfcprop, Cldprop,  &
        !pedro    Statein, tracer1, lmk, lmp, lm, lya, lyb, im, me, kd, &
        !pedro    deltaq, plvl, plyr, tlyr, qlyr, tvly,   &
        !pedro    rhly, qstl, clouds, cldsa, mtopa, mbota)

!>  - Obtain cloud information for radiation calculations
!!    (clouds,cldsa,mtopa,mbota)
!!\n   for  prognostic cloud:
!!    - For Zhao/Moorthi's prognostic cloud scheme,
!!      call module_radiation_clouds::progcld1()
!!    - For Zhao/Moorthi's prognostic cloud+pdfcld,
!!      call module_radiation_clouds::progcld3()
!!      call module_radiation_clouds::progclduni() for unified cloud and ncld=2
!>  - If cloud condensate is not computed (ntcw=0), using the legacy
!!   cloud scheme, compute cloud information based on Slingo's
!!   diagnostic cloud scheme (call module_radiation_clouds::diagcld1())

!  --- ...  obtain cloud information for radiation calculations

      if (Model%ntcw > 0) then                   ! prognostic cloud scheme
        if (Model%uni_cld .and. Model%ncld >= 2) then
          clw(:,:) = tracer1(:,1:LMK,Model%ntcw)              ! cloud water amount
          ciw(:,:) = 0.0
          do j = 2, Model%ncld
            ciw(:,:) = ciw(:,:) + tracer1(:,1:LMK,Model%ntcw+j-1)   ! cloud ice amount
          enddo

          do k = 1, LMK
            do i = 1, IM
              if ( clw(i,k) < EPSQ ) clw(i,k) = 0.0
              if ( ciw(i,k) < EPSQ ) ciw(i,k) = 0.0
            enddo
          enddo
        else
          clw(:,:) = 0.0
          do j = 1, Model%ncld
            clw(:,:) = clw(:,:) + tracer1(:,1:LMK,Model%ntcw+j-1)   ! cloud condensate amount
          enddo

          do k = 1, LMK
            do i = 1, IM
              if ( clw(i,k) < EPSQ ) clw(i,k) = 0.0
            enddo
          enddo
        endif
!
!  --- add suspended convective cloud water to grid-scale cloud water
!      only for cloud fraction & radiation computation
!      it is to enhance cloudiness due to suspended convec cloud water
!      for zhao/moorthi's (icmphys=1) &
!          ferrier's (icmphys=2) microphysics schemes
!
        if (Model%shoc_cld) then                                       ! all but MG microphys
          cldcov(:,1:LM) = Tbd%phy_f3d(:,1:LM,Model%ntot3d-2)
        elseif (Model%ncld == 2) then                                  ! MG microphys (icmphys = 1)
          cldcov(:,1:LM) = Tbd%phy_f3d(:,1:LM,1)
        else                                                           ! neither of the other two cases
          cldcov = 0
        endif

        if ((Model%num_p3d == 4) .and. (Model%npdf3d == 3)) then       ! icmphys = 3
          deltaq(:,1:LM) = Tbd%phy_f3d(:,1:LM,5)
          cnvw  (:,1:LM) = Tbd%phy_f3d(:,1:LM,6)
          cnvc  (:,1:LM) = Tbd%phy_f3d(:,1:LM,7)
        elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then  ! icmphys = 1
          deltaq(:,1:LM) = 0.
          cnvw  (:,1:LM) = Tbd%phy_f3d(:,1:LM,Model%num_p3d+1)
          cnvc  (:,1:LM) = 0.
        else                                                           !  icmphys = 1 (ncld=2)
          deltaq = 0.0
          cnvw   = 0.0
          cnvc   = 0.0
        endif

        if (lextop) then
          cldcov(:,lyb) = cldcov(:,lya)
          deltaq(:,lyb) = deltaq(:,lya)
          cnvw  (:,lyb) = cnvw  (:,lya)
          cnvc  (:,lyb) = cnvc  (:,lya)
        endif

        if (icmphys == 1) then
          clw(:,1:LMK) = clw(:,1:LMK) + cnvw(:,1:LMK)
        endif
!

        if (icmphys == 1) then           ! zhao/moorthi's prognostic cloud scheme
                                         ! or unified cloud and/or with MG microphysics

          if (Model%uni_cld .and. Model%ncld >= 2) then
            call progclduni (plyr, plvl, tlyr, tvly, clw, ciw,    &    !  ---  inputs
                             Grid%xlat, Grid%xlon, Sfcprop%slmsk, &
                             IM, LMK, LMP, cldcov(:,1:LMK),       &
                             clouds, cldsa, mtopa, mbota)              !  ---  outputs
          else
            call progcld1 (plyr ,plvl, tlyr, tvly, qlyr, qstl,    &    !  ---  inputs
                           rhly, clw, Grid%xlat,Grid%xlon,        &
                           Sfcprop%slmsk, IM, LMK, LMP,           &
                           Model%uni_cld, Model%lmfshal,          &
                           Model%lmfdeep2, cldcov(:,1:LMK),       &
                           clouds, cldsa, mtopa, mbota)                !  ---  outputs
          endif

        elseif(icmphys == 3) then      ! zhao/moorthi's prognostic cloud+pdfcld

          call progcld3 (plyr, plvl, tlyr, tvly, qlyr, qstl, rhly,&    ! ---  inputs
                         clw, cnvw, cnvc, Grid%xlat, Grid%xlon,   &
                         Sfcprop%slmsk,im, lmk, lmp, deltaq,      &
                         Model%sup, Model%kdt, me,                &
                         clouds, cldsa, mtopa, mbota)                  ! ---  outputs

        endif                            ! end if_icmphys

      else                               ! diagnostic cloud scheme

        cvt1(:) = 0.01 * Cldprop%cvt(:)
        cvb1(:) = 0.01 * Cldprop%cvb(:)

        do k = 1, LM
          k1 = k + kd
          vvel(:,k1) = 0.01 * Statein%vvl(:,k)
        enddo
        if (lextop) then
          vvel(:,lyb) = vvel(:,lya)
        endif

!  ---  compute diagnostic cloud related quantities

        call diagcld1 (plyr, plvl, tlyr, rhly, vvel, Cldprop%cv,  &    !  ---  inputs
                       cvt1, cvb1, Grid%xlat, Grid%xlon,          &
                       Sfcprop%slmsk, IM, LMK, LMP,               &
                       clouds, cldsa, mtopa, mbota)                    !  ---  outputs

      endif                                ! end_if_ntcw

      !pedro Setup surface albedo for SW calculation
      !pedro  call Set_sfc_albedo (Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr,&    !  ---  inputs:
      !pedro      Sfcprop%snoalb, Sfcprop%zorl, Radtend%coszen, tsfg, tsfa, &
      !pedro      Sfcprop%hprim, Sfcprop%alvsf, Sfcprop%alnsf, Sfcprop%alvwf, &
      !pedro      Sfcprop%alnwf, Sfcprop%facsf, Sfcprop%facwf, Sfcprop%fice,  &
      !pedro      Sfcprop%tisfc, im, Model%lsswr,  &
      !pedro      sfcalb, Radtend%sfalb)                            !  --- outputs

!  --- ...  start radiation calculations
!           remember to set heating rate unit to k/sec!
!> -# Start SW radiation calculations
      if (Model%lsswr) then

!>  - Call module_radiation_surface::setalb() to setup surface albedo.
!!  for SW radiation.

        call setalb (Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr,&    !  ---  inputs:
                     Sfcprop%snoalb, Sfcprop%zorl, Radtend%coszen,&
                     tsfg, tsfa, Sfcprop%hprim, Sfcprop%alvsf,    &
                     Sfcprop%alnsf, Sfcprop%alvwf, Sfcprop%alnwf, &
                     Sfcprop%facsf, Sfcprop%facwf, Sfcprop%fice,  &
                     Sfcprop%tisfc, IM,                           &
                     sfcalb)                                           !  ---  outputs

!> -# Approximate mean surface albedo from vis- and nir-  diffuse values.
        Radtend%sfalb(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))


      endif  ! Model%lsswr

       !zhang: should called before 
            !pedro Setup surface emissivity for LW radiation.
       call setemis (Grid%xlon, Grid%xlat, Sfcprop%slmsk, & !  --- inputs
           Sfcprop%snowd, Sfcprop%sncovr, Sfcprop%zorl,   &
           tsfg, tsfa, Sfcprop%hprim, im, Model%lslwr,    &
           Radtend%semis)                                   !  --- outputs


      end subroutine GFS_RRTMG_pre_run
   
!!\section arg_table_GFS_RRTMG_pre_finalize Argument Table
!!
      subroutine GFS_RRTMG_pre_finalize
      end subroutine GFS_RRTMG_pre_finalize

!! @}
      end module GFS_RRTMG_pre


