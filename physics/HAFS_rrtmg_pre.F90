!> \file HAFS_rrtmg_pre.f90
!! This file contains
      module HAFS_rrtmg_pre

      public HAFS_rrtmg_pre_run

      contains

!> \defgroup HAFS_rrtmg_pre HAFS RRTMG Scheme Pre
!! @{
!! \section arg_table_HAFS_rrtmg_pre_init Argument Table
!!
      subroutine HAFS_rrtmg_pre_init ()

      use module_ra_cam_support , only : oznini


!-- calculate radiation time step

    STEPRA = nint(RADT*60./DT)
    STEPRA = max(STEPRA,1)

!-- initialization

   IF(start_of_simulation)THEN
     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        RTHRATEN(i,k,j)=0.
        RTHRATENLW(i,k,j)=0.
        RTHRATENSW(i,k,j)=0.
        CLDFRA(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO

     IF( PRESENT(cldfra_dp) ) THEN
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           cldfra_dp(i,k,j)=0.
           cldfra_sh(i,k,j)=0.
        ENDDO
        ENDDO
        ENDDO
     ENDIF

     if( present(cldfra_old) ) then
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           cldfra_old(i,k,j) = 0.
        ENDDO
        ENDDO
        ENDDO
     end if
   ENDIF


      !MZ read CAM ozone data and aerosol data

!-- ww: attempt to use CAM ozone and some aerosol profiles in all rad
!schemes
!   note that CAM option will still do the same. 
!   n_ozmixm: no of months; levsiz: = 59, vertical dim
!   Read in CAM ozone data, and interpolate data to model grid
!   Interpolation of aerosols is done on domain 1 only

      IF ( config_flags%o3input .EQ. 2 ) THEN
        CALL oznini(ozmixm,pin,levsiz,n_ozmixm,XLAT,                &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte)
      ENDIF

      IF ( config_flags%aer_opt .EQ. 1 .AND. id .EQ. 1 ) THEN
        CALL aerosol_in(aerodm,pina,alevsiz,n_ozmixm-1,no_src_types,XLAT,XLONG,  &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte)
       ENDIF



      end subroutine HAFS_rrtmg_pre_init

!> \section arg_table_HAFS_rrtmg_pre_run Argument Table
!! \htmlinclude HAFS_rrtmg_pre_run.html
!!
      ! Attention - the output arguments lm, im, lmk, lmp must not be set
      ! in the CCPP version - they are defined in the interstitial_create routine
      subroutine HAFS_rrtmg_pre_run (Model, Grid, Sfcprop, Statein,   & ! input
          Tbd, Cldprop, Coupling,                                    &
          Radtend,                                                   & ! input/output
!          f_ice, f_rain, f_rimef, flgmin, cwm,                       & ! F-A mp scheme only
          lm, im, lmk, lmp,                                          & ! input
          dx,                                                        &
          kd, kt, kb, raddt, delp, dz, plvl, plyr,                   & ! output
          tlvl, tlyr, tsfg, tsfa, qlyr, olyr,                        &
          gasvmr_co2,   gasvmr_n2o,   gasvmr_ch4,   gasvmr_o2,       &
          gasvmr_co,    gasvmr_cfc11, gasvmr_cfc12,                  &
          gasvmr_cfc22, gasvmr_ccl4,  gasvmr_cfc113,                 &
          faersw1,  faersw2,  faersw3,                               &
          faerlw1, faerlw2, faerlw3, aerodp,                         &
          clouds1, clouds2, clouds3, clouds4, clouds5, clouds6,      &
          clouds7, clouds8, clouds9, cldsa,                          &
          mtopa, mbota, de_lgth, alb1d, errmsg, errflg)

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
     &,                                    rog   => con_rog          &
     &,                                    rocp  => con_rocp
      use radcons,                   only: itsfc,ltp, lextop, qmin,  &
                                           qme5, qme6, epsq, prsmin
      use funcphys,                  only: fpvs

      use module_radiation_astronomy,only: coszmn                         ! sol_init, sol_update
      use module_radiation_gases,    only: NF_VGAS, getgases, getozn      ! gas_init, gas_update,
      use module_radiation_aerosols, only: NF_AESW, NF_AELW, setaer,   &  ! aer_init, aer_update,
     &                                     NSPC1
      use module_radiation_clouds,   only: NF_CLDS,                    &  ! cld_init
     &                                     progcld1, progcld3,         &
     &                                     progcld2,                   &
     &                                     progcld4, progcld5,         &
     &                                     progclduni,                 &
     &                                     cal_cldfra3,                &  ! Thompson partial cloud scheme
     &                                     find_cloudLayers,           &
     &                                     adjust_cloudICE,            &
     &                                     adjust_cloudH2O,            &
     &                                     adjust_cloudFinal           
      use module_radsw_parameters,   only: topfsw_type, sfcfsw_type,   &
     &                                     profsw_type, NBDSW
      use module_radlw_parameters,   only: topflw_type, sfcflw_type,   &
     &                                     proflw_type, NBDLW
      use surface_perturbation,      only: cdfnor

      use module_ra_cam_support,     only: oznini
      implicit none

      type(GFS_control_type),              intent(in)    :: Model
      type(GFS_grid_type),                 intent(in)    :: Grid
      type(GFS_sfcprop_type),              intent(in)    :: Sfcprop
      type(GFS_statein_type),              intent(in)    :: Statein
      type(GFS_radtend_type),              intent(inout) :: Radtend
      type(GFS_tbd_type),                  intent(in)    :: Tbd
      type(GFS_cldprop_type),              intent(in)    :: Cldprop
      type(GFS_coupling_type),             intent(in)    :: Coupling

      integer,              intent(in)  :: im, lm, lmk, lmp
      integer,              intent(out) :: kd, kt, kb

! F-A mp scheme only
!      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP), intent(in) :: f_ice
!      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP), intent(in) :: f_rain
!      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP), intent(in) :: f_rimef
!      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP), intent(out) :: cwm
!      real(kind=kind_phys), dimension(size(Grid%xlon,1)),                intent(in)  :: flgmin
      real(kind=kind_phys), intent(out) :: raddt

      real(kind=kind_phys), dimension(size(Grid%xlon,1)),                  INTENT(IN) :: dx
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: delp
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: dz
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+1+LTP), intent(out) :: plvl
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: plyr
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+1+LTP), intent(out) :: tlvl
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: tlyr
      real(kind=kind_phys), dimension(size(Grid%xlon,1)),                  intent(out) :: tsfg
      real(kind=kind_phys), dimension(size(Grid%xlon,1)),                  intent(out) :: tsfa
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: qlyr
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: olyr

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_co2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_n2o
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_ch4
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_o2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_co
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_cfc11
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_cfc12
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_cfc22
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_ccl4
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_cfc113

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDSW), intent(out) :: faersw1
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDSW), intent(out) :: faersw2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDSW), intent(out) :: faersw3
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDLW), intent(out) :: faerlw1
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDLW), intent(out) :: faerlw2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDLW), intent(out) :: faerlw3

      real(kind=kind_phys), dimension(size(Grid%xlon,1),NSPC1),            intent(out) :: aerodp
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds1
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds3
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds4
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds5
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds6
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds7
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds8
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds9
      real(kind=kind_phys), dimension(size(Grid%xlon,1),5),                intent(out) :: cldsa
      integer,              dimension(size(Grid%xlon,1),3),                intent(out) :: mbota
      integer,              dimension(size(Grid%xlon,1),3),                intent(out) :: mtopa
      real(kind=kind_phys), dimension(size(Grid%xlon,1)),                  intent(out) :: de_lgth
      real(kind=kind_phys), dimension(size(Grid%xlon,1)),                  intent(out) :: alb1d

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! Local variables
      integer :: me, nfxr, ntrac, ntcw, ntiw, ncld, ntrw, ntsw, ntgl, ncndl

      integer :: i, j, k, k1, k2, lsk, lv, n, itop, ibtc, LP1, lla, llb, lya, lyb

      real(kind=kind_phys) :: es, qs, delt, tem0d,gridkm
      real(kind=kind_phys), dimension(size(Grid%xlon,1)) :: sm, xland

      real(kind=kind_phys), dimension(size(Grid%xlon,1)) :: cvt1, cvb1, tem1d, tskn

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP) :: &
                          htswc, htlwc, gcice, grain, grime, htsw0, htlw0, &
                          rhly, tvly,qstl, vvel, clw, ciw, prslk1, tem2da, &
                          cldcov, deltaq, cnvc, cnvw,                      &
                          effrl, effri, effrr, effrs

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP+1) :: tem2db
!     real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP+1) :: hz

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,min(4,Model%ncnd)) :: ccnd
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,2:Model%ntrac) :: tracer1
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NF_CLDS) :: clouds
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NF_VGAS) :: gasvmr
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDSW,NF_AESW)::faersw
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDLW,NF_AELW)::faerlw

      ! Local variables: dimensions used in cal_cldfra3
      integer               :: ids,ide, jds,jde, kds,kde, &
                               ims,ime, jms,jme, kms,kme, &
                               its,ite, jts,jte, kts,kte

!
!===> ...  begin here
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Set internal dimensions
      ids = 1
      ims = 1
      its = 1
      ide = size(Grid%xlon,1)
      ime = size(Grid%xlon,1)
      ite = size(Grid%xlon,1)
      jds = 1
      jms = 1
      jts = 1
      jde = 1
      jme = 1
      jte = 1
      kds = 1
      kms = 1
      kts = 1
      kde = Model%levr+LTP
      kme = Model%levr+LTP
      kte = Model%levr+LTP


      if (.not. (Model%lsswr .or. Model%lslwr)) return

      !--- set commonly used integers
      me = Model%me
      NFXR = Model%nfxr
      NTRAC = Model%ntrac        ! tracers in grrad strip off sphum - start tracer1(2:NTRAC)
      ntcw  = Model%ntcw
      ntiw  = Model%ntiw
      ncld  = Model%ncld
      ntrw  = Model%ntrw
      ntsw  = Model%ntsw
      ntgl  = Model%ntgl
      ncndl = min(Model%ncnd,4)

      LP1 = LM + 1               ! num of in/out levels


!  --- ...  set local /level/layer indexes corresponding to in/out
!  variables

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


!> -# Prepare atmospheric profiles for radiation input.
!

      lsk = 0
      if (ivflip == 0 .and. lm < Model%levs) lsk = Model%levs - lm

!           convert pressure unit from pa to mb
      do k = 1, LM
        k1 = k + kd
        k2 = k + lsk
        do i = 1, IM
          plvl(i,k1+kb) = Statein%prsi(i,k2+kb) * 0.01   ! pa to mb (hpa)
          plyr(i,k1)    = Statein%prsl(i,k2)    * 0.01   ! pa to mb (hpa)
          tlyr(i,k1)    = Statein%tgrs(i,k2)
          prslk1(i,k1)  = Statein%prslk(i,k2)

!>  - Compute relative humidity.
          es  = min( Statein%prsl(i,k2),  fpvs( Statein%tgrs(i,k2) ) )  ! fpvs and prsl in pa
          qs  = max( QMIN, eps * es / (Statein%prsl(i,k2) + epsm1*es) )
          rhly(i,k1) = max( 0.0, min( 1.0, max(QMIN, Statein%qgrs(i,k2,1))/qs ) )
          qstl(i,k1) = qs
        enddo
      enddo

      !--- recast remaining all tracers (except sphum) forcing them all to be positive
      do j = 2, NTRAC
        do k = 1, LM
          k1 = k + kd
          k2 = k + lsk
          tracer1(:,k1,j) = max(0.0, Statein%qgrs(:,k2,j))
        enddo
      enddo
!
      if (ivflip == 0) then                                ! input data from toa to sfc
        do i = 1, IM
          plvl(i,1+kd) = 0.01 * Statein%prsi(i,1)          ! pa to mb (hpa)
        enddo
        if (lsk /= 0) then
          do i = 1, IM
            plvl(i,1+kd)  = 0.5 * (plvl(i,2+kd) + plvl(i,1+kd))
          enddo
        endif
      else                                                 ! input data from sfc to top
        do i = 1, IM
          plvl(i,LP1+kd) = 0.01 * Statein%prsi(i,LP1+lsk)  ! pa to mb (hpa)
        enddo
        if (lsk /= 0) then
          do i = 1, IM
            plvl(i,LM+kd)  = 0.5 * (plvl(i,LP1+kd) + plvl(i,LM+kd))
          enddo
        endif
      endif

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


!>  - Get layer ozone mass mixing ratio (if use ozone climatology data,
!!    call getozn()).

!MZ HAFS 's ozone data is adapeted from CAM radiation
!      if (Model%ntoz > 0) then            ! interactive ozone generation
!        do k=1,lmk
!          do i=1,im
!            olyr(i,k) = max( QMIN, tracer1(i,k,Model%ntoz) )
!          enddo
!        enddo
!      else                                ! climatological ozone
!        call getozn (prslk1, Grid%xlat, IM, LMK,    &     !  ---  inputs
!                     olyr)                                !  ---  outputs
!      endif                               ! end_if_ntoz

!-- ww: attempt to use CAM ozone and some aerosol profiles in all rad schemes
!   note that CAM option will still do the same. 
!   n_ozmixm: no of months; levsiz: = 59, vertical dim
!   Read in CAM ozone data, and interpolate data to model grid
!   Interpolation of aerosols is done on domain 1 only

!MZ      IF ( config_flags%o3input .EQ. 2 ) THEN
      if(Model%hafs_o3input .EQ.2 ) then
        CALL oznini(ozmixm,pin,levsiz,n_ozmixm,XLAT,                &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte)
      endif


!>  - Call coszmn(), to compute cosine of zenith angle (only when SW is called)
      if (Model%lsswr) then
        call coszmn (Grid%xlon,Grid%sinlat,           &     !  ---  inputs
                     Grid%coslat,Model%solhr, IM, me, &
                     Radtend%coszen, Radtend%coszdg)        !  --- outputs
      endif

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
!  - gasvmr(:,:,10) -  cfc113 volumne mixing ratio

!  --- ...  set up non-prognostic gas volume mixing ratioes

      call getgases (plvl, Grid%xlon, Grid%xlat, IM, LMK,  & !  --- inputs
                     gasvmr)                                 !  --- outputs

!CCPP: re-assign gasvmr(:,:,NF_VGAS) to gasvmr_X(:,:)
      do k = 1, LMK
        do i = 1, IM
           gasvmr_co2    (i,k)  = gasvmr(i,k,1)
           gasvmr_n2o    (i,k)  = gasvmr(i,k,2)
           gasvmr_ch4    (i,k)  = gasvmr(i,k,3)
           gasvmr_o2     (i,k)  = gasvmr(i,k,4)
           gasvmr_co     (i,k)  = gasvmr(i,k,5)
           gasvmr_cfc11  (i,k)  = gasvmr(i,k,6)
           gasvmr_cfc12  (i,k)  = gasvmr(i,k,7)
           gasvmr_cfc22  (i,k)  = gasvmr(i,k,8)
           gasvmr_ccl4   (i,k)  = gasvmr(i,k,9)
           gasvmr_cfc113 (i,k)  = gasvmr(i,k,10)
         enddo
      enddo

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
          tem2db(i,1) = log( max(prsmin, plvl(i,1)) )
          tem2db(i,LMP) = log( plvl(i,LMP) )
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
            delp(i,k1) = plvl(i,k1+1) - plvl(i,k1)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
            delp(i,lyb) = plvl(i,lla) - plvl(i,llb)
          enddo
        endif

        do k = 2, LMK
          do i = 1, IM
            tlvl(i,k) = tlyr(i,k) + (tlyr(i,k-1) - tlyr(i,k))           &
     &                * (tem2db(i,k)   - tem2da(i,k))                   &
     &                / (tem2da(i,k-1) - tem2da(i,k))
          enddo
        enddo

!  ---  ...  level height and layer thickness (km)

        tem0d = 0.001 * rog
        do i = 1, IM
          do k = 1, LMK
            dz(i,k) = tem0d * (tem2db(i,k+1) - tem2db(i,k)) * tvly(i,k)
          enddo

!         hz(i,LMP) = 0.0
!         do k = LMK, 1, -1
!           hz(i,k) = hz(i,k+1) + dz(i,k)
!         enddo
        enddo

      else                               ! input data from sfc to toa

        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = log( plvl(i,1) )
          tem2db(i,LMP) = log( max(prsmin, plvl(i,LMP)) )
          tsfa  (i)   = tlyr(i,1)                    ! sfc layer air temp
          tlvl(i,1)   = tskn(i)
          tlvl(i,LMP) = tlyr(i,LMK)
        enddo

        do k = LM, 1, -1
          do i = 1, IM
            qlyr(i,k) = max( tem1d(i), Statein%qgrs(i,k,1) )
            tem1d(i)  = min( QME5, qlyr(i,k) )
            tvly(i,k) = Statein%tgrs(i,k) * (1.0 + fvirt*qlyr(i,k)) ! virtual T (K)
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
            delp(i,lyb) = plvl(i,lla) - plvl(i,llb)
          enddo
        endif

        do k = 1, LMK-1
          do i = 1, IM
            tlvl(i,k+1) = tlyr(i,k) + (tlyr(i,k+1) - tlyr(i,k))         &
     &                  * (tem2db(i,k+1) - tem2da(i,k))                 &
     &                  / (tem2da(i,k+1) - tem2da(i,k))
          enddo
        enddo

!  ---  ...  level height and layer thickness (km)

        tem0d = 0.001 * rog
        do i = 1, IM
          do k = LMK, 1, -1
            dz(i,k) = tem0d * (tem2db(i,k) - tem2db(i,k+1)) * tvly(i,k)
          enddo

!         hz(i,1) = 0.0
!         do k = 1, LMP
!           hz(i,k+1) = hz(i,k) + dz(i,k)
!         enddo
        enddo

      endif                              ! end_if_ivflip

!>  - Call module_radiation_aerosols::setaer(),to setup aerosols
!! property profile for radiation.

!check  print *,' in grrad : calling setaer '

      call setaer (plvl, plyr, prslk1, tvly, rhly, Sfcprop%slmsk, & !  ---  inputs
                   tracer1, Grid%xlon, Grid%xlat, IM, LMK, LMP,   &
                   Model%lsswr, Model%lslwr,                      &
                   faersw, faerlw, aerodp)                          !  ---  outputs

! CCPP
      do j = 1,NBDSW
        do k = 1, LMK
          do i = 1, IM
            ! NF_AESW = 3
            faersw1(i,k,j) = faersw(i,k,j,1)
            faersw2(i,k,j) = faersw(i,k,j,2)
            faersw3(i,k,j) = faersw(i,k,j,3)
          enddo
        enddo
       enddo

      do j = 1,NBDLW
        do k = 1, LMK
          do i = 1, IM
            ! NF_AELW = 3
            faerlw1(i,k,j) = faerlw(i,k,j,1)
            faerlw2(i,k,j) = faerlw(i,k,j,2)
            faerlw3(i,k,j) = faerlw(i,k,j,3)
          enddo
        enddo
       enddo

!>  - Obtain cloud information for radiation calculations
!!    (clouds,cldsa,mtopa,mbota)
!!\n   for  prognostic cloud:
!!    - For Zhao/Moorthi's prognostic cloud scheme,
!!      call module_radiation_clouds::progcld1()
!!    - For Zhao/Moorthi's prognostic cloud+pdfcld,
!!      call module_radiation_clouds::progcld3()
!!      call module_radiation_clouds::progclduni() for unified cloud and ncld=2

!  --- ...  obtain cloud information for radiation calculations

        ccnd = 0.0_kind_phys
          do k=1,LMK
            do i=1,IM
              ccnd(i,k,1) = tracer1(i,k,ntcw)                     ! liquid water
              ccnd(i,k,2) = tracer1(i,k,ntiw)                     ! ice water
              ccnd(i,k,3) = tracer1(i,k,ntrw)                     ! rain water
              ccnd(i,k,4) = tracer1(i,k,ntsw)                     ! snow water
            enddo
          enddo

        do n=1,ncndl
          do k=1,LMK
            do i=1,IM
              if (ccnd(i,k,n) < epsq) ccnd(i,k,n) = 0.0
            enddo
          enddo
        enddo

!
!  --- add suspended convective cloud water to grid-scale cloud water
!      only for cloud fraction & radiation computation
!      it is to enhance cloudiness due to suspended convec cloud water
!      for zhao/moorthi's (imp_phys=99) &
!          ferrier's (imp_phys=5) microphysics schemes

        if ((Model%num_p3d == 4) .and. (Model%npdf3d == 3)) then       ! same as Model%imp_physics = 99
          do k=1,lm
            k1 = k + kd
            do i=1,im
              deltaq(i,k1) = Tbd%phy_f3d(i,k,5)
              cnvw  (i,k1) = Tbd%phy_f3d(i,k,6)
              cnvc  (i,k1) = Tbd%phy_f3d(i,k,7)
            enddo
          enddo
        elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then ! same as MOdel%imp_physics=98
          do k=1,lm
            k1 = k + kd
            do i=1,im
              deltaq(i,k1) = 0.0
              cnvw  (i,k1) = Tbd%phy_f3d(i,k,Model%num_p3d+1)
              cnvc  (i,k1) = 0.0
            enddo
          enddo
        else                                                           ! all the rest
          do k=1,lmk
            do i=1,im
              deltaq(i,k) = 0.0
              cnvw  (i,k) = 0.0
              cnvc  (i,k) = 0.0
            enddo
          enddo
        endif

        if (lextop) then
          do i=1,im
            cldcov(i,lyb) = cldcov(i,lya)
            deltaq(i,lyb) = deltaq(i,lya)
            cnvw  (i,lyb) = cnvw  (i,lya)
            cnvc  (i,lyb) = cnvc  (i,lya)
          enddo
          if (Model%effr_in) then
            do i=1,im
              effrl(i,lyb) = effrl(i,lya)
              effri(i,lyb) = effri(i,lya)
              effrr(i,lyb) = effrr(i,lya)
              effrs(i,lyb) = effrs(i,lya)
            enddo
          endif
        endif

        if(Model%imp_physics == 15) then

!MZ*HWRF phys/module_radiation_driver.F
      ! This allows radiation schemes (mainly HWRF) to correctly
      ! interface with the convection scheme when convection is only
      ! enabled in some domains:

!      if(present(explicit_convection)) then
!         expl_conv=explicit_convection
!      else
!         expl_conv=.true. ! backward compatibility for ARW
!      endif

!HWRF         IF (PRESENT(dxkm)) then
!HWRF            gridkm = 1.414*SQRT(dxkm*dxkm + dy*0.001*dy*0.001)
!HWRF         ELSE IF (PRESENT(dx)) then
!            gridkm = SQRT(dx*0.001*dx*0.001 + dy*0.001*dy*0.001)
             gridkm = SQRT(dx(1)*0.001*dx(1)*0.001 + dx(1)*0.001*dx(1)*0.001)
!HWRF         endif

         if (Model%kdt .LE. 2) then
           WRITE ( 0 , * ) 'Grid spacing in km ',  dx(1), gridkm
         endif

! Choose how to compute cloud fraction (since 3.6)
! Initialize to zero 
      do k=1,lmk
      do i=1,im
         cldcov(i,k) = 0.
      END DO
      END DO


!MZ* ICLOUD=3 in hwrf
!MZ* adapted from module_radiation_driver.F
!+---+-----------------------------------------------------------------+ 
!..New cloud fraction scheme added by G. Thompson (2014Oct31)
!+---+-----------------------------------------------------------------+ 
!      ELSEIF (ICLOUD == 3) THEN
!        IF (PRESENT(CLDFRA) .AND. &
!                           PRESENT(F_QC) .AND. PRESENT ( F_QI ) ) THEN

!MZ* after HWRF rrtmg sw/lw, back to micro-only qc and qi
!           DO k = 1,lmk   !kts,kte
!           DO i = 1,im    !its,ite
!              qc_save(i,k) = ccnd(i,k,1)    !qc(i,k,j)
!              qi_save(i,k) = ccnd(i,k,2)    !qi(i,k,j)
!           ENDDO
!           ENDDO

!           DO k = 1,lmk    !kts,kte
!           DO i = 1,im     !its,ite
!               qs_save(i,k,j) = ccnd(i,k,4)   !qs(i,k,j)
!           ENDDO
!           ENDDO
            
          !MZ*
           DO I = 1,im
             sm(i) = 1. 
             if (Sfcprop%slmsk(i) > 0.5) sm(i)=0.
             xland(i)= sm(i)+1.
           ENDDO

          CALL cal_cldfra3(cldcov(:,1:LMK),qlyr,ccnd(:,:,1),ccnd(:,:,2),&
     &                      ccnd(:,:,4),plyr,tlyr,rhly,xland, gridkm,   &
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
     &                 its,ite, jts,jte, kts,kte)

        ELSE

         write(errmsg,'(*(a))') "Can not use Thompson cloud fraction    &
                                 option, missing QC or QI field."
         errflg = 1
         return
       
        ENDIF


!     Interpolating climatological ozone and aerosol to model time and levels
!     Adapted from camrad code

       IF ( o3input .EQ. 2 ) THEN
! interpolate to model time-step
        call ozn_time_int(julday,julian,ozmixm,ozmixt,levsiz,n_ozmixm,&
                              ids , ide , jds , jde , kds , kde ,     &
                              ims , ime , jms , jme , kms , kme ,     &
                              its , ite , jts , jte , kts , kte )

! interpolate to model model levels
        call ozn_p_int(p ,pin, levsiz, ozmixt, o3rad, &
                              ids , ide , jds , jde , kds , kde ,     &
                              ims , ime , jms , jme , kms , kme ,     &
                              its , ite , jts , jte , kts , kte )
       ELSE 
        write(errmsg,'(*(a))') "HAFS RRMTG uses CAM ozone data    &
                                 , reset o3input = 2 "
        errflg = 1

      ENDIF

     IF ( PRESENT( AEROD ) ) THEN
     IF ( aer_opt .EQ. 1 .AND. id .EQ. 1 ) THEN
        call aer_time_int(julday,julian,aerodm,aerodt,alevsiz,n_ozmixm-1,no_src_types, &
                              ids , ide , jds , jde , kds , kde ,     &
                              ims , ime , jms , jme , kms , kme ,     &
                              its , ite , jts , jte , kts , kte )

        call aer_p_int(p ,pina, alevsiz, aerodt, aerod, no_src_types, p8w, AODTOT, &
                              ids , ide , jds , jde , kds , kde ,     &
                              ims , ime , jms , jme , kms , kme ,     &
                              its , ite , jts , jte , kts , kte )
     ENDIF
     ENDIF



!!!MZ* works for GFS_RRTMG
!          if (Model%kdt == 1) then
!            Tbd%phy_f3d(:,:,Model%nleffr) = 10.
!            Tbd%phy_f3d(:,:,Model%nieffr) = 50.
!            Tbd%phy_f3d(:,:,Model%nseffr) = 250.
!          endif

!          call progcld5 (plyr,plvl,tlyr,qlyr,qstl,rhly,tracer1,     &  !  --- inputs
!                         Grid%xlat,Grid%xlon,Sfcprop%slmsk,dz,delp, &
!                         ntrac-1, ntcw-1,ntiw-1,ntrw-1,             &
!                         ntsw-1,ntgl-1,                             &
!                         im, lmk, lmp, Model%uni_cld,               &
!                         Model%lmfshal,Model%lmfdeep2,              &
!                         cldcov(:,1:LMK),Tbd%phy_f3d(:,:,1),        &
!                         Tbd%phy_f3d(:,:,2), Tbd%phy_f3d(:,:,3),    &
!                         clouds,cldsa,mtopa,mbota, de_lgth)            !  --- outputs

!! CCPP
!       do k = 1, LMK
!         do i = 1, IM
!            clouds1(i,k)  = clouds(i,k,1)
!            clouds2(i,k)  = clouds(i,k,2)
!            clouds3(i,k)  = clouds(i,k,3)
!            clouds4(i,k)  = clouds(i,k,4)
!            clouds5(i,k)  = clouds(i,k,5)
!            clouds6(i,k)  = clouds(i,k,6)
!            clouds7(i,k)  = clouds(i,k,7)
!            clouds8(i,k)  = clouds(i,k,8)
!            clouds9(i,k)  = clouds(i,k,9)
!         enddo
!       enddo
!!MZ*



!MZ      IF ( config_flags%aer_opt .EQ. 1 .AND. id .EQ. 1 ) THEN
!      CALL aerosol_in(aerodm,pina,alevsiz,n_ozmixm-1,no_src_types,  &
!                     XLAT,XLONG,                                    & 
!                     ids, ide, jds, jde, kds, kde,                  &
!                     ims, ime, jms, jme, kms, kme,                  &
!                     its, ite, jts, jte, kts, kte)

       end subroutine HAFS_rrtmg_pre_run

!> \section arg_table_HAFS_rrtmg_pre_finalize Argument Table
!!
      subroutine HAFS_rrtmg_pre_finalize ()
      end subroutine HAFS_rrtmg_pre_finalize

!! @}
      end module HAFS_rrtmg_pre
