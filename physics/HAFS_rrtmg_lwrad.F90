!>\file HAFS_rrtmg_lwrad.F90
!! This file contains HWRF RRTMG LWRAD scheme.

!>\defgroup HAFSLWRAD HAFS RRTMG LWRAD Module
!! This module contains the HAFS RRTMG LWRAD scheme.
module HAFS_rrtmg_lwrad

      use machine, only : kind_phys

      use module_ra_rrtmg_lw, only : rrtmg_lwrad, rrtmg_lwinit

      implicit none

      public :: HAFS_rrtmg_lwrad_init, HAFS_rrtmg_lwrad_run, HAFS_rrtmg_lwrad_finalize

      private

      logical :: aclwalloc = .False.

   contains

!> This subroutine is a wrapper around the actual rrtmg_lw_init() in HWRF.
#if 0
!! \section arg_table_HAFS_rrtmg_lwrad_init Argument Table
!! \htmlinclude HAFS_rrtmg_lwrad_init.html
!!
#endif
        subroutine HAFS_rrtmg_lwrad_init (ltp,mpirank,mpiroot           &
      &                                   , Grid, Model,                &
!MZ      &                                   ,allowed_to_read,             &
      &                                   ,ilw_physics                  &
      &                                   ,ilw_hafs_rrtmg               &
      &                                   ,re_cloud, re_ice, re_snow    & 
      &                                   ,has_reqc, has_reqi, has_reqs &
      &                                   ,rthraten,rthratenlw,cldfra   &
      &                                   ,cldfra_dp, cldfra_sh, cldfra_old, &
      &                                   ,errmsg, errflg)

         implicit none

         type(GFS_control_type),              intent(in)    :: Model
         type(GFS_grid_type),                 intent(in)    :: Grid

         ! Interface variables
         integer,                   intent(in)    :: mpicomm
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         integer,                   intent(in)    :: ilw_physics
         integer,                   intent(in)    :: ilw_hafs_rrtmg
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                             INTENT(INOUT) :: re_cloud, &
                                                       re_ice, re_snow
         INTEGER, INTENT(INOUT):: has_reqc, has_reqi, has_reqs
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                                     intent(inout) ::   &
                                                              RTHRATEN, & !< Theta tendency due to radiation (K/s)
                                                            RTHRATENLW, & !< Theta tendency due to long wave radiation (K/s) 
                                                            CLDFRA
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                           OPTIONAL, INTENT(INOUT) ::   &
                                                             CLDFRA_OLD

         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP) , &
                                           OPTIONAL, INTENT(INOUT)   :: & ! ckay for subgrid cloud
                                                    cldfra_dp, cldfra_sh

         integer,                   intent(in)    :: ltp
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Local variables: dimensions used in rrtmg_lwinit
         integer               :: ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte
         logical  :: allowed_to_read

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (aclwalloc) return

         if (mpirank==mpiroot) then
            write(0,*) ' --------------------------------------------------------'
            write(0,*) ' ---                   WARNING                        ---' 
            write(0,*) ' ---  the CCPP HAFS RRTMG LWRAD scheme is currently   ---'
            write(0,*) ' ---  under development, use at your own risk         ---'
            write(0,*) ' ---                   WARNING                        ---'
            write(0,*) ' --------------------------------------------------------'
         end if

         if (ilw_physics /= ilw_hafs_rrtmg) then
            write(errmsg,'(*(a))') "Logic error: namelist choice of RADLW is different from HAFS RRTMG LWRAD"
            errflg = 1
            return
         end if


!   jtf=min0(jte,jde-1)
!   ktf=min0(kte,kde-1)
!   itf=min0(ite,ide-1)
!---------------------------------------------------------------------

!-- calculate radiation time step

    !STEPRA = nint(RADT*60./DT)
    !STEPRA = max(STEPRA,1)

!-- initialization

         ! Set internal dimensions
         ids = 1
         ims = 1
         its = 1
         ide = ncol
         ime = ncol
         ite = ncol
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

!..Fill initial starting values of radiative effective radii for
!.. cloud water (2.51 microns), cloud ice (5.01 microns), and
!.. snow (10.01 microns).
    if (has_reqc.ne.0) then
          do k=kts,kte !ktf
          do i=its,ite !itf
             re_cloud(i,k) = 2.51E-6
          end do
          end do
    endif
    if (has_reqi.ne.0) then
          do k=kts,kte !ktf
          do i=its,ite !itf
             re_ice(i,k) = 5.01E-6
          end do
          end do
    endif
    if (has_reqs.ne.0) then
          do k=kts,kte !ktf
          do i=its,ite !itf
             re_snow(i,k) = 10.01E-6
          end do
          end do
    endif

!   IF(start_of_simulation)THEN
     DO k=kts,kte   !ktf
     DO i=its,ite   !itf
        RTHRATEN(i,k)=0.
        RTHRATENLW(i,k)=0.
!MZ* SW        RTHRATENSW(i,k)=0.
        CLDFRA(i,k)=0.
     ENDDO
     ENDDO
!     ENDDO

     IF( PRESENT(cldfra_dp) ) THEN
        DO k=kts,kte !ktf
        DO i=its,ite !itf
           cldfra_dp(i,k)=0.
           cldfra_sh(i,k)=0.
        ENDDO
        ENDDO
     ENDIF

     if( present(cldfra_old) ) then
        DO k=kts,kte !ktf
        DO i=its,ite !itf
           cldfra_old(i,k) = 0.
        ENDDO
        ENDDO
     end if

!> - Call rrtmg_lwinit() to calculte:
!!  - calculate nlayers = model levels + new levels
!!  - call rrtmg_lwlookuptable() to read in absorption coefficients and other
!! data
!!  - call rrtmg_lw_ini() to perform g-point reduction and other initializations

     
!MZ* 
        allowed_to_read = .true.
         CALL rrtmg_lwinit(                             &
                       ltp, allowed_to_read,            &
                       ids, ide, jds, jde, kds, kde,    &
                       ims, ime, jms, jme, kms, kme,    &
                       its, ite, jts, jte, kts, kte     )

          if (errflg /= 0) return

          aclwalloc = .true.

      end subroutine HAFS_rrtmg_lwrad_init


#if 0
!> \section arg_table_HAFS_rrtmg_lwrad_run Argument Table
!! \htmlinclude HAFS_rrtmg_lwrad_run.html
!!
#endif
!>\ingroup HAFS_lwrad_scheme
!>\section gen_HAFS_RADLW HAFS RRTMG LWRAD General Algorithm
!>@{

      subroutine HAFS_rrtmg_lwrad_run(ncol, nlay, nlp1, RTHRATEN,       &
                                      RTHRATENLW, HRLWPD,               &
                    LWUPT, LWUPTC, LWUPTCLN, LWDNT, LWDNTC, LWDNTCLN,   &
                    LWUPB, LWUPBC, LWUPBCLN, LWDNB, LWDNBC, LWDNBCLN,   &
                    GLW, RLWTOA, LWCF, EMISS, p8w, p, pi, dz8w, tsk, t, &
                    t8w, rho, r_d, g, icloud, warm_rain, CLDFRA,        &
                    cldovrlp,                                           & ! J. Henderson AER: cldovrlp namelist value
!#if (EM_CORE == 1)
!                    lradius, iradius,                                   &
!#endif
                    is_cammgmp_used,                                    &
                    F_ICE_PHY, F_RAIN_PHY,                              &
                    XLAND, XICE, SNOW,                                  &
                    QV, QC, QR,                                         &
                    QI, QS, QG,                                         &
                    O3INPUT, O3RAD,                                     &
                    F_QV, F_QC, F_QR,                                   &
                    F_QI, F_QS, F_QG,                                   &
                    re_cloud, re_ice, re_snow,                          & ! G. Thompson
                    has_reqc, has_reqi, has_reqs,                       & ! G. Thompson
!#if ( WRF_CHEM == 1 )
!                    tauaerlw1, tauaerlw2, tauaerlw3, tauaerlw4,         & ! jcb
!                    tauaerlw5, tauaerlw6, tauaerlw7, tauaerlw8,         & ! jcb
!                    tauaerlw9, tauaerlw10,tauaerlw11,tauaerlw12,        & ! jcb
!                    tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16,        & ! jcb
!                    aer_ra_feedback, progn,                             &
!#endif
                    calc_clean_atm_diag, qndrop, f_qndrop,              &
                    YR, JULIAN,                                         & ! Added for time-varying trace gases.
!mz                    ids,ide, jds,jde, kds,kde,                          &
!MZ                    ims,ime, jms,jme, kms,kme,                          &
!                    its,ite, jts,jte, kts,kte,                          &
                    LWUPFLX,LWUPFLXC,                                   & !MZ: rad_reset
                    LWDNFLX,LWDNFLXC,                                   &
                    mp_physics,                                         & 
                    mpicomm, mpirank, mpiroot,                          &
                    errmsg, errflg)

         implicit none

         type(GFS_control_type),              intent(in)    :: Model
         type(GFS_grid_type),                 intent(in)    :: Grid
         type(GFS_sfcprop_type),              intent(in)    :: Sfcprop
         type(GFS_statein_type),              intent(in)    :: Statein
         type(GFS_radtend_type),              intent(inout) :: Radtend
         type(GFS_tbd_type),                  intent(in)    :: Tbd
         type(GFS_cldprop_type),              intent(in)    :: Cldprop
         type(GFS_coupling_type),             intent(in)    :: Coupling


         ! Interface variables

         ! Dimensions and constants
         integer,                   intent(in   ) :: ncol, nlay,nlp1

         ! module_radiation_driver.F90
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                                     intent(inout) ::   &
                                                              RTHRATEN, & !< Theta tendency due to radiation (K/s)
                                                            RTHRATENLW, & !< Theta tendency due to long wave radiation (K/s) 
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                                     intent(out) ::     &
                                                              HRLWPD 

         ! TOA and surface, upward and downward, total, clear (no cloud), and clean (no aerosol) fluxes
         real(kind_phys), dimension(size(Grid%xlon,1)), OPTIONAL,       &
                                                      INTENT(INOUT) ::  &
                    LWUPT,  LWUPTC, LWUPTCLN,  LWDNT,  LWDNTC, LWDNTCLN,&
                    LWUPB,  LWUPBC, LWUPBCLN,  LWDNB,  LWDNBC, LWDNBCLN
         real(kind_phys), dimension(size(Grid%xlon,1)), intent(inout) :: glw !< downward long wave flux at ground surface (W/m^2)
         real(kind_phys), dimension(size(Grid%xlon,1)), intent(inout) :: rlwtoa !< upward long wave at top of atmosphere (W/m^2)
         real(kind_phys), dimension(size(Grid%xlon,1)), optional, intent(inout) :: LWCF
         real(kind_phys), dimension(size(Grid%xlon,1)), intent(in)    :: EMISS  !< surface emissivity (between 0 and 1)
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                                     intent(in) :: p8w, &   !< pressure at full levels (Pa)
                                                                     p, &
                                                                    pi, &
                                                                  dz8w, &   !< dz between full levels (m)
                                                                     t, &   
                                                                    t8w, &  !< temperature at full levels (K)
                                                                    rho
         real(kind_phys), intent(in) :: r_d  !< gas constant for dry air (J/kg/K)
         real(kind_phys), intent(in) :: g    !< acceleration due to gravity (m/s^2)
         integer,         intent(in) :: icloud
         logical,         intent(in) :: warm_rain
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                    optional, intent(inout) ::  CLDFRA !< cloud fraction (between 0 and 1)
         integer,         intent(in) :: cldovrlp      ! < J. Henderson AER: cldovrlp namelist value
         logical,         intent(in) :: IS_CAMMGMP_USED ! BSINGH:01/31/2013: Added for CAM5 RRTMG

         integer, intent(inout) :: has_reqc, has_reqi, has_reqs
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                                             optional,  &
                                               intent(in) :: f_ice_phy, &
                                                            f_rain_phy 

         ! Upward and downward, total and clear sky layer fluxes (W m-2)
         REAL, DIMENSION(size(Grid%xlon,1),Model%levr+LTP+2),           &
                                           OPTIONAL, INTENT(INOUT) ::   &
                                      LWUPFLX,LWUPFLXC,LWDNFLX,LWDNFLXC

         INTEGER, INTENT(IN)            ::  mp_physics

         ! Cloud effective radii
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                              intent(out) :: re_cloud,  &
                                                             re_ice,    &
                                                             re_snow

         ! Surface property
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                                  intent(in) :: xland,  &
                                                                xice,   &
                                                                tsk,    &
                                                                snow

         ! Cloud mixing ratio
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                                            optional,   &
                                                  intent(inout) :: qv,  & !< water vapor mixing ratio (kg/kg)
                                                                   qc,  & !< cloud water mixing ratio (kg/kg) 
                                                                   qr,  & !< rain water mixing ratio (kg/kg) 
                                                                   qi,  & !< cloud ice mixing ratio (kg/kg) 
                                                                   qs,  & !< snow mixing ratio (kg/kg)
                                                                   qg,  & !< graupel mixing ratio 
                                                                   qndrop !< cloud droplet number (#/kg)
         integer, intent(in)  :: o3input
         real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),  &
                                         optional, intent(in) :: o3rad

         logical, optional :: f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qndrop
         INTEGER, INTENT(IN   )  ::   calc_clean_atm_diag
         INTEGER, INTENT(IN  ),OPTIONAL ::                          YR
         REAL(kind_phys), INTENT(IN  )   ::    julian

         ! MPI information
         integer,                   intent(in)    :: mpicomm
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Local variables

         ! Air density
         real(kind_phys) :: rho(1:ncol,1:nlev)              !< kg m-3
         ! Hydrometeors
         real(kind_phys) :: qv_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qc_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qr_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qi_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qs_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         real(kind_phys) :: qg_mp(1:ncol,1:nlev)            !< kg kg-1 (dry mixing ratio)
         ! Vertical velocity and level width
         real(kind_phys) :: w(1:ncol,1:nlev)                !< m s-1
         real(kind_phys) :: dz(1:ncol,1:nlev)               !< m
         ! Rain/snow/graupel fall amounts
         real(kind_phys) :: rain_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: graupel_mp(1:ncol)              ! mm, dummy, not used
         real(kind_phys) :: ice_mp(1:ncol)                  ! mm, dummy, not used
         real(kind_phys) :: snow_mp(1:ncol)                 ! mm, dummy, not used
         real(kind_phys) :: delta_rain_mp(1:ncol)           ! mm
         real(kind_phys) :: delta_graupel_mp(1:ncol)        ! mm
         real(kind_phys) :: delta_ice_mp(1:ncol)            ! mm
         real(kind_phys) :: delta_snow_mp(1:ncol)           ! mm
         ! Radar reflectivity
         !logical         :: diagflag                        ! must be true if do_radar_ref is true, not used otherwise
         !integer         :: do_radar_ref_mp                 ! integer instead of logical do_radar_ref
         ! Effective cloud radii
         !logical         :: do_effective_radii
         !integer         :: has_reqc
         !integer         :: has_reqi
         !integer         :: has_reqs
         ! Dimensions used in mp_gt_driver
         integer         :: ids,ide, jds,jde, kds,kde, &
                            ims,ime, jms,jme, kms,kme, &
                            its,ite, jts,jte, kts,kte

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Check initialization state
         if (.not.aclwalloc) then
            write(errmsg, fmt='((a))') 'HAFS_rrtmg_lwrad_run() called before HAFS_rrtmg_lwrad_init()'
            errflg = 1
            return
         end if

         !> - Convert specific humidity/moist mixing ratios to dry mixing ratios
         !MZ* qv in FV3 is specific humidity; RRTMG needs mixing ratio
         !qv_mp = spechum/(1.0_kind_phys-spechum)
         do i = 
         do k =
         qv_mp (i,k) = qv (i,k)/(1.0_kind_phys - qv (i,k))
         qc_mp (i,k) = qc(i,k)/(1.0_kind_phys-spechum)
         qr_mp = qr/(1.0_kind_phys-spechum)
         qi_mp = qi/(1.0_kind_phys-spechum)
         qs_mp = qs/(1.0_kind_phys-spechum)
         qg_mp = qg/(1.0_kind_phys-spechum)

         end do
         end do

        ! if (is_aerosol_aware .and. .not. (present(nc)     .and. &
        !                                   present(nwfa)   .and. &
        !                                   present(nifa)   .and. &
        !                                   present(nwfa2d) .and. &
        !                                   present(nifa2d)       )) then
        !    write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_run:',  &
        !                               ' aerosol-aware microphysics require all of the', &
        !                               ' following optional arguments:', &
        !                               ' nc, nwfa, nifa, nwfa2d, nifa2d'
        !    errflg = 1
        !    return
        ! end if

         !> - Density of air in kg m-3
        ! rho = prsl/(con_rd*tgrs)

         !> - Convert omega in Pa s-1 to vertical velocity w in m s-1
        ! w = -omega/(rho*con_g)

         !> - Layer width in m from geopotential in m2 s-2
         !dz = (phii(:,2:nlev+1) - phii(:,1:nlev)) / con_g

         ! Accumulated values inside Thompson scheme, not used;
         ! only use delta and add to inout variables (different units)
         !rain_mp          = 0
         !graupel_mp       = 0
         !ice_mp           = 0
         !snow_mp          = 0
         !delta_rain_mp    = 0
         !delta_graupel_mp = 0
         !delta_ice_mp     = 0
         !delta_snow_mp    = 0

         ! Flags for calculating radar reflectivity; diagflag is redundant
         !if (do_radar_ref) then
         !    diagflag = .true.
         !    do_radar_ref_mp = 1
         !else
         !    diagflag = .false.
         !    do_radar_ref_mp = 0
         !end if

         !if (present(re_cloud) .and. present(re_ice) .and. present(re_snow)) then
         !    do_effective_radii = .true.
         !    has_reqc = 1
         !    has_reqi = 1
         !    has_reqs = 1
         !    ! Initialize to zero, intent(out) variables
         !    re_cloud = 0
         !    re_ice   = 0
         !    re_snow  = 0
         !else if (.not.present(re_cloud) .and. .not.present(re_ice) .and. .not.present(re_snow)) then
         !    do_effective_radii = .false.
         !    has_reqc = 0
         !    has_reqi = 0
         !    has_reqs = 0
         !else
         !    write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_run:',  &
         !                               ' all or none of the following optional', &
         !                               ' arguments are required: re_cloud, re_ice, re_snow'
         !    errflg = 1
         !!    return
         !end if

         ! Set internal dimensions
         ids = 1
         ims = 1
         its = 1
         ide = ncol
         ime = ncol
         ite = ncol
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

!MZ* cal_cldfra3 need to insert here-------------  

             CALL RRTMG_LWRAD(                                      &
                  RTHRATENLW=RTHRATEN,                              &
                  HRLWPD=HRLWPD,                                    & ! J. Henderson AER
                  LWUPT=LWUPT,LWUPTC=LWUPTC,LWUPTCLN=LWUPTCLN,      &
                  LWDNT=LWDNT,LWDNTC=LWDNTC,LWDNTCLN=LWDNTCLN,      &
                  LWUPB=LWUPB,LWUPBC=LWUPBC,LWUPBCLN=LWUPBCLN,      &
                  LWDNB=LWDNB,LWDNBC=LWDNBC,LWDNBCLN=LWDNBCLN,      &
                  GLW=GLW,OLR=RLWTOA,LWCF=LWCF,                     &
                  EMISS=EMISS,                                      &
                  P8W=p8w,P3D=p,PI3D=pi,DZ8W=dz8w,TSK=tsk,T3D=t,    &
                  T8W=t8w,RHO3D=rho,R=R_d,G=G,                      &
                  ICLOUD=icloud,WARM_RAIN=warm_rain,CLDFRA3D=CLDFRA,&
                  cldovrlp=cldovrlp,                                & ! J. Henderson AER: cldovrlp namelist value
#if (EM_CORE == 1)
                  LRADIUS=lradius, IRADIUS=iradius,                 &
#endif
                  IS_CAMMGMP_USED=is_cammgmp_used,                  &
!ckay
!                 CLDFRA_KF3D=cldfra_KF,QC_KF3D=qc_KF,QI_KF3D=qi_KF,&
                  F_ICE_PHY=F_ICE_PHY,F_RAIN_PHY=F_RAIN_PHY,        &
                  XLAND=XLAND,XICE=XICE,SNOW=SNOW,                  &
                  QV3D=QV,QC3D=QC,QR3D=QR,                          &
                  QI3D=QI,QS3D=QS,QG3D=QG,                          &
                  O3INPUT=O3INPUT,O33D=O3RAD,                       &
                  F_QV=F_QV,F_QC=F_QC,F_QR=F_QR,                    &
                  F_QI=F_QI,F_QS=F_QS,F_QG=F_QG,                    &
                  RE_CLOUD=re_cloud,RE_ICE=re_ice,RE_SNOW=re_snow,  & ! G. Thompson
                  has_reqc=has_reqc,has_reqi=has_reqi,has_reqs=has_reqs, & ! G. Thompson
#if ( WRF_CHEM == 1 )
                  TAUAERLW1=tauaerlw1,TAUAERLW2=tauaerlw2,          & ! jcb
                  TAUAERLW3=tauaerlw3,TAUAERLW4=tauaerlw4,          & ! jcb
                  TAUAERLW5=tauaerlw5,TAUAERLW6=tauaerlw6,          & ! jcb
                  TAUAERLW7=tauaerlw7,TAUAERLW8=tauaerlw8,          & ! jcb
                  TAUAERLW9=tauaerlw9,TAUAERLW10=tauaerlw10,        & ! jcb
                  TAUAERLW11=tauaerlw11,TAUAERLW12=tauaerlw12,      & ! jcb
                  TAUAERLW13=tauaerlw13,TAUAERLW14=tauaerlw14,      & ! jcb
                  TAUAERLW15=tauaerlw15,TAUAERLW16=tauaerlw16,      & ! jcb
                  aer_ra_feedback=aer_ra_feedback,                  &
!jdfcz            progn=progn,prescribe=prescribe,                   &
                  progn=progn,                                      &
#endif
                  calc_clean_atm_diag=calc_clean_atm_diag,          &
                  QNDROP3D=qndrop,F_QNDROP=f_qndrop,                &
!ccc Added for time-varying trace gases.
                  YR=YR,JULIAN=JULIAN,                              &
!ccc
                  IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde,&
                  IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme,&
                  ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte,&
                  LWUPFLX=LWUPFLX,LWUPFLXC=LWUPFLXC,                &
                  LWDNFLX=LWDNFLX,LWDNFLXC=LWDNFLXC,                &
                  mp_physics=mp_physics                             )


         if (errflg/=0) return

!accumulate_lw_select

!   IF(PRESENT(LWUPTC))THEN
!  NMM calls the driver every RADT time steps, EM calls every DT
!MZ#if (EM_CORE == 1)
!MZ   DTaccum = DT
!MZ#else
!MZ: radt: time for calling radiation (min)
!   DTaccum = RADT*60
!MZ#endif
   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij ,i,j,k,its,ite,jts,jte)

!   DO ij = 1 , num_tiles
!     its = i_start(ij)
!     ite = i_end(ij)
!     jts = j_start(ij)
!     jte = j_end(ij)

!        DO j=jts,jte
!        DO i=its,ite
!           ACLWUPT(I,J) = ACLWUPT(I,J) + LWUPT(I,J)*DTaccum
!           ACLWUPTC(I,J) = ACLWUPTC(I,J) + LWUPTC(I,J)*DTaccum
!           ACLWDNT(I,J) = ACLWDNT(I,J) + LWDNT(I,J)*DTaccum
!           ACLWDNTC(I,J) = ACLWDNTC(I,J) + LWDNTC(I,J)*DTaccum
!           ACLWUPB(I,J) = ACLWUPB(I,J) + LWUPB(I,J)*DTaccum
!           ACLWUPBC(I,J) = ACLWUPBC(I,J) + LWUPBC(I,J)*DTaccum
!           ACLWDNB(I,J) = ACLWDNB(I,J) + LWDNB(I,J)*DTaccum
!           ACLWDNBC(I,J) = ACLWDNBC(I,J) + LWDNBC(I,J)*DTaccum
!        ENDDO
!        ENDDO
!   ENDDO
!   !$OMP END PARALLEL DO
!   ENDIF
!     CASE DEFAULT
!     END SELECT accumulate_lw_select


!mz:leave here to investigate later: SHOULD WE DO THESE TO F-A
         !> - Convert dry mixing ratios to specific humidity/moist mixing ratios
         spechum = qv_mp/(1.0_kind_phys+qv_mp)
         qc      = qc_mp/(1.0_kind_phys+qv_mp)
         qr      = qr_mp/(1.0_kind_phys+qv_mp)
         qi      = qi_mp/(1.0_kind_phys+qv_mp)
         qs      = qs_mp/(1.0_kind_phys+qv_mp)
         qg      = qg_mp/(1.0_kind_phys+qv_mp)


      end subroutine HAFS_rrtmg_lwrad_run
!>@}

      subroutine HAFS_rrtmg_lwrad_finalize()
      end subroutine HAFS_rrtmg_lwrad_finalize

end module HAFS_rrtmg_lwrad
