!> \file module_MYNNPBL_wrapper.F90
!!  This file contains all of the code related to running the MYNN 
!! eddy-diffusivity mass-flux scheme. 

!>\ingroup gsd_mynn_edmf
!> The following references best describe the code within
!!    Olson et al. (2018, NOAA Technical Memorandum)
!!    Nakanishi and Niino (2009 ) \cite NAKANISHI_2009
      MODULE mynnedmf_wrapper

      contains

      subroutine mynnedmf_wrapper_init (lheatstrg, errmsg, errflg)
        implicit none

        logical,          intent(in)  :: lheatstrg
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (lheatstrg) then
          errmsg = 'Logic error: lheatstrg not implemented for MYNN PBL'
          errflg = 1
          return
        end if

      end subroutine mynnedmf_wrapper_init

      subroutine mynnedmf_wrapper_finalize ()
      end subroutine mynnedmf_wrapper_finalize

! \brief This scheme (1) performs pre-mynnedmf work, (2) runs the mynnedmf, and (3) performs post-mynnedmf work
!> \section arg_table_mynnedmf_wrapper_run Argument Table
!! \htmlinclude mynnedmf_wrapper_run.html
!!
SUBROUTINE mynnedmf_wrapper_run(        &
     &  im,levs,                        &
     &  flag_init,flag_restart,         &
     &  lssav, ldiag3d, qdiag3d,        &
     &  lsidea, cplflx,                 &
     &  delt,dtf,dx,zorl,               &
     &  phii,u,v,omega,t3d,             &
     &  qgrs_water_vapor,               &
     &  qgrs_liquid_cloud,              &
     &  qgrs_ice_cloud,                 &
     &  qgrs_cloud_droplet_num_conc,    &
     &  qgrs_cloud_ice_num_conc,        &
     &  qgrs_ozone,                     &
     &  qgrs_water_aer_num_conc,        &
     &  qgrs_ice_aer_num_conc,          &
     &  prsl,exner,                     &
     &  slmsk,tsurf,qsfc,ps,            &
     &  ust,ch,hflx,qflx,wspd,rb,       &
     &  dtsfc1,dqsfc1,                  &
     &  dusfc1,dvsfc1,                  &
     &  dusfci_diag,dvsfci_diag,        &
     &  dtsfci_diag,dqsfci_diag,        &
     &  dusfc_diag,dvsfc_diag,          &
     &  dtsfc_diag,dqsfc_diag,          &
     &  dusfc_cice,dvsfc_cice,          &
     &  dtsfc_cice,dqsfc_cice,          &
     &  hflx_ocn,qflx_ocn,stress_ocn,   &
     &  oceanfrac,fice,wet,icy,dry,     &
     &  dusfci_cpl,dvsfci_cpl,          &
     &  dtsfci_cpl,dqsfci_cpl,          &
     &  dusfc_cpl,dvsfc_cpl,            &
     &  dtsfc_cpl,dqsfc_cpl,            &
     &  recmol,                         &
     &  qke,qke_adv,Tsq,Qsq,Cov,        &
     &  el_pbl,sh3d,exch_h,exch_m,      &
     &  Pblh,kpbl,                      &
     &  qc_bl,qi_bl,cldfra_bl,          &
     &  edmf_a,edmf_w,edmf_qt,          &
     &  edmf_thl,edmf_ent,edmf_qc,      &
     &  sub_thl,sub_sqv,det_thl,det_sqv,&
     &  nupdraft,maxMF,ktop_plume,      &
     &  dudt, dvdt, dtdt,                                  &
     &  dqdt_water_vapor, dqdt_liquid_cloud,               &
     &  dqdt_ice_cloud, dqdt_ozone,                        &
     &  dqdt_cloud_droplet_num_conc, dqdt_ice_num_conc,    &
     &  dqdt_water_aer_num_conc, dqdt_ice_aer_num_conc,    &
     &  flag_for_pbl_generic_tend,                         &
     &  du3dt_PBL, du3dt_OGWD, dv3dt_PBL, dv3dt_OGWD,      &
     &  do3dt_PBL, dq3dt_PBL, dt3dt_PBL,                   &
     &  htrsw, htrlw, xmu,                                 &
     &  grav_settling, bl_mynn_tkebudget, bl_mynn_tkeadvect, &
     &  bl_mynn_cloudpdf, bl_mynn_mixlength,               &
     &  bl_mynn_edmf, bl_mynn_edmf_mom, bl_mynn_edmf_tke,  &
     &  bl_mynn_edmf_part, bl_mynn_cloudmix, bl_mynn_mixqt,&
     &  bl_mynn_output,                                    &
     &  icloud_bl, do_mynnsfclay,                          &
     &  imp_physics, imp_physics_gfdl,                     &
     &  imp_physics_thompson, imp_physics_wsm6,            &
     &  ltaerosol, lprnt, errmsg, errflg  )

! should be moved to inside the mynn:
      use machine , only : kind_phys
!      use funcphys, only : fpvs

      use physcons, only : cp     => con_cp,              &
     &                     g      => con_g,               &
     &                     r_d    => con_rd,              &
     &                     r_v    => con_rv,              &
     &                     cpv    => con_cvap,            &
     &                     cliq   => con_cliq,            &
     &                     Cice   => con_csol,            &
     &                     rcp    => con_rocp,            &
     &                     XLV    => con_hvap,            &
     &                     XLF    => con_hfus,            &
     &                     EP_1   => con_fvirt,           &
     &                     EP_2   => con_eps

      USE module_bl_mynn, only : mynn_bl_driver

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 
!  ---  constant parameters:
!      real(kind=kind_phys), parameter :: rvovrd  = r_v/r_d
!      real(kind=kind_phys), parameter :: karman  = 0.4
!      real(kind=kind_phys), parameter :: XLS     = 2.85E6
!      real(kind=kind_phys), parameter :: p1000mb = 100000.
      real(kind=kind_phys), parameter :: SVP1    = 0.6112
!      real(kind=kind_phys), parameter :: SVP2    = 17.67
!      real(kind=kind_phys), parameter :: SVP3    = 29.65
!      real(kind=kind_phys), parameter :: SVPT0   = 273.15

!   INTEGER , PARAMETER :: param_first_scalar = 1, &
!       &                  p_qc = 2, &
!       &                  p_qr = 0, &
!       &                  p_qi = 2, &
!       &                  p_qs = 0, &
!       &                  p_qg = 0, &
!       &                  p_qnc= 0, &
!       &                  p_qni= 0

!-------------------------------------------------------------------
!For WRF:
!-------------------------------------------------------------------
!  USE module_model_constants, only: &
!       &karman, g, p1000mb, &
!       &cp, r_d, r_v, rcp, xlv, xlf, xls, &
!      &svp1, svp2, svp3, svpt0, ep_1, ep_2, rvovrd, &
!       &cpv, cliq, cice

!  USE module_state_description, only: param_first_scalar, &
!       &p_qc, p_qr, p_qi, p_qs, p_qg, p_qnc, p_qni 

!-------------------------------------------------------------------
!For reference
!   REAL    , PARAMETER :: karman       = 0.4
!   REAL    , PARAMETER :: g            = 9.81
!   REAL    , PARAMETER :: r_d          = 287.
!   REAL    , PARAMETER :: cp           = 7.*r_d/2.
!   REAL    , PARAMETER :: r_v          = 461.6
!   REAL    , PARAMETER :: cpv          = 4.*r_v
!   REAL    , PARAMETER :: cliq         = 4190.
!   REAL    , PARAMETER :: Cice         = 2106.
!   REAL    , PARAMETER :: rcp          = r_d/cp
!   REAL    , PARAMETER :: XLS          = 2.85E6
!   REAL    , PARAMETER :: XLV          = 2.5E6
!   REAL    , PARAMETER :: XLF          = 3.50E5
!   REAL    , PARAMETER :: p1000mb      = 100000.
!   REAL    , PARAMETER :: rvovrd       = r_v/r_d
!   REAL    , PARAMETER :: SVP1         = 0.6112
!   REAL    , PARAMETER :: SVP2         = 17.67
!   REAL    , PARAMETER :: SVP3         = 29.65
!   REAL    , PARAMETER :: SVPT0        = 273.15
!   REAL    , PARAMETER :: EP_1         = R_v/R_d-1.
!   REAL    , PARAMETER :: EP_2         = R_d/R_v
!

  REAL, PARAMETER :: xlvcp=xlv/cp, xlscp=(xlv+xlf)/cp, ev=xlv, rd=r_d, &
       &rk=cp/rd, svp11=svp1*1.e3, p608=ep_1, ep_3=1.-ep_2

  REAL, PARAMETER :: tref=300.0     !< reference temperature (K)
  REAL, PARAMETER :: TKmin=253.0    !< for total water conversion, Tripoli and Cotton (1981)
  REAL, PARAMETER :: tv0=p608*tref, tv1=(1.+p608)*tref, gtr=g/tref, g_inv=1./g

  REAL, PARAMETER :: zero=0.0d0, one=1.0d0
  REAL, PARAMETER :: huge=9.9692099683868690E36 ! NetCDF float FillValue, same as in GFS_typedefs.F90

  character(len=*), intent(out) :: errmsg
  integer, intent(out) :: errflg

  LOGICAL, INTENT(IN) :: lssav, ldiag3d, lsidea, qdiag3d
  LOGICAL, INTENT(IN) :: cplflx

! NAMELIST OPTIONS (INPUT):
      LOGICAL, INTENT(IN) :: bl_mynn_tkeadvect, ltaerosol,  &
                             lprnt, do_mynnsfclay,          &
                             flag_for_pbl_generic_tend
      INTEGER, INTENT(IN) ::                                &
     &       bl_mynn_cloudpdf,                              &
     &       bl_mynn_mixlength,                             &
     &       icloud_bl,                                     &
     &       bl_mynn_edmf,                                  &
     &       bl_mynn_edmf_mom,                              &
     &       bl_mynn_edmf_tke,                              &
     &       bl_mynn_edmf_part,                             &
     &       bl_mynn_cloudmix,                              &
     &       bl_mynn_mixqt,                                 &
     &       bl_mynn_tkebudget,                             &
     &       bl_mynn_output,                                &
     &       grav_settling,                                 &
     &       imp_physics, imp_physics_wsm6,                 &
     &       imp_physics_thompson, imp_physics_gfdl

!MISC CONFIGURATION OPTIONS
      INTEGER, PARAMETER ::                                 &
     &       spp_pbl=0,                                     &
     &       bl_mynn_mixscalars=1,                          &
     &       levflag=2
      LOGICAL ::                                            &
     &       FLAG_QI, FLAG_QNI, FLAG_QC, FLAG_QNC,          &
     &       FLAG_QNWFA, FLAG_QNIFA
      ! Define locally until needed from CCPP
      LOGICAL, PARAMETER :: cycling = .false.
      INTEGER, PARAMETER :: param_first_scalar = 1
      INTEGER ::                                            &
       &      p_qc, p_qr, p_qi, p_qs, p_qg, p_qnc, p_qni

!MYNN-1D
      REAL(kind=kind_phys), intent(in) :: delt, dtf
      INTEGER, intent(in) :: im, levs
      LOGICAL, intent(in) :: flag_init, flag_restart
      INTEGER :: initflag, k, i
      INTEGER :: IDS,IDE,JDS,JDE,KDS,KDE,                                &
     &            IMS,IME,JMS,JME,KMS,KME,                               &
     &            ITS,ITE,JTS,JTE,KTS,KTE
      INTEGER :: kdvel, num_vert_mix
      INTEGER, PARAMETER :: nchem=1, ndvel=1
      REAL(kind=kind_phys) :: tem

!MYNN-3D
      real(kind=kind_phys), dimension(im,levs+1), intent(in) :: phii
      real(kind=kind_phys), dimension(im,levs  ), intent(inout) ::       &
     &        dtdt, dudt, dvdt,                                          &
     &        dqdt_water_vapor, dqdt_liquid_cloud, dqdt_ice_cloud,       &
     &        dqdt_cloud_droplet_num_conc, dqdt_ice_num_conc,            &
     &        dqdt_ozone, dqdt_water_aer_num_conc, dqdt_ice_aer_num_conc
      real(kind=kind_phys), dimension(im,levs), intent(inout) ::         &
     &        qke, qke_adv, EL_PBL, Sh3D,                                &
     &        qc_bl, qi_bl, cldfra_bl
!These 10 arrays are only allocated when bl_mynn_output > 0
      real(kind=kind_phys), dimension(:,:), intent(inout) ::             &
     &        edmf_a,edmf_w,edmf_qt,                                     &
     &        edmf_thl,edmf_ent,edmf_qc,                                 &
     &        sub_thl,sub_sqv,det_thl,det_sqv
     real(kind=kind_phys), dimension(im,levs), intent(in) ::             &
    &        u,v,omega,t3d,                                              &
    &        exner,prsl,                                                 &
    &        qgrs_water_vapor,                                           &
    &        qgrs_liquid_cloud,                                          &
    &        qgrs_ice_cloud,                                             &
    &        qgrs_cloud_droplet_num_conc,                                &
    &        qgrs_cloud_ice_num_conc,                                    &
    &        qgrs_ozone,                                                 &
    &        qgrs_water_aer_num_conc,                                    &
    &        qgrs_ice_aer_num_conc
     real(kind=kind_phys), dimension(im,levs), intent(out) ::            &
    &        Tsq, Qsq, Cov, exch_h, exch_m
     real(kind=kind_phys), dimension(:,:), intent(inout) ::              &
    &        du3dt_PBL, du3dt_OGWD, dv3dt_PBL, dv3dt_OGWD,               &
    &        do3dt_PBL, dq3dt_PBL, dt3dt_PBL
    real(kind=kind_phys), dimension(im), intent(in) :: xmu
    real(kind=kind_phys), dimension(im, levs), intent(in) :: htrsw, htrlw
     !LOCAL
      real(kind=kind_phys), dimension(im,levs) ::                        &
     &        sqv,sqc,sqi,qnc,qni,ozone,qnwfa,qnifa,                     &
     &        dz, w, p, rho, th, qv,                                     &
     &        RUBLTEN, RVBLTEN, RTHBLTEN, RQVBLTEN,                      &
     &        RQCBLTEN, RQNCBLTEN, RQIBLTEN, RQNIBLTEN,                  &
     &        RQNWFABLTEN, RQNIFABLTEN,                                  &
     &        dqke,qWT,qSHEAR,qBUOY,qDISS,                               &
     &        pattern_spp_pbl

!MYNN-CHEM arrays
      real(kind=kind_phys), dimension(im,nchem) :: chem3d
      real(kind=kind_phys), dimension(im,ndvel) :: vd3d
      REAL(kind=kind_phys), DIMENSION( levs, nchem ) :: chem1
      REAL(kind=kind_phys), DIMENSION( levs+1, nchem ) :: s_awchem1
      REAL(kind=kind_phys), DIMENSION( ndvel ) :: vd1

!MYNN-2D
      real(kind=kind_phys), dimension(im), intent(in) ::                 &
     &        dx,zorl,slmsk,tsurf,qsfc,ps,                               &
     &        hflx,qflx,ust,wspd,rb,recmol

      real(kind=kind_phys), dimension(im), intent(in) ::                 &
     &        dusfc_cice,dvsfc_cice,dtsfc_cice,dqsfc_cice,               &
     &        stress_ocn,hflx_ocn,qflx_ocn,                              &
     &        oceanfrac,fice

      logical, dimension(im), intent(in) ::                              &
     &        wet, dry, icy

      real(kind=kind_phys), dimension(im), intent(inout) ::              &
     &        pblh
      real(kind=kind_phys), dimension(im), intent(out) ::                &
     &        ch,dtsfc1,dqsfc1,dusfc1,dvsfc1,                            &
     &        dtsfci_diag,dqsfci_diag,dtsfc_diag,dqsfc_diag,             &
     &        dusfci_diag,dvsfci_diag,dusfc_diag,dvsfc_diag,             &
     &        maxMF
      integer, dimension(im), intent(inout) ::                           &
     &        kpbl,nupdraft,ktop_plume

      real(kind=kind_phys), dimension(:), intent(inout) ::              &
     &        dusfc_cpl,dvsfc_cpl,dtsfc_cpl,dqsfc_cpl
      real(kind=kind_phys), dimension(:), intent(inout) ::              &
     &        dusfci_cpl,dvsfci_cpl,dtsfci_cpl,dqsfci_cpl

     !LOCAL
      real, dimension(im) ::                                             &
     &        WSTAR,DELTA,qcg,hfx,qfx,rmol,xland,                        &
     &        uoce,voce,vdfg,znt,ts

      real, dimension(im) :: dusfci1,dvsfci1,dtsfci1,dqsfci1

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lprnt) then
         write(0,*)"=============================================="
         write(0,*)"in mynn wrapper..."
         write(0,*)"flag_init=",flag_init
         write(0,*)"flag_restart=",flag_restart
      endif

      ! DH* TODO: Use flag_restart to distinguish which fields need
      ! to be initialized and which are read from restart files
      if (flag_init) then
         initflag=1
         !print*,"in MYNN, initflag=",initflag
      else
         initflag=0
         !print*,"in MYNN, initflag=",initflag
      endif

  ! Assign variables for each microphysics scheme
        if (imp_physics == imp_physics_wsm6) then
  ! WSM6
         FLAG_QI = .true.
         FLAG_QNI= .false.
         FLAG_QC = .true.
         FLAG_QNC= .false.
         FLAG_QNWFA= .false.
         FLAG_QNIFA= .false.
         p_qc = 2
         p_qr = 0
         p_qi = 2 
         p_qs = 0 
         p_qg = 0
         p_qnc= 0 
         p_qni= 0 
         do k=1,levs
            do i=1,im
              sqv(i,k)   = qgrs_water_vapor(i,k)
              sqc(i,k)   = qgrs_liquid_cloud(i,k)
              sqi(i,k)   = qgrs_ice_cloud(i,k)
              ozone(i,k) = qgrs_ozone(i,k)
              qnc(i,k)   = 0.
              qni(i,k)   = 0.
              qnwfa(i,k) = 0.
              qnifa(i,k) = 0.
            enddo
          enddo
        elseif (imp_physics == imp_physics_thompson) then
  ! Thompson
          if(ltaerosol) then
            FLAG_QI = .true.
            FLAG_QNI= .true.
            FLAG_QC = .true.
            FLAG_QNC= .true.
            FLAG_QNWFA= .true.
            FLAG_QNIFA= .true.
            p_qc = 2
            p_qr = 0
            p_qi = 2
            p_qs = 0
            p_qg = 0
            p_qnc= 0
            p_qni= 0
            do k=1,levs
              do i=1,im
                sqv(i,k)   = qgrs_water_vapor(i,k)
                sqc(i,k)   = qgrs_liquid_cloud(i,k)
                sqi(i,k)   = qgrs_ice_cloud(i,k)
                qnc(i,k)   = qgrs_cloud_droplet_num_conc(i,k)
                qni(i,k)   = qgrs_cloud_ice_num_conc(i,k)
                ozone(i,k) = qgrs_ozone(i,k)
                qnwfa(i,k) = qgrs_water_aer_num_conc(i,k)
                qnifa(i,k) = qgrs_ice_aer_num_conc(i,k)
              enddo
            enddo
          else
            FLAG_QI = .true.
            FLAG_QNI= .true.
            FLAG_QC = .true.
            FLAG_QNC= .false.
            FLAG_QNWFA= .false.
            FLAG_QNIFA= .false.
            p_qc = 2
            p_qr = 0
            p_qi = 2
            p_qs = 0
            p_qg = 0
            p_qnc= 0
            p_qni= 0
            do k=1,levs
              do i=1,im
                sqv(i,k)   = qgrs_water_vapor(i,k)
                sqc(i,k)   = qgrs_liquid_cloud(i,k)
                sqi(i,k)   = qgrs_ice_cloud(i,k)
                qnc(i,k)   = 0.
                qni(i,k)   = qgrs_cloud_ice_num_conc(i,k)
                ozone(i,k) = qgrs_ozone(i,k)
                qnwfa(i,k) = 0.
                qnifa(i,k) = 0.
              enddo
            enddo
          endif
        elseif (imp_physics == imp_physics_gfdl) then
  ! GFDL MP
          FLAG_QI = .true.
          FLAG_QNI= .false.
          FLAG_QC = .true.
          FLAG_QNC= .false.
          FLAG_QNWFA= .false.
          FLAG_QNIFA= .false.
          p_qc = 2
          p_qr = 0
          p_qi = 2
          p_qs = 0
          p_qg = 0
          p_qnc= 0
          p_qni= 0
          do k=1,levs
            do i=1,im
                sqv(i,k)   = qgrs_water_vapor(i,k)
                sqc(i,k)   = qgrs_liquid_cloud(i,k)
                sqi(i,k)   = qgrs_ice_cloud(i,k)
                qnc(i,k)   = 0.
                qni(i,k)   = 0.
                qnwfa(i,k) = 0.
                qnifa(i,k) = 0.
                ozone(i,k) = qgrs_ozone(i,k)
            enddo
          enddo
        else
          print*,"In MYNN wrapper. Unknown microphysics scheme, imp_physics=",imp_physics
          print*,"Defaulting to qc and qv species only..."
          FLAG_QI = .false.
          FLAG_QNI= .false.
          FLAG_QC = .true.
          FLAG_QNC= .false.
          FLAG_QNWFA= .false.
          FLAG_QNIFA= .false.
          p_qc = 2
          p_qr = 0
          p_qi = 0
          p_qs = 0
          p_qg = 0
          p_qnc= 0
          p_qni= 0
          do k=1,levs
            do i=1,im
                sqv(i,k)   = qgrs_water_vapor(i,k)
                sqc(i,k)   = qgrs_liquid_cloud(i,k)
                sqi(i,k)   = 0.
                qnc(i,k)   = 0.
                qni(i,k)   = 0.
                qnwfa(i,k) = 0.
                qnifa(i,k) = 0.
                ozone(i,k) = qgrs_ozone(i,k)
            enddo
          enddo
        endif

       if (lprnt)write(0,*)"prepping MYNN-EDMF variables..."

       do k=1,levs
          do i=1,im
             dz(i,k)=(phii(i,k+1) - phii(i,k))*g_inv
             th(i,k)=t3d(i,k)/exner(i,k)
          ! keep as specific humidity
          !   qv(i,k)=qvsh(i,k)/(1.0 - qvsh(i,k))
          !   qc(i,k)=qc(i,k)/(1.0 - qvsh(i,k))
          !   qi(i,k)=qi(i,k)/(1.0 - qvsh(i,k))
             rho(i,k)=prsl(i,k)/(r_d*t3d(i,k))
             w(i,k) = -omega(i,k)/(rho(i,k)*g)
             pattern_spp_pbl(i,k)=0.0
         enddo
      enddo
      do i=1,im
         if (slmsk(i)==1. .or. slmsk(i)==2.) then !sea/land/ice mask (=0/1/2) in FV3
            xland(i)=1.0                          !but land/water = (1/2) in SFCLAY_mynn
         else
            xland(i)=2.0
         endif
         uoce(i)=0.0
         voce(i)=0.0
         vdfg(i)=0.0
         !ust(i) = sqrt(stress(i))
         ch(i)=0.0
         hfx(i)=hflx(i)*rho(i,1)*cp
         qfx(i)=qflx(i)*rho(i,1)
         wstar(i)=0.0
         delta(i)=0.0
         qcg(i)=0.0

         dtsfc1(i) = hfx(i)
         dqsfc1(i) = qfx(i)*XLV
         dusfc1(i) = -1.*rho(i,1)*ust(i)*ust(i)*u(i,1)/wspd(i)
         dvsfc1(i) = -1.*rho(i,1)*ust(i)*ust(i)*v(i,1)/wspd(i)

         !BWG: diagnostic surface fluxes for scalars & momentum
         dtsfci_diag(i) = dtsfc1(i)
         dqsfci_diag(i) = dqsfc1(i)
         dtsfc_diag(i)  = dtsfc_diag(i) + dtsfc1(i)*delt
         dqsfc_diag(i)  = dqsfc_diag(i) + dqsfc1(i)*delt
         dusfci_diag(i) = dusfc1(i)
         dvsfci_diag(i) = dvsfc1(i)
         dusfc_diag(i)  = dusfc_diag(i) + dusfci_diag(i)*delt
         dvsfc_diag(i)  = dvsfc_diag(i) + dvsfci_diag(i)*delt

         znt(i)=zorl(i)*0.01 !cm -> m?
         if (do_mynnsfclay) then
           rmol(i)=recmol(i)
         else
           if (hfx(i) .ge. 0.)then
             rmol(i)=-hfx(i)/(200.*dz(i,1)*0.5)
           else
             rmol(i)=ABS(rb(i))*1./(dz(i,1)*0.5)
           endif
           !if (rb(i) .ge. 0.)then
           !  rmol(i)=rb(i)*8./(dz(i,1)*0.5)
           !else
           !  rmol(i)=MAX(rb(i)*5.,-10.)/(dz(i,1)*0.5)
           !endif
         endif
         ts(i)=tsurf(i)/exner(i,1)  !theta
!        qsfc(i)=qss(i)
!        ps(i)=pgr(i)
!        wspd(i)=wind(i)
      enddo

      ! BWG: Coupling insertion
      if (cplflx) then
        do i=1,im
          if (oceanfrac(i) > zero) then ! Ocean only, NO LAKES
            if ( .not. wet(i)) then ! no open water, use results from CICE
              dusfci_cpl(i) = dusfc_cice(i)
              dvsfci_cpl(i) = dvsfc_cice(i)
              dtsfci_cpl(i) = dtsfc_cice(i)
              dqsfci_cpl(i) = dqsfc_cice(i)
            elseif (icy(i) .or. dry(i)) then ! use stress_ocean for opw component at mixed point
              if (wspd(i) > zero) then
                dusfci_cpl(i) = -1.*rho(i,1)*stress_ocn(i)*u(i,1)/wspd(i)   ! U-momentum flux
                dvsfci_cpl(i) = -1.*rho(i,1)*stress_ocn(i)*v(i,1)/wspd(i)   ! V-momentum flux
              else
                dusfci_cpl(i) = zero
                dvsfci_cpl(i) = zero
              endif
              dtsfci_cpl(i) =  cp*rho(i,1)*hflx_ocn(i) ! sensible heat flux over open ocean
              dqsfci_cpl(i) = XLV*rho(i,1)*qflx_ocn(i) ! latent heat flux over open ocean
            else                                       ! use results from this scheme for 100% open ocean
              dusfci_cpl(i) = dusfci_diag(i)
              dvsfci_cpl(i) = dvsfci_diag(i)
              dtsfci_cpl(i) = dtsfci_diag(i)
              dqsfci_cpl(i) = dqsfci_diag(i)
            endif
!
            dusfc_cpl (i) = dusfc_cpl(i) + dusfci_cpl(i) * delt
            dvsfc_cpl (i) = dvsfc_cpl(i) + dvsfci_cpl(i) * delt
            dtsfc_cpl (i) = dtsfc_cpl(i) + dtsfci_cpl(i) * delt
            dqsfc_cpl (i) = dqsfc_cpl(i) + dqsfci_cpl(i) * delt
          else ! If no ocean
            dusfc_cpl(i) = huge
            dvsfc_cpl(i) = huge
            dtsfc_cpl(i) = huge
            dqsfc_cpl(i) = huge
          endif ! Ocean only, NO LAKES
        enddo
      endif ! End coupling insertion

      if (lprnt) then
         print*
         write(0,*)"===CALLING mynn_bl_driver; input:"
         print*,"bl_mynn_tkebudget=",bl_mynn_tkebudget," bl_mynn_tkeadvect=",bl_mynn_tkeadvect
         print*,"bl_mynn_cloudpdf=",bl_mynn_cloudpdf," bl_mynn_mixlength=",bl_mynn_mixlength
         print*,"bl_mynn_edmf=",bl_mynn_edmf," bl_mynn_edmf_mom=",bl_mynn_edmf_mom
         print*,"bl_mynn_edmf_tke=",bl_mynn_edmf_tke," bl_mynn_edmf_part=",bl_mynn_edmf_part
         print*,"bl_mynn_cloudmix=",bl_mynn_cloudmix," bl_mynn_mixqt=",bl_mynn_mixqt
         print*,"icloud_bl=",icloud_bl
         print*,"T:",t3d(1,1),t3d(1,2),t3d(1,levs)
         print*,"TH:",th(1,1),th(1,2),th(1,levs)
         print*,"rho:",rho(1,1),rho(1,2),rho(1,levs)
         print*,"exner:",exner(1,1),exner(1,2),exner(1,levs)
         print*,"prsl:",prsl(1,1),prsl(1,2),prsl(1,levs)
         print*,"dz:",dz(1,1),dz(1,2),dz(1,levs)
         print*,"u:",u(1,1),u(1,2),u(1,levs)
         print*,"v:",v(1,1),v(1,2),v(1,levs)
         print*,"sqv:",sqv(1,1),sqv(1,2),sqv(1,levs)
         print*,"sqc:",sqc(1,1),sqc(1,2),sqc(1,levs)
         print*,"sqi:",sqi(1,1),sqi(1,2),sqi(1,levs)
         print*,"rmol:",rmol(1)," ust:",ust(1)
         print*," dx=",dx(1),"initflag=",initflag
         print*,"Tsurf:",tsurf(1)," Thetasurf:",ts(1)
         print*,"HFX:",hfx(1)," qfx",qfx(1)
         print*,"qsfc:",qsfc(1)," ps:",ps(1)
         print*,"wspd:",wspd(1)," rb=",rb(1)
         print*,"znt:",znt(1)," delt=",delt
         print*,"im=",im," levs=",levs
         print*,"PBLH=",pblh(1)," KPBL=",KPBL(1)," xland=",xland(1)
         print*,"vdfg=",vdfg(1)," ch=",ch(1)
         !print*,"TKE:",TKE_PBL(1,1),TKE_PBL(1,2),TKE_PBL(1,levs)
         print*,"qke:",qke(1,1),qke(1,2),qke(1,levs)
         print*,"el_pbl:",el_pbl(1,1),el_pbl(1,2),el_pbl(1,levs)
         print*,"Sh3d:",Sh3d(1,1),sh3d(1,2),sh3d(1,levs)
         !print*,"exch_h:",exch_h(1,1),exch_h(1,2),exch_h(1,levs) ! - intent(out)
         !print*,"exch_m:",exch_m(1,1),exch_m(1,2),exch_m(1,levs) ! - intent(out)
         print*,"max cf_bl:",maxval(cldfra_bl(1,:))
      endif


              CALL  mynn_bl_driver(                                    &
     &             initflag=initflag,restart=flag_restart,             &
     &             cycling=cycling,                                    &
     &             grav_settling=grav_settling,                        &
     &             delt=delt,dz=dz,dx=dx,znt=znt,                      &
     &             u=u,v=v,w=w,th=th,sqv3D=sqv,sqc3D=sqc,              &
     &             sqi3D=sqi,qni=qni,qnc=qnc,                          &
     &             qnwfa=qnwfa,qnifa=qnifa,ozone=ozone,                &
     &             p=prsl,exner=exner,rho=rho,T3D=t3d,                 &
     &             xland=xland,ts=ts,qsfc=qsfc,qcg=qcg,ps=ps,          &
     &             ust=ust,ch=ch,hfx=hfx,qfx=qfx,rmol=rmol,            &
     &             wspd=wspd,uoce=uoce,voce=voce,vdfg=vdfg,            & !input
     &             qke=QKE,sh3d=Sh3d,                                  & !output
     &             qke_adv=qke_adv,bl_mynn_tkeadvect=bl_mynn_tkeadvect,&
#if (WRF_CHEM == 1)
     &             chem3d=chem,vd3d=vd,nchem=nchem,kdvel=kdvel,        &
     &             ndvel=ndvel,num_vert_mix=num_vert_mix,              &
#endif
     &             Tsq=tsq,Qsq=qsq,Cov=cov,                            & !output
     &             RUBLTEN=RUBLTEN,RVBLTEN=RVBLTEN,RTHBLTEN=RTHBLTEN,  & !output
     &             RQVBLTEN=RQVBLTEN,RQCBLTEN=rqcblten,                &
     &             RQIBLTEN=rqiblten,RQNCBLTEN=rqncblten,              & !output
     &             RQNIBLTEN=rqniblten,RQNWFABLTEN=RQNWFABLTEN,        & !output
     &             RQNIFABLTEN=RQNIFABLTEN,dozone=dqdt_ozone,          & !output
     &             EXCH_H=exch_h,EXCH_M=exch_m,                        & !output
     &             pblh=pblh,KPBL=KPBL                                 & !output
     &             ,el_pbl=el_pbl                                      & !output
     &             ,dqke=dqke                                          & !output
     &             ,qWT=qWT,qSHEAR=qSHEAR,qBUOY=qBUOY,qDISS=qDISS      & !output
     &             ,WSTAR=wstar,DELTA=delta                            & !unused input
     &             ,bl_mynn_tkebudget=bl_mynn_tkebudget                & !input parameter
     &             ,bl_mynn_cloudpdf=bl_mynn_cloudpdf                  & !input parameter
     &             ,bl_mynn_mixlength=bl_mynn_mixlength                & !input parameter
     &             ,icloud_bl=icloud_bl                                & !input parameter
     &             ,qc_bl=qc_bl,qi_bl=qi_bl,cldfra_bl=cldfra_bl        & !output
     &             ,levflag=levflag,bl_mynn_edmf=bl_mynn_edmf          & !input parameter
     &             ,bl_mynn_edmf_mom=bl_mynn_edmf_mom                  & !input parameter
     &             ,bl_mynn_edmf_tke=bl_mynn_edmf_tke                  & !input parameter
     &             ,bl_mynn_mixscalars=bl_mynn_mixscalars              & !input parameter
     &             ,bl_mynn_output=bl_mynn_output                      & !input parameter
     &             ,bl_mynn_cloudmix=bl_mynn_cloudmix                  & !input parameter
     &             ,bl_mynn_mixqt=bl_mynn_mixqt                        & !input parameter
     &             ,edmf_a=edmf_a,edmf_w=edmf_w,edmf_qt=edmf_qt        & !output
     &             ,edmf_thl=edmf_thl,edmf_ent=edmf_ent,edmf_qc=edmf_qc &!output
     &             ,sub_thl3D=sub_thl,sub_sqv3D=sub_sqv                &
     &             ,det_thl3D=det_thl,det_sqv3D=det_sqv                &
     &             ,nupdraft=nupdraft,maxMF=maxMF                      & !output
     &             ,ktop_plume=ktop_plume                              & !output
     &             ,spp_pbl=spp_pbl,pattern_spp_pbl=pattern_spp_pbl    & !input
     &             ,RTHRATEN=htrlw                                     & !input
     &             ,FLAG_QI=flag_qi,FLAG_QNI=flag_qni                  & !input
     &             ,FLAG_QC=flag_qc,FLAG_QNC=flag_qnc                  & !input
     &             ,FLAG_QNWFA=FLAG_QNWFA,FLAG_QNIFA=FLAG_QNIFA        & !input
     &             ,IDS=1,IDE=im,JDS=1,JDE=1,KDS=1,KDE=levs            & !input
     &             ,IMS=1,IME=im,JMS=1,JME=1,KMS=1,KME=levs            & !input
     &             ,ITS=1,ITE=im,JTS=1,JTE=1,KTS=1,KTE=levs)             !input


     ! POST MYNN (INTERSTITIAL) WORK:
        !update/save MYNN-only variables
        !do k=1,levs
        !   do i=1,im
        !      gq0(i,k,4)=qke(i,k,1)      !tke*2
        !   enddo
        !enddo
        !For MYNN, convert TH-tend to T-tend
        do k = 1, levs
           do i = 1, im
              dtdt(i,k) = dtdt(i,k) + RTHBLTEN(i,k)*exner(i,k)
              dudt(i,k) = dudt(i,k) + RUBLTEN(i,k)
              dvdt(i,k) = dvdt(i,k) + RVBLTEN(i,k)
           enddo
        enddo
        accum_duvt3dt: if(lssav) then
          if(ldiag3d .and. .not. flag_for_pbl_generic_tend) then
            do k = 1, levs
              do i = 1, im
                du3dt_PBL(i,k) = du3dt_PBL(i,k) + RUBLTEN(i,k)*dtf
                dv3dt_PBL(i,k) = dv3dt_PBL(i,k) + RVBLTEN(i,k)*dtf
              enddo
            enddo
          endif
          
          if (lsidea .or. (ldiag3d .and. .not. flag_for_pbl_generic_tend)) then
            do k = 1, levs
               do i = 1, im
                 dt3dt_PBL(i,k) = dt3dt_PBL(i,k) + RTHBLTEN(i,k)*exner(i,k)*dtf
               enddo
            enddo   
          endif
        endif accum_duvt3dt
        !Update T, U and V:
        !do k = 1, levs
        !   do i = 1, im
        !      T3D(i,k) = T3D(i,k) + RTHBLTEN(i,k)*exner(i,k)*delt
        !      u(i,k)   = u(i,k) + RUBLTEN(i,k)*delt
        !      v(i,k)   = v(i,k) + RVBLTEN(i,k)*delt
        !   enddo
        !enddo

        !DO moist/scalar/tracer tendencies:
        if (imp_physics == imp_physics_wsm6) then
           ! WSM6
           do k=1,levs
             do i=1,im
               dqdt_water_vapor(i,k)  = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_liquid_cloud(i,k) = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_ice_cloud(i,k)    = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
               !dqdt_ozone(i,k)        = 0.0
             enddo
           enddo
           !Update moist species:
           !do k=1,levs
           !  do i=1,im
           !    qgrs_water_vapor(i,k)  = qgrs_water_vapor(i,k)  + (RQVBLTEN(i,k)/(1.0+RQVBLTEN(i,k)))*delt
           !    qgrs_liquid_cloud(i,k) = qgrs_liquid_cloud(i,k) + RQCBLTEN(i,k)*delt
           !    qgrs_ice_cloud(i,k)    = qgrs_ice_cloud(i,k)    + RQIBLTEN(i,k)*delt
           !    !dqdt_ozone(i,k)        = 0.0
           !  enddo
           !enddo
        elseif (imp_physics == imp_physics_thompson) then
           ! Thompson-Aerosol
           if(ltaerosol) then
             do k=1,levs
               do i=1,im
                 dqdt_water_vapor(i,k)             = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_liquid_cloud(i,k)            = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_cloud_droplet_num_conc(i,k)  = RQNCBLTEN(i,k)
                 dqdt_ice_cloud(i,k)               = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_ice_num_conc(i,k)            = RQNIBLTEN(i,k)
                 !dqdt_ozone(i,k)                   = 0.0
                 dqdt_water_aer_num_conc(i,k)      = RQNWFABLTEN(i,k)
                 dqdt_ice_aer_num_conc(i,k)        = RQNIFABLTEN(i,k)
               enddo
             enddo
             !do k=1,levs
             !  do i=1,im
             !    qgrs_water_vapor(i,k)            = qgrs_water_vapor(i,k)    + (RQVBLTEN(i,k)/(1.0+RQVBLTEN(i,k)))*delt
             !    qgrs_liquid_cloud(i,k)           = qgrs_liquid_cloud(i,k)   + RQCBLTEN(i,k)*delt
             !    qgrs_ice_cloud(i,k)              = qgrs_ice_cloud(i,k)      + RQIBLTEN(i,k)*delt
             !    qgrs_cloud_droplet_num_conc(i,k) = qgrs_cloud_droplet_num_conc(i,k) + RQNCBLTEN(i,k)*delt
             !    qgrs_cloud_ice_num_conc(i,k)     = qgrs_cloud_ice_num_conc(i,k)     + RQNIBLTEN(i,k)*delt
             !    !dqdt_ozone(i,k)        = 0.0
             !    !qgrs_water_aer_num_conc(i,k)     = qgrs_water_aer_num_conc(i,k)     + RQNWFABLTEN(i,k)*delt
             !    !qgrs_ice_aer_num_conc(i,k)       = qgrs_ice_aer_num_conc(i,k)       + RQNIFABLTEN(i,k)*delt
             !  enddo
             !enddo
           else
             !Thompson (2008)
             do k=1,levs
               do i=1,im
                 dqdt_water_vapor(i,k)   = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_liquid_cloud(i,k)  = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_ice_cloud(i,k)     = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_ice_num_conc(i,k)  = RQNIBLTEN(i,k)
                 !dqdt_ozone(i,k)         = 0.0
               enddo
             enddo
             !do k=1,levs
             !  do i=1,im
             !    qgrs_water_vapor(i,k)            = qgrs_water_vapor(i,k)    + (RQVBLTEN(i,k)/(1.0+RQVBLTEN(i,k)))*delt
             !    qgrs_liquid_cloud(i,k)           = qgrs_liquid_cloud(i,k)   + RQCBLTEN(i,k)*delt
             !    qgrs_ice_cloud(i,k)              = qgrs_ice_cloud(i,k)      + RQIBLTEN(i,k)*delt
             !    qgrs_cloud_ice_num_conc(i,k)     = qgrs_cloud_ice_num_conc(i,k)     + RQNIBLTEN(i,k)*delt
             !    !dqdt_ozone(i,k)        = 0.0
             !  enddo
             !enddo
           endif !end thompson choice
        elseif (imp_physics == imp_physics_gfdl) then
           ! GFDL MP
           do k=1,levs
             do i=1,im
               dqdt_water_vapor(i,k)   = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_liquid_cloud(i,k)  = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_ice_cloud(i,k)     = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
               !dqdt_rain(i,k)          = 0.0
               !dqdt_snow(i,k)          = 0.0
               !dqdt_graupel(i,k)       = 0.0
               !dqdt_ozone(i,k)         = 0.0
             enddo
           enddo
           !do k=1,levs
           !  do i=1,im
           !    qgrs_water_vapor(i,k)            = qgrs_water_vapor(i,k)    + (RQVBLTEN(i,k)/(1.0+RQVBLTEN(i,k)))*delt
           !    qgrs_liquid_cloud(i,k)           = qgrs_liquid_cloud(i,k)   + RQCBLTEN(i,k)*delt
           !    qgrs_ice_cloud(i,k)              = qgrs_ice_cloud(i,k)      + RQIBLTEN(i,k)*delt
           !    !dqdt_ozone(i,k)        = 0.0
           !  enddo
           !enddo
       else
!          print*,"In MYNN wrapper. Unknown microphysics scheme, imp_physics=",imp_physics
           do k=1,levs
             do i=1,im
               dqdt_water_vapor(i,k)   = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_liquid_cloud(i,k)  = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_ice_cloud(i,k)     = 0.0
               !dqdt_rain(i,k)          = 0.0
               !dqdt_snow(i,k)          = 0.0
               !dqdt_graupel(i,k)       = 0.0
               !dqdt_ozone(i,k)         = 0.0
             enddo
           enddo
       endif
       
       if(lssav .and. (ldiag3d .and. qdiag3d .and. .not. flag_for_pbl_generic_tend)) then
         do k=1,levs
           do i=1,im
             dq3dt_PBL(i,k)  = dq3dt_PBL(i,k) + dqdt_water_vapor(i,k)*dtf
           enddo
         enddo
       endif

       if (lprnt) then
          print*
          print*,"===Finished with mynn_bl_driver; output:"
          print*,"T:",t3d(1,1),t3d(1,2),t3d(1,levs)
          print*,"TH:",th(1,1),th(1,2),th(1,levs)
          print*,"rho:",rho(1,1),rho(1,2),rho(1,levs)
          print*,"exner:",exner(1,1),exner(1,2),exner(1,levs)
          print*,"prsl:",prsl(1,1),prsl(1,2),prsl(1,levs)
          print*,"dz:",dz(1,1),dz(1,2),dz(1,levs)
          print*,"u:",u(1,1),u(1,2),u(1,levs)
          print*,"v:",v(1,1),v(1,2),v(1,levs)
          print*,"sqv:",sqv(1,1),sqv(1,2),sqv(1,levs)
          print*,"sqc:",sqc(1,1),sqc(1,2),sqc(1,levs)
          print*,"sqi:",sqi(1,1),sqi(1,2),sqi(1,levs)
          print*,"rmol:",rmol(1)," ust:",ust(1)
          print*,"dx(1)=",dx(1),"initflag=",initflag
          print*,"Tsurf:",tsurf(1)," Thetasurf:",ts(1)
          print*,"HFX:",hfx(1)," qfx",qfx(1)
          print*,"qsfc:",qsfc(1)," ps:",ps(1)
          print*,"wspd:",wspd(1)," rb=",rb(1)
          print*,"znt:",znt(1)," delt=",delt
          print*,"im=",im," levs=",levs
          print*,"PBLH=",pblh(1)," KPBL=",KPBL(1)," xland=",xland(1)
          print*,"vdfg=",vdfg(1)," ch=",ch(1)
          !print*,"TKE:",TKE_PBL(1,1),TKE_PBL(1,2),TKE_PBL(1,levs)
          print*,"qke:",qke(1,1),qke(1,2),qke(1,levs)
          print*,"el_pbl:",el_pbl(1,1),el_pbl(1,2),el_pbl(1,levs)
          print*,"Sh3d:",Sh3d(1,1),sh3d(1,2),sh3d(1,levs)
          print*,"exch_h:",exch_h(1,1),exch_h(1,2),exch_h(1,levs)
          print*,"exch_m:",exch_m(1,1),exch_m(1,2),exch_m(1,levs)
          print*,"max cf_bl:",maxval(cldfra_bl(1,:))
          print*,"max qc_bl:",maxval(qc_bl(1,:))
          print*,"dtdt:",dtdt(1,1),dtdt(1,2),dtdt(1,levs)
          print*,"dudt:",dudt(1,1),dudt(1,2),dudt(1,levs)
          print*,"dvdt:",dvdt(1,1),dvdt(1,2),dvdt(1,levs)
          print*,"dqdt:",dqdt_water_vapor(1,1),dqdt_water_vapor(1,2),dqdt_water_vapor(1,levs)
          print*,"ktop_plume:",ktop_plume(1)," maxmf:",maxmf(1)
          print*,"nup:",nupdraft(1)
          print*
       endif


  END SUBROUTINE mynnedmf_wrapper_run

!###=================================================================

END MODULE mynnedmf_wrapper
