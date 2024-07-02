!> \file mynnedmf_wrapper.F90
!!  This file contains all of the code related to running the MYNN 
!! eddy-diffusivity mass-flux scheme. 

!> The following references best describe the code within
!!    Olson et al. (2019, NOAA Technical Memorandum)
!!    Nakanishi and Niino (2009) \cite NAKANISHI_2009
      MODULE mynnedmf_wrapper

      contains

!> \section arg_table_mynnedmf_wrapper_init Argument Table
!! \htmlinclude mynnedmf_wrapper_init.html
!!
      subroutine mynnedmf_wrapper_init (          &
        &  con_cp, con_grav, con_rd, con_rv,      &
        &  con_cpv, con_cliq, con_cice, con_rcp,  &
        &  con_XLV, con_XLF, con_p608, con_ep2,   &
        &  con_karman, con_t0c,                   &
        &  do_mynnedmf,                           &
        &  errmsg, errflg                         )

        use machine,  only : kind_phys
        use bl_mynn_common

        implicit none

        logical,        intent(in)  :: do_mynnedmf
        character(len=*),intent(out):: errmsg
        integer,        intent(out) :: errflg

        real(kind_phys),intent(in)  :: con_xlv
        real(kind_phys),intent(in)  :: con_xlf
        real(kind_phys),intent(in)  :: con_rv
        real(kind_phys),intent(in)  :: con_rd
        real(kind_phys),intent(in)  :: con_ep2
        real(kind_phys),intent(in)  :: con_grav
        real(kind_phys),intent(in)  :: con_cp
        real(kind_phys),intent(in)  :: con_cpv
        real(kind_phys),intent(in)  :: con_rcp
        real(kind_phys),intent(in)  :: con_p608
        real(kind_phys),intent(in)  :: con_cliq
        real(kind_phys),intent(in)  :: con_cice
        real(kind_phys),intent(in)  :: con_karman
        real(kind_phys),intent(in)  :: con_t0c

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        xlv    = con_xlv
        xlf    = con_xlf
        r_v    = con_rv
        r_d    = con_rd
        ep_2   = con_ep2
        grav   = con_grav
        cp     = con_cp
        cpv    = con_cpv
        rcp    = con_rcp
        p608   = con_p608
        cliq   = con_cliq
        cice   = con_cice
        karman = con_karman
        t0c    = con_t0c
       
        xls    = xlv+xlf      != 2.85E6 (J/kg) sublimation                                      
        rvovrd = r_v/r_d      != 1.608
        ep_3   = 1.-ep_2      != 0.378
        gtr    = grav/tref
        rk     = cp/r_d
        tv0    = p608*tref
        tv1    = (1.+p608)*tref
        xlscp  = (xlv+xlf)/cp
        xlvcp  = xlv/cp
        g_inv  = 1./grav

        ! Consistency checks
        if (.not. do_mynnedmf) then
          errmsg = 'Logic error: do_mynnedmf = .false.'
          errflg = 1
            return
         end if

      end subroutine mynnedmf_wrapper_init

!>\defgroup gp_mynnedmf MYNN-EDMF PBL and Shallow Convection Module  
!> This scheme (1) performs pre-mynnedmf work, (2) runs the mynnedmf, and (3) performs post-mynnedmf work
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
     &  qgrs_ice,                       &
     &  qgrs_snow,                      &
     &  qgrs_cloud_droplet_num_conc,    &
     &  qgrs_cloud_ice_num_conc,        &
     &  qgrs_ozone,                     &
     &  qgrs_water_aer_num_conc,        &
     &  qgrs_ice_aer_num_conc,          &
     &  qgrs_cccn,                      &
     &  prsl,prsi,exner,                &
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
     &  hflx_wat,qflx_wat,stress_wat,   &
     &  oceanfrac,fice,wet,icy,dry,     &
     &  dusfci_cpl,dvsfci_cpl,          &
     &  dtsfci_cpl,dqsfci_cpl,          &
     &  dusfc_cpl,dvsfc_cpl,            &
     &  dtsfc_cpl,dqsfc_cpl,            &
     &  recmol,                         &
     &  qke,qke_adv,Tsq,Qsq,Cov,        &
     &  el_pbl,sh3d,sm3d,exch_h,exch_m, &
     &  dqke,qwt,qshear,qbuoy,qdiss,    &
     &  Pblh,kpbl,                      &
     &  qc_bl,qi_bl,cldfra_bl,          &
     &  edmf_a,edmf_w,edmf_qt,          &
     &  edmf_thl,edmf_ent,edmf_qc,      &
     &  sub_thl,sub_sqv,det_thl,det_sqv,&
     &  maxwidth,maxMF,ztop_plume,      &
     &  ktop_plume,                     &
     &  dudt, dvdt, dtdt,                                  &
     &  dqdt_water_vapor,            dqdt_liquid_cloud,    & ! <=== ntqv, ntcw
     &  dqdt_ice,                    dqdt_snow,            & ! <=== ntiw, ntsw
     &  dqdt_ozone,                                        & ! <=== ntoz
     &  dqdt_cloud_droplet_num_conc, dqdt_ice_num_conc,    & ! <=== ntlnc, ntinc
     &  dqdt_water_aer_num_conc,     dqdt_ice_aer_num_conc,& ! <=== ntwa, ntia
     &  dqdt_cccn,                                         & ! <=== ntccn
     &  flag_for_pbl_generic_tend,                         &
     &  dtend, dtidx, index_of_temperature,                &
     &  index_of_x_wind, index_of_y_wind, ntke,            &
     &  ntqv, ntcw, ntiw, ntsw,                            &
     &  ntoz, ntlnc, ntinc, ntwa, ntia,                    &
     &  index_of_process_pbl, htrsw, htrlw, xmu,           &
     &  tke_budget,            bl_mynn_tkeadvect,          &
     &  bl_mynn_cloudpdf,      bl_mynn_mixlength,          &
     &  bl_mynn_edmf,                                      &
     &  bl_mynn_edmf_mom,      bl_mynn_edmf_tke,           &
     &  bl_mynn_cloudmix,      bl_mynn_mixqt,              &
     &  bl_mynn_output,        bl_mynn_closure,            &
     &  icloud_bl, do_mynnsfclay,                          &
     &  imp_physics, imp_physics_gfdl,                     &
     &  imp_physics_thompson, imp_physics_wsm6,            &
     &  imp_physics_fa,                                    &
     &  chem3d, frp, mix_chem, rrfs_sd, enh_mix,           &
     &  nchem, ndvel, vdep, smoke_dbg,                     &
     &  imp_physics_nssl, nssl_ccn_on,                     &
     &  ltaerosol, mraerosol, spp_wts_pbl, spp_pbl,        &
     &  lprnt, huge, errmsg, errflg                        )

! should be moved to inside the mynn:
     use machine,        only: kind_phys
     use bl_mynn_common, only: cp, r_d, grav, g_inv, zero, &
         xlv, xlvcp, xlscp, p608
     use module_bl_mynn, only: mynn_bl_driver

!------------------------------------------------------------------- 
     implicit none
!------------------------------------------------------------------- 

     real(kind_phys),  intent(in)  :: huge
     character(len=*), intent(out) :: errmsg
     integer, intent(out)          :: errflg

     logical, intent(in) :: lssav, ldiag3d, lsidea, qdiag3d
     logical, intent(in) :: cplflx

     !smoke/chem
     integer, intent(in) :: nchem, ndvel
     integer, parameter  :: kdvel=1
     logical, intent(in) :: smoke_dbg

! NAMELIST OPTIONS (INPUT):
     logical, intent(in) ::                                 &
     &       bl_mynn_tkeadvect,                             &
     &       ltaerosol,                                     &
     &       mraerosol,                                     &
     &       lprnt,                                         &
     &       do_mynnsfclay,                                 &
     &       flag_for_pbl_generic_tend,                     &
     &       nssl_ccn_on
      integer, intent(in) ::                                &
     &       bl_mynn_cloudpdf,                              &
     &       bl_mynn_mixlength,                             &
     &       icloud_bl,                                     &
     &       bl_mynn_edmf,                                  &
     &       bl_mynn_edmf_mom,                              &
     &       bl_mynn_edmf_tke,                              &
     &       bl_mynn_cloudmix,                              &
     &       bl_mynn_mixqt,                                 &
     &       bl_mynn_output,                                &
     &       imp_physics, imp_physics_wsm6,                 &
     &       imp_physics_thompson, imp_physics_gfdl,        &
     &       imp_physics_nssl, imp_physics_fa,              &
     &       spp_pbl,                                       &
     &       tke_budget
      real(kind_phys), intent(in) ::                        &
     &       bl_mynn_closure

!TENDENCY DIAGNOSTICS
      real(kind_phys), intent(inout), optional :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:)
      integer, intent(in) :: index_of_temperature, index_of_x_wind
      integer, intent(in) :: index_of_y_wind, index_of_process_pbl
      integer, intent(in) :: ntoz, ntqv, ntcw, ntiw, ntsw, ntlnc
      integer, intent(in) :: ntinc, ntwa, ntia, ntke

!MISC CONFIGURATION OPTIONS
      INTEGER, PARAMETER ::                                              &
     &       bl_mynn_mixscalars=1
      LOGICAL ::                                                         &
     &       FLAG_QI, FLAG_QNI, FLAG_QC, FLAG_QS, FLAG_QNC,              &
     &       FLAG_QNWFA, FLAG_QNIFA, FLAG_QNBCA, FLAG_OZONE
      ! Define locally until needed from CCPP
      LOGICAL, PARAMETER :: cycling = .false.

!MYNN-1D
      REAL(kind_phys), intent(in) :: delt, dtf
      INTEGER, intent(in) :: im, levs
      LOGICAL, intent(in) :: flag_init, flag_restart
      INTEGER :: initflag, k, i
      INTEGER :: IDS,IDE,JDS,JDE,KDS,KDE,                                &
     &           IMS,IME,JMS,JME,KMS,KME,                                &
     &           ITS,ITE,JTS,JTE,KTS,KTE

      REAL(kind_phys) :: tem

!MYNN-3D
      real(kind_phys), dimension(:,:), intent(in)    :: phii
      real(kind_phys), dimension(:,:), intent(inout) ::                  &
     &        dtdt, dudt, dvdt,                                          &
     &        dqdt_water_vapor, dqdt_liquid_cloud, dqdt_ice,             &
     &        dqdt_snow, dqdt_ice_num_conc, dqdt_ozone
      real(kind_phys), dimension(:,:), intent(inout), optional ::        &
     &        dqdt_cloud_droplet_num_conc, dqdt_water_aer_num_conc,      &
     &        dqdt_ice_aer_num_conc
      real(kind_phys), dimension(:,:), intent(inout), optional :: qke,   &
     &        EL_PBL, Sh3D, Sm3D, qc_bl, qi_bl, cldfra_bl, dqdt_cccn
      real(kind_phys), dimension(:,:), intent(inout) ::                  &
     &        qke_adv
     !These 10 arrays are only allocated when bl_mynn_output > 0
      real(kind_phys), dimension(:,:), intent(inout), optional ::        &
     &        edmf_a,edmf_w,edmf_qt,                                     &
     &        edmf_thl,edmf_ent,edmf_qc,                                 &
     &        sub_thl,sub_sqv,det_thl,det_sqv
      real(kind_phys), dimension(:,:), intent(inout) ::                  &
     &        t3d,qgrs_water_vapor,qgrs_liquid_cloud,qgrs_ice,           &
     &        qgrs_snow
      real(kind_phys), dimension(:,:), intent(in) ::                     &
     &        qgrs_cloud_ice_num_conc,                                   &
     &        u,v,omega,                                                 &
     &        exner,prsl,prsi,                                           &
     &        qgrs_ozone
      real(kind_phys), dimension(:,:), intent(in), optional  ::          &
     &        qgrs_water_aer_num_conc,                                   &
     &        qgrs_cloud_droplet_num_conc,                               &
     &        qgrs_ice_aer_num_conc
      real(kind_phys), dimension(:,:), intent(in), optional :: qgrs_cccn
      real(kind_phys), dimension(:,:), intent(out), optional ::          &
     &        Tsq, Qsq, Cov, exch_h, exch_m, dqke, qWT, qSHEAR, qBUOY,   &
     &        qDISS
      real(kind_phys), dimension(:), intent(in) :: xmu
      real(kind_phys), dimension(:,:), intent(in) :: htrsw, htrlw
      ! spp_wts_pbl only allocated if spp_pbl == 1
      real(kind_phys), dimension(:,:), intent(in), optional :: spp_wts_pbl

     !LOCAL
      real(kind_phys), dimension(im,levs) ::                             &
     &        sqv,sqc,sqi,sqs,qnc,qni,ozone,qnwfa,qnifa,qnbca,           &
     &        dz, w, p, rho, th, qv, delp,                               &
     &        RUBLTEN, RVBLTEN, RTHBLTEN, RQVBLTEN,                      &
     &        RQCBLTEN, RQNCBLTEN, RQIBLTEN, RQNIBLTEN, RQSBLTEN,        &
     &        RQNWFABLTEN, RQNIFABLTEN, RQNBCABLTEN
      real(kind_phys), allocatable :: old_ozone(:,:)

!smoke/chem arrays
      real(kind_phys), dimension(:), intent(inout), optional :: frp
      logical, intent(in) :: mix_chem, enh_mix, rrfs_sd
      real(kind_phys), dimension(:,:,:), intent(inout), optional :: chem3d
      real(kind_phys), dimension(:,:  ), intent(in), optional :: vdep
      real(kind_phys), dimension(im)   :: emis_ant_no

!MYNN-2D
      real(kind_phys), dimension(:), intent(in) ::                       &
     &        dx,zorl,slmsk,tsurf,qsfc,ps,                               &
     &        hflx,qflx,ust,wspd,rb,recmol

      real(kind_phys), dimension(:), intent(in), optional ::             &
     &        dusfc_cice,dvsfc_cice,dtsfc_cice,dqsfc_cice
      real(kind_phys), dimension(:), intent(in) ::                       &
     &        stress_wat,hflx_wat,qflx_wat,                              &
     &        oceanfrac,fice

      logical, dimension(:), intent(in) ::                               &
     &        wet, dry, icy

      real(kind_phys), dimension(:), intent(inout) ::                    &
     &        pblh,dusfc_diag,dvsfc_diag,dtsfc_diag,dqsfc_diag
      real(kind_phys), dimension(:), intent(out) ::                      &
     &        ch,dtsfc1,dqsfc1,dusfc1,dvsfc1,                            &
     &        dtsfci_diag,dqsfci_diag,dusfci_diag,dvsfci_diag
     real(kind_phys), dimension(:), intent(out), optional ::             &
     &        maxMF,maxwidth,ztop_plume
      integer, dimension(:), intent(inout) ::                            &
     &        kpbl
      integer, dimension(:), intent(inout), optional ::                  &
     &        ktop_plume

      real(kind_phys), dimension(:), intent(inout), optional ::          &
     &        dusfc_cpl,dvsfc_cpl,dtsfc_cpl,dqsfc_cpl
      real(kind_phys), dimension(:), intent(inout), optional ::          &
     &        dusfci_cpl,dvsfci_cpl,dtsfci_cpl,dqsfci_cpl

     !LOCAL
      real(kind_phys), dimension(im) ::                                  &
     &        hfx,qfx,rmol,xland,uoce,voce,znt,ts
      integer :: idtend
      real(kind_phys), dimension(im) :: dusfci1,dvsfci1,dtsfci1,dqsfci1
      real(kind_phys), allocatable :: save_qke_adv(:,:)
      real(kind_phys), dimension(levs) :: kzero

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lprnt) then
         write(0,*)"=============================================="
         write(0,*)"in mynn wrapper..."
         write(0,*)"flag_init=",flag_init
         write(0,*)"flag_restart=",flag_restart
      endif

      if (.not. flag_for_pbl_generic_tend .and. ldiag3d) then
         idtend = dtidx(ntke+100,index_of_process_pbl)
         if (idtend>=1) then
            allocate(save_qke_adv(im,levs))
            save_qke_adv=qke_adv
         endif
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

      kzero = zero !generic zero array
      !initialize arrays for test
      EMIS_ANT_NO = 0.

      FLAG_OZONE = ntoz>0

  ! Assign variables for each microphysics scheme
        if (imp_physics == imp_physics_wsm6 .or. imp_physics == imp_physics_fa) then
  ! WSM6 or Ferrier-Aligo
         FLAG_QI = .true.
         FLAG_QNI= .false.
         FLAG_QC = .true.
         FLAG_QNC= .false.
         FLAG_QS = .false.
         FLAG_QNWFA= .false.
         FLAG_QNIFA= .false.
         FLAG_QNBCA= .false.
         do k=1,levs
            do i=1,im
              sqv(i,k)   = qgrs_water_vapor(i,k)
              sqc(i,k)   = qgrs_liquid_cloud(i,k)
              sqi(i,k)   = qgrs_ice(i,k)
              sqs(i,k)   = 0.
              ozone(i,k) = qgrs_ozone(i,k)
              qnc(i,k)   = 0.
              qni(i,k)   = 0.
              qnwfa(i,k) = 0.
              qnifa(i,k) = 0.
              qnbca(i,k) = 0.
            enddo
          enddo
        elseif (imp_physics == imp_physics_nssl ) then
  ! NSSL
         FLAG_QI = .true.
         FLAG_QNI= .true.
         FLAG_QC = .true.
         FLAG_QNC= .true.
         FLAG_QS = .true.
         FLAG_QNWFA= nssl_ccn_on ! ERM: Perhaps could use this field for CCN field?
         FLAG_QNIFA= .false.
         FLAG_QNBCA= .false.
         do k=1,levs
            do i=1,im
              sqv(i,k)   = qgrs_water_vapor(i,k)
              sqc(i,k)   = qgrs_liquid_cloud(i,k)
              sqi(i,k)   = qgrs_ice(i,k)
              sqs(i,k)   = qgrs_snow(i,k)
              ozone(i,k) = qgrs_ozone(i,k)
              qnc(i,k)   = qgrs_cloud_droplet_num_conc(i,k)
              qni(i,k)   = qgrs_cloud_ice_num_conc(i,k)
              qnwfa(i,k) = 0.
              IF ( nssl_ccn_on ) THEN
                qnwfa(i,k) = qgrs_cccn(i,k)
              ENDIF
              qnifa(i,k) = 0.
              qnbca(i,k) = 0.
            enddo
          enddo
        elseif (imp_physics == imp_physics_thompson) then
  ! Thompson
          if(ltaerosol) then
            FLAG_QI = .true.
            FLAG_QNI= .true.
            FLAG_QC = .true.
            FLAG_QS = .true. !pipe it in, but do not mix
            FLAG_QNC= .true.
            FLAG_QNWFA= .true.
            FLAG_QNIFA= .true.
            FLAG_QNBCA= .false.
            do k=1,levs
              do i=1,im
                sqv(i,k)   = qgrs_water_vapor(i,k)
                sqc(i,k)   = qgrs_liquid_cloud(i,k)
                sqi(i,k)   = qgrs_ice(i,k)
                sqs(i,k)   = qgrs_snow(i,k)
                qnc(i,k)   = qgrs_cloud_droplet_num_conc(i,k)
                qni(i,k)   = qgrs_cloud_ice_num_conc(i,k)
                ozone(i,k) = qgrs_ozone(i,k)
                qnwfa(i,k) = qgrs_water_aer_num_conc(i,k)
                qnifa(i,k) = qgrs_ice_aer_num_conc(i,k)
                qnbca(i,k) = 0.
              enddo
            enddo
          else if(mraerosol) then
            FLAG_QI = .true.
            FLAG_QNI= .true.
            FLAG_QC = .true.
            FLAG_QS = .true.
            FLAG_QNC= .true.
            FLAG_QNWFA= .false.
            FLAG_QNIFA= .false.
            FLAG_QNBCA= .false.
            do k=1,levs
              do i=1,im
                sqv(i,k)   = qgrs_water_vapor(i,k)
                sqc(i,k)   = qgrs_liquid_cloud(i,k)
                sqi(i,k)   = qgrs_ice(i,k)
                sqs(i,k)   = qgrs_snow(i,k)
                qnc(i,k)   = qgrs_cloud_droplet_num_conc(i,k)
                qni(i,k)   = qgrs_cloud_ice_num_conc(i,k)
                ozone(i,k) = qgrs_ozone(i,k)
                qnwfa(i,k) = 0.
                qnifa(i,k) = 0.
                qnbca(i,k) = 0.
              enddo
            enddo
          else
            FLAG_QI = .true.
            FLAG_QNI= .true.
            FLAG_QC = .true.
            FLAG_QS = .true.
            FLAG_QNC= .false.
            FLAG_QNWFA= .false.
            FLAG_QNIFA= .false.
            FLAG_QNBCA= .false.
            do k=1,levs
              do i=1,im
                sqv(i,k)   = qgrs_water_vapor(i,k)
                sqc(i,k)   = qgrs_liquid_cloud(i,k)
                sqi(i,k)   = qgrs_ice(i,k)
                sqs(i,k)   = qgrs_snow(i,k)
                qnc(i,k)   = 0.
                qni(i,k)   = qgrs_cloud_ice_num_conc(i,k)
                ozone(i,k) = qgrs_ozone(i,k)
                qnwfa(i,k) = 0.
                qnifa(i,k) = 0.
                qnbca(i,k) = 0.
              enddo
            enddo
          endif
        elseif (imp_physics == imp_physics_gfdl) then
  ! GFDL MP
          FLAG_QI = .true.
          FLAG_QNI= .false.
          FLAG_QC = .true.
          FLAG_QNC= .false.
          FLAG_QS = .false.
          FLAG_QNWFA= .false.
          FLAG_QNIFA= .false.
          FLAG_QNBCA= .false.
          do k=1,levs
            do i=1,im
                sqv(i,k)   = qgrs_water_vapor(i,k)
                sqc(i,k)   = qgrs_liquid_cloud(i,k)
                sqi(i,k)   = qgrs_ice(i,k)
                qnc(i,k)   = 0.
                qni(i,k)   = 0.
                sqs(i,k)   = 0.
                qnwfa(i,k) = 0.
                qnifa(i,k) = 0.
                qnbca(i,k) = 0.
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
          FLAG_QS = .false.
          FLAG_QNWFA= .false.
          FLAG_QNIFA= .false.
          FLAG_QNBCA= .false.
          do k=1,levs
            do i=1,im
                sqv(i,k)   = qgrs_water_vapor(i,k)
                sqc(i,k)   = qgrs_liquid_cloud(i,k)
                sqi(i,k)   = 0.
                sqs(i,k)   = 0.
                qnc(i,k)   = 0.
                qni(i,k)   = 0.
                qnwfa(i,k) = 0.
                qnifa(i,k) = 0.
                qnbca(i,k) = 0.
                ozone(i,k) = qgrs_ozone(i,k)
            enddo
          enddo
        endif
       if(ldiag3d .and. dtidx(100+ntoz,index_of_process_pbl)>1) then
         allocate(old_ozone(im,levs))
         old_ozone = ozone
       endif

       do k=1,levs
          do i=1,im
             th(i,k)=t3d(i,k)/exner(i,k)
             rho(i,k)=prsl(i,k)/(r_d*t3d(i,k)*(1.+p608*max(sqv(i,k),1e-8)))
             w(i,k) = -omega(i,k)/(rho(i,k)*grav)
          enddo
       enddo

  ! Check incoming moist species to ensure non-negative values
  ! First, create height difference (dz)
      do k=1,levs
         do i=1,im
            dz(i,k)=(phii(i,k+1) - phii(i,k))*g_inv
         enddo
      enddo

      do i=1,im
         do k=1,levs
            delp(i,k) = prsi(i,k) - prsi(i,k+1)
         enddo
      enddo

      do i=1,im
         call moisture_check2(levs, delt,            &
                              delp(i,:), exner(i,:), &
                              sqv(i,:),  sqc(i,:),   &
                              sqi(i,:),  kzero(:),   &
                              t3d(i,:)               )
      enddo

      !intialize more variables
      do i=1,im
         if (slmsk(i)==1. .or. slmsk(i)==2.) then !sea/land/ice mask (=0/1/2) in FV3
            xland(i)=1.0                          !but land/water = (1/2) in SFCLAY_mynn
         else
            xland(i)=2.0
         endif
         uoce(i)=0.0
         voce(i)=0.0
         !ust(i) = sqrt(stress(i))
         ch(i)=0.0
         hfx(i)=hflx(i)*rho(i,1)*cp
         qfx(i)=qflx(i)*rho(i,1)
         !filter bad incoming fluxes
         if (hfx(i) > 1200.)hfx(i) = 1200.
         if (hfx(i) < -500.)hfx(i) = -500.
         if (qfx(i) > .0005)qfx(i) = 0.0005
         if (qfx(i) < -.0002)qfx(i) = -0.0002

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
                dusfci_cpl(i) = -1.*rho(i,1)*stress_wat(i)*u(i,1)/wspd(i)   ! U-momentum flux
                dvsfci_cpl(i) = -1.*rho(i,1)*stress_wat(i)*v(i,1)/wspd(i)   ! V-momentum flux
              else
                dusfci_cpl(i) = zero
                dvsfci_cpl(i) = zero
              endif
              dtsfci_cpl(i) =  cp*rho(i,1)*hflx_wat(i) ! sensible heat flux over open ocean
              dqsfci_cpl(i) = XLV*rho(i,1)*qflx_wat(i) ! latent heat flux over open ocean
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
         print*,"tke_budget=",tke_budget," bl_mynn_tkeadvect=",bl_mynn_tkeadvect
         print*,"bl_mynn_cloudpdf=",bl_mynn_cloudpdf," bl_mynn_mixlength=",bl_mynn_mixlength
         print*,"bl_mynn_edmf=",bl_mynn_edmf," bl_mynn_edmf_mom=",bl_mynn_edmf_mom
         print*,"bl_mynn_edmf_tke=",bl_mynn_edmf_tke
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
         print*,"ch=",ch(1)
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
     &             delt=delt,dz=dz,dx=dx,znt=znt,                      &
     &             u=u,v=v,w=w,th=th,sqv3D=sqv,sqc3D=sqc,              &
     &             sqi3D=sqi,sqs3D=sqs,qnc=qnc,qni=qni,                &
     &             qnwfa=qnwfa,qnifa=qnifa,qnbca=qnbca,ozone=ozone,    &
     &             p=prsl,exner=exner,rho=rho,T3D=t3d,                 &
     &             xland=xland,ts=ts,qsfc=qsfc,ps=ps,                  &
     &             ust=ust,ch=ch,hfx=hfx,qfx=qfx,rmol=rmol,            &
     &             wspd=wspd,uoce=uoce,voce=voce,                      & !input
     &             qke=QKE,qke_adv=qke_adv,                            & !output
     &             sh3d=Sh3d,sm3d=Sm3d,                                &
!chem/smoke
     &             nchem=nchem,kdvel=kdvel,ndvel=ndvel,                &
     &             Chem3d=chem3d,Vdep=vdep,smoke_dbg=smoke_dbg,        &
     &             FRP=frp,EMIS_ANT_NO=emis_ant_no,                    &
     &             mix_chem=mix_chem,enh_mix=enh_mix,                  &
     &             rrfs_sd=rrfs_sd,                                    &
!-----
     &             Tsq=tsq,Qsq=qsq,Cov=cov,                            & !output
     &             RUBLTEN=RUBLTEN,RVBLTEN=RVBLTEN,RTHBLTEN=RTHBLTEN,  & !output
     &             RQVBLTEN=RQVBLTEN,RQCBLTEN=rqcblten,                &
     &             RQIBLTEN=rqiblten,RQNCBLTEN=rqncblten,              & !output
     &             RQSBLTEN=rqsblten,                                  & !output
     &             RQNIBLTEN=rqniblten,RQNWFABLTEN=RQNWFABLTEN,        & !output
     &             RQNIFABLTEN=RQNIFABLTEN,RQNBCABLTEN=RQNBCABLTEN,    & !output
     &             dozone=dqdt_ozone,                                  & !output
     &             EXCH_H=exch_h,EXCH_M=exch_m,                        & !output
     &             pblh=pblh,KPBL=KPBL,                                & !output
     &             el_pbl=el_pbl,                                      & !output
     &             dqke=dqke,                                          & !output
     &             qWT=qWT,qSHEAR=qSHEAR,qBUOY=qBUOY,qDISS=qDISS,      & !output
     &             bl_mynn_tkeadvect=bl_mynn_tkeadvect,                &
     &             tke_budget=tke_budget,                              & !input parameter
     &             bl_mynn_cloudpdf=bl_mynn_cloudpdf,                  & !input parameter
     &             bl_mynn_mixlength=bl_mynn_mixlength,                & !input parameter
     &             icloud_bl=icloud_bl,                                & !input parameter
     &             qc_bl=qc_bl,qi_bl=qi_bl,cldfra_bl=cldfra_bl,        & !output
     &             closure=bl_mynn_closure,bl_mynn_edmf=bl_mynn_edmf,  & !input parameter
     &             bl_mynn_edmf_mom=bl_mynn_edmf_mom,                  & !input parameter
     &             bl_mynn_edmf_tke=bl_mynn_edmf_tke,                  & !input parameter
     &             bl_mynn_mixscalars=bl_mynn_mixscalars,              & !input parameter
     &             bl_mynn_output=bl_mynn_output,                      & !input parameter
     &             bl_mynn_cloudmix=bl_mynn_cloudmix,                  & !input parameter
     &             bl_mynn_mixqt=bl_mynn_mixqt,                        & !input parameter
     &             edmf_a=edmf_a,edmf_w=edmf_w,edmf_qt=edmf_qt,        & !output
     &             edmf_thl=edmf_thl,edmf_ent=edmf_ent,edmf_qc=edmf_qc,& !output
     &             sub_thl3D=sub_thl,sub_sqv3D=sub_sqv,                &
     &             det_thl3D=det_thl,det_sqv3D=det_sqv,                &
     &             maxwidth=maxwidth,maxMF=maxMF,ztop_plume=ztop_plume,& !output
     &             ktop_plume=ktop_plume,                              & !output
     &             spp_pbl=spp_pbl,pattern_spp_pbl=spp_wts_pbl,        & !input
     &             RTHRATEN=htrlw,                                     & !input
     &             FLAG_QI=flag_qi,FLAG_QNI=flag_qni,                  & !input
     &             FLAG_QC=flag_qc,FLAG_QNC=flag_qnc,FLAG_QS=flag_qs,  & !input
     &             FLAG_QNWFA=FLAG_QNWFA,FLAG_QNIFA=FLAG_QNIFA,        & !input
     &             FLAG_QNBCA=FLAG_QNBCA,FLAG_OZONE=FLAG_OZONE,        & !input
     &             IDS=1,IDE=im,JDS=1,JDE=1,KDS=1,KDE=levs,            & !input
     &             IMS=1,IME=im,JMS=1,JME=1,KMS=1,KME=levs,            & !input
     &             ITS=1,ITE=im,JTS=1,JTE=1,KTS=1,KTE=levs             ) !input


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
        accum_duvt3dt: if(ldiag3d .or. lsidea) then
          call dtend_helper(index_of_x_wind,RUBLTEN)
          call dtend_helper(index_of_y_wind,RVBLTEN)
          call dtend_helper(index_of_temperature,RTHBLTEN,exner)
          if(ldiag3d) then
            call dtend_helper(100+ntoz,dqdt_ozone)
            ! idtend = dtidx(100+ntoz,index_of_process_pbl)
            ! if(idtend>=1) then
            !   dtend(:,:,idtend) = dtend(:,:,idtend) + (ozone-old_ozone)
            !   deallocate(old_ozone)
            ! endif
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
        if (imp_physics == imp_physics_wsm6 .or. imp_physics == imp_physics_fa) then
           ! WSM6
           do k=1,levs
             do i=1,im
               dqdt_water_vapor(i,k)  = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_liquid_cloud(i,k) = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_ice(i,k)          = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_snow(i,k)         = RQSBLTEN(i,k) !/(1.0 + qv(i,k))
               !dqdt_ozone(i,k)        = 0.0
             enddo
           enddo
           if(ldiag3d .and. .not. flag_for_pbl_generic_tend) then
             call dtend_helper(100+ntqv,RQVBLTEN)
             call dtend_helper(100+ntcw,RQCBLTEN)
             call dtend_helper(100+ntiw,RQIBLTEN)
           endif
           !Update moist species:
           !do k=1,levs
           !  do i=1,im
           !    qgrs_water_vapor(i,k)  = qgrs_water_vapor(i,k)  + (RQVBLTEN(i,k)/(1.0+RQVBLTEN(i,k)))*delt
           !    qgrs_liquid_cloud(i,k) = qgrs_liquid_cloud(i,k) + RQCBLTEN(i,k)*delt
           !    qgrs_ice(i,k)          = qgrs_ice(i,k)          + RQIBLTEN(i,k)*delt
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
                 dqdt_ice(i,k)                     = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_ice_num_conc(i,k)            = RQNIBLTEN(i,k)
                 dqdt_snow(i,k)                    = RQSBLTEN(i,k) !/(1.0 + qv(i,k))
                 !dqdt_ozone(i,k)                   = 0.0
                 dqdt_water_aer_num_conc(i,k)      = RQNWFABLTEN(i,k)
                 dqdt_ice_aer_num_conc(i,k)        = RQNIFABLTEN(i,k)
               enddo
             enddo
             if(ldiag3d .and. .not. flag_for_pbl_generic_tend) then
               call dtend_helper(100+ntqv,RQVBLTEN)
               call dtend_helper(100+ntcw,RQCBLTEN)
               call dtend_helper(100+ntlnc,RQNCBLTEN)
               call dtend_helper(100+ntiw,RQIBLTEN)
               call dtend_helper(100+ntinc,RQNIBLTEN)
               call dtend_helper(100+ntwa,RQNWFABLTEN)
               call dtend_helper(100+ntia,RQNIFABLTEN)
             endif
             !do k=1,levs
             !  do i=1,im
             !    qgrs_water_vapor(i,k)            = qgrs_water_vapor(i,k)    + (RQVBLTEN(i,k)/(1.0+RQVBLTEN(i,k)))*delt
             !    qgrs_liquid_cloud(i,k)           = qgrs_liquid_cloud(i,k)   + RQCBLTEN(i,k)*delt
             !    qgrs_ice(i,k)                    = qgrs_ice(i,k)            + RQIBLTEN(i,k)*delt
             !    qgrs_cloud_droplet_num_conc(i,k) = qgrs_cloud_droplet_num_conc(i,k) + RQNCBLTEN(i,k)*delt
             !    qgrs_cloud_ice_num_conc(i,k)     = qgrs_cloud_ice_num_conc(i,k)     + RQNIBLTEN(i,k)*delt
             !    !dqdt_ozone(i,k)        = 0.0
             !    !qgrs_water_aer_num_conc(i,k)     = qgrs_water_aer_num_conc(i,k)     + RQNWFABLTEN(i,k)*delt
             !    !qgrs_ice_aer_num_conc(i,k)       = qgrs_ice_aer_num_conc(i,k)       + RQNIFABLTEN(i,k)*delt
             !  enddo
             !enddo
           else if(mraerosol) then
             do k=1,levs
               do i=1,im
                 dqdt_water_vapor(i,k)             = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_liquid_cloud(i,k)            = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_cloud_droplet_num_conc(i,k)  = RQNCBLTEN(i,k)
                 dqdt_ice(i,k)                     = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_ice_num_conc(i,k)            = RQNIBLTEN(i,k)
                 dqdt_snow(i,k)                    = RQSBLTEN(i,k) !/(1.0 + qv(i,k))
               enddo
             enddo
             if(ldiag3d .and. .not. flag_for_pbl_generic_tend) then
               call dtend_helper(100+ntqv,RQVBLTEN)
               call dtend_helper(100+ntcw,RQCBLTEN)
               call dtend_helper(100+ntlnc,RQNCBLTEN)
               call dtend_helper(100+ntiw,RQIBLTEN)
               call dtend_helper(100+ntinc,RQNIBLTEN)
             endif
           else
             !Thompson (2008)
             do k=1,levs
               do i=1,im
                 dqdt_water_vapor(i,k)   = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_liquid_cloud(i,k)  = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_ice(i,k)           = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_ice_num_conc(i,k)  = RQNIBLTEN(i,k)
                 dqdt_snow(i,k)          = RQSBLTEN(i,k) !/(1.0 + qv(i,k))
                 !dqdt_ozone(i,k)         = 0.0
               enddo
             enddo
             if(ldiag3d .and. .not. flag_for_pbl_generic_tend) then
               call dtend_helper(100+ntqv,RQVBLTEN)
               call dtend_helper(100+ntcw,RQCBLTEN)
               call dtend_helper(100+ntiw,RQIBLTEN)
               call dtend_helper(100+ntinc,RQNIBLTEN)
               !call dtend_helper(100+ntsw,RQSBLTEN)
             endif
             !do k=1,levs
             !  do i=1,im
             !    qgrs_water_vapor(i,k)            = qgrs_water_vapor(i,k)    + (RQVBLTEN(i,k)/(1.0+RQVBLTEN(i,k)))*delt
             !    qgrs_liquid_cloud(i,k)           = qgrs_liquid_cloud(i,k)   + RQCBLTEN(i,k)*delt
             !    qgrs_ice(i,k)                    = qgrs_ice(i,k)            + RQIBLTEN(i,k)*delt
             !    qgrs_cloud_ice_num_conc(i,k)     = qgrs_cloud_ice_num_conc(i,k)     + RQNIBLTEN(i,k)*delt
             !    !dqdt_ozone(i,k)        = 0.0
             !  enddo
             !enddo
           endif !end thompson choice
        elseif (imp_physics == imp_physics_nssl) then
           ! NSSL
             do k=1,levs
               do i=1,im
                 dqdt_water_vapor(i,k)             = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_liquid_cloud(i,k)            = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_cloud_droplet_num_conc(i,k)  = RQNCBLTEN(i,k)
                 dqdt_ice(i,k)                     = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
                 dqdt_ice_num_conc(i,k)            = RQNIBLTEN(i,k)
                 dqdt_snow(i,k)                    = RQSBLTEN(i,k) !/(1.0 + qv(i,k))
                 IF ( nssl_ccn_on ) THEN ! 
                   dqdt_cccn(i,k)      = RQNWFABLTEN(i,k)
                 ENDIF
               enddo
             enddo

        elseif (imp_physics == imp_physics_gfdl) then
           ! GFDL MP
           do k=1,levs
             do i=1,im
               dqdt_water_vapor(i,k)   = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_liquid_cloud(i,k)  = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_ice(i,k)           = RQIBLTEN(i,k) !/(1.0 + qv(i,k))
               !dqdt_rain(i,k)          = 0.0
               !dqdt_snow(i,k)          = 0.0
               !dqdt_graupel(i,k)       = 0.0
               !dqdt_ozone(i,k)         = 0.0
             enddo
           enddo
           if(ldiag3d .and. .not. flag_for_pbl_generic_tend) then
             call dtend_helper(100+ntqv,RQVBLTEN)
             call dtend_helper(100+ntcw,RQCBLTEN)
             call dtend_helper(100+ntiw,RQIBLTEN)
           endif
           !do k=1,levs
           !  do i=1,im
           !    qgrs_water_vapor(i,k)            = qgrs_water_vapor(i,k)    + (RQVBLTEN(i,k)/(1.0+RQVBLTEN(i,k)))*delt
           !    qgrs_liquid_cloud(i,k)           = qgrs_liquid_cloud(i,k)   + RQCBLTEN(i,k)*delt
           !    qgrs_ice(i,k)                    = qgrs_ice(i,k)            + RQIBLTEN(i,k)*delt
           !    !dqdt_ozone(i,k)        = 0.0
           !  enddo
           !enddo
       else
!          print*,"In MYNN wrapper. Unknown microphysics scheme, imp_physics=",imp_physics
           do k=1,levs
             do i=1,im
               dqdt_water_vapor(i,k)   = RQVBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_liquid_cloud(i,k)  = RQCBLTEN(i,k) !/(1.0 + qv(i,k))
               dqdt_ice(i,k)           = 0.0
               !dqdt_rain(i,k)          = 0.0
               !dqdt_snow(i,k)          = 0.0
               !dqdt_graupel(i,k)       = 0.0
               !dqdt_ozone(i,k)         = 0.0
             enddo
           enddo
           if(ldiag3d .and. .not. flag_for_pbl_generic_tend) then
             call dtend_helper(100+ntqv,RQVBLTEN)
             call dtend_helper(100+ntcw,RQCBLTEN)
             call dtend_helper(100+ntiw,RQIBLTEN)
           endif
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
          print*,"ch=",ch(1)
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
          print*,"ztop_plume:",ztop_plume(1)," maxmf:",maxmf(1)
          print*,"maxwidth:",maxwidth(1)
          print*
       endif

       if(allocated(save_qke_adv)) then
         if(ldiag3d .and. .not. flag_for_pbl_generic_tend) then
           idtend = dtidx(100+ntke,index_of_process_pbl)
           if(idtend>=1) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + qke_adv-save_qke_adv
           endif
         endif
         deallocate(save_qke_adv)
       endif

  CONTAINS

    SUBROUTINE dtend_helper(itracer,field,mult)
      real(kind_phys), intent(in) :: field(im,levs)
      real(kind_phys), intent(in), optional :: mult(im,levs)
      integer, intent(in) :: itracer
      integer :: idtend
      
      idtend=dtidx(itracer,index_of_process_pbl)
      if(idtend>=1) then
        if(present(mult)) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + field*dtf*mult
        else
          dtend(:,:,idtend) = dtend(:,:,idtend) + field*dtf
        endif
      endif
    END SUBROUTINE dtend_helper

! ==================================================================
  SUBROUTINE moisture_check2(kte, delt, dp, exner, &
                             qv, qc, qi, qs, th    )
  !
  ! If qc < qcmin, qi < qimin, or qv < qvmin happens in any layer,
  ! force them to be larger than minimum value by (1) condensating 
  ! water vapor into liquid or ice, and (2) by transporting water vapor 
  ! from the very lower layer.
  ! 
  ! We then update the final state variables and tendencies associated
  ! with this correction. If any condensation happens, update theta/temperature too.
  ! Note that (qv,qc,qi,th) are the final state variables after
  ! applying corresponding input tendencies and corrective tendencies.

    implicit none
    integer,  intent(in)     :: kte
    real(kind_phys), intent(in)     :: delt
    real(kind_phys), dimension(kte), intent(in)     :: dp, exner
    real(kind_phys), dimension(kte), intent(inout)  :: qv, qc, qi, qs, th
    integer   k
    real ::  dqc2, dqi2, dqs2, dqv2, sum, aa, dum
    real, parameter :: qvmin1= 1e-8,    & !min at k=1
                       qvmin = 1e-20,   & !min above k=1
                       qcmin = 0.0,     &
                       qimin = 0.0

    do k = kte, 1, -1  ! From the top to the surface
       dqc2 = max(0.0, qcmin-qc(k)) !qc deficit (>=0)
       dqi2 = max(0.0, qimin-qi(k)) !qi deficit (>=0)
       dqs2 = max(0.0, qimin-qs(k)) !qs deficit (>=0)

       !update species
       qc(k)  = qc(k)  +  dqc2
       qi(k)  = qi(k)  +  dqi2
       qs(k)  = qs(k)  +  dqs2
       qv(k)  = qv(k)  -  dqc2 - dqi2 - dqs2
       !for theta
       !th(k)  = th(k)  +  xlvcp/exner(k)*dqc2 + &
       !                   xlscp/exner(k)*dqi2
       !for temperature
       th(k)  = th(k)  +  xlvcp*dqc2 + &
                          xlscp*(dqi2+dqs2)

       !then fix qv if lending qv made it negative
       if (k .eq. 1) then
          dqv2   = max(0.0, qvmin1-qv(k)) !qv deficit (>=0)
          qv(k)  = qv(k)  + dqv2
          qv(k)  = max(qv(k),qvmin1)
          dqv2   = 0.0
       else
          dqv2   = max(0.0, qvmin-qv(k)) !qv deficit (>=0)
          qv(k)  = qv(k)  + dqv2
          qv(k-1)= qv(k-1)  - dqv2*dp(k)/dp(k-1)
          qv(k)  = max(qv(k),qvmin)
       endif
       qc(k) = max(qc(k),qcmin)
       qi(k) = max(qi(k),qimin)
       qs(k) = max(qs(k),qimin)
    end do

    ! Extra moisture used to satisfy 'qv(1)>=qvmin' is proportionally
    ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
    ! preserves column moisture.
    if( dqv2 .gt. 1.e-20 ) then
        sum = 0.0
        do k = 1, kte
           if( qv(k) .gt. 2.0*qvmin ) sum = sum + qv(k)*dp(k)
        enddo
        aa = dqv2*dp(1)/max(1.e-20,sum)
        if( aa .lt. 0.5 ) then
            do k = 1, kte
               if( qv(k) .gt. 2.0*qvmin ) then
                   dum    = aa*qv(k)
                   qv(k)  = qv(k) - dum
               endif
            enddo
        else
        ! For testing purposes only (not yet found in any output):
        !    write(*,*) 'Full moisture conservation is impossible'
        endif
    endif

    return

  END SUBROUTINE moisture_check2

  END SUBROUTINE mynnedmf_wrapper_run

!###=================================================================

END MODULE mynnedmf_wrapper
