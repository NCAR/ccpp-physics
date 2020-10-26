!> \file m_micro.F90
!! This file contains the subroutine that call Morrison-Gettelman microphysics (MG1, MG2 and MG3)
!! MG1 forecasts cloud ice and cloud liquid and their number
!! MG2 forecasts cloud ice, cloud liquid, rain, snow,  and their number
!! MG3 forecasts cloud ice, cloud liquid, rain, snow, Graupel/Hail  and their number

!> This module contains the CCPP-compliant Morrison-Gettelman microphysics (MG1, MG2 and MG3) scheme.
module m_micro

  implicit none
  public :: m_micro_init, m_micro_run, m_micro_finalize
  private
  logical :: is_initialized = .False.

contains

!>\ingroup mg_driver
!! This subroutine is the MG initialization.
!> \section arg_table_m_micro_init  Argument Table
!! \htmlinclude m_micro_init.html
!!
subroutine m_micro_init(imp_physics, imp_physics_mg, fprcp, gravit, rair, rh2o, cpair,&
                        tmelt, latvap, latice, mg_dcs, mg_qcvar, mg_ts_auto_ice,      &
                        mg_rhmini, microp_uniform, do_cldice, hetfrz_classnuc,        &
                        mg_precip_frac_method, mg_berg_eff_factor, sed_supersat,      &
                        do_sb_physics, mg_do_hail,  mg_do_graupel, mg_nccons,         &
                        mg_nicons, mg_ngcons, mg_ncnst, mg_ninst, mg_ngnst,           &
                        mg_do_ice_gmao, mg_do_liq_liu, errmsg, errflg)

    use machine,            only: kind_phys
    use cldwat2m_micro,     only: ini_micro
    use micro_mg2_0,        only: micro_mg_init2_0 => micro_mg_init
    use micro_mg3_0,        only: micro_mg_init3_0 => micro_mg_init
    use aer_cloud,          only: aer_cloud_init

    integer,              intent(in) :: imp_physics, imp_physics_mg, fprcp
    logical,              intent(in) :: microp_uniform, do_cldice, hetfrz_classnuc,     &
                                        sed_supersat, do_sb_physics, mg_do_hail,        &
                                        mg_do_graupel, mg_nccons, mg_nicons, mg_ngcons, &
                                        mg_do_ice_gmao, mg_do_liq_liu
    real(kind=kind_phys), intent(in) :: gravit, rair, rh2o, cpair, tmelt, latvap, latice
    real(kind=kind_phys), intent(in) :: mg_dcs, mg_qcvar, mg_ts_auto_ice(2), mg_rhmini, &
                                        mg_berg_eff_factor, mg_ncnst, mg_ninst, mg_ngnst
    character(len=16),    intent(in) :: mg_precip_frac_method
    character(len=*),     intent(out) :: errmsg
    integer,              intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if (is_initialized) return

    if (imp_physics/=imp_physics_mg) then
       write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Morrison-Gettelman MP"
       errflg = 1
       return
    end if

    if (fprcp <= 0) then
      call ini_micro (mg_dcs, mg_qcvar, mg_ts_auto_ice(1))
    elseif (fprcp == 1) then
      call micro_mg_init2_0(kind_phys, gravit, rair, rh2o, cpair, &
                            tmelt, latvap, latice, mg_rhmini,     &
                            mg_dcs, mg_ts_auto_ice,               &
                            mg_qcvar,                             &
                            microp_uniform, do_cldice,            &
                            hetfrz_classnuc,                      &
                            mg_precip_frac_method,                &
                            mg_berg_eff_factor,                   &
                            sed_supersat, do_sb_physics,          &
                            mg_do_ice_gmao, mg_do_liq_liu,        &
                            mg_nccons, mg_nicons,                 &
                            mg_ncnst, mg_ninst)
    elseif (fprcp == 2) then
      call micro_mg_init3_0(kind_phys, gravit, rair, rh2o, cpair, &
                            tmelt, latvap, latice, mg_rhmini,     &
                            mg_dcs, mg_ts_auto_ice,               &
                            mg_qcvar,                             &
                            mg_do_hail,       mg_do_graupel,      &
                            microp_uniform,   do_cldice,          &
                            hetfrz_classnuc,                      &
                            mg_precip_frac_method,                &
                            mg_berg_eff_factor,                   &
                            sed_supersat, do_sb_physics,          &
                            mg_do_ice_gmao, mg_do_liq_liu,        &
                            mg_nccons,    mg_nicons,              &
                            mg_ncnst,     mg_ninst,               &
                            mg_ngcons,    mg_ngnst)
    else
      write(0,*)' fprcp = ',fprcp,' is not a valid option - aborting'
      stop
    endif
    call aer_cloud_init ()

    is_initialized = .true.

end subroutine m_micro_init

! \brief Brief description of the subroutine
!
!> \section arg_table_m_micro_finalize  Argument Table
!!
       subroutine m_micro_finalize
       end subroutine m_micro_finalize

!> \defgroup mg2mg3 Morrison-Gettelman MP scheme Module
!! This module contains the the entity of MG2 and MG3 schemes. 
!> @{
!> \defgroup mg_driver Morrison-Gettelman MP Driver Module
!! \brief This subroutine is the Morrison-Gettelman MP driver, which computes 
!! grid-scale condensation and evaporation of cloud condensate.

#if 0

!> \section arg_table_m_micro_run Argument Table
!! \htmlinclude m_micro_run.html
!!
#endif
!>\ingroup mg_driver
!>\section detail_m_micro_run MG m_micro_run Detailed Algorithm
!> @{
      subroutine m_micro_run(   im,       ix,     lm,     flipv, dt_i   &
     &,                         prsl_i,   prsi_i, phil,   phii          &
     &,                         omega_i,  QLLS_i, QLCN_i, QILS_i, QICN_i&
     &,                         lwheat_i, swheat_i, w_upi, cf_upi       &
     &,                         FRLAND,   ZPBL, CNV_MFD_i               &
     &,                         CNV_DQLDT_i, CLCN_i, u_i, v_i           &
     &,                         TAUGWX,   TAUGWY                        &
     &,                         TAUOROX,  TAUOROY, CNV_FICE_i           &
     &,                         CNV_NDROP_i,CNV_NICE_i, q_io, lwm_o     &
     &,                         qi_o,     t_io,    rn_o, sr_o           &
     &,                         ncpl_io,  ncpi_io, fprcp, rnw_io, snw_io&
     &,                         qgl_io,   ncpr_io, ncps_io, ncgl_io     &
     &,                         CLLS_io,  KCBL                          &
     &,                         CLDREFFL, CLDREFFI, CLDREFFR, CLDREFFS  &
     &,                         CLDREFFG, aerfld_i                      &
     &,                         aero_in,  naai_i, npccn_i, iccn         &
     &,                         skip_macro                              &
     &,                         lprnt, alf_fac, qc_min, pdfflag         &
     &,                         ipr, kdt, xlat, xlon, rhc_i,            &
     &                          errmsg, errflg)

       use machine ,      only: kind_phys
       use physcons,           grav   => con_g,    pi     => con_pi,    &
     &                         rgas   => con_rd,   cp     => con_cp,    &
     &                         hvap   => con_hvap, hfus   => con_hfus,  &
     &                         ttp    => con_ttp,  tice   => con_t0c,   &
     &                         eps    => con_eps,  epsm1  => con_epsm1, &
     &                         VIREPS => con_fvirt,                     &
     &                         latvap => con_hvap, latice => con_hfus

!      use funcphys,      only: fpvs                !< saturation vapor pressure for water-ice mixed
!      use funcphys,      only: fpvsl, fpvsi, fpvs  !< saturation vapor pressure for water,ice & mixed
       use aer_cloud,     only: AerProps, getINsubset,init_aer,         &
     &                          aerosol_activate,AerConversion1
       use cldmacro,      only: macro_cloud,meltfrz_inst,update_cld,    &
     &                          meltfrz_inst, fix_up_clouds_2M
       use cldwat2m_micro,only: mmicro_pcond
       use micro_mg2_0,   only: micro_mg_tend2_0 => micro_mg_tend, qcvar2 => qcvar
       use micro_mg3_0,   only: micro_mg_tend3_0 => micro_mg_tend, qcvar3 => qcvar
       ! DH* TODO - make this an input argument, no cross-import!
       use aerclm_def,    only: ntrcaer

!      use wv_saturation, only: aqsat

       implicit none
!   Anning Cheng July 2015 writing the interface for GSM. Based on GMAO version of M-2M,
!    and Donifan's nuclei activation, notice the vertical coordinate is top-down
!    opposite to the GSM dynamic core, much work is still needed to consistently
!    treat other parts of the model
!    Anning Cheng 9/29/2017 implemented the MG2 from NCAR
!                 alphar8 for qc_var scaled from climatology value
!
!   Feb 2018 : S. Moorthi Updated for MG3 with graupel as prognostic variable
!------------------------------------
!   input
!      real,   parameter  :: r_air = 3.47d-3
       real,   parameter  :: one=1.0, oneb3=one/3.0, onebcp=one/cp,      &
     &                       kapa=rgas*onebcp,  cpbg=cp/grav,            &
     &                       lvbcp=hvap*onebcp, lsbcp=(hvap+hfus)*onebcp,&
     &                       qsmall=1.e-14, rainmin = 1.0e-13,           &
     &                       fourb3=4.0/3.0, RL_cub=1.0e-15, nmin=1.0

       integer, parameter :: ncolmicro = 1
       integer,intent(in) :: im, ix,lm, ipr, kdt, fprcp, pdfflag
       logical,intent(in) :: flipv, aero_in, skip_macro, lprnt, iccn
       real (kind=kind_phys), intent(in):: dt_i, alf_fac, qc_min(2)

       real (kind=kind_phys), dimension(ix,lm),intent(in)  ::           &
     &                prsl_i,u_i,v_i,phil,   omega_i, QLLS_i,QILS_i,    &
     &                                       lwheat_i,swheat_i
       real (kind=kind_phys), dimension(ix,0:lm),intent(in):: prsi_i,   &
     &                                                        phii
! GJF* These variables are conditionally allocated depending on whether the
!     Morrison-Gettelman microphysics is used, so they must be declared 
!     using assumed shape.
       real (kind=kind_phys), dimension(:,:),  intent(in)  ::           &
     &       CNV_DQLDT_i, CLCN_i,     QLCN_i, QICN_i,                   &
     &       CNV_MFD_i,               cf_upi, CNV_FICE_i, CNV_NDROP_i,  &
     &       CNV_NICE_i,  w_upi
! *GJF
       real (kind=kind_phys), dimension(im,lm),intent(in)  ::           &
     &       rhc_i, naai_i, npccn_i
       real (kind=kind_phys), dimension(im,lm,ntrcaer),intent(in) ::    &
     &       aerfld_i
       real (kind=kind_phys),dimension(im),intent(in):: TAUGWX,         &
     &       TAUGWY, TAUOROX, TAUOROY, FRLAND,ZPBL,xlat,xlon
!    &       TAUGWY, TAUX, TAUY, TAUOROX, TAUOROY,ps_i,FRLAND,ZPBL
!    &       CNVPRCP

!   output
       real (kind=kind_phys),dimension(ix,lm), intent(out) :: lwm_o, qi_o,  &
                        cldreffl, cldreffi, cldreffr, cldreffs, cldreffg
       real (kind=kind_phys),dimension(im), intent(out) :: rn_o,  sr_o
       character(len=*),                    intent(out) :: errmsg
       integer,                             intent(out) :: errflg

!   input and output
!      Anning Cheng 10/24/2016 twat for total water, diagnostic purpose
       integer, dimension(IM), intent(inout):: KCBL
       real (kind=kind_phys),dimension(ix,lm),intent(inout):: q_io, t_io,   &
     &                                             ncpl_io,ncpi_io,CLLS_io
! GJF* These variables are conditionally allocated depending on whether the
!     Morrison-Gettelman microphysics is used, so they must be declared 
!     using assumed shape.
       real (kind=kind_phys),dimension(:,:),intent(inout):: rnw_io,snw_io,&
     &                                             ncpr_io, ncps_io,        &
     &                                             qgl_io,  ncgl_io
! *GJF
!Moo   real (kind=kind_phys),dimension(im,lm),intent(inout):: CLLS_io


!   Local variables
       integer kcldtopcvn,i,k,ll, kbmin, NAUX, nbincontactdust,l
       integer, dimension(im) :: kct
       real (kind=kind_phys) T_ICE_ALL, USE_AV_V,BKGTAU,LCCIRRUS,       &
     &    NPRE_FRAC, Nct, Wct, fcn, ksa1, tauxr8, DT_Moist, dt_r8,      &
     &    TMAXLL, USURF,LTS_UP, LTS_LOW, MIN_EXP, fracover, c2_gw, est3

       real(kind=kind_phys), allocatable, dimension(:,:) ::             &
     &            CNV_MFD,         CNV_FICE,CNV_NDROP,CNV_NICE
!    &            CNV_MFD,CNV_PRC3,CNV_FICE,CNV_NDROP,CNV_NICE

       real(kind=kind_phys), dimension(IM,LM)::ncpl,ncpi,omega,SC_ICE,  &
     & RAD_CF, radheat,Q1,U1,V1,     PLO, ZLO,    temp,                 &
     & QLLS, QLCN, QILS,QICN,        CNV_CVW,CNV_UPDF,                  &
!    & QLLS, QLCN, QILS,QICN,        CNV_CVW,CNV_UPDF,SMAXL,SMAXI,      &
!    & NHET_NUC, NLIM_NUC, CDNC_NUC,INC_NUC,CNN01,CNN04,CNN1,DNHET_IMM, &
     & NHET_NUC, NLIM_NUC, CDNC_NUC,INC_NUC,                 DNHET_IMM, &
     & NHET_IMM,NHET_DEP,NHET_DHF,DUST_IMM,DUST_DEP, DUST_DHF,WSUB,     &
     & SIGW_GW,SIGW_CNV,SIGW_TURB,                                      &
!    & SIGW_GW,SIGW_CNV,SIGW_TURB,SIGW_RC,REV_CN_X,REV_LS_X,            &
     &                                     rnw,snw,ncpr,ncps,qgl,ncgl,  &
!    & RSU_LS_X, ALPHT_X, DLPDF_X, DIPDF_X,rnw,snw,ncpr,ncps,qgl,ncgl,  &
!    & ACLL_CN_X,ACIL_CN_X, PFRZ, FQA,QCNTOT,QTOT,QL_TOT,qi_tot,blk_l,rhc
     &                            FQA,QL_TOT,qi_tot,blk_l,rhc

       real(kind=kind_phys) :: QCNTOT, QTOT

       real(kind=kind_phys), dimension(IM,LM):: CNV_DQLDT, CLCN, CLLS

!      real(kind=kind_phys), dimension(IM,LM):: DQRL_X,                 &
!      real(kind=kind_phys), dimension(IM,LM):: CNV_DQLDT, CLCN,CLLS,   &
!    &                                          CCN01,CCN04,CCN1

!      real(kind=kind_phys), allocatable, dimension(:,:) :: RHX_X       &
!    &,        CFPDF_X, VFALLSN_CN_X, QSNOW_CN, VFALLRN_CN_X, QRAIN_CN  &

       real(kind=kind_phys), allocatable, dimension(:,:) ::             &
     &                                               ALPHT_X, PFRZ
!    &                           QSNOW_CN, QRAIN_CN, ALPHT_X, PFRZ

!      real(kind=kind_phys), allocatable, dimension(:,:) ::             &
!    &                                QSNOW_CN,               QRAIN_CN  &
!!   &,        CFPDF_X,               QSNOW_CN,               QRAIN_CN  &
!    &,                                              ALPHT_X, PFRZ
!!   &,        REV_CN_X, RSU_CN_X, DLPDF_X, DIPDF_X, ALPHT_X, PFRZ      &
!!   &,        ACLL_CN_X, ACIL_CN_X, DQRL_X                             &
!!   &,        PFI_CN_X,  PFL_CN_X, QST3, DZET, QDDF3
!!     real(kind=kind_phys), allocatable, dimension(:) :: vmip

!      real(kind=kind_phys), dimension(IM,LM) :: QDDF3
!      real(kind=kind_phys), dimension(IM,LM):: QST3, DZET, QDDF3
!    &                                          MASS, RHX_X,  CFPDF_X,  &
!    &                                          VFALLSN_CN_X, QSNOW_CN, &
!    &                                          VFALLRN_CN_X, QRAIN_CN
!    &                                          VFALLRN_CN_X, QRAIN_CN, dum

       real(kind=kind_phys), dimension(IM,LM+1) :: ZET
       real(kind=kind_phys), dimension(IM,0:LM) :: PLE, kh

!      real(kind=kind_phys), dimension(IM,0:LM) :: PLE, PKE, kh
!    &,                                            PFI_CN_X, PFL_CN_X

       real(kind=kind_phys),dimension(LM)      :: rhdfdar8, rhu00r8,    &
     &         ttendr8,qtendr8, cwtendr8,npre8, npccninr8,ter8,         &
     &         plevr8,ndropr8,qir8,qcr8,wparc_turb,qvr8, nir8,ncr8,     &
     &         nimmr8,nsootr8,rnsootr8,omegr8,qrr8,qsr8,nrr8,nsr8,      &
     &         qgr8,  ngr8

       real(kind=kind_phys), dimension(1:LM,10) :: rndstr8,naconr8

!      real(kind=kind_phys), dimension(IM)      :: CN_PRC2,CN_SNR,CN_ARFX,&
!    &                                             LS_SNR, LS_PRC2, TPREC
       real(kind=kind_phys), dimension(IM)      :: LS_SNR, LS_PRC2
!    &                                             VMIP, twat

       real(kind=kind_phys), dimension (LM) :: uwind_gw,vwind_gw,       &
     &   tm_gw, pm_gw, nm_gw, h_gw, rho_gw, khaux, qcaux,               &
     &   dummyW , wparc_cgw,             cfaux,            dpre8,       &
     &   wparc_ls,wparc_gw, swparc,smaxliq,smaxicer8,nheticer8,         &
     &   nhet_immr8,dnhet_immr8,nhet_depr8,nhet_dhfr8,sc_icer8,         &
     &   dust_immr8,dust_depr8,dust_dhfr8,nlimicer8,cldfr8,liqcldfr8,   &
     &   icecldfr8,cldor8, pdelr8,                                      &
     &   rpdelr8,lc_turb,zmr8,ficer8,rate1ord_cw2pr, tlatr8, qvlatr8,   &
     &   qctendr8, qitendr8, nctendr8, nitendr8, effcr8, effc_fnr8,     &
     &   effir8,  nevaprr8, evapsnowr8, prainr8,                        &
     &   prodsnowr8, cmeoutr8, deffir8, pgamradr8, lamcradr8,qsoutr8,   &
     &   qroutr8,droutr8, qcsevapr8,qisevapr8, qvresr8,                 &
     &   cmeioutr8, dsoutr8, qcsinksum_rate1ord,qrtend,nrtend,          &
     &   qstend,    nstend,  alphar8, rhr8,                             &

     &   qgtend, ngtend, qgoutr8, ngoutr8, dgoutr8

       real(kind=kind_phys),  dimension(1)      :: prectr8, precir8

       real(kind=kind_phys), dimension (LM)    :: vtrmcr8,vtrmir8,      &
     &   qcsedtenr8,qisedtenr8, praor8,prcor8,mnucccor8, mnucctor8,     &
     &   msacwior8,psacwsor8, bergsor8,bergor8,meltor8, homoor8,        &
     &   qcresor8,                                                      &
     &   prcior8, praior8,qiresor8, mnuccror8,pracsor8, meltsdtr8,      &
     &   frzrdtr8,                                                      &
     &   ncalr8, ncair8, mnuccdor8, nnucctor8, nsoutr8, nroutr8,        &
     &   nnuccdor8, nnucccor8,naair8,                                   &
     &   nsacwior8, nsubior8, nprcior8, npraior8, npccnor8, npsacwsor8, &
     &   nsubcor8, npraor8, nprc1or8, tlatauxr8,pfrz_inc_r8,sadice,     &
     &   sadsnow,  am_evp_st, reff_rain, reff_snow,                     &
     &   umr,ums,qrsedten,qssedten,refl,arefl,areflz,frefl,csrfl,       &
     &   acsrfl,fcsrfl,rercld,qrout2,qsout2,nrout2,nsout2,drout2,       &
     &   dsout2,freqs,freqr,nfice,qcrat,prer_evap,                      &
!    graupel related
     &              reff_grau, umg,      qgsedtenr8, mnuccrior8,        &
     &   pracgr8,   psacwgr8,  pgsacwr8, pgracsr8,   prdgr8,   qmultgr8,&
     &   qmultrgr8, psacrr8,   npracgr8, nscngr8,    ngracsr8, nmultgr8,&
     &   nmultrgr8, npsacwgr8, qgout2,   ngout2,     dgout2,   freqg

       real(kind=kind_phys), dimension (0:LM) :: pi_gw, rhoi_gw,        &
     &                                           ni_gw, ti_gw

       real(kind=kind_phys), dimension(LM+1)    :: pintr8, kkvhr8
       real(kind=kind_phys), dimension(2:LM+1)  :: lflx, iflx, rflx,    &
                                                   sflx, gflx

!      real (kind=kind_phys), parameter :: disp_liu=2., ui_scale=1.0    &
!    &,                                    dcrit=20.0e-6                &
       real (kind=kind_phys), parameter :: disp_liu=1.0, ui_scale=1.0   &
     &,                                    dcrit=1.0e-6                 &
!    &,                                    ts_autice=1800.0             &
!    &,                                    ts_autice=3600.0             & !time scale
     &,                                    ninstr8 = 0.1e6              &
     &,                                    ncnstr8 = 100.0e6

       real(kind=kind_phys):: k_gw, maxkh, tausurf_gw, overscale, tx1, rh1_r8
       real(kind=kind_phys)::  t_ice_denom

       integer, dimension(1)           :: lev_sed_strt      ! sedimentation start level
       real(kind=kind_phys), parameter :: sig_sed_strt=0.05 ! normalized pressure at sedimentation start

       real(kind=kind_phys),dimension(3)  :: ccn_diag
       real(kind=kind_phys),dimension(58) :: cloudparams

       integer, parameter              :: CCN_PARAM=2, IN_PARAM=5

       real(kind=kind_phys), parameter ::fdust_drop=1.0,   fsoot_drop=0.1 &
     &,                                  sigma_nuc_r8=0.28,SCLMFDFR=0.03
!    &,                                  sigma_nuc_r8=0.28,SCLMFDFR=0.1

       type (AerProps), dimension (IM,LM)  :: AeroProps
       type (AerProps)                     :: AeroAux, AeroAux_b
       real, allocatable, dimension(:,:,:) :: AERMASSMIX

       logical :: use_average_v, ltrue, lprint

!==================================
!====2-moment Microhysics=
!================== Start Stratiform cloud processes==========================================
!set up initial values

       data USE_AV_V/1./, BKGTAU/0.015/, LCCIRRUS/500./, NPRE_FRAC/1./, &
     &      TMAXLL/296./, fracover/1./,  LTS_LOW/12./,   LTS_UP/24./,   &
     &      MIN_EXP/0.5/

       data cloudparams/                                                &
     &  10.0, 4.0  , 4.0   , 1.0   , 2.e-3, 8.e-4, 2.0  , 1.0   , -1.0  &
     &, 0.0 , 1.3  , 1.0e-9, 3.3e-4, 20.0 , 4.8  , 4.8  , 230.0 ,  1.0  &
     &, 1.0 , 230.0, 14400., 50.0  , 0.01 , 0.1  , 200.0, 0.0   ,  0.0  &
     &, 0.5 , 0.5  , 2000.0, 0.8   , 0.5  , -40.0, 1.0  , 4.0   ,  0.0  &
     &, 0.0 , 0.0  , 1.0e-3, 8.0e-4, 1.0  , 0.95 , 1.0  , 0.0   ,  900.0&
!    &, 0.0 , 0.0  , 1.0e-3, 8.0e-4, 1.0  , 0.95 , 1.0  , 0.0   ,  880.0&
!    &, 0.0 , 0.0  , 1.0e-3, 8.0e-4, 1.0  , 0.95 , 1.0  , 0.0   ,  980.0&
     &, 1.0 , 1.0  , 1.0   , 0.0   , 0.0  , 1.e-5, 2.e-5, 2.1e-5,  4.e-5&
!    &, 3e-5, 0.1  , 4.0   , 250./          ! Annings version
     &, 3e-5, 0.1  , 4.0   , 150./          ! Annings version
!    &, 3e-5, 0.1  , 1.0   , 150./

! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

!      rhr8 = 1.0
       if(flipv) then
         DO K=1, LM
           ll = lm-k+1
           DO I = 1,IM
             Q1(i,k)        = q_io(i,ll)
             U1(i,k)        = u_i(i,ll)
             V1(i,k)        = v_i(i,ll)
             omega(i,k)     = omega_i(i,ll)
             ncpl(i,k)      = ncpl_io(i,ll)
             ncpi(i,k)      = ncpi_io(i,ll)
             rnw(i,k)       = rnw_io(i,ll)
             snw(i,k)       = snw_io(i,ll)
             qgl(i,k)       = qgl_io(i,ll)
             ncpr(i,k)      = ncpr_io(i,ll)
             ncps(i,k)      = ncps_io(i,ll)
             ncgl(i,k)      = ncgl_io(i,ll)
!                                                  QLLS is the total cloud water
             QLLS(i,k)      = QLLS_i(i,ll)-QLCN_i(i,ll)
             QLCN(i,k)      = QLCN_i(i,ll)
             QILS(i,k)      = QILS_i(i,ll)-QICN_i(i,ll)
             QICN(i,k)      = QICN_i(i,ll)
             CNV_CVW(i,k)   = w_upi(i,ll)
             CNV_UPDF(i,k)  = cf_upi(i,ll)
             CNV_DQLDT(I,K) = CNV_DQLDT_i(I,ll)
             CLCN(I,k)      = CLCN_i(I,ll)
             CLLS(I,k)      = max(CLLS_io(I,ll)-CLCN_i(I,ll),0.0)
             PLO(i,k)       = prsl_i(i,ll)*0.01
             zlo(i,k)       = phil(i,ll) * (1.0/grav)
             temp(i,k)      = t_io(i,ll)
             radheat(i,k)   = lwheat_i(i,ll) + swheat_i(i,ll)
             rhc(i,k)       = rhc_i(i,ll)
             if (iccn) then
               CDNC_NUC(i,k) = npccn_i(i,ll)
               INC_NUC(i,k)  = naai_i (i,ll)
             endif

           END DO
         END DO
         DO K=0, LM
           ll = lm-k
           DO I = 1,IM
             PLE(i,k)   = prsi_i(i,ll) *.01      ! interface pressure in hPa
             zet(i,k+1) = phii(i,ll) * (1.0/grav)
           END DO
         END DO
         if (.not. skip_macro) then
!          allocate(CNV_MFD(im,lm),   CNV_PRC3(im,lm), CNV_FICE(im,lm)  &
           allocate(CNV_MFD(im,lm),                    CNV_FICE(im,lm)  &
     &,             CNV_NDROP(im,lm), CNV_NICE(im,lm))
           DO K=1, LM
             ll = lm-k+1
             DO I = 1,IM
               CNV_MFD(i,k)   = CNV_MFD_i(i,ll)
!              CNV_PRC3(i,k)  = CNV_PRC3_i(i,ll)
               CNV_FICE(i,k)  = CNV_FICE_i(i,ll)
               CNV_NDROP(i,k) = CNV_NDROP_i(i,ll)
               CNV_NICE(i,k)  = CNV_NICE_i(i,ll)
             enddo
           enddo
         endif

       else
         DO K=1, LM
           DO I = 1,IM
             Q1(i,k)        = q_io(i,k)
             U1(i,k)        = u_i(i,k)
             V1(i,k)        = v_i(i,k)
             omega(i,k)     = omega_i(i,k)
             ncpl(i,k)      = ncpl_io(i,k)
             ncpi(i,k)      = ncpi_io(i,k)
             ncpi(i,k)      = ncpi_io(i,k)
             rnw(i,k)       = rnw_io(i,k)
             snw(i,k)       = snw_io(i,k)
             qgl(i,k)       = qgl_io(i,k)
             ncpr(i,k)      = ncpr_io(i,k)
             ncps(i,k)      = ncps_io(i,k)
             ncgl(i,k)      = ncgl_io(i,k)
!                                                  QLLS is the total cloud water
             QLLS(i,k)      = QLLS_i(i,k)-QLCN_i(i,k)
             QLCN(i,k)      = QLCN_i(i,k)
             QILS(i,k)      = QILS_i(i,k)-QICN_i(i,k)
             QICN(i,k)      = QICN_i(i,k)
             CNV_CVW(i,k)   = w_upi(i,k)
             CNV_UPDF(i,k)  = cf_upi(i,k)
             CNV_DQLDT(I,K) = CNV_DQLDT_i(I,k)
             CLCN(I,k)      = CLCN_i(I,k)
             CLLS(I,k)      = max(CLLS_io(I,k)-CLCN_i(I,k),0.0)
             PLO(i,k)       = prsl_i(i,k)*0.01
             zlo(i,k)       = phil(i,k) * (1.0/grav)
             temp(i,k)      = t_io(i,k)
             radheat(i,k)   = lwheat_i(i,k) + swheat_i(i,k)
             rhc(i,k)       = rhc_i(i,k)
             if (iccn) then
               CDNC_NUC(i,k) = npccn_i(i,k)
               INC_NUC(i,k)  = naai_i (i,k)
             endif

           END DO
         END DO
         DO K=0, LM
           DO I = 1,IM
             PLE(i,k)   = prsi_i(i,k) *.01      ! interface pressure in hPa
             zet(i,k+1) = phii(i,k) * (1.0/grav)
           END DO
         END DO
         if (.not. skip_macro) then
!          allocate(CNV_MFD(im,lm),   CNV_PRC3(im,lm), CNV_FICE(im,lm)  &
           allocate(CNV_MFD(im,lm),                    CNV_FICE(im,lm)  &
     &,             CNV_NDROP(im,lm), CNV_NICE(im,lm))
           DO K=1, LM
             DO I = 1,IM
               CNV_MFD(i,k)   = CNV_MFD_i(i,k)
!              CNV_PRC3(i,k)  = CNV_PRC3_i(i,k)
               CNV_FICE(i,k)  = CNV_FICE_i(i,k)
               CNV_NDROP(i,k) = CNV_NDROP_i(i,k)
               CNV_NICE(i,k)  = CNV_NICE_i(i,k)
             enddo
           enddo
         endif
       endif
!
       DT_MOIST = dt_i
       dt_r8    = dt_i

       if (kdt == 1) then
         DO K=1, LM
           DO I = 1,IM
             CALL fix_up_clouds_2M(Q1(I,K),   TEMP(i,k), QLLS(I,K),     &
     &                            QILS(I,K), CLLS(I,K), QLCN(I,K),      &
     &                            QICN(I,K), CLCN(I,K), NCPL(I,K),      &
     &                            NCPI(I,K), qc_min)
             if (rnw(i,k) <= qc_min(1)) then
               ncpl(i,k) = 0.0
             elseif (ncpl(i,k) <= nmin) then ! make sure NL > 0 if Q >0
               ncpl(i,k) = max(rnw(i,k) / (fourb3 * PI *RL_cub*997.0), nmin)
             endif
             if (snw(i,k) <= qc_min(2)) then
               ncpl(i,k) = 0.0
             elseif (ncps(i,k) <= nmin) then
               ncps(i,k) = max(snw(i,k) / (fourb3 * PI *RL_cub*500.0), nmin)
             endif
             if (qgl(i,k) <= qc_min(2)) then
               ncgl(i,k) = 0.0
             elseif (ncgl(i,k) <= nmin) then
               ncgl(i,k) = max(qgl(i,k) / (fourb3 * PI *RL_cub*500.0), nmin)
             endif

           enddo
         enddo
       endif
       do i=1,im
         KCBL(i)     = max(LM-KCBL(i),10)
         KCT(i)      = 10
       enddo

       DO I=1, IM
         DO K = LM-2, 10, -1
           If ((CNV_DQLDT(I,K)  <= 1.0e-9) .and.                        &
     &         (CNV_DQLDT(I,K+1) > 1.0e-9)) then
             KCT(I) = K+1
             exit
           end if
         END DO
       END DO

!      do L=LM,1,-1
!        do i=1,im
!          DZET(i,L)   = ZET(i,L)   - ZET(i,L+1)
!          tx1       = plo(i,l)*100.0
!          est3      = min(tx1, fpvs(temp(i,l)))
!          qst3(i,l) = min(eps*est3/max(tx1+epsm1*est3,1.0e-10),1.0)
!          MASS(i,l) = (ple(i,l) - ple(i,l-1)) * (100.0/grav)
!        enddo
!      enddo
!------------------------------------------------------------------------------
!      call aqsat(temp,plo*100.,est3,qst3,im,im,lm,1,lm)
!      do k=1,lm
!        do i=1,im
!          DZET(i,k) = TH1(i,k) * (pke(i,k)-pke(i,k-1))           &
!    &                          * cpbg * (1.0 + vireps*q1(i,k))
!          MASS(i,k) = (ple(i,k) - ple(i,k-1)) * (100.0/grav)
!        end do
!      end do

!      do k=1,lm
!        do i=1,im
!          temp(i,k) = th1(i,k) * PK(i,k)
!          est3      = fpvs(temp(i,k))
!          qst3(i,k) = min(eps*est3/max(plo(i,k)*100.0+epsm1*est3,1.0e-10),1.0)
!        enddo
!      enddo
!      call aqsat(temp,plo*100.,est3,qst3,im,im,lm,1,lm)
!      do k=1,lm
!        do i=1,im
!          DZET(i,k) = TH1(i,k) * (pke(i,k)-pke(i,k-1))           &
!    &                          * cpbg * (1.0 + vireps*q1(i,k))
!          MASS(i,k) = (ple(i,k) - ple(i,k-1)) * (100.0/grav)
!        enddo
!      enddo

!      do i=1,im
!        ZET(i,LM+1) = 0.0
!        vmip(i)     = 0.0
!      enddo
!------------------------------------------------------------------------------

!      if (.not. skip_macro) then
!        allocate(qddf3(im,lm))
!        allocate(vmip(im))
!        do i=1,im
!          vmip(i) = 0.0
!        enddo
!        DO K = LM, 1, -1
!          do i=1,im
!            if (zet(i,k) < 3000.0) then
!!             qddf3(i,k) = - (zet(i,k) - 3000.0) * zet(i,k) * mass(i,k)
!              qddf3(i,k) = - (zet(i,k) - 3000.0) * zet(i,k)            &
!    &                    * (ple(i,k) - ple(i,k-1)) * (100.0/grav)
!            else
!              qddf3(i,k) = 0.0
!            endif
!            vmip(i) = vmip(i) + qddf3(i,k)
!          enddo
!        END DO
!        do i=1,im
!          if (vmip(i) /= 0.0) vmip(i) = 1.0 / vmip(i)
!        enddo
!        DO K = 1,LM
!          do i=1,im
!            QDDF3(i,K) = QDDF3(i,K) * VMIP(i)
!          enddo
!        END DO
!        deallocate (vmip)
!      endif


       do l=lm-1,1,-1
         do i=1,im
           tx1     = 0.5 * (temp(i,l+1) + temp(i,l))
           kh(i,l) = 3.55e-7*tx1**2.5*(rgas*0.01) / ple(i,l)   !kh molecule diff only needing refinement
         enddo
       end do
       do i=1,im
         kh(i,0)  = kh(i,1)
         kh(i,lm) = kh(i,lm-1)
       enddo
       do L=LM,1,-1
         do i=1,im
           blk_l(i,l)  = 1.0 / ( 1.0/max(0.15*ZPBL(i),0.4*zlo(i,lm-1))&
     &                 + 1.0/(zlo(i,l)*.4) )

           SC_ICE(i,l) = 1.0
           NCPL(i,l)   = MAX( NCPL(i,l), 0.)
           NCPI(i,l)   = MAX( NCPI(i,l), 0.)
           RAD_CF(i,l) = max(0.0, min(CLLS(i,l)+CLCN(i,l), 1.0))
           if (.not. iccn) then
             CDNC_NUC(i,l) = 0.0
             INC_NUC(i,l)  = 0.0
           endif

         enddo
       end do
!      T_ICE_ALL = TICE - 40.0
       T_ICE_ALL = CLOUDPARAMS(33) + TICE
       t_ice_denom = 1.0 / (tice - t_ice_all)


       do l=1,lm
         rhdfdar8(l)  = 1.e-8
         rhu00r8(l)   = 0.95

         ttendr8(l)   = 0.
         qtendr8(l)   = 0.
         cwtendr8(l)  = 0.

         npccninr8(l) = 0.
       enddo
       do k=1,10
         do l=1,lm
           rndstr8(l,k) = 2.0e-7
         enddo
       enddo

!need an estimate of convective area
!=======================================================================================================================
!=======================================================================================================================
!> -# Nucleation of cloud droplets and ice crystals 
!! Aerosol cloud interactions. Calculate maxCCN tendency using Fountoukis and Nenes (2005) or Abdul Razzak and Ghan (2002)
!! liquid Activation Parameterization
!! Ice activation follows the Barahona & Nenes ice activation scheme, ACP, (2008, 2009).
!! Written by Donifan Barahona and described in Barahona et al. (2013)
!=======================================================================================================================
!=======================================================================================================================
!=======================================================================================================================
!      if(aero_in) then
!        allocate(AERMASSMIX  (IM,LM, 15))
!        AERMASSMIX = 1.e-15
!        call AerConversion1 (AERMASSMIX,  AeroProps)
!        deallocate(AERMASSMIX)
!      end if

!
!>  - Call init_aer()
       do k=1,lm
         do i=1,im
           call init_Aer(AeroProps(I, K))
         enddo
       enddo
!

       allocate(AERMASSMIX(IM,LM,15))
       if ( aero_in ) then
         AERMASSMIX(:,:,1:ntrcaer) = aerfld_i(:,:,1:ntrcaer)
       else
         AERMASSMIX(:,:,1:5) = 1.e-6
         AERMASSMIX(:,:,6:15) = 2.e-14
       end if
!>  - Call aerconversion1()
       call AerConversion1 (AERMASSMIX,  AeroProps)
       deallocate(AERMASSMIX)

       use_average_v = .false.
       if (USE_AV_V > 0.0) then
         use_average_v = .true.
       end if

       k_gw = (pi+pi) / LCCIRRUS

!-------------------------------------------------------------------------------
       do I=1,IM                         ! beginning of first big I loop

         kcldtopcvn = KCT(I)

         tausurf_gw = min(0.5*SQRT(TAUOROX(I)*TAUOROX(I)                &
     &                           + TAUOROY(I)*TAUOROY(I)), 10.0)
         do k=1,lm

           uwind_gw(k)   = min(0.5*SQRT( U1(I,k)*U1(I,k)                &
     &                                 + V1(I,k)*V1(I,k)), 50.0)

!  tausurf_gw   =tausurf_gw  + max (tausurf_gw,  min(0.5*SQRT(TAUX(I , J)**2+TAUY(I , J)**2), 10.0)*BKGTAU) !adds a minimum value from unresolved sources


           pm_gw(k)    = 100.0*PLO(I,k)
           tm_gw(k)    = TEMP(I,k)

           nm_gw(k)    = 0.0
           rho_gw(k)   = pm_gw(k) /(RGAS*tm_gw(k))

           ter8(k)    = TEMP(I,k)
           plevr8(k)  = 100.*PLO(I,k)
           ndropr8(k) = NCPL(I,k)
           qir8(k)    = QILS(I,k) + QICN(I,k)
           qcr8(k)    = QLLS(I,k) + QLCN(I,k)
           qcaux(k)   = qcr8(k)

           npccninr8(k) = 0.0
           naair8(k)    = 0.0

           npre8(k)     = 0.0

           if (RAD_CF(I,k) > 0.01 .and. qir8(k) > 0.0) then
             npre8(k) = NPRE_FRAC*NCPI(I,k)
           else
             npre8(k) = 0.0
           endif

           omegr8(k)      = OMEGA(I,k)
           lc_turb(k)     = max(blk_l(I,k), 50.0)
!          rad_cooling(k) = RADheat(I,k)

           if (npre8(k) > 0.0 .and. qir8(k) > 0.0) then
             dpre8(k) = ( qir8(k)/(6.0*npre8(k)*900.0*PI))**(1.0/3.0)
           else
             dpre8(k) = 1.0e-9
           endif
           wparc_ls(k) = -omegr8(k) / (rho_gw(k)*GRAV)                  &
     &                 +  cpbg * radheat(i,k)
!    &                 +  cpbg * rad_cooling(k)
         enddo
         do k=0,lm
           pi_gw(k)   = 100.0*PLE(I,k)
           rhoi_gw(k) = 0.0
           ni_gw(k)   = 0.0
           ti_gw(k)   = 0.0
         enddo


! ====================================================================
!> -# Call gw_prof() to calculate subgrid scale distribution in vertical velocity
! ====================================================================


         call gw_prof (1, LM, 1, tm_gw, pm_gw, pi_gw, rhoi_gw, ni_gw,   &
     &                 ti_gw, nm_gw, q1(i,:))

         do k=1,lm
           nm_gw(k)    = max(nm_gw(k), 0.005)
           h_gw(k)     = k_gw*rho_gw(k)*uwind_gw(k)*nm_gw(k)
           if (h_gw(K) > 0.0) then
             h_gw(K)    = sqrt(2.0*tausurf_gw/h_gw(K))
           end if

           wparc_gw(k)  = k_gw*uwind_gw(k)*h_gw(k)*0.133

           wparc_cgw(k) = 0.0
         end do

!>  - Subgrid variability from convective sources according to Barahona et al. 2014 (in preparation)

         if (kcldtopcvn > 20) then

           ksa1 = 1.0
           Nct  = nm_gw(kcldtopcvn)
           Wct  = max(CNV_CVW(I,kcldtopcvn), 0.0)

           fcn  = maxval(CNV_UPDF(I,kcldtopcvn:LM))

           do k=1,kcldtopcvn
             c2_gw    = (nm_gw(k) + Nct) / Nct
             wparc_cgw(k) = sqrt(ksa1*fcn*fcn*12.56*              &
     &                      1.806*c2_gw*c2_gw)*Wct*0.133
           enddo

         end if

         do k=1,lm
           dummyW(k) = 0.133*k_gw*uwind_gw(k)/nm_gw(k)
         enddo

         do K=1, LM-5, 1
           if (wparc_cgw(K)+wparc_gw(K) > dummyW(K)) then
             exit
           end if
         end do

         do l=1,min(k,lm-5)
           wparc_cgw(l) = 0.0
           wparc_gw(l)  = 0.0
         enddo



         kbmin = KCBL(I)
         kbmin = min(kbmin, LM-1) - 4
         do K = 1, LM
           wparc_turb(k) = KH(I,k) / lc_turb(k)
           dummyW(k)     = 10.0
         enddo

         if (FRLAND(I)  < 0.1   .and. ZPBL(I)    < 800.0 .and.          &
     &       TEMP(I,LM) < 298.0 .and. TEMP(I,LM) > 274.0 ) then
           do K = 1, LM
             dummyW(k) = max(min((ZET(I,k+1)-ZPBL(I))*0.01,10.0),-10.0)
             dummyW(k) = 1.0 / (1.0+exp(dummyW(k)))
           enddo
           maxkh = max(maxval(KH(I,kbmin:LM-1)*nm_gw(kbmin:LM-1)/       &
     &                                                  0.17), 0.3)
           do K = 1, LM
             wparc_turb(k) = (1.0-dummyW(k))*wparc_turb(k)              &
     &                     + dummyW(k)*maxkh
           enddo

         end if

         wparc_turb(kbmin:LM) = max(wparc_turb(kbmin:LM), 0.2)



!>  - Compute total variance

         do K = 1, LM
           swparc(k) = sqrt(wparc_gw(k)   * wparc_gw(k)                 &
     &                    + wparc_turb(k) * wparc_turb(k)               &
     &                    + wparc_cgw(k)  * wparc_cgw(k))
         enddo


! ==========================================================================================
! ========================Activate the aerosols ============================================

         do K = 1, LM

           if (plevr8(K) > 70.0) then

             ccn_diag(1) = 0.001
             ccn_diag(2) = 0.004
             ccn_diag(3) = 0.01

             if (K > 2 .and. K <= LM-2) then
               tauxr8 = (ter8(K-1) + ter8(K+1) + ter8(K)) * oneb3
             else
               tauxr8 = ter8(K)
             endif

!            if(aero_in) then
               AeroAux = AeroProps(I, K)
!            else
!              call init_Aer(AeroAux)
!              call init_Aer(AeroAux_b)
!            endif

             pfrz_inc_r8(k) = 0.0
             rh1_r8         = 0.0 !related to cnv_dql_dt, needed to changed soon

!     if (lprnt) write(0,*)' bef aero npccninr8=',npccninr8(k),' k=',k  &
!    &,' ccn_param=',ccn_param,' in_param=',in_param                    &
!    &,' AeroAux%kap=',AeroAux%kap

!> -# Call aerosol_activate() to activate the aerosols
             call aerosol_activate(tauxr8, plevr8(K), swparc(K),        &
     &            wparc_ls(K), AeroAux, npre8(k), dpre8(k), ccn_diag,   &
     &            ndropr8(k),          npccninr8(K), smaxliq(K),        &
!    &            ndropr8(k), qcr8(K), npccninr8(K), smaxliq(K),        &
     &            naair8(K), smaxicer8(K), nheticer8(K), nhet_immr8(K), &
     &            dnhet_immr8(K), nhet_depr8(k), nhet_dhfr8(k),         &
     &            sc_icer8(k),    dust_immr8(K), dust_depr8(k),         &
     &            dust_dhfr8(k),  nlimicer8(k),  use_average_v,         &
     &            CCN_PARAM,      IN_PARAM,      fdust_drop,            &
     &            fsoot_drop,pfrz_inc_r8(K),sigma_nuc_r8, rh1_r8,       &
     &            size(ccn_diag))
!    &            size(ccn_diag), lprnt)
!     if (lprnt) write(0,*)' aft aero npccninr8=',npccninr8(k),' k=',k

             if (npccninr8(k) < 1.0e-12) npccninr8(k) = 0.0

!            CCN01(I,K) = max(ccn_diag(1)*1e-6, 0.0)
!            CCN04(I,K) = max(ccn_diag(2)*1e-6, 0.0)
!            CCN1 (I,K) = max(ccn_diag(3)*1e-6, 0.0)



           else
             ccn_diag(:)    = 0.0
             smaxliq(K)     = 0.0
             swparc(K)      = 0.0
             smaxicer8(K)   = 0.0
             nheticer8(K)   = 0.0
             sc_icer8(K)    = 2.0
!            sc_icer8(K)    = 1.0
             naair8(K)      = 0.0
             npccninr8(K)   = 0.0
             nlimicer8(K)   = 0.0
             nhet_immr8(K)  = 0.0
             dnhet_immr8(K) = 0.0
             nhet_depr8(K)  = 0.0
             nhet_dhfr8(K)  = 0.0
             dust_immr8(K)  = 0.0
             dust_depr8(K)  = 0.0
             dust_dhfr8(K)  = 0.0

           end if

!          SMAXL(I,k)      = smaxliq(k)   * 100.0
!          SMAXI(I,k)      = smaxicer8(k) * 100.0
           NHET_NUC(I,k)   = nheticer8(k) * 1e-6
           NLIM_NUC(I,k)   = nlimicer8(k) * 1e-6
           SC_ICE(I,k)     = min(max(sc_icer8(k),1.0),2.0)
!          SC_ICE(I,k)     = min(max(sc_icer8(k),1.0),1.2)
!          if(temp(i,k) < T_ICE_ALL) SC_ICE(i,k) = max(SC_ICE(I,k), 1.2)
!          if(temp(i,k) < T_ICE_ALL) SC_ICE(i,k) = max(SC_ICE(I,k), 1.5)
!          if(temp(i,k) > T_ICE_ALL) SC_ICE(i,k) = 1.0
!          if(temp(i,k) > TICE)      SC_ICE(i,k) = rhc(i,k)
!
           if(temp(i,k) < T_ICE_ALL) then
!            SC_ICE(i,k) = max(SC_ICE(I,k), 1.2)
             SC_ICE(i,k) = max(SC_ICE(I,k), 1.5)
           elseif(temp(i,k) > TICE) then
             SC_ICE(i,k) = rhc(i,k)
           else
!            SC_ICE(i,k) = 1.0
!            tx1 = max(SC_ICE(I,k), 1.2)
             tx1 = max(SC_ICE(I,k), 1.5)
             SC_ICE(i,k) = ((tice-temp(i,k))*tx1 + (temp(i,k)-t_ice_all)*rhc(i,k)) &
                         * t_ice_denom
           endif
           if (.not. iccn) then
             CDNC_NUC(I,k) = npccninr8(k)
             INC_NUC (I,k) = naair8(k)
           endif
           NHET_IMM(I,k)   = max(nhet_immr8(k), 0.0)
           DNHET_IMM(I,k)  = max(dnhet_immr8(k), 0.0)
           NHET_DEP(I,k)   = nhet_depr8(k) * 1e-6
           NHET_DHF(I,k)   = nhet_dhfr8(k) * 1e-6
           DUST_IMM(I,k)   = max(dust_immr8(k), 0.0)*1e-6
           DUST_DEP(I,k)   = max(dust_depr8(k), 0.0)*1e-6
           DUST_DHF(I,k)   = max(dust_dhfr8(k), 0.0)*1e-6
           WSUB (I,k)      = wparc_ls(k)   + swparc(k)*0.8
           SIGW_GW (I,k)   = wparc_gw(k)   * wparc_gw(k)
           SIGW_CNV (I,k)  = wparc_cgw(k)  * wparc_cgw(k)
           SIGW_TURB (I,k) = wparc_turb(k) * wparc_turb(k)

         enddo                     ! end of K loop
       enddo                       ! end of first big I loop
!-------------------------------------------------------------------------------

!      SC_ICE=MIN(MAX(SC_ICE, 1.0), 2.0)
!      WHERE (TEMP .gt. T_ICE_ALL)
!      SC_ICE=1.0
!      END WHERE

!===========================End cloud particle nucleation=======================
!                           -----------------------------
!
!> -# Begin cloud macrophysics

!     do k=1,lm
!       do i=1,im
!         REV_CN_X(i,k) = 0.0
!         REV_LS_X(i,k) = 0.0
!         RSU_CN_X(i,k) = 0.0
!         RSU_LS_X(i,k) = 0.0
!         CFX(i,k) = INC_NUC(i,k) + NHET_IMM(i,k)
!       enddo
!     enddo
!     do k=0,lm
!       do i=1,im
!         PFI_CN_X(i,k) = 0.0
!         PFL_CN_X(i,k) = 0.0
!       enddo
!     enddo

!     if(lprnt) write(0,*)' skip_macro=',skip_macro

      if (.not. skip_macro) then

!      if (lprnt) write(0,*) ' in micro qicn2=',qicn(ipr,25),' kdt=',kdt&
!    &,' qils=',qils(ipr,25)
!       if(lprnt) write(0,*)' bef macro_cloud clcn=',clcn(ipr,:)
!       if(lprnt) write(0,*)' bef macro_cloud clls=',clls(ipr,:)

!       allocate(RHX_X(im,lm),     CFPDF_X(im,lm), VFALLSN_CN_X(im,lm), &
        allocate(                                                       &
!    &           QSNOW_CN(im,lm),  VFALLRN_CN_X(im,lm), QRAIN_CN(im,lm),&
!    &           QSNOW_CN(im,lm),                       QRAIN_CN(im,lm),&
!    &           REV_CN_X(im,lm),  RSU_CN_X(im,lm), DLPDF_X(im,lm),     &
!    &           DIPDF_X(im,lm),   ALPHT_X(im,lm),  PFRZ(im,lm),        &
     &                             ALPHT_X(im,lm),  PFRZ(im,lm))
!    &           ACLL_CN_X(im,lm), ACIL_CN_X(im,lm), DQRL_X(im,lm)
!    &           ACLL_CN_X(im,lm), ACIL_CN_X(im,lm), DQRL_X(im,lm),     &
!    &           DZET(im,lm))
!    &           DZET(im,lm),      qst3(im,lm))
!       allocate (PFI_CN_X(im,0:lm), PFL_CN_X(im,0:lm))

!       do L=LM,1,-1
!         do i=1,im
!           DZET(i,L)   = ZET(i,L)   - ZET(i,L+1)
!           tx1       = plo(i,l)*100.0
!           est3      = min(tx1, fpvs(temp(i,l)))
!           qst3(i,l) = min(eps*est3/max(tx1+epsm1*est3,1.0e-10),1.0)
!!          MASS(i,l) = (ple(i,l) - ple(i,l-1)) * (100.0/grav)
!         enddo
!       enddo
!       do k=1,lm
!         do i=1,im
!           REV_CN_X(i,k) = 0.0
!           RSU_CN_X(i,k) = 0.0
!         enddo
!       enddo
!       do k=0,lm
!         do i=1,im
!           PFI_CN_X(i,k) = 0.0
!           PFL_CN_X(i,k) = 0.0
!         enddo
!       enddo

!       call macro_cloud (IM, LM, DT_MOIST, PLO, PLE, PK, FRLAND,       &
!       call macro_cloud (IM, LM, DT_MOIST, PLO, PLE,     FRLAND,       &
!>  - Call macro_cloud() for cloud macrophysics
        call macro_cloud (IM, LM, DT_MOIST, alf_fac, PLO, PLE,          &
     &                             CNV_DQLDT,                           &
!    &                    CNV_MFD, CNV_DQLDT,                           &
!    &                    CNV_MFD, CNV_DQLDT, CNV_PRC3, CNV_UPDF,       &
!    &                    U1, V1, temp, Q1, QLLS,  QLCN, QILS, QICN,    &
     &                            temp, Q1, QLLS,  QLCN, QILS, QICN,    &
!    &                    U1, V1, TH1, Q1, QLLS,  QLCN, QILS, QICN,     &
     &                    CLCN, CLLS,                                   &
!    &                    CLCN, CLLS, CN_PRC2, CN_ARFX, CN_SNR,         &
     &                    CLOUDPARAMS, SCLMFDFR,                        &
!    &                    CLOUDPARAMS, SCLMFDFR, QST3, DZET, QDDF3,     &
!    &                    RHX_X, REV_CN_X, RSU_CN_X,                    &
!    &                    ACLL_CN_X, ACIL_CN_X, PFL_CN_X,               &
!    &                    PFI_CN_X, DLPDF_X, DIPDF_X,                   &
!    &                    ALPHT_X,  CFPDF_X, DQRL_X, VFALLSN_CN_X,      &
!    &                    VFALLRN_CN_X, CNV_FICE, CNV_NDROP, CNV_NICE,  &
     &                    ALPHT_X,                                      &
     &                                  CNV_FICE, CNV_NDROP, CNV_NICE,  &
     &                    SC_ICE,   NCPL,     NCPI,  PFRZ,              &
     &                              lprnt, ipr, rhc, pdfflag, qc_min)
!    &                    QRAIN_CN, QSNOW_CN, KCBL,  lprnt, ipr, rhc)


!       if (lprnt) write(0,*) ' in micro qicn3=',qicn(ipr,25)
!       if(lprnt) write(0,*)' aft macro_cloud clcn=',clcn(ipr,:)
!       if(lprnt) write(0,*)' aft macro_cloud clls=',clls(ipr,:)
!       if(lprnt) write(0,*)' aft macro_cloud q1=',q1(ipr,:)
!       if(lprnt) write(0,*)' aft macro_cloud qils=',qils(ipr,:)

        do k=1,lm
          do i=1,im
            if (CNV_MFD(i,k) > 1.0e-6) then
              tx1 = 1.0 / CNV_MFD(i,k)
              CNV_NDROP(i,k) = CNV_NDROP(i,k) * tx1
              CNV_NICE(i,k)  = CNV_NICE(i,k)  * tx1
            else
              CNV_NDROP(i,k) = 0.0
              CNV_NICE(i,k)  = 0.0
            endif
!           temp(i,k)   = th1(i,k) * PK(i,k)
            RAD_CF(i,k) = min(CLLS(i,k)+CLCN(i,k), 1.0)
!
            if (.not. iccn) then
              if (PFRZ(i,k) > 0.0) then
                INC_NUC(i,k)  = INC_NUC(i,k)  * PFRZ(i,k)
                NHET_NUC(i,k) = NHET_NUC(i,k) * PFRZ(i,k)
              else
                INC_NUC(i,k)  = 0.0
                NHET_NUC(i,k) = 0.0
              endif
            endif

          enddo
        enddo


!make sure QI , NI stay within T limits
!        call meltfrz_inst(IM, LM, TEMP, QLLS, QLCN, QILS, QICN, NCPL, NCPI)

!============ a little treatment of cloud before micorphysics
!        call update_cld(im,lm,DT_MOIST, ALPHT_X, qc_min                &
!    &,                 pdfflag,       PLO , Q1, QLLS                   &
!    &,                 QLCN, QILS,    QICN,   TEMP                     &
!    &,                 CLLS, CLCN,    SC_ICE, NCPI                     &
!    &,                 NCPL)
!!   &,                 NCPL, INC_NUC)
!============ Put cloud fraction back in contact with the PDF (Barahona et al., GMD, 2014)============

!make sure QI , NI stay within T limits
!>  - Call meltfrz_inst() to calculate instantaneous freezing or condensate
         call meltfrz_inst(IM, LM, TEMP, QLLS, QLCN, QILS, QICN, NCPL, NCPI)


!       deallocate(RHX_X,    CFPDF_X,      VFALLSN_CN_X,                &
        deallocate(                                                     &
!    &             QSNOW_CN, VFALLRN_CN_X, QRAIN_CN, REV_CN_X, RSU_CN_X,&
!    &             QSNOW_CN,               QRAIN_CN,                    &
     &                               PFRZ)
!    &             DLPDF_X, DIPDF_X, PFRZ, ACLL_CN_X, ACIL_CN_X, DQRL_X,&
!    &             PFI_CN_X, PFL_CN_X)
!    &             PFI_CN_X, PFL_CN_X,  DZET, qst3, qddf3)

      else
!       do i=1,im
!         CN_PRC2(i) = 0.0
!         CN_SNR(i)  = 0.0
!       enddo


      endif    ! .not. skip_macro


!===========================End of  Cloud Macrophysics ========================
!                           --------------------------
!


!TVQX1    = SUM( (  Q1 +  QLCN +  QICN )*DM, 3)

      do k=1,lm
        do i=1,im
          QCNTOT      = QLCN(i,k)   + QICN(i,k)
          QL_TOT(i,k) = QLCN(i,k)   + QLLS(i,k)
          QI_TOT(i,k) = QICN(i,k)   + QILS(i,k)
!         Anning if negative, borrow water and ice from vapor 11/23/2016
          if (QL_TOT(i,k) < 0.0) then
            Q1(i,k)     = Q1(i,k)   + QL_TOT(i,k)
            TEMP(i,k)   = TEMP(i,k) - lvbcp*QL_TOT(i,k)
            QL_TOT(i,k) = 0.0
          endif
          if (QI_TOT(i,k) < 0.0) then
            Q1(i,k)     = Q1(i,k)   + QI_TOT(i,k)
            TEMP(i,k)   = TEMP(i,k) - lsbcp*QI_TOT(i,k)
            QI_TOT(i,k) = 0.0
          endif
          QTOT = QL_TOT(i,k) + QI_TOT(i,k)
          if (QTOT > 0.0) then
            FQA(i,k) = min(max(QCNTOT / QTOT, 0.0), 1.0)
          else
            FQA(i,k) = 0.0
          endif
        enddo
      enddo

!=============================================================================================
!===========================Two-moment stratiform microphysics ===============================
!===========This is the implementation of the Morrison and Gettelman (2008) microphysics =====
!=============================================================================================
!> -# Two-moment stratiform microphysics: this is the implementation of the Morrison and 
!! Gettelman (2008) microphysics \cite Morrison_2008

      do I=1,IM
        LS_SNR(i)  = 0.0
        LS_PRC2(i) = 0.0

        nbincontactdust = 1

        do l=1,10
          do k=1,lm
            naconr8(k,l) = 0.0
            rndstr8(k,l) = 2.0e-7
          enddo
        enddo
        do k=1,lm
          npccninr8(k) = 0.0
          naair8(k)    = 0.0
          omegr8(k)    = 0.0

!         tx1 = MIN(CLLS(I,k) + CLCN(I,k), 0.99)
          tx1 = MIN(CLLS(I,k) + CLCN(I,k), 1.00)
          if (tx1 > 0.0) then
            cldfr8(k) = min(max(tx1, 0.00001), 1.0)
          else
            cldfr8(k) = 0.0
          endif

          if (temp(i,k) > tice) then
            liqcldfr8(k) = cldfr8(k)
            icecldfr8(k) = 0.0
          elseif (temp(i,k) <= t_ice_all) then
            liqcldfr8(k) = 0.0
            icecldfr8(k) = cldfr8(k)
          else
            icecldfr8(k) = cldfr8(k) * (tice - temp(i,k))/(tice-t_ice_all)
            liqcldfr8(k) = cldfr8(k) - icecldfr8(k)
          endif


          cldor8(k) = cldfr8(k)
          ter8(k)   = TEMP(I,k)
          qvr8(k)   = Q1(I,k)

          qcr8(k)   = QL_TOT(I,k)
          qir8(k)   = QI_TOT(I,k)
          ncr8(k)   = MAX(NCPL(I,k), 0.0)
          nir8(k)   = MAX(NCPI(I,k), 0.0)
          qrr8(k)   = rnw(I,k)
          qsr8(k)   = snw(I,k)
          qgr8(k)   = qgl(I,k)
          nrr8(k)   = MAX(NCPR(I,k), 0.0)
          nsr8(k)   = MAX(NCPS(I,k), 0.0)
          ngr8(k)   = MAX(ncgl(I,k), 0.0)


          naair8(k)    = INC_NUC(I,k)
          npccninr8(k) = CDNC_NUC(I,k)

          if (cldfr8(k) >= 0.001) then
            nimmr8(k) = min(DNHET_IMM(I,k),ncr8(k)/(cldfr8(k)*DT_MOIST))
          else
            nimmr8(k) = 0.0
          endif


!         if(aero_in) then
            AeroAux = AeroProps(I, K)
!         else
!           call init_Aer(AeroAux)
!         end if
!>  - Call getinsubset() to extract dust properties 
          call getINsubset(1, AeroAux, AeroAux_b)
          naux = AeroAux_b%nmods
          if (nbincontactdust < naux) then
            nbincontactdust = naux
          endif
          naconr8(K, 1:naux) = AeroAux_b%num(1:naux)
          rndstr8(K, 1:naux) = AeroAux_b%dpg(1:naux) * 0.5

! The following moved inside of if(fprcp <= 0) then loop
! Get black carbon properties for contact ice nucleation
!         call getINsubset(2, AeroAux, AeroAux_b)
!         nsootr8 (K)  = sum(AeroAux_b%num)
!         naux         = AeroAux_b%nmods
!         rnsootr8 (K) = sum(AeroAux_b%dpg(1:naux))/naux

          pdelr8(k)  = (PLE(I,k) - PLE(I,k-1)) * 100.0
          rpdelr8(k) = 1. / pdelr8(k)
          plevr8(k)  = 100. * PLO(I,k)
          zmr8(k)    = ZLO(I,k)
          ficer8(k)  = qir8(k) / (qcr8(k)+qir8(k) + 1.e-10)
          omegr8(k)  = WSUB(I,k)
!         alphar8(k) = max(alpht_x(i,k)/maxval(alpht_x(i,:))*8.,0.5)
!         alphar8(k) = qcvar2
          rhr8(k)    = rhc(i,k)

        END DO
        do k=1,lm+1
          pintr8(k)  = PLE(I,k-1) * 100.0
          kkvhr8(k)  = KH(I,k-1)
        END DO

        lev_sed_strt = 0
        tx1 = 1.0 / pintr8(lm+1)
        do k=1,lm
          if (plevr8(k)*tx1 < sig_sed_strt) then
            lev_sed_strt(1) = k
          endif
        enddo
        lev_sed_strt(1) = max(lm/6, min(lm/3,lev_sed_strt(1)))
!       if (kdt == 1)   &
!       write(0,*)' lev_sed_strt=',lev_sed_strt,' plevr8=',plevr8(lev_sed_strt), &
!                 ' pintr8=',pintr8(lm+1),' sig_sed_strt=',sig_sed_strt
!
!       do k=1,lm
!         if (cldfr8(k) <= 0.2 ) then
!           alphar8(k) = 0.5
!         elseif (cldfr8(k) <= 0.999) then
!!          tx1 = 0.0284 * exp(4.4*cldfr8(k))
!!          alphar8(k) = tx1 / (cldfr8(k) - tx1*(one-cldfr8(k)))
!!          alphar8(k) = 0.5 + (7.5/0.799)*(cldfr8(k)-0.2)
!           alphar8(k) = 0.5 + (7.5/0.799)*(cldfr8(k)-0.2)
!         else
!           alphar8(k) = 8.0
!         endif
!         alphar8(k) = min(8.0, max(alphar8(k), 0.5))
!       enddo

        kbmin = KCBL(I)

!!!Call to MG microphysics. Lives in cldwat2m_micro.f
!  ttendr8, qtendr8,cwtendr8, not used so far Anning noted August 2015

        if (fprcp <= 0) then  ! if fprcp = -1, then Anning's code for MG2 will be used
                              ! if fprcp = 0,  then MG1 is used

! Get black carbon properties for contact ice nucleation
          do k=1,lm
            call getINsubset(2, AeroAux, AeroAux_b)
            nsootr8 (K)  = sum(AeroAux_b%num)
            naux         = AeroAux_b%nmods
            rnsootr8 (K) = sum(AeroAux_b%dpg(1:naux))/naux
          enddo

          call mmicro_pcond ( ncolmicro,  ncolmicro,                                 &
     &                        dt_r8,      ter8, ttendr8,                             &
     &                        ncolmicro,  LM , qvr8,                                 &
     &                        qtendr8,    cwtendr8,   qcr8,     qir8,    ncr8, nir8, &
     &                        abs(fprcp), qrr8,       qsr8,     nrr8,    nsr8,       &
     &                        plevr8,     pdelr8,     cldfr8,   liqcldfr8,           &
     &                        icecldfr8,  cldor8,     pintr8,                        &
     &                        rpdelr8,    zmr8,       rate1ord_cw2pr,                &
     &                        naair8,     npccninr8,                                 &
     &                        rndstr8,    naconr8,    rhdfdar8, rhu00r8, ficer8,     &
     &                        tlatr8,     qvlatr8,    qctendr8,                      &
     &                        qitendr8,   nctendr8,   nitendr8, effcr8,              &
     &                        effc_fnr8,  effir8,     prectr8,  precir8,             &
     &                        nevaprr8,   evapsnowr8,                                &
     &                        prainr8,    prodsnowr8, cmeoutr8,                      &
     &                        deffir8,    pgamradr8,  lamcradr8,                     &
     &                        qsoutr8,    dsoutr8,    qroutr8,  droutr8,             &
     &                        qcsevapr8,  qisevapr8,  qvresr8,                       &
     &                        cmeioutr8,  vtrmcr8,    vtrmir8,                       &
     &                        qcsedtenr8, qisedtenr8, praor8,   prcor8,              &
     &                        mnucccor8,  mnucctor8,  msacwior8,                     &
     &                        psacwsor8,  bergsor8,   bergor8,                       &
     &                        meltor8,    homoor8,    qcresor8, prcior8,             &
     &                        praior8,    qiresor8,   mnuccror8,                     &
     &                        pracsor8,   meltsdtr8,  frzrdtr8, ncalr8,              &
     &                        ncair8,     mnuccdor8,                                 &
     &                        nnucctor8,  nsoutr8,    nroutr8,                       &
     &                        ncnstr8,    ninstr8,    nimmr8,   disp_liu,            &
     &                        nsootr8,    rnsootr8,   ui_scale, dcrit,               &
     &                        nnuccdor8,  nnucccor8,                                 &
     &                        nsacwior8,  nsubior8,   nprcior8,                      &
     &                        npraior8,   npccnor8,   npsacwsor8,                    &
     &                        nsubcor8,   npraor8,    nprc1or8,                      &
     &                        tlatauxr8,  nbincontactdust,                           &
     &                        lprnt,      xlat(i),    xlon(i),  rhr8)

!         if (lprint) write(0,*)' prectr8=',prectr8(1), &
!    &     ' precir8=',precir8(1)

          LS_PRC2(I) = max(1000.*(prectr8(1)-precir8(1)), 0.0)
          LS_SNR(I)  = max(1000.*precir8(1), 0.0)


          do k=1,lm
            QL_TOT(I,k) = QL_TOT(I,k) + qctendr8(k)*DT_R8
            QI_TOT(I,k) = QI_TOT(I,k) + qitendr8(k)*DT_R8
            Q1(I,k)     = Q1(I,k)     + qvlatr8(k)*DT_R8
!     if(lprnt .and. i == ipr) write(0,*)' k=',k,' q1aftm=',q1(i,k)     &
!    &,' qvlatr8=',qvlatr8(k)
            TEMP(I,k)   = TEMP(I,k)   + tlatr8(k)*DT_R8*onebcp

            NCPL(I,k)   = MAX(NCPL(I,k)   + nctendr8(k) * DT_R8, 0.0)
            NCPI(I,k)   = MAX(NCPI(I,k)   + nitendr8(k) * DT_R8, 0.0)
            rnw(I,k)    = qrr8(k)
            snw(I,k)    = qsr8(k)
            NCPR(I,k)   = nrr8(k)
            NCPS(I,k)   = nsr8(k)

            CLDREFFL(I,k) = min(max(effcr8(k), 10.), 150.)
            CLDREFFI(I,k) = min(max(effir8(k), 20.), 150.)
            CLDREFFR(I,k) = max(droutr8(k)*0.5*1.e6, 150.)
            CLDREFFS(I,k) = max(0.192*dsoutr8(k)*0.5*1.e6, 250.)

          enddo       ! K loop

        elseif (fprcp == 1) then           ! Call MG2
!                                            --------
!     if (lprnt .and. i == ipr) then
!       write(0,*)' bef micro_mg_tend ter8= ', ter8(:)
!       write(0,*)' bef micro_mg_tend qvr8= ', qvr8(:),'dt_r8=',dt_r8
!       write(0,*)' bef micro_mg_tend rhr8= ', rhr8(:)
!     endif
          lprint = lprnt .and. i == ipr
          ltrue = any(qcr8 >= qsmall) .or. any(qir8 >= qsmall)          &
             .or. any(qsr8 >= qsmall) .or. any(qrr8 >= qsmall)
          if (ltrue) then
            alphar8(:) = qcvar2

!           if(lprint) then
!             write(0,*)' calling micro_mg_tend2_0 qcvar2=',qcvar2
!             write(0,*)' qcr8=',qcr8(:)
!             write(0,*)' ncr8=',ncr8(:)
!             write(0,*)' npccninr8=',npccninr8(:)
!             write(0,*)' plevr8=',plevr8(:)
!             write(0,*)' ter8=',ter8(:)
!           endif

            call micro_mg_tend2_0 (                                     &
     &         ncolmicro,          lm,                 dt_r8,           &
     &         ter8,                         qvr8,                      &
     &         qcr8,                         qir8,                      &
     &         ncr8,                         nir8,                      &
     &         qrr8,                         qsr8,                      &
     &         nrr8,                         nsr8,                      &
     &         alphar8,                      1.,                        &
     &         plevr8,                       pdelr8,                    &
!    &         cldfr8,  liqcldfr8,      icecldfr8,     rhc,             &
     &         cldfr8,  liqcldfr8,      icecldfr8,     rhr8,            &
     &         qcsinksum_rate1ord,                                      &
     &         naair8,                       npccninr8,                 &
     &         rndstr8,                      naconr8,                   &
     &         tlatr8,                       qvlatr8,                   &
     &         qctendr8,                     qitendr8,                  &
     &         nctendr8,                     nitendr8,                  &
     &         qrtend,                       qstend,                    &
     &         nrtend,                       nstend,                    &
     &         effcr8,             effc_fnr8,          effir8,          &
     &         sadice,                       sadsnow,                   &
     &         prectr8,                      precir8,                   &
     &         nevaprr8,                     evapsnowr8,                &
     &         am_evp_st,                                               &
     &         prainr8,                      prodsnowr8,                &
     &         cmeoutr8,                     deffir8,                   &
     &         pgamradr8,                    lamcradr8,                 &
     &         qsoutr8,                      dsoutr8,                   &
     &         lflx,               iflx,                                &
     &         rflx,               sflx,               qroutr8,         &
     &         reff_rain,                    reff_snow,                 &
     &         qcsevapr8,          qisevapr8,          qvresr8,         &
     &         cmeioutr8,          vtrmcr8,            vtrmir8,         &
     &         umr,                          ums,                       &
     &         qcsedtenr8,                   qisedtenr8,                &
     &         qrsedten,                     qssedten,                  &
     &         praor8,                       prcor8,                    &
     &         mnucccor8,          mnucctor8,          msacwior8,       &
     &         psacwsor8,          bergsor8,           bergor8,         &
     &         meltor8,                      homoor8,                   &
     &         qcresor8,           prcior8,            praior8,         &
     &         qiresor8,           mnuccror8,          pracsor8,        &
     &         meltsdtr8,          frzrdtr8,           mnuccdor8,       &
     &         nroutr8,                      nsoutr8,                   &
     &         refl,               arefl,              areflz,          &
     &         frefl,              csrfl,              acsrfl,          &
     &         fcsrfl,                       rercld,                    &
     &         ncair8,                       ncalr8,                    &
     &         qrout2,                       qsout2,                    &
     &         nrout2,                       nsout2,                    &
     &         drout2,                       dsout2,                    &
     &         freqs,                        freqr,                     &
     &         nfice,                        qcrat,                     &
     &         prer_evap, xlat(i), xlon(i), lprint, iccn, aero_in,      &
     &         lev_sed_strt)
!
            LS_PRC2(I) = max(1000.*(prectr8(1)-precir8(1)), 0.0)
            LS_SNR(I)  = max(1000.*precir8(1), 0.0)
            do k=1,lm
              QL_TOT(I,k)   = QL_TOT(I,k) + qctendr8(k)*DT_R8
              QI_TOT(I,k)   = QI_TOT(I,k) + qitendr8(k)*DT_R8
              Q1(I,k)       = Q1(I,k)     + qvlatr8(k)*DT_R8
              TEMP(I,k)     = TEMP(I,k)   + tlatr8(k)*DT_R8*onebcp
              rnw(I,k)      = rnw(I,k)    + qrtend(k)*dt_r8
              snw(I,k)      = snw(I,k)    + qstend(k)*dt_r8

              NCPL(I,k)     = MAX(NCPL(I,k) + nctendr8(k)*DT_R8, 0.0)
              NCPI(I,k)     = MAX(NCPI(I,k) + nitendr8(k)*DT_R8, 0.0)
              NCPR(I,k)     = max(NCPR(I,k) + nrtend(k)*dt_r8,   0.0)
              NCPS(I,k)     = max(NCPS(I,k) + nstend(k)*dt_r8,   0.0)

              CLDREFFL(I,k) = min(max(effcr8(k), 10.),150.)
              CLDREFFI(I,k) = min(max(effir8(k), 20.),150.)
              CLDREFFR(I,k) = max(reff_rain(k),150.)
              CLDREFFS(I,k) = max(reff_snow(k),250.)
            enddo       ! K loop
!     if (lprint) then
!       write(0,*)' aft micro_mg_tend temp= ', temp(i,:)
!       write(0,*)' aft micro_mg_tend q1= ', q1(i,:)
!       write(0,*)' aft micro_mg_tend LS_PRC2= ', LS_PRC2(i),' ls_snr=',ls_snr(i)
!     endif
          else
            LS_PRC2(I) = 0.
            LS_SNR(I)  = 0.
            do k=1,lm
              CLDREFFL(I,k) = 10.
              CLDREFFI(I,k) = 50.
              CLDREFFR(I,k) = 1000.
              CLDREFFS(I,k) = 250.
            enddo       ! K loop
          endif
!
        else                               ! Call MG3
!                                            --------
          ltrue = any(qcr8 >= qsmall) .or. any(qir8 >= qsmall)          &
             .or. any(qsr8 >= qsmall) .or. any(qrr8 >= qsmall)          &
             .or. any(qgr8 >= qsmall)
          lprint = lprnt .and. i == ipr
          if (ltrue) then
            alphar8(:) = qcvar3

!           if(lprint) then
!             write(0,*)' calling micro_mg_tend3_0 qcvar3=',qcvar3,' i=',i
!             write(0,*)' qcr8=',qcr8(:)
!             write(0,*)' ncr8=',ncr8(:)
!             write(0,*)' npccninr8=',npccninr8(:)
!             write(0,*)' plevr8=',plevr8(:)
!             write(0,*)' ter8=',ter8(:)
!           endif
!>  - Call micro_mg3_0::micro_mg_tend(), which is the main microphysics routine to
!! calculate microphysical processes and other utilities.
            call micro_mg_tend3_0 (                                     &
     &         ncolmicro,          lm,                 dt_r8,           &
     &         ter8,                         qvr8,                      &
     &         qcr8,                         qir8,                      &
     &         ncr8,                         nir8,                      &
     &         qrr8,                         qsr8,                      &
     &         nrr8,                         nsr8,                      &
     &         qgr8,                         ngr8,                      &
     &         alphar8,                      1.,                        &
     &         plevr8,                       pdelr8,                    &
!    &         cldfr8,  liqcldfr8,      icecldfr8,     rhc,             &
     &         cldfr8,  liqcldfr8,      icecldfr8,     rhr8,            &
     &         qcsinksum_rate1ord,                                      &
     &         naair8,                       npccninr8,                 &
     &         rndstr8,                      naconr8,                   &
     &         tlatr8,                       qvlatr8,                   &
     &         qctendr8,                     qitendr8,                  &
     &         nctendr8,                     nitendr8,                  &
     &         qrtend,                       qstend,                    &
     &         nrtend,                       nstend,                    &
!
     &         qgtend,                       ngtend,                    &
!
     &         effcr8,             effc_fnr8,          effir8,          &
     &         sadice,                       sadsnow,                   &
     &         prectr8,                      precir8,                   &
     &         nevaprr8,                     evapsnowr8,                &
     &         am_evp_st,                                               &
     &         prainr8,                      prodsnowr8,                &
     &         cmeoutr8,                     deffir8,                   &
     &         pgamradr8,                    lamcradr8,                 &
     &         qsoutr8,                      dsoutr8,                   &
!
     &         qgoutr8,            ngoutr8,  dgoutr8,                   &
!
     &         lflx,               iflx,     gflx,                      &
!
     &         rflx,               sflx,               qroutr8,         &
!
     &         reff_rain,          reff_snow, reff_grau,                &
!
     &         qcsevapr8,          qisevapr8,          qvresr8,         &
     &         cmeioutr8,          vtrmcr8,            vtrmir8,         &
     &         umr,                          ums,                       &
!
     &         umg,                          qgsedtenr8,                &
!
     &         qcsedtenr8,                   qisedtenr8,                &
     &         qrsedten,                     qssedten,                  &
     &         praor8,                       prcor8,                    &
     &         mnucccor8,          mnucctor8,          msacwior8,       &
     &         psacwsor8,          bergsor8,           bergor8,         &
     &         meltor8,                      homoor8,                   &
     &         qcresor8,           prcior8,            praior8,         &
!
     &         qiresor8,           mnuccror8, mnuccrior8, pracsor8,     &
!
     &         meltsdtr8,          frzrdtr8,           mnuccdor8,       &
!
     &         pracgr8,            psacwgr8,           pgsacwr8,        &
     &         pgracsr8,           prdgr8,                              &
     &         qmultgr8,           qmultrgr8,          psacrr8,         &
     &         npracgr8,           nscngr8,            ngracsr8,        &
     &         nmultgr8,           nmultrgr8,          npsacwgr8,       &
!
     &         nroutr8,                      nsoutr8,                   &
     &         refl,               arefl,              areflz,          &
     &         frefl,              csrfl,              acsrfl,          &
     &         fcsrfl,                       rercld,                    &
     &         ncair8,                       ncalr8,                    &
     &         qrout2,                       qsout2,                    &
     &         nrout2,                       nsout2,                    &
     &         drout2,                       dsout2,                    &
!
     &         qgout2,             ngout2,   dgout2, freqg,             &
     &         freqs,                        freqr,                     &
     &         nfice,                        qcrat,                     &
     &         prer_evap, xlat(i), xlon(i), lprint, iccn, aero_in,      &
     &         lev_sed_strt)

            LS_PRC2(I) = max(1000.*(prectr8(1)-precir8(1)), 0.0)
            LS_SNR(I)  = max(1000.*precir8(1), 0.0)
            do k=1,lm
              QL_TOT(I,k)   = QL_TOT(I,k) + qctendr8(k)*DT_R8
              QI_TOT(I,k)   = QI_TOT(I,k) + qitendr8(k)*DT_R8
              Q1(I,k)       = Q1(I,k)     + qvlatr8(k)*DT_R8
              TEMP(I,k)     = TEMP(I,k)   + tlatr8(k)*DT_R8*onebcp
              rnw(I,k)      = rnw(I,k)    + qrtend(k)*dt_r8
              snw(I,k)      = snw(I,k)    + qstend(k)*dt_r8
              qgl(I,k)      = qgl(I,k)    + qgtend(k)*dt_r8

              NCPL(I,k)     = MAX(NCPL(I,k) + nctendr8(k)*DT_R8, 0.0)
              NCPI(I,k)     = MAX(NCPI(I,k) + nitendr8(k)*DT_R8, 0.0)
              NCPR(I,k)     = max(NCPR(I,k) + nrtend(k)*dt_r8,   0.0)
              NCPS(I,k)     = max(NCPS(I,k) + nstend(k)*dt_r8,   0.0)
              NCGL(I,k)     = max(NCGL(I,k) + ngtend(k)*dt_r8,   0.0)

              CLDREFFL(I,k) = min(max(effcr8(k), 10.),150.)
              CLDREFFI(I,k) = min(max(effir8(k), 20.),150.)
              CLDREFFR(I,k) = max(reff_rain(k),150.)
              CLDREFFS(I,k) = max(reff_snow(k),250.)
              CLDREFFG(I,k) = max(reff_grau(k),250.)
            enddo       ! K loop
!     if (lprint) then
!       write(0,*)' aft micro_mg_tend temp= ', temp(i,:)
!       write(0,*)' aft micro_mg_tend q1= ', q1(i,:)
!       write(0,*)' aft micro_mg_tend LS_PRC2= ', LS_PRC2(i),' ls_snr=',ls_snr(i)
!     endif
          else
            LS_PRC2(I) = 0.
            LS_SNR(I)  = 0.
            do k=1,lm
              CLDREFFL(I,k) = 10.
              CLDREFFI(I,k) = 50.
              CLDREFFR(I,k) = 1000.
              CLDREFFS(I,k) = 250.
              CLDREFFG(I,k) = 250.
            enddo       ! K loop
          endif
        endif

      enddo        ! I loop
!============================================Finish 2-moment micro implementation===========================

!TVQX1    = SUM( (  Q1 +  QL_TOT + QI_TOT(1:im,:,:))*DM, 3) &


      if (skip_macro) then
        do k=1,lm
          do i=1,im
            CALL fix_up_clouds_2M(Q1(I,K),   TEMP(i,k), QLLS(I,K),      &
     &                            QILS(I,K), CLLS(I,K), QLCN(I,K),      &
     &                            QICN(I,K), CLCN(I,K), NCPL(I,K),      &
     &                            NCPI(I,K), qc_min)
            if (rnw(i,k) <= qc_min(1)) then
              ncpl(i,k) = 0.0
            elseif (ncpl(i,k) <= nmin) then ! make sure NL > 0 if Q >0
              ncpl(i,k) = max(rnw(i,k) / (fourb3 * PI *RL_cub*997.0), nmin)
            endif
            if (snw(i,k) <= qc_min(2)) then
              ncpl(i,k) = 0.0
            elseif (ncps(i,k) <= nmin) then
              ncps(i,k) = max(snw(i,k) / (fourb3 * PI *RL_cub*500.0), nmin)
            endif
            if (qgl(i,k) <= qc_min(2)) then
              ncgl(i,k) = 0.0
            elseif (ncgl(i,k) <= nmin) then
              ncgl(i,k) = max(qgl(i,k) / (fourb3 * PI *RL_cub*500.0), nmin)
            endif
          enddo
        enddo
      else
        do k=1,lm
          do i=1,im
            QLCN(i,k) = QL_TOT(i,k) * FQA(i,k)
            QLLS(i,k) = QL_TOT(i,k) - QLCN(i,k)
            QICN(i,k) = QI_TOT(i,k) * FQA(i,k)
            QILS(i,k) = QI_TOT(i,k) - QICN(i,k)
          enddo
        enddo

!>  - Call update_cld()
        call update_cld(im, lm,  DT_MOIST,   ALPHT_X, qc_min            &
     &,                 pdfflag, PLO,  Q1,   QLLS,    QLCN              &
     &,                 QILS,    QICN, TEMP, CLLS,    CLCN              &
     &,                 SC_ICE,  NCPI, NCPL)

!       if(lprnt) write(0,*)' aft update_cloud clls=',clls(ipr,:)

        do k=1,lm
          do i=1,im
            QL_TOT(I,K) = QLLS(I,K) + QLCN(I,K)
            QI_TOT(I,K) = QILS(I,K) + QICN(I,K)
!
            if (rnw(i,k) <= qc_min(1)) then
              ncpl(i,k) = 0.0
            elseif (ncpl(i,k) <= nmin) then ! make sure NL > 0 if Q >0
              ncpl(i,k) = max(rnw(i,k) / (fourb3 * PI *RL_cub*997.0), nmin)
            endif
            if (snw(i,k) <= qc_min(2)) then
              ncpl(i,k) = 0.0
            elseif (ncps(i,k) <= nmin) then
              ncps(i,k) = max(snw(i,k) / (fourb3 * PI *RL_cub*500.0), nmin)
            endif
            if (qgl(i,k) <= qc_min(2)) then
              ncgl(i,k) = 0.0
            elseif (ncgl(i,k) <= nmin) then
              ncgl(i,k) = max(qgl(i,k) / (fourb3 * PI *RL_cub*500.0), nmin)
            endif
          enddo
        enddo
        deallocate(CNV_MFD,CNV_FICE,CNV_NDROP,CNV_NICE)
!       deallocate(CNV_MFD,CNV_PRC3,CNV_FICE,CNV_NDROP,CNV_NICE)
      endif

!     do I=1,IM
!       TPREC(i) = CN_PRC2(i) + LS_PRC2(i) +  CN_SNR(i) + LS_SNR(i)
!     enddo

      do K= 1, LM
        do I=1,IM
          if (QI_TOT(i,k) <= 0.0) NCPI(i,k) = 0.0
          if (QL_TOT(i,k) <= 0.0) NCPL(i,k) = 0.0
        end do
      end do


!=============================================End Stratiform cloud processes==========================================
!======================================================================================================================
!===========================Clean stuff and send it to radiation ======================================================
!======================================================================================================================
! outputs
       if(flipv) then
         DO K=1, LM
           ll = lm-k+1
           DO I = 1,IM
             t_io(i,k)    = TEMP(i,ll)
             q_io(i,k)    = Q1(i,ll)
             ncpi_io(i,k) = NCPI(i,ll)
             ncpl_io(i,k) = NCPL(i,ll)
             rnw_io(i,k)  = rnw(i,ll)
             snw_io(i,k)  = snw(i,ll)
             qgl_io(i,k)  = qgl(i,ll)
             ncpr_io(i,k) = NCPR(i,ll)
             ncps_io(i,k) = NCPS(i,ll)
             ncgl_io(i,k) = NCGL(i,ll)
             lwm_o(i,k)   = QL_TOT(i,ll)
             qi_o(i,k)    = QI_TOT(i,ll)
           END DO
         END DO
         if (skip_macro) then
           DO K=1, LM
             ll = lm-k+1
             DO I = 1,IM
               CLLS_io(i,k) = max(0.0, min(CLLS(i,ll)+CLCN(i,ll),1.0))
             enddo
           enddo
         else
           DO K=1, LM
             ll = lm-k+1
             DO I = 1,IM
               CLLS_io(i,k) = CLLS(i,ll)
             enddo
           enddo
         endif
       else
         DO K=1, LM
           DO I = 1,IM
             t_io(i,k)    = TEMP(i,k)
             q_io(i,k)    = Q1(i,k)
             ncpi_io(i,k) = NCPI(i,k)
             ncpl_io(i,k) = NCPL(i,k)
             rnw_io(i,k)  = rnw(i,k)
             snw_io(i,k)  = snw(i,k)
             qgl_io(i,k)  = qgl(i,k)
             ncpr_io(i,k) = NCPR(i,k)
             ncps_io(i,k) = NCPS(i,k)
             ncgl_io(i,k) = NCGL(i,k)
             lwm_o(i,k)   = QL_TOT(i,k)
             qi_o(i,k)    = QI_TOT(i,k)
           END DO
         END DO
         if (skip_macro) then
           DO K=1, LM
             DO I = 1,IM
               CLLS_io(i,k) = max(0.0, min(CLLS(i,k)+CLCN(i,k),1.0))
             enddo
           enddo
         else
           DO K=1, LM
             DO I = 1,IM
               CLLS_io(i,k) = CLLS(i,k)
             enddo
           enddo
         endif
       endif       ! end of flipv if

       DO I = 1,IM
         tx1     = LS_PRC2(i) + LS_SNR(i)
         rn_o(i) = tx1 * dt_i * 0.001

         if (rn_o(i) < rainmin) then
           sr_o(i) = 0.
         else
           sr_o(i) = LS_SNR(i) / tx1
         endif
       ENDDO

       if (allocated(ALPHT_X)) deallocate (ALPHT_X)

!     if (lprnt) then
!       write(0,*)' rn_o=',rn_o(ipr),' ls_prc2=',ls_prc2(ipr),' ls_snr=',ls_snr(ipr)
!       write(0,*)' end micro_mg_tend t_io= ', t_io(ipr,:)
!       write(0,*)' end micro_mg_tend clls_io= ', clls_io(ipr,:)
!     endif
!      do k=1,lm
!        do i=1,im
!          dum(i,k) = clls_io(i,k)
!        enddo
!      enddo
!      do k=2,lm-1
!        do i=1,im
!          clls_io(i,k) = 0.25*dum(i,k-1) + 0.5*dum(i,k)+0.25*dum(i,k+1)
!        enddo
!      enddo
!      do i=1,im
!        clls_io(i,lm) = 0.5 * (dum(i,lm-1) + dum(i,lm))
!      enddo



!=======================================================================

       end subroutine m_micro_run
!> @}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!DONIF Calculate the Brunt_Vaisala frequency

!===============================================================================
!>\ingroup mg_driver
!> This subroutine computes profiles of background state quantities for 
!! the multiple gravity wave drag parameterization.
!!\section gw_prof_gen MG gw_prof General Algorithm
!> @{
       subroutine gw_prof (pcols, pver, ncol, t, pm, pi, rhoi, ni, ti,  &
                           nm, sph)
       use machine , only : kind_phys
       use physcons, grav => con_g, cp => con_cp, rgas => con_rd,       &
                     fv   => con_fvirt
       implicit none
!-----------------------------------------------------------------------
! Compute profiles of background state quantities for the multiple
! gravity wave drag parameterization.
!
! The parameterization is assumed to operate only where water vapor
! concentrations are negligible in determining the density.
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
       integer, intent(in) :: ncol, pcols, pver



       real(kind=kind_phys), intent(in) :: t(pcols,pver)
       real(kind=kind_phys), intent(in) :: pm(pcols,pver)
       real(kind=kind_phys), intent(in) :: pi(pcols,0:pver)
       real(kind=kind_phys), intent(in) :: sph(pcols,pver)

       real(kind=kind_phys), intent(out) :: rhoi(pcols,0:pver)
       real(kind=kind_phys), intent(out) :: ni(pcols,0:pver)
       real(kind=kind_phys), intent(out) :: ti(pcols,0:pver)
       real(kind=kind_phys), intent(out) :: nm(pcols,pver)

       real(kind=kind_phys), parameter :: r=rgas, cpair=cp, g=grav, &
                                          oneocp=1.0/cp, n2min=1.e-8

!---------------------------Local storage-------------------------------
       integer :: ix,kx

       real :: dtdp, n2

!-----------------------------------------------------------------------------
!> -# Determine the interface densities and Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------

! The top interface values are calculated assuming an isothermal atmosphere
! above the top level.
       kx = 0
       do ix = 1, ncol
         ti(ix,kx)   = t(ix,kx+1)
         rhoi(ix,kx) = pi(ix,kx) / (r*(ti(ix,kx)*(1.0+fv*sph(ix,kx+1))))
         ni(ix,kx)   = sqrt (g*g / (cpair*ti(ix,kx)))
       end do

! Interior points use centered differences
       do kx = 1, pver-1
         do ix = 1, ncol
           ti(ix,kx)   = 0.5 * (t(ix,kx) + t(ix,kx+1))
           rhoi(ix,kx) = pi(ix,kx) / (r*ti(ix,kx)*(1.0+0.5*fv*(sph(ix,kx)+sph(ix,kx+1))))
           dtdp = (t(ix,kx+1)-t(ix,kx)) / (pm(ix,kx+1)-pm(ix,kx))
           n2   = g*g/ti(ix,kx) * (oneocp - rhoi(ix,kx)*dtdp)
           ni(ix,kx) = sqrt (max (n2min, n2))
         end do
       end do

! Bottom interface uses bottom level temperature, density; next interface
! B-V frequency.
       kx = pver
       do ix = 1, ncol
         ti(ix,kx)   = t(ix,kx)
         rhoi(ix,kx) = pi(ix,kx) / (r*ti(ix,kx)*(1.0+fv*sph(ix,kx)))
         ni(ix,kx)   = ni(ix,kx-1)
       end do

!-----------------------------------------------------------------------------
!> -# Determine the midpoint Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------
       do kx=1,pver
         do ix=1,ncol
           nm(ix,kx) = 0.5 * (ni(ix,kx-1) + ni(ix,kx))
         end do
       end do

       return
       end subroutine gw_prof
!> @}

!>\ingroup mg_driver
!! This subroutine is to find cloud top based on cloud fraction.
      subroutine find_cldtop(ncol, pver, cf, kcldtop)
       implicit none

       integer, intent(in)  :: pver , ncol
       real,    intent(in)  :: cf(ncol,pver)
       integer, intent(out) :: kcldtop
       integer              :: kuppest, ibot, k
       real                 :: stab, cfcrit, cf00, cfp1


       ibot    = pver-1
       kcldtop = ibot+1
       kuppest = 20
       cfcrit  = 1e-2


       do k = kuppest , ibot
         cfp1 = cf(ncol, k+1)

         if ( ( cfp1 >= cfcrit ) ) then
           kcldtop = k +1
           exit
         end if
       end do

       if (kcldtop >= ibot) then
         kcldtop = pver
         return
       endif


      end subroutine find_cldtop
!> @}

end module m_micro
