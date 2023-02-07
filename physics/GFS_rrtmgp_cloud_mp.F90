!> \file GFS_rrtmgp_cloud_mp.F90
!!
!> \defgroup GFS_rrtmgp_cloud_mp GFS_rrtmgp_cloud_mp.F90
!!
!! \brief This module contains the interface for ALL cloud microphysics assumptions and 
!! the RRTMGP radiation scheme. Specific details below in subroutines.
!!
module GFS_rrtmgp_cloud_mp
  use machine,      only: kind_phys
  use radiation_tools,   only: check_error_msg
  use module_radiation_clouds, only: progcld_thompson
  use rrtmgp_lw_cloud_optics, only: &
       radliq_lwr => radliq_lwrLW, radliq_upr => radliq_uprLW,&
       radice_lwr => radice_lwrLW, radice_upr => radice_uprLW  
  use module_mp_thompson, only: calc_effectRad, Nt_c_l, Nt_c_o, re_qc_min, re_qc_max,  &
       re_qi_min, re_qi_max, re_qs_min, re_qs_max
  use module_mp_thompson_make_number_concentrations, only: make_IceNumber,             &
       make_DropletNumber, make_RainNumber
  
  real (kind_phys), parameter :: &
       cld_limit_lower = 0.001, &
       cld_limit_ovcst = 1.0 - 1.0e-8, &
       reliq_def  = 10.0 ,      & ! Default liq radius to 10 micron (used when effr_in=F)
       reice_def  = 50.0,       & ! Default ice radius to 50 micron (used when effr_in=F)
       rerain_def = 1000.0,     & ! Default rain radius to 1000 micron (used when effr_in=F)
       resnow_def = 250.0,      & ! Default snow radius to 250 micron (used when effr_in=F)
       reice_min  = 10.0,       & ! Minimum ice size allowed by GFDL MP scheme
       reice_max  = 150.0         ! Maximum ice size allowed by GFDL MP scheme  
  
  public GFS_rrtmgp_cloud_mp_init, GFS_rrtmgp_cloud_mp_run, GFS_rrtmgp_cloud_mp_finalize

contains  

!>\defgroup gfs_rrtmgp_cloud_mp_mod GFS RRTMGP Cloud MP Module
!! \section arg_table_GFS_rrtmgp_cloud_mp_run
!! \htmlinclude GFS_rrtmgp_cloud_mp_run_html
!!
!> \ingroup GFS_rrtmgp_cloud_mp
!!
!! Here the cloud-radiative properties (optical-path, particle-size and sometimes cloud-
!! fraction) are computed for cloud producing physics schemes (e.g GFDL-MP, Thompson-MP,
!! MYNN-EDMF-pbl, GF-convective, and SAMF-convective clouds).
!!
!! \section GFS_rrtmgp_cloud_mp_run
  subroutine GFS_rrtmgp_cloud_mp_run(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice,     &
       i_cldrain, i_cldsnow, i_cldgrpl, i_cldtot, i_cldliq_nc, i_cldice_nc, i_twa, kdt,  &
       imfdeepcnv, imfdeepcnv_gf, imfdeepcnv_samf, doSWrad, doLWrad, effr_in, lmfshal,   &
       ltaerosol,mraerosol, icloud, imp_physics, imp_physics_thompson, imp_physics_gfdl, &
       lgfdlmprad, do_mynnedmf, uni_cld, lmfdeep2, p_lev, p_lay, t_lay, qs_lay, q_lay,   &
       relhum, lsmask, xlon, xlat, dx, tv_lay, effrin_cldliq, effrin_cldice,             &
       effrin_cldrain, effrin_cldsnow, tracer, cnv_mixratio, cld_cnv_frac, qci_conv,     &
       deltaZ, deltaZc, deltaP, qc_mynn, qi_mynn, cld_pbl_frac, con_g, con_rd, con_eps,  &
       con_ttp, doGP_cldoptics_PADE, doGP_cldoptics_LUT, doGP_smearclds,                 &
       cld_frac, cld_lwp, cld_reliq,                                                     &
       cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp, cld_rerain, precip_frac,        &
       cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp, cld_cnv_reice, cld_pbl_lwp,              &
       cld_pbl_reliq, cld_pbl_iwp, cld_pbl_reice, lwp_ex, iwp_ex, lwp_fc, iwp_fc,        &
       cldfra2d, errmsg, errflg)
    implicit none

    ! Inputs   
    integer, intent(in)    :: &
         nCol,                      & ! Number of horizontal grid points
         nLev,                      & ! Number of vertical layers
         ncnd,                      & ! Number of cloud condensation types.
         nTracers,                  & ! Number of tracers from model. 
         i_cldliq,                  & ! Index into tracer array for cloud liquid. 
         i_cldice,                  & ! Index into tracer array for cloud ice.
         i_cldrain,                 & ! Index into tracer array for cloud rain.
         i_cldsnow,                 & ! Index into tracer array for cloud snow.
         i_cldgrpl,                 & ! Index into tracer array for cloud groupel.
         i_cldtot,                  & ! Index into tracer array for cloud total amount.
         i_cldliq_nc,               & !                             cloud liquid number concentration.
         i_cldice_nc,               & !                             cloud ice number concentration.
         i_twa,                     & !                             water friendly aerosol.
         imfdeepcnv,                & ! Choice of mass-flux deep convection scheme
         imfdeepcnv_gf,             & ! Flag for Grell-Freitas deep convection scheme
         imfdeepcnv_samf,           & ! Flag for scale awware mass flux convection scheme
         kdt,                       & ! Current forecast iteration
         imp_physics,               & ! Choice of microphysics scheme
         imp_physics_thompson,      & ! Choice of Thompson
         imp_physics_gfdl,          & ! Choice of GFDL
         icloud                       ! Control for cloud are fraction option
    logical, intent(in) :: &
         doSWrad,                   & ! Call SW radiation?
         doLWrad,                   & ! Call LW radiation?
         effr_in,                   & ! Provide hydrometeor radii from macrophysics?
         lmfshal,                   & ! Flag for mass-flux shallow convection scheme used by Xu-Randall
         ltaerosol,                 & ! Flag for aerosol option
         mraerosol,                 & ! Flag for aerosol option
         lgfdlmprad,                & ! Flag for GFDLMP radiation interaction
         do_mynnedmf,               & ! Flag to activate MYNN-EDMF 
         uni_cld,                   & ! Flag for unified cloud scheme
         lmfdeep2,                  & ! Flag for mass flux deep convection 
         doGP_cldoptics_LUT,        & ! Flag to do GP cloud-optics (LUTs)
         doGP_cldoptics_PADE,       & !                            (PADE approximation)
         doGP_smearclds               ! If true, add sgs clouds to gridmean clouds
    real(kind_phys), intent(in) :: &
         con_g,                     & ! Physical constant: gravitational constant
         con_rd,                    & ! Physical constant: gas-constant for dry air
         con_ttp,                   & ! Triple point temperature of water (K)  
         con_eps                      ! Physical constant: gas constant air / gas constant H2O
    real(kind_phys), dimension(:), intent(in) :: &
         lsmask,                    & ! Land/Sea mask
         xlon,                      & ! Longitude
         xlat,                      & ! Latitude 
         dx                           ! Characteristic grid lengthscale (m)
    real(kind_phys), dimension(:,:), intent(in) :: &         
         tv_lay,                    & ! Virtual temperature (K)
         t_lay,                     & ! Temperature (K)
         qs_lay,                    & ! Saturation vapor pressure (Pa)
         q_lay,                     & ! water-vapor mixing ratio (kg/kg)
         relhum,                    & ! Relative humidity
         p_lay,                     & ! Pressure at model-layers (Pa)
         cnv_mixratio,              & ! Convective cloud mixing-ratio (kg/kg)
         qci_conv,                  & ! Convective cloud condesate after rainout (kg/kg)
         deltaZ,                    & ! Layer-thickness (m)
         deltaZc,                   & ! Layer-thickness, from layer centers (m)
         deltaP,                    & ! Layer-thickness (Pa)
         qc_mynn,                   & !
         qi_mynn,                   & !
         cld_pbl_frac                 !
    real(kind_phys), dimension(:,:), intent(inout) :: &
         effrin_cldliq,             & ! Effective radius for stratiform liquid cloud-particles (microns)
         effrin_cldice,             & ! Effective radius for stratiform ice cloud-particles (microns)
         effrin_cldsnow               ! Effective radius for stratiform snow cloud-particles (microns)
    real(kind_phys), dimension(:,:), intent(in) :: &
         effrin_cldrain               ! Effective radius for stratiform rain cloud-particles (microns)
    real(kind_phys), dimension(:,:), intent(in) :: &
         p_lev                        ! Pressure at model-level interfaces (Pa)
    real(kind_phys), dimension(:,:,:),intent(in) :: &
         tracer                       ! Cloud condensate amount in layer by type ()

    ! Outputs
    real(kind_phys), dimension(:), intent(inout) :: &
         lwp_ex,                    & ! Total liquid water path from explicit microphysics
         iwp_ex,                    & ! Total ice    water path from explicit microphysics
         lwp_fc,                    & ! Total liquid water path from cloud fraction scheme
         iwp_fc                       ! Total ice    water path from cloud fraction scheme
    real(kind_phys), dimension(:), intent(out) :: &
         cldfra2d                     ! Instantaneous 2D (max-in-column) cloud fraction
    real(kind_phys), dimension(:,:),intent(inout) :: &
         cld_frac,                  & ! Cloud-fraction for   stratiform   clouds
         cld_lwp,                   & ! Water path for       stratiform   liquid cloud-particles
         cld_reliq,                 & ! Effective radius for stratiform   liquid cloud-particles
         cld_iwp,                   & ! Water path for       stratiform   ice    cloud-particles
         cld_reice,                 & ! Effective radius for stratiform   ice    cloud-particles
         cld_swp,                   & ! Water path for                    snow   hydrometeors
         cld_resnow,                & ! Effective radius for              snow   hydrometeors
         cld_rwp,                   & ! Water path for                    rain   hydrometeors
         cld_rerain,                & ! Effective radius for              rain   hydrometeors
         precip_frac,               & ! Precipitation fraction
         cld_cnv_frac,              & ! Cloud-fraction for   convective clouds
         cld_cnv_lwp,               & ! Water path for       convective   liquid cloud-particles
         cld_cnv_reliq,             & ! Effective radius for convective   liquid cloud-particles
         cld_cnv_iwp,               & ! Water path for       convective   ice    cloud-particles
         cld_cnv_reice,             & ! Effective radius for convective   ice    cloud-particles
         cld_pbl_lwp,               & ! Water path for       SGS PBL liquid cloud-particles
         cld_pbl_reliq,             & ! Effective radius for SGS PBL liquid cloud-particles
         cld_pbl_iwp,               & ! Water path for       SGS PBL ice    cloud-particles
         cld_pbl_reice                ! Effective radius for SGS PBL ice    cloud-particles
    character(len=*), intent(out) :: &
         errmsg                       ! Error message
    integer, intent(out) :: &  
         errflg                       ! Error flag

    ! Local
    integer :: iCol, iLay
    real(kind_phys) :: alpha0
    real(kind_phys), dimension(nCol,nLev) :: cldcov, cldtot, cldcnv

    if (.not. (doSWrad .or. doLWrad)) return

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! ###################################################################################
    ! GFDL Microphysics
    ! ("Implicit" SGS cloud-coupling to the radiation)
    ! ###################################################################################
    if (imp_physics == imp_physics_gfdl) then
       ! GFDL-Lin
       if (.not. lgfdlmprad) then
          errflg = 1
          errmsg = "ERROR: MP choice not available with RRTMGP"
          return
       ! GFDL-EMC
       else

          ! "cld_frac" is modified prior to include subgrid scale cloudiness, see 
          ! module_SGSCloud_RadPre.F90.
          do iLay = 1, nLev
             do iCol = 1, nCol
                ! 
                ! SGS clouds present, use cloud-fraction modified to include sgs clouds.
                !
                if ((imfdeepcnv==imfdeepcnv_gf .or. do_mynnedmf) .and. kdt>1) then
                   ! MYNN sub-grid cloud fraction.
                   if (do_mynnedmf) then
                      ! If rain/snow present, use GFDL MP cloud-fraction...
                      if (tracer(iCol,iLay,i_cldrain)>1.0e-7 .OR. tracer(iCol,iLay,i_cldsnow)>1.0e-7) then
                         cld_frac(iCol,iLay) = tracer(iCol,iLay,i_cldtot)
                      endif
                   ! GF sub-grid cloud fraction.
                   else
                      ! If no convective cloud condensate present, use GFDL MP cloud-fraction....
                      if (qci_conv(iCol,iLay) <= 0.) then
                         cld_frac(iCol,iLay) = tracer(iCol,iLay,i_cldtot)
                      endif
                   endif
                !
                ! No SGS clouds, use GFDL MP cloud-fraction...
                !
                else
                   cld_frac(iCol,iLay) = tracer(iCol,iLay,i_cldtot)
                endif
             enddo
          enddo

          call cloud_mp_uni(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice, i_cldrain,  &
               i_cldsnow, i_cldgrpl, i_cldtot,  effr_in, kdt, lsmask, p_lev, p_lay,     &
               t_lay, tv_lay, effrin_cldliq, effrin_cldice, effrin_cldsnow, tracer,     &
               con_g, con_rd, con_ttp, cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice,&
               cld_swp, cld_resnow, cld_rwp, cld_rerain, effrin_cldrain=effrin_cldrain)
       end if
    endif

    ! ###################################################################################
    ! Thompson Microphysics
    ! ("Explicit" SGS cloud-coupling to the radiation)
    ! ###################################################################################
    if (imp_physics == imp_physics_thompson) then

       ! MYNN-EDMF PBL clouds?
       if(do_mynnedmf) then
          call cloud_mp_MYNN(nCol, nLev, lsmask, t_lay, p_lev, p_lay, qs_lay, relhum,   &
               qc_mynn, qi_mynn, con_ttp, con_g,                                        &
               cld_pbl_lwp, cld_pbl_reliq, cld_pbl_iwp, cld_pbl_reice, cld_pbl_frac)
       endif

       ! Grell-Freitas convective clouds?
       if (imfdeepcnv == imfdeepcnv_gf) then
          alpha0 = 100.
          call cloud_mp_GF(nCol, nLev, lsmask, t_lay, p_lev, p_lay, qs_lay, relhum,     &
               qci_conv, con_ttp, con_g, alpha0,                                        &
               cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp, cld_cnv_reice, cld_cnv_frac)
       endif

       ! SAMF scale & aerosol-aware mass-flux convective clouds?
       if (imfdeepcnv == imfdeepcnv_samf) then
          alpha0 = 200.
          call cloud_mp_SAMF(nCol, nLev, t_lay, p_lev, p_lay, qs_lay, relhum,           &
               cnv_mixratio, con_ttp, con_g, alpha0,                                    &
               cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp, cld_cnv_reice, cld_cnv_frac)
       endif

       ! Update particle size using modified mixing-ratios from Thompson.
       call cmp_reff_Thompson(nLev, nCol, i_cldliq, i_cldice, i_cldsnow, i_cldice_nc,   &
            i_cldliq_nc, i_twa, q_lay, p_lay, t_lay, tracer, con_eps, con_rd, ltaerosol,&
            mraerosol, lsmask,  effrin_cldliq, effrin_cldice, effrin_cldsnow)
       cld_reliq  = effrin_cldliq
       cld_reice  = effrin_cldice
       cld_resnow = effrin_cldsnow

       ! Thomson MP using modified Xu-Randall cloud-fraction (additionally conditioned on RH)
       alpha0 = 2000.
       if (lmfshal) then
          alpha0 = 100.
          if (lmfdeep2) alpha0 = 200.
       endif
       call cloud_mp_thompson(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice, i_cldrain,&
            i_cldsnow, i_cldgrpl, p_lev, p_lay, tv_lay, t_lay, tracer, qs_lay, q_lay,   &
            relhum, con_ttp, con_g, con_rd, con_eps, alpha0, cnv_mixratio, lwp_ex,      &
            iwp_ex, lwp_fc, iwp_fc, cld_frac, cld_lwp, cld_iwp, cld_swp, cld_rwp,       &
            cond_cfrac_onRH = .true., doGP_smearclds = doGP_smearclds)
    endif

    ! Bound effective radii for RRTMGP, LUT's for cloud-optics go from
    !   2.5 - 21.5 microns for liquid clouds,
    !   10  - 180  microns for ice-clouds
    if (doGP_cldoptics_PADE .or. doGP_cldoptics_LUT) then
       where(cld_reliq .lt. radliq_lwr) cld_reliq = radliq_lwr
       where(cld_reliq .gt. radliq_upr) cld_reliq = radliq_upr
       where(cld_reice .lt. radice_lwr) cld_reice = radice_lwr
       where(cld_reice .gt. radice_upr) cld_reice = radice_upr
       if (imfdeepcnv == imfdeepcnv_samf .or. imfdeepcnv == imfdeepcnv_gf) then
          where(cld_cnv_reliq .lt. radliq_lwr) cld_cnv_reliq = radliq_lwr
          where(cld_cnv_reliq .gt. radliq_upr) cld_cnv_reliq = radliq_upr
          where(cld_cnv_reice .lt. radice_lwr) cld_cnv_reice = radice_lwr
          where(cld_cnv_reice .gt. radice_upr) cld_cnv_reice = radice_upr
       endif
       if (do_mynnedmf) then
          where(cld_pbl_reliq .lt. radliq_lwr) cld_pbl_reliq = radliq_lwr
          where(cld_pbl_reliq .gt. radliq_upr) cld_pbl_reliq = radliq_upr
          where(cld_pbl_reice .lt. radice_lwr) cld_pbl_reice = radice_lwr
          where(cld_pbl_reice .gt. radice_upr) cld_pbl_reice = radice_upr
       endif
    endif

    ! Instantaneous 2D (max-in-column) cloud fraction
    do iCol = 1, nCol
       cldfra2d(iCol) = 0._kind_phys
       do iLay = 1, nLev-1
          cldfra2d(iCol) = max(cldfra2d(iCol), cld_frac(iCol,iLay))
       enddo
    enddo

    precip_frac(1:nCol,1:nLev) = cld_frac(1:nCol,1:nLev)

  end subroutine GFS_rrtmgp_cloud_mp_run

!> \ingroup GFS_rrtmgp_cloud_mp
!! Compute cloud radiative properties for Grell-Freitas convective cloud scheme.
!!                 (Adopted from module_SGSCloud_RadPre)
!!  
!! - The total convective cloud condensate is partitoned by phase, using temperature, into
!!     liquid/ice convective cloud mixing-ratios. Compute convective cloud LWP and IWP's.
!!
!! - The liquid and ice cloud effective particle sizes are assigned reference values*.
!!   *TODO* Find references, include DOIs, parameterize magic numbers, etc...
!!
!! - The convective cloud-fraction is computed using Xu-Randall (1996).
!!   (DJS asks: Does the GF scheme produce a cloud-fraction? If so, maybe use instead of 
!!              Xu-Randall? Xu-Randall is consistent with the Thompson MP scheme, but 
!!              not GFDL-EMC)
!!
!! \section cloud_mp_GF_gen General Algorithm
  subroutine cloud_mp_GF(nCol, nLev, lsmask, t_lay, p_lev, p_lay, qs_lay, relhum,        &
       qci_conv, con_ttp, con_g, alpha0, cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp,        &
       cld_cnv_reice, cld_cnv_frac)
    implicit none

    ! Inputs
    integer, intent(in)    :: &
         nCol,          & ! Number of horizontal grid points
         nLev             ! Number of vertical layers
    real(kind_phys), dimension(:), intent(in) :: &
         lsmask           ! Land/Sea mask
    real(kind_phys), intent(in) :: &
         con_g,         & ! Physical constant: gravitational constant 
         con_ttp,       & ! Triple point temperature of water (K)
         alpha0           !
    real(kind_phys), dimension(:,:),intent(in) :: &
         t_lay,         & ! Temperature at layer centers (K)
         p_lev,         & ! Pressure at layer interfaces (Pa)
         p_lay,         & !
         qs_lay,        & !
         relhum,        & !
         qci_conv         !
    ! Outputs
    real(kind_phys), dimension(:,:),intent(inout) :: &
         cld_cnv_lwp,   & ! Convective cloud liquid water path
         cld_cnv_reliq, & ! Convective cloud liquid effective radius
         cld_cnv_iwp,   & ! Convective cloud ice water path
         cld_cnv_reice, & ! Convective cloud ice effecive radius
         cld_cnv_frac     ! Convective cloud-fraction (1)
    ! Local
    integer :: iCol, iLay
    real(kind_phys) :: tem1, deltaP, clwc, qc, qi

    tem1 = 1.0e5/con_g
    do iLay = 1, nLev
       do iCol = 1, nCol
          if (qci_conv(iCol,iLay) > 0.) then
             ! Partition the convective clouds by phase.
             qc = qci_conv(iCol,iLay)*(     min(1., max(0., (t_lay(iCol,iLay)-244.)*0.04)))
             qi = qci_conv(iCol,iLay)*(1. - min(1., max(0., (t_lay(iCol,iLay)-244.)*0.04)))

             ! Compute LWP/IWP
             deltaP = abs(p_lev(iCol,iLay+1)-p_lev(iCol,iLay))*0.01
             cld_cnv_lwp(iCol,iLay) = max(0., qc * tem1 * deltaP)
             cld_cnv_iwp(iCol,iLay) = max(0., qi * tem1 * deltaP)

             ! Particle sizes
             if (nint(lsmask(iCol)) == 1) then !land
                if(qc > 1.E-8) cld_cnv_reliq(iCol,iLay) = 5.4
             else
                !eff radius cloud water (microns), from Miles et al. 
                if(qc > 1.E-8) cld_cnv_reliq(iCol,iLay) = 9.6
             endif
             !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 6b) 
             if(qi > 1.E-8) cld_cnv_reice(iCol,iLay) = max(173.45 + 2.14*(t_lay(iCol,iLay)-273.15), 20.)
      
             ! Xu-Randall (1996) cloud-fraction.
             cld_cnv_frac(iCol,iLay) = cld_frac_XuRandall(p_lay(iCol,iLay),            &
                  qs_lay(iCol,iLay), relhum(iCol,iLay), qc+qi, alpha0)
          endif
       enddo
    enddo
  end subroutine cloud_mp_GF

!> \ingroup GFS_rrtmgp_cloud_mp 
!! Compute cloud radiative properties for MYNN-EDMF PBL cloud scheme.
!!                    (Adopted from module_SGSCloud_RadPre)
!!
!! - Cloud-fraction, liquid, and ice condensate mixing-ratios from MYNN-EDMF cloud scheme
!!   are provided as inputs. Cloud LWP and IWP are computed.
!!
!! - The liquid and ice cloud effective particle sizes are assigned reference values*.
!!   *TODO* Find references, include DOIs, parameterize magic numbers, etc...
!!
!! \section cloud_mp_MYNN_gen General Algorithm
  subroutine cloud_mp_MYNN(nCol, nLev, lsmask, t_lay, p_lev, p_lay, qs_lay, relhum,      &
       qc_mynn, qi_mynn, con_ttp, con_g, cld_pbl_lwp, cld_pbl_reliq, cld_pbl_iwp,     &
       cld_pbl_reice, cld_pbl_frac)
    implicit none

    ! Inputs
    integer, intent(in)    :: &
         nCol,          & ! Number of horizontal grid points
         nLev             ! Number of vertical layers
    real(kind_phys), dimension(:), intent(in) :: &
         lsmask           ! Land/Sea mask
    real(kind_phys), intent(in) :: &
         con_g,         & ! Physical constant: gravitational constant 
         con_ttp          ! Triple point temperature of water (K)
    real(kind_phys), dimension(:,:),intent(in) :: &
         t_lay,         & ! Temperature at layer centers (K)
         p_lev,         & ! Pressure at layer interfaces (Pa)
         p_lay,         & !
         qs_lay,        & !
         relhum,        & !
         qc_mynn,       & ! Liquid cloud mixing-ratio (MYNN PBL cloud)
         qi_mynn,       & ! Ice cloud mixing-ratio (MYNN PBL cloud)
         cld_pbl_frac    ! Cloud-fraction (MYNN PBL cloud)
    ! Outputs
    real(kind_phys), dimension(:,:),intent(inout) :: &
         cld_pbl_lwp,   & ! Convective cloud liquid water path
         cld_pbl_reliq, & ! Convective cloud liquid effective radius
         cld_pbl_iwp,   & ! Convective cloud ice water path
         cld_pbl_reice    ! Convective cloud ice effecive radius
    
    ! Local
    integer :: iCol, iLay
    real(kind_phys) :: tem1, qc, qi, deltaP

    tem1 = 1.0e5/con_g
    do iLay = 1, nLev
       do iCol = 1, nCol
          if (cld_pbl_frac(iCol,iLay) > cld_limit_lower) then
             ! Cloud mixing-ratios (DJS asks: Why is this done?)
             qc = qc_mynn(iCol,iLay)*cld_pbl_frac(iCol,iLay)
             qi = qi_mynn(iCol,iLay)*cld_pbl_frac(iCol,iLay)

             ! LWP/IWP
             deltaP = abs(p_lev(iCol,iLay+1)-p_lev(iCol,iLay))
             cld_pbl_lwp(iCol,iLay) = max(0., qc * tem1 * deltaP)
             cld_pbl_iwp(iCol,iLay) = max(0., qi * tem1 * deltaP)

             ! Particle sizes
             if (nint(lsmask(iCol)) == 1) then
                if(qc > 1.E-8) cld_pbl_reliq(iCol,iLay) = 5.4
             else
                ! Cloud water (microns), from Miles et al.
                if(qc > 1.E-8) cld_pbl_reliq(iCol,iLay) = 9.6
             endif
             ! Cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 6b) 
             if(qi > 1.E-8) cld_pbl_reice(iCol,iLay) = max(173.45 + 2.14*(t_lay(iCol,iLay)-273.15), 20.)
          endif
       enddo
    enddo
  end subroutine cloud_mp_MYNN


!> \ingroup GFS_rrtmgp_cloud_mp 
!! Compute cloud radiative properties for SAMF convective cloud scheme.
!!
!! - The total-cloud convective mixing-ratio is partitioned by phase into liquid/ice 
!!   cloud properties. LWP and IWP are computed.
!!
!! - The liquid and ice cloud effective particle sizes are assigned reference values.
!!
!! - The convective cloud-fraction is computed using Xu-Randall (1996).
!!   (DJS asks: Does the SAMF scheme produce a cloud-fraction?)
!!
!! \section cloud_mp_SAMF_gen General Algorithm
  subroutine cloud_mp_SAMF(nCol, nLev, t_lay, p_lev, p_lay, qs_lay, relhum,              &
       cnv_mixratio, con_ttp, con_g, alpha0, cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp,    &
       cld_cnv_reice, cld_cnv_frac)
    implicit none

    ! Inputs
    integer, intent(in)    :: &
         nCol,          & ! Number of horizontal grid points
         nLev             ! Number of vertical layers
    real(kind_phys), intent(in) :: &
         con_g,         & ! Physical constant: gravity         (m s-2)
         con_ttp,       & ! Triple point temperature of water  (K)
         alpha0           !
    real(kind_phys), dimension(:,:),intent(in) :: &
         t_lay,         & ! Temperature at layer-centers       (K)
         p_lev,         & ! Pressure at layer-interfaces       (Pa)
         p_lay,         & ! Presure at layer-centers           (Pa)
         qs_lay,        & ! Specific-humidity at layer-centers (kg/kg)
         relhum,        & ! Relative-humidity                  (1)
         cnv_mixratio     ! Convective cloud mixing-ratio      (kg/kg)
    ! Outputs
    real(kind_phys), dimension(:,:),intent(inout) :: &
         cld_cnv_lwp,   & ! Convective cloud liquid water path
         cld_cnv_reliq, & ! Convective cloud liquid effective radius
         cld_cnv_iwp,   & ! Convective cloud ice water path
         cld_cnv_reice, & ! Convective cloud ice effecive radius
         cld_cnv_frac     ! Convective cloud-fraction
    ! Local
    integer :: iCol, iLay
    real(kind_phys) :: tem0, tem1, deltaP, clwc

    tem0 = 1.0e5/con_g
    do iLay = 1, nLev
       do iCol = 1, nCol
          if (cnv_mixratio(iCol,iLay) > 0._kind_phys) then
             tem1   = min(1.0, max(0.0, (con_ttp-t_lay(iCol,iLay))*0.05))
             deltaP = abs(p_lev(iCol,iLay+1)-p_lev(iCol,iLay))*0.01
             clwc   = max(0.0, cnv_mixratio(iCol,iLay)) * tem0 * deltaP
             cld_cnv_iwp(iCol,iLay)   = clwc * tem1
             cld_cnv_lwp(iCol,iLay)   = clwc - cld_cnv_iwp(iCol,iLay)
             cld_cnv_reliq(iCol,iLay) = reliq_def
             cld_cnv_reice(iCol,iLay) = reice_def

             ! Xu-Randall (1996) cloud-fraction.
             cld_cnv_frac(iCol,iLay) = cld_frac_XuRandall(p_lay(iCol,iLay),              &
               qs_lay(iCol,iLay), relhum(iCol,iLay), cnv_mixratio(iCol,iLay), alpha0)
          endif
       enddo
    enddo

  end subroutine cloud_mp_SAMF

!> \ingroup GFS_rrtmgp_cloud_mp 
!! This routine computes the cloud radiative properties for a "unified cloud".
!! - "unified cloud" implies that the cloud-fraction is PROVIDED.
!! - The cloud water path is computed for all provided cloud mixing-ratios and hydrometeors.
!! - If particle sizes are provided, they are used. If not, default values are assigned.
!! \section cloud_mp_uni_gen General Algorithm
  subroutine cloud_mp_uni(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice, i_cldrain,     &
       i_cldsnow, i_cldgrpl, i_cldtot, effr_in, kdt, lsmask, p_lev, p_lay, t_lay, tv_lay,&
       effrin_cldliq, effrin_cldice, effrin_cldsnow, tracer, con_g, con_rd, con_ttp,     &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,   &
       cld_rerain, effrin_cldrain)
    implicit none
    
    ! Inputs   
    integer, intent(in)    :: &
         nCol,                 & ! Number of horizontal grid points
         nLev,                 & ! Number of vertical layers
         ncnd,                 & ! Number of cloud condensation types.
         nTracers,             & ! Number of tracers from model. 
         i_cldliq,             & ! Index into tracer array for cloud liquid. 
         i_cldice,             & ! Index into tracer array for cloud ice.
         i_cldrain,            & ! Index into tracer array for cloud rain.
         i_cldsnow,            & ! Index into tracer array for cloud snow.
         i_cldgrpl,            & ! Index into tracer array for cloud groupel.
         i_cldtot,             & ! Index into tracer array for cloud total amount.
         kdt
    logical, intent(in) :: &
    	 effr_in                 ! Provide hydrometeor radii from macrophysics?
    real(kind_phys), intent(in) :: &
         con_g,                & ! Physical constant: gravitational constant
         con_ttp,              & ! Triple point temperature of water (K)  
         con_rd                  ! Physical constant: gas-constant for dry air
    real(kind_phys), dimension(:), intent(in) :: &
         lsmask
    real(kind_phys), dimension(:,:), intent(in) :: &         
         t_lay,                & ! Temperature at model-layers (K)
         tv_lay,               & ! Virtual temperature (K)
         p_lay,                & ! Pressure at model-layers (Pa)
         cld_frac,             & ! Total cloud fraction 
         effrin_cldliq,        & ! Effective radius for liquid cloud-particles (microns)
         effrin_cldice,        & ! Effective radius for ice cloud-particles (microns)
         effrin_cldsnow          ! Effective radius for snow cloud-particles (microns)
    real(kind_phys), dimension(:,:), intent(in) ,optional :: &
         effrin_cldrain          ! Effective radius for rain cloud-particles (microns) 
    real(kind_phys), dimension(:,:), intent(in) :: &
         p_lev                   ! Pressure at model-level interfaces (Pa)
    real(kind_phys), dimension(:,:,:),intent(in) :: &
         tracer                  ! Cloud condensate amount in layer by type ()         
    
    ! Outputs
    real(kind_phys), dimension(:,:),intent(inout) :: &
         cld_lwp,              & ! Cloud liquid water path
         cld_reliq,            & ! Cloud liquid effective radius
         cld_iwp,              & ! Cloud ice water path
         cld_reice,            & ! Cloud ice effecive radius
         cld_swp,              & ! Cloud snow water path
         cld_resnow,           & ! Cloud snow effective radius
         cld_rwp,              & ! Cloud rain water path
         cld_rerain              ! Cloud rain effective radius       

    ! Local variables
    real(kind_phys) :: tem1,tem2,tem3,pfac,deltaP
    real(kind_phys), dimension(nCol, nLev, min(4,ncnd)) :: cld_condensate
    integer :: iCol,iLay,l,ncndl

    ! Cloud condensate
    cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,i_cldliq)        ! -liquid water
    cld_condensate(1:nCol,1:nLev,2) = tracer(1:nCol,1:nLev,i_cldice)        ! -ice water
    if (ncnd > 2) then
       cld_condensate(1:nCol,1:nLev,3) = tracer(1:nCol,1:nLev,i_cldrain)    ! -rain water
       cld_condensate(1:nCol,1:nLev,4) = tracer(1:nCol,1:nLev,i_cldsnow) + &! -snow + grapuel
                                         tracer(1:nCol,1:nLev,i_cldgrpl) 
    endif

    ! Cloud water path (g/m2)
    tem1 = 1.0e5/con_g
    do iLay = 1, nLev
       do iCol = 1, nCol
          ! Compute liquid/ice condensate path from mixing ratios (kg/kg)->(g/m2)   
          if (cld_frac(iCol,iLay) > cld_limit_lower) then
             deltaP = abs(p_lev(iCol,iLay+1)-p_lev(iCol,iLay))*0.01
             cld_lwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,1) * tem1 * deltaP)
             cld_iwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,2) * tem1 * deltaP)
             if (ncnd > 2) then
                cld_rwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,3) * tem1 * deltaP)
                cld_swp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,4) * tem1 * deltaP)
             endif
          endif
       enddo
    enddo

    ! Particle size
    do iLay = 1, nLev
       do iCol = 1, nCol
          ! Use radii provided from the macrophysics
          if (effr_in) then
             cld_reliq(iCol,iLay)  = effrin_cldliq(iCol,iLay)
             cld_reice(iCol,iLay)  = max(reice_min, min(reice_max,effrin_cldice(iCol,iLay)))
             cld_resnow(iCol,iLay) = effrin_cldsnow(iCol,iLay)
             if (present(effrin_cldrain)) then
                cld_rerain(iCol,iLay) = effrin_cldrain(iCol,iLay)
             else
                cld_rerain(iCol,iLay) = rerain_def
             endif
          else
             ! Compute effective liquid cloud droplet radius over land.
             if (nint(lsmask(iCol)) == 1) then
                cld_reliq(iCol,iLay) = 5.0 + 5.0 * min(1.0, max(0.0, (con_ttp-t_lay(iCol,iLay))*0.05))
             endif
             ! Compute effective ice cloud droplet radius following Heymsfield
             ! and McFarquhar (1996) \cite heymsfield_and_mcfarquhar_1996.
             tem2 = t_lay(iCol,iLay) - con_ttp
             if (cld_iwp(iCol,iLay) > 0.0) then
                deltaP = abs(p_lev(iCol,iLay+1)-p_lev(iCol,iLay))*0.01
                tem3 = (con_g/con_rd ) * cld_iwp(iCol,iLay) * (0.01*p_lay(iCol,iLay)) / (deltaP*tv_lay(iCol,iLay))
                if (tem2 < -50.0) then
                   cld_reice(iCol,iLay) = (1250.0/9.917) * tem3 ** 0.109
                elseif (tem2 < -40.0) then
                   cld_reice(iCol,iLay) = (1250.0/9.337) * tem3 ** 0.08
                elseif (tem2 < -30.0) then
                   cld_reice(iCol,iLay) = (1250.0/9.208) * tem3 ** 0.055
                else
                   cld_reice(iCol,iLay) = (1250.0/9.387) * tem3 ** 0.031
                endif
                cld_reice(iCol,iLay)   = max(10.0, min(cld_reice(iCol,iLay), 150.0))
             endif
          endif ! effr_in
       enddo    ! nCol
    enddo       ! nLev

  end subroutine cloud_mp_uni
!> \ingroup GFS_rrtmgp_cloud_mp 
!! This routine computes the cloud radiative properties for the Thompson cloud micro-
!! physics scheme.
!!
!! - The cloud water path is computed for all provided cloud mixing-ratios and hydrometeors.
!!
!! - There are no assumptions about particle size applied here. Effective particle sizes 
!!   are updated prior to this routine, see cmp_reff_Thompson().
!!
!! - The cloud-fraction is computed using Xu-Randall** (1996).
!!   **Additionally, Conditioned on relative-humidity**
!!
!! \section cloud_mp_thompson_gen General Algorithm
  subroutine cloud_mp_thompson(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice, i_cldrain,&
       i_cldsnow, i_cldgrpl, p_lev, p_lay, tv_lay, t_lay, tracer, qs_lay, q_lay, relhum, &
       con_ttp, con_g, con_rd, con_eps, alpha0, cnv_mixratio, lwp_ex, iwp_ex, lwp_fc,    &
       iwp_fc, cld_frac, cld_lwp, cld_iwp, cld_swp, cld_rwp, cond_cfrac_onRH, doGP_smearclds)
    implicit none

    ! Inputs
    logical, intent(in), optional :: &
         cond_cfrac_onRH,   & ! If true, cloud-fracion set to unity when rh>99%
         doGP_smearclds       ! If true, add sgs clouds to gridmean clouds
    integer, intent(in)    :: &
         nCol,              & ! Number of horizontal grid points
         nLev,              & ! Number of vertical layers
         ncnd,              & ! Number of cloud condensation types.
         nTracers,          & ! Number of tracers from model.
         i_cldliq,          & ! Index into tracer array for cloud liquid amount.
         i_cldice,          & !                             cloud ice amount.
         i_cldrain,         & !                             cloud rain amount.
         i_cldsnow,         & !                             cloud snow amount.
         i_cldgrpl            !                             cloud groupel amount.
    real(kind_phys), intent(in) :: &
         con_ttp,           & ! Triple point temperature of water (K)  
         con_g,             & ! Physical constant: gravitational constant
         con_rd,            & ! Physical constant: gas-constant for dry air
         con_eps,           & ! Physical constant: gas constant air / gas constant H2O
         alpha0               !
    real(kind_phys), dimension(:,:), intent(in) :: &
         tv_lay,            & ! Virtual temperature (K)
         t_lay,             & ! Temperature (K)
         qs_lay,            & ! Saturation vapor pressure (Pa)
         q_lay,             & ! water-vapor mixing ratio (kg/kg)
         relhum,            & ! Relative humidity
         p_lay,             & ! Pressure at model-layers (Pa)
         cnv_mixratio         ! Convective cloud mixing-ratio (kg/kg) 
    real(kind_phys), dimension(:,:), intent(in) :: &
         p_lev                ! Pressure at model-level interfaces (Pa)
    real(kind_phys), dimension(:,:,:),intent(in) :: &
         tracer               ! Cloud condensate amount in layer by type ()

    ! In/Outs
    real(kind_phys), dimension(:), intent(inout) :: &
         lwp_ex,            & ! total liquid water path from explicit microphysics 
         iwp_ex,            & ! total ice    water path from explicit microphysics 
         lwp_fc,            & ! total liquid water path from cloud fraction scheme 
         iwp_fc               ! total ice    water path from cloud fraction scheme 
    real(kind_phys), dimension(:,:), intent(inout) :: &
         cld_frac,          & ! Total cloud fraction
         cld_lwp,           & ! Cloud liquid water path
         cld_iwp,           & ! Cloud ice water path
         cld_swp,           & ! Cloud snow water path
         cld_rwp              ! Cloud rain water path

    ! Local variables
    real(kind_phys) :: tem1, pfac, cld_mr, deltaP, tem2
    real(kind_phys), dimension(nCol, nLev, min(4,ncnd)) :: cld_condensate
    integer :: iCol,iLay,l

    ! Cloud condensate
    cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,i_cldliq)     ! -liquid water
    cld_condensate(1:nCol,1:nLev,2) = tracer(1:nCol,1:nLev,i_cldice)     ! -ice water
    cld_condensate(1:nCol,1:nLev,3) = tracer(1:nCol,1:nLev,i_cldrain)    ! -rain hydrometeors
    cld_condensate(1:nCol,1:nLev,4) = tracer(1:nCol,1:nLev,i_cldsnow)    ! -snow hydrometeors

    cld_lwp(:,:) = 0.0
    cld_iwp(:,:) = 0.0
    cld_rwp(:,:) = 0.0
    cld_swp(:,:) = 0.0
    cld_frac(:,:) = 0.0
    tem1 = 1.0e5/con_g
    do iLay = 1, nLev-1
       do iCol = 1, nCol
          ! Add convective cloud to gridmean cloud?
          if (doGP_smearclds) then
             tem2 = min(1.0, max(0.0, (con_ttp-t_lay(iCol,iLay))*0.05))
             cld_condensate(iCol,iLay,1) = cld_condensate(iCol,iLay,1) + cnv_mixratio(iCol,iLay)*(1._kind_phys - tem2)
             cld_condensate(iCol,iLay,2) = cld_condensate(iCol,iLay,2) + cnv_mixratio(iCol,iLay)*tem2
          endif
          ! Compute liquid/ice condensate path from mixing ratios (kg/kg)->(g/m2)
          deltaP              = abs(p_lev(iCol,iLay+1)-p_lev(iCol,iLay))*0.01
          cld_lwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,1) * tem1 * deltaP)
          cld_iwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,2) * tem1 * deltaP)
          cld_rwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,3) * tem1 * deltaP)
          cld_swp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,4) * tem1 * deltaP)
       
          ! Xu-Randall (1996) cloud-fraction. **Additionally, Conditioned on relative-humidity**
          if (present(cond_cfrac_onRH) .and. relhum(iCol,iLay) > 0.99) then
             cld_frac(iCol,iLay) = 1._kind_phys
          else
             cld_mr = cld_condensate(iCol,iLay,1) + cld_condensate(iCol,iLay,2) +  &
                  cld_condensate(iCol,iLay,3) + cld_condensate(iCol,iLay,4)
             cld_frac(iCol,iLay) = cld_frac_XuRandall(p_lay(iCol,iLay),            &
                  qs_lay(iCol,iLay), relhum(iCol,iLay), cld_mr, alpha0)
          endif
       enddo
    enddo

    ! Sum the liquid water and ice paths that come from explicit micro 
    ! What portion of water and ice contents is associated with the partly cloudy boxes?
    do iCol = 1, nCol
       lwp_ex(iCol) = 0.0
       iwp_ex(iCol) = 0.0
       lwp_fc(iCol) = 0.0
       iwp_fc(iCol) = 0.0
       do iLay = 1, nLev-1
          lwp_ex(iCol) = lwp_ex(iCol) + cld_lwp(iCol,iLay)
          iwp_ex(iCol) = iwp_ex(iCol) + cld_iwp(iCol,iLay) + cld_swp(iCol,iLay)
          if (cld_frac(iCol,iLay) .ge. cld_limit_lower .and. &
              cld_frac(iCol,iLay) .lt. cld_limit_ovcst) then
             lwp_fc(iCol) = lwp_fc(iCol) + cld_lwp(iCol,iLay)
             iwp_fc(iCol) = iwp_fc(iCol) + cld_iwp(iCol,iLay) + cld_swp(iCol,iLay)
          endif
       enddo
       lwp_fc(iCol) = lwp_fc(iCol)*1.E-3
       iwp_fc(iCol) = iwp_fc(iCol)*1.E-3
       lwp_ex(iCol) = lwp_ex(iCol)*1.E-3
       iwp_ex(iCol) = iwp_ex(iCol)*1.E-3
    enddo

  end subroutine cloud_mp_thompson

!> \ingroup GFS_rrtmgp_cloud_mp 
!! This function computes the cloud-fraction following.
!! Xu-Randall(1996) A Semiempirical Cloudiness Parameterization for Use in Climate Models
!! https://doi.org/10.1175/1520-0469(1996)053<3084:ASCPFU>2.0.CO;2
!!
!! cld_frac = {1-exp[-alpha*cld_mr/((1-relhum)*qs_lay)**lambda]}*relhum**P
!!
!! \section cld_frac_XuRandall_gen General Algorithm
  function cld_frac_XuRandall(p_lay, qs_lay, relhum, cld_mr, alpha)
    implicit none
    ! Inputs
    real(kind_phys), intent(in) :: &
       p_lay,    & ! Pressure (Pa)
       qs_lay,   & ! Saturation vapor-pressure (Pa)
       relhum,   & ! Relative humidity
       cld_mr,   & ! Total cloud mixing ratio
       alpha       ! Scheme parameter (default=100)

    ! Outputs
    real(kind_phys) :: cld_frac_XuRandall

    ! Locals
    real(kind_phys) :: clwt, clwm, onemrh, tem1, tem2, tem3

    ! Parameters
    real(kind_phys) :: &
       lambda = 0.50, & !
       P      = 0.25

    clwt = 1.0e-6 * (p_lay*0.001)
    if (cld_mr > clwt) then
       onemrh = max(1.e-10, 1.0 - relhum)
       tem1   = alpha / min(max((onemrh*qs_lay)**lambda,0.0001),1.0)
       tem2   = max(min(tem1*(cld_mr - clwt), 50.0 ), 0.0 )
       tem3   = sqrt(sqrt(relhum)) ! This assumes "p" = 0.25. Identical, but cheaper than relhum**p
       !
       cld_frac_XuRandall = max( tem3*(1.0-exp(-tem2)), 0.0 )
    else
       cld_frac_XuRandall = 0.0
    endif

    return
  end function

  ! ######################################################################################
  ! This routine is a wrapper to update the Thompson effective particle sizes used by the
  ! RRTMGP radiation scheme.
  !
  ! ######################################################################################
  subroutine cmp_reff_Thompson(nLev, nCol, i_cldliq, i_cldice, i_cldsnow, i_cldice_nc,   &
       i_cldliq_nc, i_twa, q_lay, p_lay, t_lay, tracer, con_eps, con_rd, ltaerosol,      &
       mraerosol, lsmask, effrin_cldliq, effrin_cldice, effrin_cldsnow)
    implicit none

    ! Inputs
    integer, intent(in) :: nLev, nCol, i_cldliq, i_cldice, i_cldsnow, i_cldice_nc,       &
         i_cldliq_nc, i_twa
    logical, intent(in) :: ltaerosol, mraerosol
    real(kind_phys), intent(in) :: con_eps,con_rd
    real(kind_phys), dimension(:,:),intent(in) :: q_lay, p_lay, t_lay
    real(kind_phys), dimension(:,:,:),intent(in) :: tracer
    real(kind_phys), dimension(:), intent(in) :: lsmask

    ! Outputs
    real(kind_phys), dimension(:,:), intent(inout) :: effrin_cldliq, effrin_cldice,      &
         effrin_cldsnow

    ! Local
    integer :: iCol, iLay
    real(kind_phys) :: rho, orho
    real(kind_phys),dimension(nCol,nLev) :: qv_mp, qc_mp, qi_mp, qs_mp, ni_mp, nc_mp,    &
         nwfa, re_cloud, re_ice, re_snow
    integer :: ilsmask 

    ! Prepare cloud mixing-ratios and number concentrations for calc_effectRa
    do iLay = 1, nLev
       do iCol = 1, nCol
          qv_mp(iCol,iLay) = q_lay(iCol,iLay)/(1.-q_lay(iCol,iLay))
          rho = con_eps*p_lay(iCol,iLay)/(con_rd*t_lay(iCol,iLay)*(qv_mp(iCol,iLay)+con_eps))
          orho = 1./rho
          qc_mp(iCol,iLay) = tracer(iCol,iLay,i_cldliq)    / (1.-q_lay(iCol,iLay))
          qi_mp(iCol,iLay) = tracer(iCol,iLay,i_cldice)    / (1.-q_lay(iCol,iLay))
          qs_mp(iCol,iLay) = tracer(iCol,iLay,i_cldsnow)   / (1.-q_lay(iCol,iLay))
          ni_mp(iCol,iLay) = tracer(iCol,iLay,i_cldice_nc) / (1.-q_lay(iCol,iLay))
          if (ltaerosol) then
             nc_mp(iCol,iLay) = tracer(iCol,iLay,i_cldliq_nc) / (1.-q_lay(iCol,iLay))
             nwfa(iCol,iLay)  = tracer(iCol,iLay,i_twa)
             if (qc_mp(iCol,iLay) > 1.e-12 .and. nc_mp(iCol,iLay) < 100.) then
               nc_mp(iCol,iLay) = make_DropletNumber(qc_mp(iCol,iLay)*rho, nwfa(iCol,iLay)*rho) * orho
             endif
          elseif (mraerosol) then
             nc_mp(iCol,iLay) = tracer(iCol,iLay,i_cldliq_nc) / (1.-q_lay(iCol,iLay))
             if (qc_mp(iCol,iLay) > 1.e-12 .and. nc_mp(iCol,iLay) < 100.) then
               nc_mp(iCol,iLay) = make_DropletNumber(qc_mp(iCol,iLay)*rho, nwfa(iCol,iLay)*rho) * orho
             endif
          else
             if (nint(lsmask(iCol)) == 1) then !land
                nc_mp(iCol,iLay) = nt_c_l*orho
             else 
                nc_mp(iCol,iLay) = nt_c_o*orho
             endif 
          endif
          if (qi_mp(iCol,iLay) > 1.e-12 .and. ni_mp(iCol,iLay) < 100.) then
             ni_mp(iCol,iLay) = make_IceNumber(qi_mp(iCol,iLay)*rho, t_lay(iCol,iLay)) * orho
          endif
       enddo
    enddo

    ! Compute effective radii for liquid/ice/snow.
    do iCol=1,nCol
       ilsmask = nint(lsmask(iCol))
       call calc_effectRad (t_lay(iCol,:), p_lay(iCol,:), qv_mp(iCol,:), qc_mp(iCol,:),  &
                            nc_mp(iCol,:), qi_mp(iCol,:), ni_mp(iCol,:), qs_mp(iCol,:),  &
                            re_cloud(iCol,:), re_ice(iCol,:), re_snow(iCol,:), ilsmask,  & 
                            1, nLev )
       do iLay = 1, nLev
          re_cloud(iCol,iLay) = MAX(re_qc_min, MIN(re_cloud(iCol,iLay), re_qc_max))
          re_ice(iCol,iLay)   = MAX(re_qi_min, MIN(re_ice(iCol,iLay),   re_qi_max))
          re_snow(iCol,iLay)  = MAX(re_qs_min, MIN(re_snow(iCol,iLay),  re_qs_max))
       enddo
    enddo

    ! Scale to microns.
    do iLay = 1, nLev
       do iCol = 1, nCol
          effrin_cldliq(iCol,iLay)  = re_cloud(iCol,iLay)*1.e6
          effrin_cldice(iCol,iLay)  = re_ice(iCol,iLay)*1.e6
          effrin_cldsnow(iCol,iLay) = re_snow(iCol,iLay)*1.e6
       enddo
    enddo

  end subroutine cmp_reff_Thompson

end module GFS_rrtmgp_cloud_mp
