! ###########update_#############################################################################
! ########################################################################################
module GFS_rrtmgp_cloud_mp
  use machine,      only: kind_phys
  use radiation_tools,   only: check_error_msg
  use rrtmgp_lw_cloud_optics, only: &
       radliq_lwr => radliq_lwrLW, radliq_upr => radliq_uprLW,&
       radice_lwr => radice_lwrLW, radice_upr => radice_uprLW  
  use module_mp_thompson, only: calc_effectRad, Nt_c, re_qc_min, re_qc_max, re_qi_min, &
       re_qi_max, re_qs_min, re_qs_max
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
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_cloud_mp_init()
  end subroutine GFS_rrtmgp_cloud_mp_init

!! \section arg_table_GFS_rrtmgp_cloud_mp_run
!! \htmlinclude GFS_rrtmgp_cloud_mp_run_html
!!
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_cloud_mp_run(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice,     &
       i_cldrain, i_cldsnow, i_cldgrpl, i_cldtot, i_cldliq_nc, i_cldice_nc, i_twa, kdt,  &
       imfdeepcnv, imfdeepcnv_gf, doSWrad, doLWrad, effr_in, lmfshal, ltaerosol, icloud, &
       imp_physics, imp_physics_thompson, imp_physics_gfdl, imp_physics_zhao_carr,       &
       imp_physics_zhao_carr_pdf, imp_physics_mg, imp_physics_wsm6, lgfdlmprad,          &
       imp_physics_fer_hires, do_mynnedmf, uni_cld, lmfdeep2, doGP_convcld, p_lev,       &
       p_lay, t_lay, qs_lay, q_lay, relhum, lsmask, tv_lay, effrin_cldliq, effrin_cldice,&
       effrin_cldrain, effrin_cldsnow, tracer, cnv_mixratio, cnv_cldfrac, qci_conv,      &
       con_g, con_rd, con_eps, con_ttp, doGP_cldoptics_PADE, doGP_cldoptics_LUT,         &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,   &
       cld_rerain, precip_frac, cnv_cld_lwp, cnv_cld_reliq, cnv_cld_iwp, cnv_cld_reice,  &
       lwp_ex, iwp_ex, lwp_fc, iwp_fc, errmsg, errflg)

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
         kdt,                       & ! Current forecast iteration
         imp_physics,               & ! Choice of microphysics scheme
         imp_physics_thompson,      & ! Choice of Thompson
         imp_physics_gfdl,          & ! Choice of GFDL
         imp_physics_zhao_carr,     & ! Choice of Zhao-Carr
         imp_physics_zhao_carr_pdf, & ! Choice of Zhao-Carr + PDF clouds
         imp_physics_mg,            & ! Choice of Morrison-Gettelman
         imp_physics_wsm6,          & ! Choice of WSM6
         imp_physics_fer_hires,     & ! Choice of Ferrier-Aligo
         icloud                       ! Control for cloud are fraction option
    logical, intent(in) :: &
         doSWrad,                   & ! Call SW radiation?
         doLWrad,                   & ! Call LW radiation?
         effr_in,                   & ! Provide hydrometeor radii from macrophysics?
         lmfshal,                   & ! Flag for mass-flux shallow convection scheme used by Xu-Randall
         ltaerosol,                 & ! Flag for aerosol option
         lgfdlmprad,                & ! Flag for GFDLMP radiation interaction
         do_mynnedmf,               & ! Flag to activate MYNN-EDMF 
         uni_cld,                   & ! Flag for unified cloud scheme
         lmfdeep2,                  & ! Flag for mass flux deep convection 
         doGP_convcld,              & ! Treat convective clouds seperately?
         doGP_cldoptics_LUT,        & ! Flag to do GP cloud-optics (LUTs)
         doGP_cldoptics_PADE          !                            (PADE approximation)
    real(kind_phys), intent(in) :: &
         con_g,                     & ! Physical constant: gravitational constant
         con_rd,                    & ! Physical constant: gas-constant for dry air
         con_ttp,                   & ! Triple point temperature of water (K)  
         con_eps                      ! Physical constant: gas constant air / gas constant H2O
    real(kind_phys), dimension(:), intent(in) :: &
         lsmask                       ! Land/Sea mask
    real(kind_phys), dimension(:,:), intent(in) :: &         
         tv_lay,                    & ! Virtual temperature (K)
         t_lay,                     & ! Temperature (K)
         qs_lay,                    & ! Saturation vapor pressure (Pa)
         q_lay,                     & ! water-vapor mixing ratio (kg/kg)
         relhum,                    & ! Relative humidity
         p_lay,                     & ! Pressure at model-layers (Pa)
         cnv_mixratio,              & ! Convective cloud mixing-ratio (kg/kg)
         cnv_cldfrac,               & ! Convective cloud-fraction (1)
         qci_conv                     !
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
    real(kind_phys), dimension(:,:),intent(inout) :: &
         cld_frac,                  & ! Total cloud fraction
         cld_lwp,                   & ! Cloud liquid water path
         cld_reliq,                 & ! Cloud liquid effective radius
         cld_iwp,                   & ! Cloud ice water path
         cld_reice,                 & ! Cloud ice effecive radius
         cld_swp,                   & ! Cloud snow water path
         cld_resnow,                & ! Cloud snow effective radius
         cld_rwp,                   & ! Cloud rain water path
         cld_rerain,                & ! Cloud rain effective radius
         precip_frac,               & ! Precipitation fraction
         cnv_cld_lwp,               & ! Water path for       convective liquid cloud-particles (microns) 
         cnv_cld_reliq,             & ! Effective radius for convective liquid cloud-particles (microns)
         cnv_cld_iwp,               & ! Water path for       convective ice cloud-particles (microns)
         cnv_cld_reice                ! Effective radius for convective ice cloud-particles (microns) 
    character(len=*), intent(out) :: &
         errmsg                       ! Error message
    integer, intent(out) :: &  
         errflg                       ! Error flag

    ! Local
    integer :: iCol, iLay

    if (.not. (doSWrad .or. doLWrad)) return

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! ###################################################################################
    ! GFDL Microphysics
    ! ###################################################################################
    if (imp_physics == imp_physics_gfdl) then
       if (.not. lgfdlmprad) then
          ! Call progcld_gfdl_lin
       else

          ! The cloud-fraction used for the radiation is conditional on other mp choices.
          do iLay = 1, nLev
             do iCol = 1, nCol
                if ((imfdeepcnv==imfdeepcnv_gf .or. do_mynnedmf) .and. kdt>1) then
                   if (do_mynnedmf) then
                      if (tracer(iCol,iLay,i_cldrain)>1.0e-7 .OR. tracer(iCol,iLay,i_cldsnow)>1.0e-7) then
                         cld_frac(iCol,iLay) = tracer(iCol,iLay,i_cldtot)
                      endif
                   else
                      if (qci_conv(iCol,iLay) <= 0.) then
                         cld_frac(iCol,iLay) = tracer(iCol,iLay,i_cldtot)
                      endif
                   endif
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
    ! ###################################################################################
    if (imp_physics == imp_physics_thompson) then
       ! Update particle size using modified mixing-ratios.
       call cmp_reff_Thompson(nLev, nCol, i_cldliq, i_cldice, i_cldsnow, i_cldice_nc,   &
            i_cldliq_nc, i_twa, q_lay, p_lay, t_lay, tracer, con_eps, con_rd, ltaerosol,&
            effrin_cldliq, effrin_cldice, effrin_cldsnow)
       cld_reliq  = effrin_cldliq
       cld_reice  = effrin_cldice
       cld_resnow = effrin_cldsnow

       if(do_mynnedmf .or. imfdeepcnv == imfdeepcnv_gf ) then
          if (icloud == 3) then
             ! Call progcld_thompson
          else
             call cloud_mp_uni(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice,          &
                  i_cldrain, i_cldsnow, i_cldgrpl, i_cldtot,  effr_in, kdt, lsmask,     &
                  p_lev, p_lay, t_lay, tv_lay, effrin_cldliq, effrin_cldice,            &
                  effrin_cldsnow, tracer, con_g, con_rd, con_ttp, cld_frac, cld_lwp,    &
                  cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,          &
                  cld_rerain)
          endif
       else
          if (icloud == 3) then
             ! Call progcld_thompson
          else
             !
             if (doGP_convcld) then
                call cloud_mp_convective(nCol, nLev, t_lay, p_lev, cnv_mixratio,        &
                     cnv_cldfrac, con_ttp, con_g, cnv_cld_lwp, cnv_cld_reliq,           &
                     cnv_cld_iwp, cnv_cld_reice)
             endif
             !
             call cloud_mp_thompson(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice,     &
                  i_cldrain, i_cldsnow, i_cldgrpl, i_cldtot, i_cldliq_nc, i_cldice_nc,  &
                  i_twa, p_lev, p_lay, tv_lay, t_lay, tracer, qs_lay, q_lay, relhum,    &
                  con_g, con_rd, con_eps, lmfshal, ltaerosol, imfdeepcnv, imfdeepcnv_gf,&
                  uni_cld, lmfdeep2, lwp_ex, iwp_ex, lwp_fc, iwp_fc, cld_frac, cld_lwp, &
                  cld_iwp, cld_swp, cld_rwp)
          endif
       endif
    endif

    ! Bound effective radii for RRTMGP, LUT's for cloud-optics go from
    !   2.5 - 21.5 microns for liquid clouds,
    !   10  - 180  microns for ice-clouds
    if (doGP_cldoptics_PADE .or. doGP_cldoptics_LUT) then
       where(cld_reliq .lt. radliq_lwr) cld_reliq = radliq_lwr
       where(cld_reliq .gt. radliq_upr) cld_reliq = radliq_upr
       where(cld_reice .lt. radice_lwr) cld_reice = radice_lwr
       where(cld_reice .gt. radice_upr) cld_reice = radice_upr
    endif

    precip_frac(1:nCol,1:nLev) = cld_frac(1:nCol,1:nLev)

  end subroutine GFS_rrtmgp_cloud_mp_run

  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_cloud_mp_finalize()
  end subroutine GFS_rrtmgp_cloud_mp_finalize

  ! ######################################################################################
  ! ######################################################################################
  subroutine cloud_mp_convective(nCol, nLev, t_lay, p_lev, cnv_mixratio, cnv_cldfrac,    &
       con_ttp, con_g, cnv_cld_lwp, cnv_cld_reliq, cnv_cld_iwp, cnv_cld_reice)
    ! Inputs
    integer, intent(in)    :: &
         nCol,          & ! Number of horizontal grid points
         nLev             ! Number of vertical layers
    real(kind_phys), intent(in) :: &
         con_g,         & ! Physical constant: gravitational constant 
         con_ttp          ! Triple point temperature of water (K)
    real(kind_phys), dimension(:,:),intent(in) :: &
         t_lay,         & ! Temperature at layer centers (K)
         p_lev,         & ! Pressure at layer interfaces (Pa)
         cnv_mixratio,  & ! Convective cloud mixing-ratio (kg/kg)
         cnv_cldfrac      ! Convective cloud-fraction (1)
    ! Outputs
    real(kind_phys), dimension(:,:),intent(inout) :: &
         cnv_cld_lwp,   & ! Convective cloud liquid water path
         cnv_cld_reliq, & ! Convective cloud liquid effective radius
         cnv_cld_iwp,   & ! Convective cloud ice water path
         cnv_cld_reice    ! Convective cloud ice effecive radius
    ! Local
    integer :: iCol, iLay
    real(kind_phys) :: tem1, deltaP, clwc

    do iLay = 1, nLev
       do iCol = 1, nCol
          if (cnv_cldfrac(iCol,iLay) > cld_limit_lower) then
             tem1   = min(1.0, max(0.0, (con_ttp-t_lay(iCol,iLay))*0.05))
             deltaP = abs(p_lev(iCol,iLay+1)-p_lev(iCol,iLay))/100.
             clwc   = max(0.0, cnv_mixratio(iCol,iLay)) * con_g * deltaP
             cnv_cld_iwp(iCol,iLay) = clwc * tem1
             cnv_cld_lwp(iCol,iLay) = clwc - cnv_cld_iwp(iCol,iLay)
             cnv_cld_reliq(iCol,iLay) = reliq_def
             cnv_cld_reice(iCol,iLay) = reice_def
        else
             cnv_cld_iwp(iCol,iLay)   = 0._kind_phys
             cnv_cld_lwp(iCol,iLay)   = 0._kind_phys
             cnv_cld_reliq(iCol,iLay) = 0._kind_phys
             cnv_cld_reice(iCol,iLay) = 0._kind_phys
          endif
       enddo
    enddo

  end subroutine cloud_mp_convective

  ! ######################################################################################
  ! ######################################################################################
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
    real(kind_phys) :: tem1,tem2,tem3,pfac
    real(kind_phys), dimension(nCol, nLev, min(4,ncnd)) :: cld_condensate
    integer :: iCol,iLay,l,ncndl
    real(kind_phys), dimension(nCol,nLev) :: deltaP

    ! Cloud condensate
    cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,i_cldliq)        ! -liquid water
    cld_condensate(1:nCol,1:nLev,2) = tracer(1:nCol,1:nLev,i_cldice)        ! -ice water
    if (ncnd > 2) then
       cld_condensate(1:nCol,1:nLev,3) = tracer(1:nCol,1:nLev,i_cldrain)    ! -rain water
       cld_condensate(1:nCol,1:nLev,4) = tracer(1:nCol,1:nLev,i_cldsnow) + &! -snow + grapuel
                                         tracer(1:nCol,1:nLev,i_cldgrpl) 
    endif

    ! Cloud water path (g/m2)
    deltaP = abs(p_lev(:,2:nLev+1)-p_lev(:,1:nLev))/100.  
    do iLay = 1, nLev
       do iCol = 1, nCol
          ! Compute liquid/ice condensate path from mixing ratios (kg/kg)->(g/m2)   
          if (cld_frac(iCol,iLay) > cld_limit_lower) then
             tem1                = (1.0e5/con_g) * deltaP(iCol,iLay)
             cld_lwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,1) * tem1)
             cld_iwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,2) * tem1)
             if (ncnd > 2) then
                cld_rwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,3) * tem1)
                cld_swp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,4) * tem1) 
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
                tem3 = (con_g/con_rd ) * cld_iwp(iCol,iLay) * (0.01*p_lay(iCol,iLay)) / (deltaP(iCol,iLay)*tv_lay(iCol,iLay))
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
  ! ######################################################################################
  ! ######################################################################################
  subroutine cloud_mp_thompson(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice, i_cldrain,&
       i_cldsnow, i_cldgrpl, i_cldtot, i_cldliq_nc, i_cldice_nc, i_twa, p_lev,  &
       p_lay, tv_lay, t_lay, tracer,       &
       qs_lay, q_lay, relhum, con_g, con_rd, con_eps, lmfshal, ltaerosol, imfdeepcnv,    &
       imfdeepcnv_gf, uni_cld, lmfdeep2,                                                 &
       lwp_ex, iwp_ex, lwp_fc, iwp_fc, cld_frac, cld_lwp, cld_iwp, cld_swp, cld_rwp)
    implicit none

    ! Inputs
    integer, intent(in)    :: &
         nCol,              & ! Number of horizontal grid points
         nLev,              & ! Number of vertical layers
         ncnd,              & ! Number of cloud condensation types.
         nTracers,          & ! Number of tracers from model.
         i_cldliq,          & ! Index into tracer array for cloud liquid amount.
         i_cldice,          & !                             cloud ice amount.
         i_cldrain,         & !                             cloud rain amount.
         i_cldsnow,         & !                             cloud snow amount.
         i_cldgrpl,         & !                             cloud groupel amount.
         i_cldtot,          & !                             cloud total amount.
         i_cldliq_nc,       & !                             cloud liquid number concentration.
         i_cldice_nc,       & !                             cloud ice number concentration.
         i_twa,             & !                             water friendly aerosol.
         imfdeepcnv,        & ! Choice of mass-flux deep convection scheme
         imfdeepcnv_gf        ! Flag for Grell-Freitas deep convection scheme
    logical, intent(in) :: &
         uni_cld,           & ! Flag for unified cloud scheme
         lmfshal,           & ! Flag for mass-flux shallow convection scheme used by Xu-Randall
         ltaerosol,         & ! Flag for aerosol option
         lmfdeep2             ! Flag for mass flux deep convection 
    real(kind_phys), intent(in) :: &
         con_g,             & ! Physical constant: gravitational constant
         con_rd,            & ! Physical constant: gas-constant for dry air
         con_eps              ! Physical constant: gas constant air / gas constant H2O

    real(kind_phys), dimension(:,:), intent(in) :: &
         tv_lay,            & ! Virtual temperature (K)
         t_lay,             & ! Temperature (K)
         qs_lay,            & ! Saturation vapor pressure (Pa)
         q_lay,             & ! water-vapor mixing ratio (kg/kg)
         relhum,            & ! Relative humidity
         p_lay                ! Pressure at model-layers (Pa)
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
    real(kind_phys) :: alpha0, pfac, tem1, cld_mr
    real(kind_phys), dimension(nCol, nLev, min(4,ncnd)) :: cld_condensate
    integer :: iCol,iLay,l
    real(kind_phys), dimension(nCol,nLev) :: deltaP

    ! Cloud condensate
    cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,i_cldliq)     ! -liquid water
    cld_condensate(1:nCol,1:nLev,2) = tracer(1:nCol,1:nLev,i_cldice)     ! -ice water
    cld_condensate(1:nCol,1:nLev,3) = tracer(1:nCol,1:nLev,i_cldrain)    ! -rain water
    cld_condensate(1:nCol,1:nLev,4) = tracer(1:nCol,1:nLev,i_cldsnow) + &! -snow + grapuel
                                      tracer(1:nCol,1:nLev,i_cldgrpl)

    ! Cloud water path (g/m2)
    deltaP = abs(p_lev(:,2:nLev+1)-p_lev(:,1:nLev))/100.
    do iLay = 1, nLev
       do iCol = 1, nCol
          ! Compute liquid/ice condensate path from mixing ratios (kg/kg)->(g/m2)
          tem1                = (1.0e5/con_g) * deltaP(iCol,iLay)
          cld_lwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,1) * tem1)
          cld_iwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,2) * tem1)
          cld_rwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,3) * tem1)
          cld_swp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,4) * tem1)
       enddo
    enddo

    ! Compute cloud-fraction. Only if not pre-computed
    if(.not. uni_cld) then
       ! Cloud-fraction
       if(.not. lmfshal) then
          alpha0 = 2000. ! Default (from GATE simulations)
       else
          if (lmfdeep2) then
             alpha0 = 200
          else 
             alpha0 = 100
          endif
       endif
       
       ! Xu-Randall (1996) cloud-fraction. Conditioned on relative-humidity
       do iLay = 1, nLev
          do iCol = 1, nCol
             if (relhum(iCol,iLay) > 0.99) then
                cld_frac(iCol,iLay) = 1._kind_phys
             else
                cld_mr = cld_condensate(iCol,iLay,1) + cld_condensate(iCol,iLay,2) +  &
                     cld_condensate(iCol,iLay,4)
                cld_frac(iCol,iLay) = cld_frac_XuRandall(p_lay(iCol,iLay),            &
                     qs_lay(iCol,iLay), relhum(iCol,iLay), cld_mr, alpha0)
             endif
          enddo
       enddo
    else
       cld_frac = tracer(:,:,i_cldtot)
    endif

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

  ! ######################################################################################
  ! This function computes the cloud-fraction following.
  ! Xu-Randall(1996) A Semiempirical Cloudiness Parameterization for Use in Climate Models
  ! https://doi.org/10.1175/1520-0469(1996)053<3084:ASCPFU>2.0.CO;2
  !
  ! cld_frac = {1-exp[-alpha*cld_mr/((1-relhum)*qs_lay)**lambda]}*relhum**P
  !
  ! ######################################################################################
  function cld_frac_XuRandall(p_lay, qs_lay, relhum, cld_mr, alpha)

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
  ! ######################################################################################
  subroutine cmp_reff_Thompson(nLev, nCol, i_cldliq, i_cldice, i_cldsnow, i_cldice_nc,   &
       i_cldliq_nc, i_twa, q_lay, p_lay, t_lay, tracer, con_eps, con_rd, ltaerosol,      &
       effrin_cldliq, effrin_cldice, effrin_cldsnow)

    implicit none

    ! Inputs
    integer, intent(in) :: nLev, nCol, i_cldliq, i_cldice, i_cldsnow, i_cldice_nc,       &
         i_cldliq_nc, i_twa
    logical, intent(in) :: ltaerosol
    real(kind_phys), intent(in) :: con_eps,con_rd
    real(kind_phys), dimension(:,:),intent(in) :: q_lay, p_lay, t_lay
    real(kind_phys), dimension(:,:,:),intent(in) :: tracer

    ! Outputs
    real(kind_phys), dimension(:,:), intent(inout) :: effrin_cldliq, effrin_cldice,      &
         effrin_cldsnow

    ! Local
    integer :: iCol, iLay
    real(kind_phys) :: rho, orho
    real(kind_phys),dimension(nCol,nLev) :: qv_mp, qc_mp, qi_mp, qs_mp, ni_mp, nc_mp,    &
         nwfa, re_cloud, re_ice, re_snow

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
          else
             nc_mp(iCol,iLay) = nt_c*orho
          endif
          if (qi_mp(iCol,iLay) > 1.e-12 .and. ni_mp(iCol,iLay) < 100.) then
             ni_mp(iCol,iLay) = make_IceNumber(qi_mp(iCol,iLay)*rho, t_lay(iCol,iLay)) * orho
          endif
       enddo
    enddo

    ! Compute effective radii for liquid/ice/snow.
    do iCol=1,nCol
       call calc_effectRad (t_lay(iCol,:), p_lay(iCol,:), qv_mp(iCol,:), qc_mp(iCol,:),  &
                            nc_mp(iCol,:), qi_mp(iCol,:), ni_mp(iCol,:), qs_mp(iCol,:),  &
                            re_cloud(iCol,:), re_ice(iCol,:), re_snow(iCol,:), 1, nLev )
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
