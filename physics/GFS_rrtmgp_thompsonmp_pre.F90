! ########################################################################################
! This module contains the interface between the THOMPSON macrophysics and the RRTMGP radiation
! schemes. Only compatable with Model%imp_physics = Model%imp_physics_thompson
! ########################################################################################
module GFS_rrtmgp_thompsonmp_pre
  use machine, only: &
       kind_phys
  use rrtmgp_aux, only: &
       check_error_msg
  use module_mp_thompson, only: &
       calc_effectRad, Nt_c,    &
       re_qc_min, re_qc_max,    &
       re_qi_min, re_qi_max,    &
       re_qs_min, re_qs_max
  use module_mp_thompson_make_number_concentrations, only: &
       make_IceNumber,      &
       make_DropletNumber, &
       make_RainNumber
  use rrtmgp_lw_cloud_optics, only: radliq_lwr => radliq_lwrLW, radliq_upr => radliq_uprLW,&
                                    radice_lwr => radice_lwrLW, radice_upr => radice_uprLW
  implicit none

  ! Parameters specific to THOMPSON MP scheme.
  real(kind_phys), parameter :: &
       rerain_def = 1000.0 ! Default rain radius to 1000 microns

  public GFS_rrtmgp_thompsonmp_pre_init, GFS_rrtmgp_thompsonmp_pre_run, GFS_rrtmgp_thompsonmp_pre_finalize

contains
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_thompsonmp_pre_init()
  end subroutine GFS_rrtmgp_thompsonmp_pre_init

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_GFS_rrtmgp_thompsonmp_pre_run
!! \htmlinclude GFS_rrtmgp_thompsonmp_pre_run.html
!!
  subroutine GFS_rrtmgp_thompsonmp_pre_run(nCol, nLev, nTracers, ncnd, doSWrad, doLWrad, &
       i_cldliq, i_cldice, i_cldrain, i_cldsnow, i_cldgrpl, i_cldtot, i_cldliq_nc,       &
       i_cldice_nc, i_twa, effr_in, p_lev, p_lay, tv_lay, t_lay, effrin_cldliq,          &
       effrin_cldice, effrin_cldsnow, tracer, qs_lay, q_lay, relhum, cld_frac_mg, con_g, &
       con_rd, con_eps, uni_cld, lmfshal, lmfdeep2, ltaerosol, do_mynnedmf, imfdeepcnv,  &
       imfdeepcnv_gf, doGP_cldoptics_PADE, doGP_cldoptics_LUT,                           &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,   &
       cld_rerain, precip_frac, errmsg, errflg)

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
         doSWrad,           & ! Call SW radiation?
         doLWrad,           & ! Call LW radiation
         effr_in,           & ! Use cloud effective radii provided by model?
         uni_cld,           & ! Use provided cloud-fraction?
         lmfshal,           & ! Flag for mass-flux shallow convection scheme used by Xu-Randall
         lmfdeep2,          & ! Flag for some scale-aware mass-flux convection scheme active
         ltaerosol,         & ! Flag for aerosol option
         do_mynnedmf,       & ! Flag to activate MYNN-EDMF
         doGP_cldoptics_LUT,& ! Flag to do GP cloud-optics (LUTs)
         doGP_cldoptics_PADE  !                            (PADE approximation)
    real(kind_phys), intent(in) :: &
         con_g,             & ! Physical constant: gravitational constant
         con_rd,            & ! Physical constant: gas-constant for dry air
         con_eps              ! Physical constant: gas constant air / gas constant H2O

    real(kind_phys), dimension(nCol,nLev), intent(in) :: &
         tv_lay,            & ! Virtual temperature (K)
         t_lay,             & ! Temperature (K)
         qs_lay,            & ! Saturation vapor pressure (Pa)
         q_lay,             & ! water-vapor mixing ratio (kg/kg)
         relhum,            & ! Relative humidity
         p_lay,             & ! Pressure at model-layers (Pa)
         cld_frac_mg          ! Cloud-fraction from MG scheme. WTF?????
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &
         p_lev                ! Pressure at model-level interfaces (Pa)
    real(kind_phys), dimension(nCol, nLev, nTracers),intent(in) :: &
         tracer               ! Cloud condensate amount in layer by type ()

    ! In/Outs
    real(kind_phys), dimension(nCol,nLev), intent(inout) :: &
         cld_frac,          & ! Total cloud fraction
         cld_lwp,           & ! Cloud liquid water path
         cld_reliq,         & ! Cloud liquid effective radius
         cld_iwp,           & ! Cloud ice water path
         cld_reice,         & ! Cloud ice effecive radius
         cld_swp,           & ! Cloud snow water path
         cld_resnow,        & ! Cloud snow effective radius
         cld_rwp,           & ! Cloud rain water path
         cld_rerain,        & ! Cloud rain effective radius
         precip_frac,       & ! Precipitation fraction
         effrin_cldliq,     & ! Effective radius for liquid cloud-particles (microns)
         effrin_cldice,     & ! Effective radius for ice cloud-particles (microns)
         effrin_cldsnow       ! Effective radius for snow cloud-particles (microns)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &
         errflg               ! Error flag

    ! Local variables
    real(kind_phys) :: alpha0, pfac, tem1, cld_mr
    real(kind_phys), dimension(nCol, nLev, min(4,ncnd)) :: cld_condensate
    integer :: iCol,iLay,l
    real(kind_phys) :: rho, orho
    real(kind_phys), dimension(nCol,nLev) :: deltaP, deltaZ, re_cloud, re_ice,&
         re_snow, qv_mp, qc_mp, qi_mp, qs_mp, nc_mp, ni_mp, nwfa
    logical :: top_at_1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. (doSWrad .or. doLWrad)) return

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

    ! Cloud particle sizes and number concentrations...

    ! Prepare cloud mixing-ratios and number concentrations for calc_effectRad,
    ! and update number concentrations, consistent with sub-grid clouds
    do iLay = 1, nLev
       do iCol = 1, nCol
          qv_mp(iCol,iLay) = q_lay(iCol,iLay)/(1.-q_lay(iCol,iLay))
          rho = con_eps*p_lay(iCol,iLay)/(con_rd*t_lay(iCol,iLay)*(qv_mp(iCol,iLay)+con_eps))
          orho = 1./rho
          qc_mp(iCol,iLay) = tracer(iCol,iLay,i_cldliq)    / (1.-q_lay(iCol,iLay))
          qi_mp(iCol,iLay) = tracer(iCol,iLay,i_cldice)    / (1.-q_lay(iCol,iLay))
          qs_mp(iCol,iLay) = tracer(iCol,iLay,i_cldsnow)   / (1.-q_lay(iCol,iLay))
          nc_mp(iCol,iLay) = tracer(iCol,iLay,i_cldliq_nc) / (1.-q_lay(iCol,iLay))
          ni_mp(iCol,iLay) = tracer(iCol,iLay,i_cldice_nc) / (1.-q_lay(iCol,iLay))
          if (ltaerosol) then
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

    ! Compute effective radii for liquid/ice/snow using subgrid scale clouds
    ! Call Thompson's subroutine to compute effective radii
    do iCol=1,nCol
       call calc_effectRad (t_lay(iCol,:), p_lay(iCol,:), qv_mp(iCol,:), qc_mp(iCol,:),  &
                            nc_mp(iCol,:), qi_mp(iCol,:), ni_mp(iCol,:), qs_mp(iCol,:),  &
                            re_cloud(iCol,:), re_ice(iCol,:), re_snow(iCol,:), 1, nLev )
       do iLay = 1, nLev
          re_cloud(iCol,iLay) = MAX(re_qc_min, MIN(re_cloud(iCol,iLay), re_qc_max))
          re_ice(iCol,iLay)   = MAX(re_qi_min, MIN(re_ice(iCol,iLay),   re_qi_max))
          re_snow(iCol,iLay)  = MAX(re_qs_min, MIN(re_snow(iCol,iLay),  re_qs_max))
       end do
    enddo

    ! Scale Thompson's effective radii from meter to micron
    effrin_cldliq(1:nCol,1:nLev)  = re_cloud(1:nCol,1:nLev)*1.e6
    effrin_cldice(1:nCol,1:nLev)  = re_ice(1:nCol,1:nLev)*1.e6
    effrin_cldsnow(1:nCol,1:nLev) = re_snow(1:nCol,1:nLev)*1.e6

    ! Bound effective radii for RRTMGP, LUT's for cloud-optics go from
    !   2.5 - 21.5 microns for liquid clouds,
    !   10  - 180  microns for ice-clouds
    if (doGP_cldoptics_PADE .or. doGP_cldoptics_LUT) then
       where(effrin_cldliq .lt. radliq_lwr) effrin_cldliq = radliq_lwr
       where(effrin_cldliq .gt. radliq_upr) effrin_cldliq = radliq_upr
       where(effrin_cldice .lt. radice_lwr) effrin_cldice = radice_lwr
       where(effrin_cldice .gt. radice_upr) effrin_cldice = radice_upr
    endif

    ! Update global effective radii arrays.
    cld_reliq(1:nCol,1:nLev)      = effrin_cldliq(1:nCol,1:nLev)
    cld_reice(1:nCol,1:nLev)      = effrin_cldice(1:nCol,1:nLev)
    cld_resnow(1:nCol,1:nLev)     = effrin_cldsnow(1:nCol,1:nLev)
    cld_rerain(1:nCol,1:nLev)     = rerain_def

    ! Compute cloud-fraction. Else, use value provided
    if(.not. do_mynnedmf .or. imfdeepcnv .ne. imfdeepcnv_gf ) then ! MYNN PBL or GF conv
       ! Cloud-fraction
       if (uni_cld) then
          cld_frac(1:nCol,1:nLev) = cld_frac_mg(1:nCol,1:nLev)
       else
          if(      lmfshal) alpha0 = 100. ! Default (from GATE simulations)
          if(.not. lmfshal) alpha0 = 2000.
          ! Xu-Randall (1996) cloud-fraction
          do iLay = 1, nLev
             do iCol = 1, nCol
                cld_mr = cld_condensate(iCol,iLay,1) + cld_condensate(iCol,iLay,2) +  &
                         cld_condensate(iCol,iLay,4)
                cld_frac(iCol,iLay) = cld_frac_XuRandall(p_lay(iCol,iLay),            &
                   qs_lay(iCol,iLay), relhum(iCol,iLay), cld_mr, alpha0)
             enddo
          enddo
       endif
    endif

    ! Precipitation fraction (Hack. For now use cloud-fraction)
    precip_frac(1:nCol,1:nLev) = cld_frac(1:nCol,1:nLev)

  end subroutine GFS_rrtmgp_thompsonmp_pre_run

  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_thompsonmp_pre_finalize()
  end subroutine GFS_rrtmgp_thompsonmp_pre_finalize

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
end module GFS_rrtmgp_thompsonmp_pre
