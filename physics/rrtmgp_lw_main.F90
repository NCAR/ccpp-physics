! ###########################################################################################
!> \file rrtmgp_lw_main.F90
!!
!> \defgroup rrtmgp_lw_main rrtmgp_lw_main.F90
!!
!! \brief This module contains the longwave RRTMGP radiation scheme.
!!
! ###########################################################################################
module rrtmgp_lw_main
  use machine,                only: kind_phys, kind_dbl_prec
  use mo_optical_props,       only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_cloud_optics,        only: ty_cloud_optics
  use mo_rte_lw,              only: rte_lw
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,  only: ty_gas_concs
  use mo_fluxes_byband,       only: ty_fluxes_byband
  use mo_source_functions,    only: ty_source_func_lw
  use radiation_tools,        only: check_error_msg
  use rrtmgp_lw_gas_optics,   only: lw_gas_props,rrtmgp_lw_gas_optics_init
  use rrtmgp_lw_cloud_optics, only: lw_cloud_props, rrtmgp_lw_cloud_optics_init, abssnow0,   &
                                    abssnow1, absrain
  use module_radiation_gases, only: NF_VGAS, getgases, getozn
  use GFS_rrtmgp_pre,         only: iStr_h2o, iStr_co2, iStr_o3, iStr_n2o, iStr_ch4,         &
                                    iStr_o2, iStr_ccl4, iStr_cfc11, iStr_cfc12, iStr_cfc22,  &
                                    eps, oneminus, ftiny
  use mersenne_twister,       only: random_setseed, random_number, random_stat 
  use rrtmgp_sampling,        only: sampled_mask, draw_samples
  implicit none

  type(ty_gas_concs)          :: gas_concs
  type(ty_optical_props_1scl) :: lw_optical_props_clrsky, lw_optical_props_aerosol_local
  type(ty_optical_props_2str) :: lw_optical_props_clouds, lw_optical_props_cloudsByBand,    &
       lw_optical_props_cnvcloudsByBand, lw_optical_props_pblcloudsByBand,                  &
       lw_optical_props_precipByBand
  type(ty_source_func_lw)     :: sources 

  public rrtmgp_lw_main_init, rrtmgp_lw_main_run
contains
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_main_init
!! \htmlinclude rrtmgp_lw_main_int.html
!!
!> \ingroup rrtmgp_lw_main
!!
!! \brief 
!!
!! \section rrtmgp_lw_main_init
!> @{
  ! #########################################################################################
  subroutine rrtmgp_lw_main_init(rrtmgp_root_dir, rrtmgp_lw_file_gas, rrtmgp_lw_file_clouds,&
       active_gases_array, doGP_cldoptics_PADE, doGP_cldoptics_LUT, doGP_sgs_pbl,           &
       doGP_sgs_cnv, nrghice, mpicomm, mpirank, mpiroot, nLay, rrtmgp_phys_blksz,           &
       errmsg, errflg)

    ! Inputs
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,       & ! RTE-RRTMGP root directory
         rrtmgp_lw_file_clouds, & ! RRTMGP file containing coefficients used to compute
                                  ! clouds optical properties
         rrtmgp_lw_file_gas       ! RRTMGP file containing coefficients used to compute
                                  ! gaseous optical properties
    character(len=*), dimension(:), intent(in) :: &
         active_gases_array ! List of active gases from namelist as array)
    logical, intent(in) :: &
         doGP_cldoptics_PADE,   & ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT,    & ! Use RRTMGP cloud-optics: LUTs?
         doGP_sgs_pbl,          & ! Flag to include sgs PBL clouds
         doGP_sgs_cnv             ! Flag to include sgs convective clouds 
    integer, intent(inout) :: &
         nrghice                  ! Number of ice-roughness categories
    integer,intent(in) :: &
         mpicomm,               & ! MPI communicator
         mpirank,               & ! Current MPI rank
         mpiroot,               & ! Master MPI rank
         rrtmgp_phys_blksz,     & ! Number of horizontal points to process at once.
         nLay

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                   ! CCPP error message
    integer,          intent(out) :: &
         errflg                   ! CCPP error code

    ! Initialize CCPP error handling variables 
    errmsg = ''
    errflg = 0

    ! RRTMGP longwave gas-optics (k-distribution) initialization
    call rrtmgp_lw_gas_optics_init(rrtmgp_root_dir, rrtmgp_lw_file_gas,                  &
         active_gases_array, mpicomm, mpirank, mpiroot, errmsg, errflg)

    ! RRTMGP longwave cloud-optics initialization
    call rrtmgp_lw_cloud_optics_init(rrtmgp_root_dir, rrtmgp_lw_file_clouds,             &
         doGP_cldoptics_PADE, doGP_cldoptics_LUT, nrghice, mpicomm, mpirank, mpiroot,    &
         errmsg, errflg)

    ! DDTs
    
    ! ty_gas_concs
    call check_error_msg('rrtmgp_lw_main_gas_concs_init',gas_concs%init(active_gases_array))

    ! ty_optical_props
    call check_error_msg('rrtmgp_lw_main_gas_optics_init',&
         lw_optical_props_clrsky%alloc_1scl(rrtmgp_phys_blksz, nLay, lw_gas_props))
    call check_error_msg('rrtmgp_lw_main_sources_init',&
         sources%alloc(rrtmgp_phys_blksz, nLay, lw_gas_props))
    call check_error_msg('rrtmgp_lw_main_cloud_optics_init',&
         lw_optical_props_cloudsByBand%alloc_2str(rrtmgp_phys_blksz, nLay, lw_gas_props%get_band_lims_wavenumber()))
    call check_error_msg('rrtmgp_lw_main_precip_optics_init',&
         lw_optical_props_precipByBand%alloc_2str(rrtmgp_phys_blksz, nLay, lw_gas_props%get_band_lims_wavenumber()))
    call check_error_msg('rrtmgp_lw_mian_cloud_sampling_init', &
         lw_optical_props_clouds%alloc_2str(rrtmgp_phys_blksz, nLay, lw_gas_props))
    call check_error_msg('rrtmgp_lw_main_aerosol_optics_init',&
         lw_optical_props_aerosol_local%alloc_1scl(rrtmgp_phys_blksz, nLay, lw_gas_props%get_band_lims_wavenumber()))
    if (doGP_sgs_cnv) then
       call check_error_msg('rrtmgp_lw_main_cnv_cloud_optics_init',&
            lw_optical_props_cnvcloudsByBand%alloc_2str(rrtmgp_phys_blksz, nLay, lw_gas_props%get_band_lims_wavenumber()))
    endif
    if (doGP_sgs_pbl) then
       call check_error_msg('rrtmgp_lw_main_pbl_cloud_optics_init',&
            lw_optical_props_pblcloudsByBand%alloc_2str(rrtmgp_phys_blksz, nLay, lw_gas_props%get_band_lims_wavenumber()))
    endif

  end subroutine rrtmgp_lw_main_init
!> @}
  ! ######################################################################################
!! \section arg_table_rrtmgp_lw_main_run
!! \htmlinclude rrtmgp_lw_main_run.html
!!
!> \ingroup rrtmgp_lw_main
!!
!! \brief
!!
!! \section rrtmgp_lw_main_run
!> @{
  ! ######################################################################################
  subroutine rrtmgp_lw_main_run(doLWrad, doLWclrsky, top_at_1, doGP_lwscat,              &
       use_LW_jacobian, doGP_sgs_cnv, doGP_sgs_pbl, nCol, nLay, nGases,rrtmgp_phys_blksz,&
       nGauss_angles, icseed_lw, iovr, iovr_convcld, iovr_max, iovr_maxrand, iovr_rand,  &
       iovr_dcorr, iovr_exp, iovr_exprand, isubc_lw, semis, tsfg, p_lay, p_lev, t_lay,   &
       t_lev,  vmr_o2, vmr_h2o, vmr_o3, vmr_ch4, vmr_n2o, vmr_co2,                       &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow,            &
       cld_rwp, cld_rerain, precip_frac, cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp,        &
       cld_cnv_reice, cld_pbl_lwp, cld_pbl_reliq, cld_pbl_iwp, cld_pbl_reice,            &
       cloud_overlap_param, active_gases_array, aerlw_tau, aerlw_ssa, aerlw_g,           &
       fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky, fluxlwDOWN_clrsky,           &
       fluxlwUP_jac, fluxlwUP_radtime, fluxlwDOWN_radtime, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doLWrad,            & ! Flag to perform longwave calculation
         doLWclrsky,         & ! Flag to compute clear-sky fluxes
         top_at_1,           & ! Flag for vertical ordering convention
         use_LW_jacobian,    & ! Flag to compute Jacobian of longwave surface flux
         doGP_sgs_pbl,       & ! Flag to include sgs PBL clouds
         doGP_sgs_cnv,       & ! Flag to include sgs convective clouds
         doGP_lwscat           ! Flag to include scattering in clouds
    integer,intent(in) :: &
         nCol,               & ! Number of horizontal points
         nLay,               & ! Number of vertical grid points.
         nGases,             & ! Number of active gases
         rrtmgp_phys_blksz,  & ! Number of horizontal points to process at once.
         nGauss_angles,      & ! Number of gaussian quadrature angles used
         iovr,               & ! Choice of cloud-overlap method
         iovr_convcld,       & ! Choice of convective cloud-overlap
         iovr_max,           & ! Flag for maximum cloud overlap method
         iovr_maxrand,       & ! Flag for maximum-random cloud overlap method
         iovr_rand,          & ! Flag for random cloud overlap method
         iovr_dcorr,         & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,           & ! Flag for exponential cloud overlap method
         iovr_exprand,       & ! Flag for exponential-random cloud overlap method
         isubc_lw              ! Flag for cloud-seeding (rng) for cloud-sampling
    integer,intent(in),dimension(:) :: &
         icseed_lw             ! Seed for random number generation for longwave radiation
    real(kind_phys), dimension(:), intent(in) :: &
         semis,              & ! Surface-emissivity (1)
         tsfg                  ! Skin temperature (K)
    real(kind_phys), dimension(:,:), intent(in) :: &
         p_lay,               & ! Pressure @ model layer-centers (Pa)
         t_lay,               & ! Temperature (K)
         p_lev,               & ! Pressure @ model layer-interfaces (Pa)
         t_lev,               & ! Temperature @ model levels (K)
         vmr_o2,              & ! Molar-mixing ratio oxygen
         vmr_h2o,             & ! Molar-mixing ratio water vapor
         vmr_o3,              & ! Molar-mixing ratio ozone
         vmr_ch4,             & ! Molar-mixing ratio methane
         vmr_n2o,             & ! Molar-mixing ratio nitrous oxide
         vmr_co2,             & ! Molar-mixing ratio carbon dioxide
         cld_frac,            & ! Cloud-fraction for   stratiform   clouds
         cld_lwp,             & ! Water path for       stratiform   liquid cloud-particles
         cld_reliq,           & ! Effective radius for stratiform   liquid cloud-particles
         cld_iwp,             & ! Water path for       stratiform   ice    cloud-particles
         cld_reice,           & ! Effective radius for stratiform   ice    cloud-particles
         cld_swp,             & ! Water path for                    snow   hydrometeors
         cld_resnow,          & ! Effective radius for              snow   hydrometeors
         cld_rwp,             & ! Water path for                    rain   hydrometeors
         cld_rerain,          & ! Effective radius for              rain   hydrometeors
         precip_frac,         & ! Precipitation fraction (not active, currently precipitation optics uses cloud-fraction)
         cld_cnv_lwp,         & ! Water path for       convective   liquid cloud-particles
         cld_cnv_reliq,       & ! Effective radius for convective   liquid cloud-particles
         cld_cnv_iwp,         & ! Water path for       convective   ice    cloud-particles
         cld_cnv_reice,       & ! Effective radius for convective   ice    cloud-particles
         cld_pbl_lwp,         & ! Water path for       PBL          liquid cloud-particles
         cld_pbl_reliq,       & ! Effective radius for PBL          liquid cloud-particles
         cld_pbl_iwp,         & ! Water path for       PBL          ice    cloud-particles
         cld_pbl_reice,       & ! Effective radius for PBL          ice    cloud-particles
         cloud_overlap_param    ! Cloud overlap parameter
    real(kind_phys), dimension(:,:,:), intent(in) :: &
          aerlw_tau,          & ! Aerosol optical depth
          aerlw_ssa,          & ! Aerosol single scattering albedo
          aerlw_g               ! Aerosol asymmetry paramter
    character(len=*), dimension(:), intent(in) :: &
         active_gases_array     ! List of active gases from namelist as array

    ! Outputs
    real(kind_phys), dimension(:,:), intent(inout) :: &
         fluxlwUP_jac,        & ! Jacobian of upwelling LW surface radiation (W/m2/K) 
         fluxlwUP_allsky,     & ! All-sky flux (W/m2)
         fluxlwDOWN_allsky,   & ! All-sky flux (W/m2)
         fluxlwUP_clrsky,     & ! Clear-sky flux (W/m2)
         fluxlwDOWN_clrsky,   & ! All-sky flux (W/m2)
         fluxlwUP_radtime,    & ! Copy of fluxes (Used for coupling)
         fluxlwDOWN_radtime     !
    character(len=*), intent(out) :: & 
         errmsg                 ! CCPP error message
    integer, intent(out) :: & 
         errflg                 ! CCPP error flag

    ! Local variables
    type(ty_fluxes_byband) :: flux_allsky, flux_clrsky
    integer :: iCol, iLay, iGas, iBand, iCol2, ix, iblck
    integer, dimension(rrtmgp_phys_blksz) :: ipseed_lw
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(rrtmgp_phys_blksz) :: zcf0, zcf1
    logical, dimension(rrtmgp_phys_blksz,nLay,lw_gas_props%get_ngpt()) :: maskMCICA
    real(kind_phys), dimension(rrtmgp_phys_blksz) :: tau_rain, tau_snow
    real(kind_dbl_prec), dimension(lw_gas_props%get_ngpt()) :: rng1D
    real(kind_dbl_prec), dimension(lw_gas_props%get_ngpt(),nLay,rrtmgp_phys_blksz) :: rng3D,rng3D2
    real(kind_dbl_prec), dimension(lw_gas_props%get_ngpt()*nLay) :: rng2D
    real(kind_phys), dimension(rrtmgp_phys_blksz,nLay+1,lw_gas_props%get_nband()),target :: &
         fluxLW_up_allsky, fluxLW_up_clrsky, fluxLW_dn_allsky, fluxLW_dn_clrsky
    real(kind_phys), dimension(rrtmgp_phys_blksz,lw_gas_props%get_ngpt()) :: lw_Ds
    real(kind_phys), dimension(lw_gas_props%get_nband(),rrtmgp_phys_blksz) :: sfc_emiss_byband

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return

    ! ######################################################################################
    !
    ! Loop over all columns...
    !
    ! ###################################################################################### 
    do iCol=1,nCol,rrtmgp_phys_blksz
       iCol2 = iCol + rrtmgp_phys_blksz - 1

       ! Initialize/reset

       ! ty_optical_props
       lw_optical_props_clrsky%tau       = 0._kind_phys
       lw_optical_props_precipByBand%tau = 0._kind_phys
       lw_optical_props_precipByBand%ssa = 0._kind_phys
       lw_optical_props_precipByBand%g   = 0._kind_phys
       lw_optical_props_cloudsByBand%tau = 0._kind_phys
       lw_optical_props_cloudsByBand%ssa = 0._kind_phys
       lw_optical_props_cloudsByBand%g   = 0._kind_phys
       lw_optical_props_clouds%tau       = 0._kind_phys
       lw_optical_props_clouds%ssa       = 0._kind_phys
       lw_optical_props_clouds%g         = 0._kind_phys
       sources%sfc_source                = 0._kind_phys
       sources%lay_source                = 0._kind_phys
       sources%lev_source_inc            = 0._kind_phys
       sources%lev_source_dec            = 0._kind_phys
       sources%sfc_source_Jac            = 0._kind_phys
       fluxLW_up_allsky                  = 0._kind_phys
       fluxLW_dn_allsky                  = 0._kind_phys
       fluxLW_up_clrsky                  = 0._kind_phys
       fluxLW_dn_clrsky                  = 0._kind_phys
       if (doGP_sgs_cnv) lw_optical_props_cnvcloudsByBand%tau = 0._kind_phys
       if (doGP_sgs_pbl) lw_optical_props_pblcloudsByBand%tau = 0._kind_phys

       ! ty_fluxes_byband
       fluxLW_up_allsky        = 0._kind_phys
       fluxLW_dn_allsky        = 0._kind_phys
       fluxLW_up_clrsky        = 0._kind_phys
       fluxLW_dn_clrsky        = 0._kind_phys
       flux_allsky%bnd_flux_up => fluxLW_up_allsky
       flux_allsky%bnd_flux_dn => fluxLW_dn_allsky
       flux_clrsky%bnd_flux_up => fluxLW_up_clrsky
       flux_clrsky%bnd_flux_dn => fluxLW_dn_clrsky

       ! ###################################################################################
       !
       ! Set gas-concentrations
       !
       ! ###################################################################################
       call check_error_msg('rrtmgp_lw_main_set_vmr_o2',  &
            gas_concs%set_vmr(trim(active_gases_array(istr_o2)), vmr_o2(iCol:iCol2,:)))
       call check_error_msg('rrtmgp_lw_main_set_vmr_co2', &
            gas_concs%set_vmr(trim(active_gases_array(istr_co2)),vmr_co2(iCol:iCol2,:)))
       call check_error_msg('rrtmgp_lw_main_set_vmr_ch4', &
            gas_concs%set_vmr(trim(active_gases_array(istr_ch4)),vmr_ch4(iCol:iCol2,:)))
       call check_error_msg('rrtmgp_lw_main_set_vmr_n2o', &
            gas_concs%set_vmr(trim(active_gases_array(istr_n2o)),vmr_n2o(iCol:iCol2,:)))
       call check_error_msg('rrtmgp_lw_main_set_vmr_h2o', &
            gas_concs%set_vmr(trim(active_gases_array(istr_h2o)),vmr_h2o(iCol:iCol2,:)))
       call check_error_msg('rrtmgp_lw_main_set_vmr_o3',  &
            gas_concs%set_vmr(trim(active_gases_array(istr_o3)), vmr_o3(iCol:iCol2,:)))

       ! ###################################################################################
       !
       ! Surface emissity in each band
       !
       ! ###################################################################################
       ! Assign same emissivity to all band
       do iblck=1,rrtmgp_phys_blksz
          if (semis(iCol+iblck-1) > eps .and. semis(iCol+iblck-1) <= 1._kind_phys) then
             do iBand=1,lw_gas_props%get_nband()
                sfc_emiss_byband(iBand,iblck) = semis(iCol+iblck-1)
             enddo
          else
             sfc_emiss_byband(1:lw_gas_props%get_nband(),iblck) = 1.0
          endif
       enddo

       ! ###################################################################################
       !
       ! Compute gas-optics...
       !
       ! ###################################################################################
       call check_error_msg('rrtmgp_lw_main_gas_optics',lw_gas_props%gas_optics(&
            p_lay(iCol:iCol2,:),              & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(iCol:iCol2,:),              & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(iCol:iCol2,:),              & ! IN  - Temperature @ layer-centers (K)
            tsfg(iCol:iCol2),                 & ! IN  - Skin-temperature (K)
            gas_concs,                        & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
            lw_optical_props_clrsky,          & ! OUT - RRTMGP DDT: longwave optical properties
            sources,                          & ! OUT - RRTMGP DDT: source functions
            tlev=t_lev(iCol:iCol2,:)))          ! IN  - Temperature @ layer-interfaces (K) (optional)

       ! ###################################################################################
       !
       ! Compute cloud-optics...
       !
       ! ###################################################################################
       ! Create clear/cloudy indicator
       zcf0(:) = 1._kind_phys
       zcf1(:) = 1._kind_phys
       do iblck = 1, rrtmgp_phys_blksz
          do iLay=1,nLay
             zcf0(iblck) = min(zcf0(iblck), 1._kind_phys - cld_frac(iCol+iblck-1,iLay))
          enddo
          if (zcf0(iblck) <= ftiny)   zcf0(iblck) = 0._kind_phys
          if (zcf0(iblck) > oneminus) zcf0(iblck) = 1._kind_phys
          zcf1(iblck) = 1._kind_phys - zcf0(iblck)
       enddo

       if (any(zcf1 .gt. eps)) then
          ! Microphysical (gridmean) cloud optics
          call check_error_msg('rrtmgp_lw_main_cloud_optics',lw_cloud_props%cloud_optics(&
               cld_lwp(iCol:iCol2,:),                & ! IN  - Cloud liquid water path (g/m2)
               cld_iwp(iCol:iCol2,:),                & ! IN  - Cloud ice water path (g/m2)
               cld_reliq(iCol:iCol2,:),              & ! IN  - Cloud liquid effective radius (microns)
               cld_reice(iCol:iCol2,:),              & ! IN  - Cloud ice effective radius (microns)
               lw_optical_props_cloudsByBand))         ! OUT - RRTMGP DDT containing cloud radiative properties
                                                       !       in each band
          ! Include convective (subgrid scale) clouds?
          if (doGP_sgs_cnv) then
             ! Compute
             call check_error_msg('rrtmgp_lw_main_cnv_cloud_optics',lw_cloud_props%cloud_optics(&
                  cld_cnv_lwp(iCol:iCol2,:),         & ! IN  - Convective cloud liquid water path (g/m2)
                  cld_cnv_iwp(iCol:iCol2,:),         & ! IN  - Convective cloud ice water path (g/m2)
                  cld_cnv_reliq(iCol:iCol2,:),       & ! IN  - Convective cloud liquid effective radius (microns)
                  cld_cnv_reice(iCol:iCol2,:),       & ! IN  - Convective cloud ice effective radius (microns)
                  lw_optical_props_cnvcloudsByBand))   ! OUT - RRTMGP DDT containing convective cloud radiative properties
                                                       !       in each band
             ! Increment
             call check_error_msg('rrtmgp_lw_main_increment_cnvclouds_to_clouds',&
                  lw_optical_props_cnvcloudsByBand%increment(lw_optical_props_cloudsByBand))
          endif

          ! Include PBL (subgrid scale) clouds?
          if (doGP_sgs_pbl) then
             ! Compute
             call check_error_msg('rrtmgp_lw_main_pbl_cloud_optics',lw_cloud_props%cloud_optics(&
                  cld_pbl_lwp(iCol:iCol2,:),         & ! IN  - PBL cloud liquid water path (g/m2)
                  cld_pbl_iwp(iCol:iCol2,:),         & ! IN  - PBL cloud ice water path (g/m2)
                  cld_pbl_reliq(iCol:iCol2,:),       & ! IN  - PBL cloud liquid effective radius (microns)
                  cld_pbl_reice(iCol:iCol2,:),       & ! IN  - PBL cloud ice effective radius (microns)
                  lw_optical_props_pblcloudsByBand))   ! OUT - RRTMGP DDT containing PBL cloud radiative properties
                                                       !       in each band
             ! Increment
             call check_error_msg('rrtmgp_lw_main_increment_pblclouds_to_clouds',&
                  lw_optical_props_pblcloudsByBand%increment(lw_optical_props_cloudsByBand))
          endif
       endif

       ! ###################################################################################
       !
       ! Cloud precipitation optics: rain and snow(+groupel)
       !
       ! ###################################################################################
       tau_rain(:) = 0._kind_phys
       tau_snow(:) = 0._kind_phys
       do ix=1,rrtmgp_phys_blksz
          do iLay=1,nLay
             if (cld_frac(iCol+ix-1,iLay) .gt. eps) then
                ! Rain optical-depth (No band dependence)
                tau_rain(ix) = absrain*cld_rwp(iCol+ix-1,iLay)
                
                ! Snow (+groupel) optical-depth (No band dependence)
                if (cld_swp(iCol+ix-1,iLay) .gt. 0. .and. cld_resnow(iCol+ix-1,iLay) .gt. 10._kind_phys) then
                   tau_snow(ix) = abssnow0*1.05756*cld_swp(iCol+ix-1,iLay)/cld_resnow(iCol+ix-1,iLay)
                else
                   tau_snow(ix) = 0.0
                endif
                do iBand=1,lw_gas_props%get_nband()
                   lw_optical_props_precipByBand%tau(ix,iLay,iBand) = tau_rain(ix) + tau_snow(ix)
                enddo
             endif
          enddo
       enddo
       ! Increment
       call check_error_msg('rrtmgp_lw_main_increment_precip_to_clouds',&
            lw_optical_props_precipByBand%increment(lw_optical_props_cloudsByBand))

       ! ###################################################################################
       !
       ! Cloud-sampling
       ! *Note* All of the included cloud-types are sampled together, not independently.
       !
       ! ###################################################################################
       if (any(zcf1 .gt. eps)) then
          ! Change random number seed value for each radiation invocation (isubc_lw =1 or 2).
          if(isubc_lw == 1) then      ! advance prescribed permutation seed
             do ix=1,rrtmgp_phys_blksz
                ipseed_lw(ix) = lw_gas_props%get_ngpt() + iCol + ix - 1
             enddo
          elseif (isubc_lw == 2) then ! use input array of permutaion seeds
             do ix=1,rrtmgp_phys_blksz
                ipseed_lw(ix) = icseed_lw(iCol+ix-1)
             enddo
          endif

          ! Call RNG
          do ix=1,rrtmgp_phys_blksz
             call random_setseed(ipseed_lw(ix),rng_stat)
             ! Use same rng for each layer
             if (iovr == iovr_max) then
                call random_number(rng1D,rng_stat)
                do iLay=1,nLay
                   rng3D(:,iLay,ix) = rng1D
                enddo
             else
                do iLay=1,nLay
                   call random_number(rng1D,rng_stat)
                   rng3D(:,iLay,ix) = rng1D
                enddo
             endif
          enddo

          ! Cloud-overlap.
          ! Maximum-random, random or maximum.
          if (iovr == iovr_maxrand .or. iovr == iovr_rand .or. iovr == iovr_max) then
             call sampled_mask(real(rng3D,kind=kind_phys), cld_frac(iCol:iCol2,:), maskMCICA)
          endif
          ! Exponential decorrelation length overlap
          if (iovr == iovr_dcorr) then
             do ix=1,rrtmgp_phys_blksz
                ! Generate second RNG
                call random_setseed(ipseed_lw(ix),rng_stat)
                call random_number(rng2D,rng_stat)
                rng3D2(:,:,ix) = reshape(source = rng2D,shape=[lw_gas_props%get_ngpt(),nLay])
             enddo
             !
             call sampled_mask(real(rng3D,kind=kind_phys), cld_frac(iCol:iCol2,:), maskMCICA,                    &
                  overlap_param = cloud_overlap_param(iCol:iCol2,1:nLay-1), randoms2 = real(rng3D2, kind=kind_phys))
          endif
          ! Exponential or Exponential-random
          if (iovr == iovr_exp .or. iovr == iovr_exprand) then
             call sampled_mask(real(rng3D,kind=kind_phys), cld_frac(iCol:iCol2,:), maskMCICA,  &
                  overlap_param = cloud_overlap_param(iCol:iCol2,1:nLay-1))
          endif
          ! Sampling. Map band optical depth to each g-point using McICA
          call check_error_msg('rrtmgp_lw_main_cloud_sampling',&
               draw_samples(maskMCICA, .true., &
               lw_optical_props_cloudsByBand, lw_optical_props_clouds))
       endif

       ! ###################################################################################
       !
       ! Compute clear-sky fluxes (gaseous+aerosol) (optional)
       !
       ! ###################################################################################
       ! Increment
       lw_optical_props_aerosol_local%tau = aerlw_tau(iCol:iCol2,:,:)
       call check_error_msg('rrtmgp_lw_main_increment_aerosol_to_clrsky',&
            lw_optical_props_aerosol_local%increment(lw_optical_props_clrsky))

       ! Call RTE solver
       if (doLWclrsky) then
          call check_error_msg('rrtmgp_lw_main_opt_angle',&
               lw_gas_props%compute_optimal_angles(lw_optical_props_clrsky,lw_Ds))
          if (nGauss_angles .gt. 1) then
             call check_error_msg('rrtmgp_lw_main_lw_rte_clrsky',rte_lw(           &
                  lw_optical_props_clrsky,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources,                         & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_clrsky,                     & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles))   ! IN  - Number of angles in Gaussian quadrature
          else
             call check_error_msg('rrtmgp_lw_main_lw_rte_clrsky',rte_lw(           &
                  lw_optical_props_clrsky,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources,                         & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_clrsky,                     & ! OUT - Fluxes
                  lw_Ds = lw_Ds))
          endif
          
          ! Store fluxes
          fluxlwUP_clrsky(iCol:iCol2,:)   = sum(flux_clrsky%bnd_flux_up, dim=3)
          fluxlwDOWN_clrsky(iCol:iCol2,:) = sum(flux_clrsky%bnd_flux_dn, dim=3)
       else
          fluxlwUP_clrsky(iCol:iCol2,:)   = 0.0
          fluxlwDOWN_clrsky(iCol:iCol2,:) = 0.0   
       endif

       ! ###################################################################################
       !
       ! All-sky fluxes (clear-sky + clouds + precipitation)
       ! *Note* CCPP does not allow for polymorphic types, they are ambiguous to the CCPP
       ! framework. rte-rrtmgp uses polymorphic types extensively, for example, querying the
       ! type to determine physics configuration/pathway/etc...
       !
       ! The logic in the code below is to satisfy the polymorphishm in the rte-rrtmgp code.
       ! The rte-rrtmgp "increment" procedures are utilized to provide the correct type to the
       ! rte solver (rte_lw). Rte_lw quieries the type determine if scattering is to be 
       ! included in the calculation. The increment procedures are called so that the correct
       ! optical properties are inherited. ugh...
       ! 
       ! ###################################################################################

       ! Include LW cloud-scattering?
       if (doGP_lwscat) then 
          ! Increment
          call check_error_msg('rrtmgp_lw_main_increment_clrsky_to_clouds',&
               lw_optical_props_clrsky%increment(lw_optical_props_clouds))
          
          if (use_LW_jacobian) then
             ! Compute LW Jacobians
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky',rte_lw(           &
                  lw_optical_props_clouds,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources,                         & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky,                     & ! OUT - Flxues 
                  n_gauss_angles = nGauss_angles,  & ! IN  - Number of angles in Gaussian quadrature
                  flux_up_Jac    = fluxlwUP_jac))    ! OUT - surface temperature flux (upward) Jacobian (W/m2/K)
          else
             call check_error_msg('rrtmgp_lw_main_lw_rte_allsky',rte_lw(           &
                  lw_optical_props_clouds,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources,                         & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky,                     & ! OUT - Flxues 
                  n_gauss_angles = nGauss_angles))   ! IN  - Number of angles in Gaussian quadrature    
          end if
       ! No scattering in LW clouds.   
       else
          ! Increment
          call check_error_msg('rrtmgp_lw_main_increment_clouds_to_clrsky', &
               lw_optical_props_clouds%increment(lw_optical_props_clrsky))
          
          if (use_LW_jacobian) then
             ! Compute LW Jacobians
             call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
                  lw_optical_props_clrsky,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources,                         & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky,                     & ! OUT - Flxues 
                  n_gauss_angles = nGauss_angles,  & ! IN  - Number of angles in Gaussian quadrature
                  flux_up_Jac    = fluxlwUP_jac))    ! OUT - surface temperature flux (upward) Jacobian (W/m2/K)
          else
             call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
                  lw_optical_props_clrsky,         & ! IN  - optical-properties
                  top_at_1,                        & ! IN  - veritcal ordering flag
                  sources,                         & ! IN  - source function
                  sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
                  flux_allsky,                     & ! OUT - Flxues 
                  n_gauss_angles = nGauss_angles))   ! IN  - Number of angles in Gaussian quadrature    
          end if
       endif
       
       ! Store fluxes
       fluxlwUP_allsky(iCol:iCol2,:)   = sum(flux_allsky%bnd_flux_up, dim=3)
       fluxlwDOWN_allsky(iCol:iCol2,:) = sum(flux_allsky%bnd_flux_dn, dim=3)
       
       ! Save fluxes for coupling
       fluxlwUP_radtime(iCol:iCol2,:)   = fluxlwUP_allsky(iCol:iCol2,:)
       fluxlwDOWN_radtime(iCol:iCol2,:) = fluxlwDOWN_allsky(iCol:iCol2,:)

    enddo

  end subroutine rrtmgp_lw_main_run
!> @}
end module rrtmgp_lw_main
