! ###########################################################################################
! ###########################################################################################
module rrtmgp_sw_main
  use machine,                only: kind_phys, kind_dbl_prec
  use mo_optical_props,       only: ty_optical_props_2str
  use mo_cloud_optics,        only: ty_cloud_optics
  use module_radsw_parameters, only: cmpfsw_type
  use mo_rte_sw,              only: rte_sw
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,  only: ty_gas_concs
  use mo_fluxes_byband,       only: ty_fluxes_byband
  use radiation_tools,        only: check_error_msg
  use rrtmgp_sw_gas_optics,   only: sw_gas_props,rrtmgp_sw_gas_optics_init
  use rrtmgp_sw_cloud_optics, only: sw_cloud_props, rrtmgp_sw_cloud_optics_init, a0r, a0s,  &
                                    a1s, b0r, b0s, b1s, c0r, c0s
  use GFS_rrtmgp_pre,         only: iStr_h2o, iStr_co2, iStr_o3, iStr_n2o, iStr_ch4,        &
                                    iStr_o2, iStr_ccl4, iStr_cfc11, iStr_cfc12, iStr_cfc22, &
                                    eps, oneminus, ftiny
  use mersenne_twister,       only: random_setseed, random_number, random_stat
  use rrtmgp_sampling,        only: sampled_mask, draw_samples
  implicit none

  type(ty_gas_concs)          :: gas_concs
  type(ty_optical_props_2str) :: sw_optical_props_accum, sw_optical_props_aerosol_local,    &
       sw_optical_props_cloudsByBand, sw_optical_props_cnvcloudsByBand,                     &
       sw_optical_props_pblcloudsByBand, sw_optical_props_precipByBand,                     &
       sw_optical_props_clouds

  public rrtmgp_sw_main_init, rrtmgp_sw_main_run

contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_main_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_main_init
!! \htmlinclude rrtmgp_sw_main_init.html
!!
  subroutine rrtmgp_sw_main_init(rrtmgp_root_dir, rrtmgp_sw_file_gas, rrtmgp_sw_file_clouds,&
       active_gases_array, doGP_cldoptics_PADE, doGP_cldoptics_LUT, doGP_sgs_pbl,           &
       doGP_sgs_cnv, nrghice, mpicomm, mpirank, mpiroot, nLay, rrtmgp_phys_blksz,           &
       errmsg, errflg)

    ! Inputs
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,       & ! RTE-RRTMGP root directory
         rrtmgp_sw_file_clouds, & ! RRTMGP file containing K-distribution data
         rrtmgp_sw_file_gas       ! RRTMGP file containing cloud-optics data
    character(len=*), dimension(:), intent(in) :: &
         active_gases_array       ! List of active gases from namelist as array)
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

    ! RRTMGP shortwave gas-optics (k-distribution) initialization
    call rrtmgp_sw_gas_optics_init(rrtmgp_root_dir, rrtmgp_sw_file_gas, active_gases_array,&
         mpicomm, mpirank, mpiroot, errmsg, errflg)

    ! RRTMGP shortwave cloud-optics initialization
    call rrtmgp_sw_cloud_optics_init(rrtmgp_root_dir, rrtmgp_sw_file_clouds,               &
         doGP_cldoptics_PADE, doGP_cldoptics_LUT, nrghice, mpicomm, mpirank, mpiroot,      &
         errmsg, errflg)

    ! DDTs

    ! ty_gas_concs
    call check_error_msg('rrtmgp_sw_main_gas_concs_init',gas_concs%init(active_gases_array))

    ! ty_optical_props
    call check_error_msg('rrtmgp_sw_main_accumulated_optics_init',&
         sw_optical_props_accum%alloc_2str(rrtmgp_phys_blksz, nLay, sw_gas_props))
    call check_error_msg('rrtmgp_sw_main_cloud_optics_init',&
         sw_optical_props_cloudsByBand%alloc_2str(rrtmgp_phys_blksz, nLay, sw_gas_props%get_band_lims_wavenumber()))
    call check_error_msg('rrtmgp_sw_main_precip_optics_init',&
         sw_optical_props_precipByBand%alloc_2str(rrtmgp_phys_blksz, nLay, sw_gas_props%get_band_lims_wavenumber()))
    call check_error_msg('rrtmgp_sw_mian_cloud_sampling_init', &
         sw_optical_props_clouds%alloc_2str(rrtmgp_phys_blksz, nLay, sw_gas_props))
    call check_error_msg('rrtmgp_sw_main_aerosol_optics_init',&
         sw_optical_props_aerosol_local%alloc_2str(rrtmgp_phys_blksz, nLay, sw_gas_props%get_band_lims_wavenumber()))
    if (doGP_sgs_cnv) then
       call check_error_msg('rrtmgp_sw_main_cnv_cloud_optics_init',&
            sw_optical_props_cnvcloudsByBand%alloc_2str(rrtmgp_phys_blksz, nLay, sw_gas_props%get_band_lims_wavenumber()))
    endif
    if (doGP_sgs_pbl) then
       call check_error_msg('rrtmgp_sw_main_pbl_cloud_optics_init',&
            sw_optical_props_pblcloudsByBand%alloc_2str(rrtmgp_phys_blksz, nLay, sw_gas_props%get_band_lims_wavenumber()))
    endif
  end subroutine rrtmgp_sw_main_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_main_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_main_run
!! \htmlinclude rrtmgp_sw_main_run.html
!!
  subroutine rrtmgp_sw_main_run(doSWrad, doSWclrsky, top_at_1, doGP_sgs_cnv, doGP_sgs_pbl,  &
       nCol, nDay, nLay, nGases, rrtmgp_phys_blksz, idx, icseed_sw, iovr, iovr_convcld,     &
       iovr_max, iovr_maxrand, iovr_rand, iovr_dcorr, iovr_exp, iovr_exprand, isubc_sw,     &
       iSFC, sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, coszen,&
       p_lay, p_lev, t_lay, t_lev, vmr_o2, vmr_h2o, vmr_o3, vmr_ch4, vmr_n2o, vmr_co2,      &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,      &
       cld_rerain, precip_frac, cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp, cld_cnv_reice,     &
       cld_pbl_lwp, cld_pbl_reliq, cld_pbl_iwp, cld_pbl_reice, cloud_overlap_param,         &
       active_gases_array, aersw_tau, aersw_ssa, aersw_g, solcon, scmpsw,                   &
       fluxswUP_allsky, fluxswDOWN_allsky, fluxswUP_clrsky, fluxswDOWN_clrsky, cldtausw,    &
       errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doSWrad,             & ! Flag to perform shortwave calculation
         doSWclrsky,          & ! Flag to compute clear-sky fluxes
         top_at_1,            & ! Flag for vertical ordering convention
         doGP_sgs_pbl,        & ! Flag to include sgs PBL clouds
         doGP_sgs_cnv           ! Flag to include sgs convective clouds
    integer,intent(in) :: &
         nCol,                & ! Number of horizontal points
         nDay,                & ! Number of daytime points
         nLay,                & ! Number of vertical grid points.
         nGases,              & ! Number of active gases
         rrtmgp_phys_blksz,   & ! Number of horizontal points to process at once.
         iovr,                & ! Choice of cloud-overlap method
         iovr_convcld,        & ! Choice of convective cloud-overlap
         iovr_max,            & ! Flag for maximum cloud overlap method
         iovr_maxrand,        & ! Flag for maximum-random cloud overlap method
         iovr_rand,           & ! Flag for random cloud overlap method
         iovr_dcorr,          & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,            & ! Flag for exponential cloud overlap method
         iovr_exprand,        & ! Flag for exponential-random cloud overlap method
         isubc_sw,            & !
         iSFC
    integer,intent(in),dimension(:) :: &
         idx,                 & ! Index array for daytime points
         icseed_sw              ! Seed for random number generation for shortwave radiation
    real(kind_phys), dimension(:), intent(in) :: &
         sfc_alb_nir_dir,     & ! Surface albedo (direct)
         sfc_alb_nir_dif,     & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir,   & ! Surface albedo (direct)
         sfc_alb_uvvis_dif,   & ! Surface albedo (diffuse)
         coszen                 ! Cosize of SZA
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
         precip_frac,         & ! Precipitation fraction
         cld_cnv_lwp,         & ! Water path for       convective   liquid cloud-particles
         cld_cnv_reliq,       & ! Effective radius for convective   liquid cloud-particles
         cld_cnv_iwp,         & ! Water path for       convective   ice    cloud-particles
         cld_cnv_reice,       & ! Effective radius for convective   ice    cloud-particles
         cld_pbl_lwp,         & ! Water path for       PBL          liquid cloud-particles
         cld_pbl_reliq,       & ! Effective radius for PBL          liquid cloud-particles
         cld_pbl_iwp,         & ! Water path for       PBL          ice    cloud-particles
         cld_pbl_reice,       & ! Effective radius for PBL          ice    cloud-particles
         cloud_overlap_param    !
    real(kind_phys), dimension(:,:,:), intent(in) :: &
          aersw_tau,          & ! Aerosol optical depth
          aersw_ssa,          & ! Aerosol single scattering albedo
          aersw_g               ! Aerosol asymmetry paramter
    character(len=*), dimension(:), intent(in) :: &
         active_gases_array     ! List of active gases from namelist as array
    real(kind_phys), intent(in) :: &
         solcon                 ! Solar constant

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                ! CCPP error message
    integer, intent(out) :: &
         errflg                ! CCPP error flag
    real(kind_phys), dimension(:,:), intent(inout) :: &
         cldtausw              ! Approx 10.mu band layer cloud optical depth  
    real(kind_phys), dimension(:,:), intent(inout) :: &
         fluxswUP_allsky,    & ! RRTMGP upward all-sky flux profiles (W/m2)
         fluxswDOWN_allsky,  & ! RRTMGP downward all-sky flux profiles (W/m2)
         fluxswUP_clrsky,    & ! RRTMGP upward clear-sky flux profiles (W/m2)
         fluxswDOWN_clrsky     ! RRTMGP downward clear-sky flux profiles (W/m2)
    type(cmpfsw_type), dimension(:), intent(inout) :: &
         scmpsw                ! 2D surface fluxes, components:
                               ! uvbfc - total sky downward uv-b flux (W/m2)
                               ! uvbf0 - clear sky downward uv-b flux (W/m2)
                               ! nirbm - downward nir direct beam flux (W/m2)
                               ! nirdf - downward nir diffused flux (W/m2)
                               ! visbm - downward uv+vis direct beam flux (W/m2)
                               ! visdf - downward uv+vis diffused flux (W/m2)

    ! Local variables
    type(cmpfsw_type), dimension(rrtmgp_phys_blksz) :: scmpsw_clrsky, scmpsw_allsky
    type(ty_fluxes_byband)      :: flux_allsky, flux_clrsky
    real(kind_phys) :: tau_rain, tau_snow, ssa_rain, ssa_snow, asy_rain, asy_snow, &
         tau_prec, asy_prec, ssa_prec, asyw, ssaw, za1, za2, flux_dir, flux_dif
    real(kind_phys), dimension(rrtmgp_phys_blksz) :: zcf0, zcf1
    real(kind_dbl_prec), dimension(sw_gas_props%get_ngpt()) :: rng1D
    real(kind_dbl_prec), dimension(sw_gas_props%get_ngpt(),nLay,rrtmgp_phys_blksz) :: rng3D,rng3D2
    real(kind_dbl_prec), dimension(sw_gas_props%get_ngpt()*nLay) :: rng2D
    logical, dimension(rrtmgp_phys_blksz,nLay,sw_gas_props%get_ngpt()) :: maskMCICA
    logical :: cloudy_column, clear_column
    real(kind_phys), dimension(sw_gas_props%get_nband(),rrtmgp_phys_blksz) :: &
         sfc_alb_dir, sfc_alb_dif
    real(kind_phys), dimension(rrtmgp_phys_blksz,nLay+1,sw_gas_props%get_nband()),target :: &
         fluxSW_up_allsky, fluxSW_up_clrsky, fluxSW_dn_dir_clrsky, fluxSW_dn_allsky, &
         fluxSW_dn_clrsky, fluxSW_dn_dir_allsky
    integer :: iBand, ibd, ibd_uv, iCol, iGas, iLay, ix, ix2, iblck
    integer, dimension(rrtmgp_phys_blksz) :: ipseed_sw, iCols
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(2,sw_gas_props%get_nband()) :: bandlimits
    real(kind_phys), dimension(2), parameter :: &
         nIR_uvvis_bnd = (/12850,16000/), &
         uvb_bnd       = (/29000,38000/)
    real(kind_phys), dimension(rrtmgp_phys_blksz,sw_gas_props%get_ngpt()) :: toa_src_sw

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doSWrad) return

    if (nDay .gt. 0) then

       bandlimits = sw_gas_props%get_band_lims_wavenumber()
       ! ######################################################################################
       !
       ! Loop over all (daylit) columns...
       !
       ! ######################################################################################
       do iCol=1,nDay,rrtmgp_phys_blksz
          !ix  = idx(iCol)
          !ix2 = idx(iCol + rrtmgp_phys_blksz - 1)
          iCols = idx(iCol:iCol + rrtmgp_phys_blksz - 1)

          ! Create clear/cloudy indicator
          zcf0(:) = 1._kind_phys
          zcf1(:) = 1._kind_phys
          do iblck = 1, rrtmgp_phys_blksz
             do iLay=1,nLay
                zcf0(iblck) = min(zcf0(iblck), 1._kind_phys - cld_frac(iCols(iblck),iLay))
             enddo
             if (zcf0(iblck) <= ftiny)   zcf0(iblck) = 0._kind_phys
             if (zcf0(iblck) > oneminus) zcf0(iblck) = 1._kind_phys
             zcf1(iblck) = 1._kind_phys - zcf0(iblck)
          enddo
          cloudy_column = any(zcf1 .gt. eps)
          clear_column  = .true.
          if (cloudy_column) clear_column = .false.

          ! ###################################################################################
          !
          ! Initialize/reset
          !
          ! ###################################################################################
          sw_optical_props_clouds%tau             = 0._kind_phys
          sw_optical_props_clouds%ssa             = 0._kind_phys
          sw_optical_props_clouds%g               = 0._kind_phys
          sw_optical_props_accum%tau              = 0._kind_phys
          sw_optical_props_accum%ssa              = 0._kind_phys
          sw_optical_props_accum%g                = 0._kind_phys
          sw_optical_props_cloudsByBand%tau       = 0._kind_phys
          sw_optical_props_cloudsByBand%ssa       = 0._kind_phys
          sw_optical_props_cloudsByBand%g         = 0._kind_phys
          sw_optical_props_precipByBand%tau       = 0._kind_phys
          sw_optical_props_precipByBand%ssa       = 0._kind_phys
          sw_optical_props_precipByBand%g         = 0._kind_phys
          if (doGP_sgs_cnv) then
             sw_optical_props_cnvcloudsByBand%tau = 0._kind_phys
             sw_optical_props_cnvcloudsByBand%ssa = 0._kind_phys
             sw_optical_props_cnvcloudsByBand%g   = 0._kind_phys
          endif
          if (doGP_sgs_pbl) then
             sw_optical_props_pblcloudsByBand%tau = 0._kind_phys
             sw_optical_props_pblcloudsByBand%ssa = 0._kind_phys
             sw_optical_props_pblcloudsByBand%g   = 0._kind_phys
          endif
          scmpsw_clrsky= cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
          scmpsw_allsky= cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
          cldtausw = 0._kind_phys

          ! ty_fluxes_byband
          fluxSW_up_allsky            = 0._kind_phys
          fluxSW_dn_allsky            = 0._kind_phys
          fluxSW_dn_dir_allsky        = 0._kind_phys
          fluxSW_up_clrsky            = 0._kind_phys
          fluxSW_dn_clrsky            = 0._kind_phys
          flux_allsky%bnd_flux_up     => fluxSW_up_allsky
          flux_allsky%bnd_flux_dn     => fluxSW_dn_allsky
          flux_allsky%bnd_flux_dn_dir => fluxSW_dn_dir_allsky
          flux_clrsky%bnd_flux_up     => fluxSW_up_clrsky
          flux_clrsky%bnd_flux_dn     => fluxSW_dn_clrsky

          ! ###################################################################################
          !
          ! Set gas-concentrations
          !
          ! ###################################################################################
          call check_error_msg('rrtmgp_sw_main_set_vmr_o2',  &
               gas_concs%set_vmr(trim(active_gases_array(istr_o2)), vmr_o2(iCols,:)))
          call check_error_msg('rrtmgp_sw_main_set_vmr_co2', &
               gas_concs%set_vmr(trim(active_gases_array(istr_co2)),vmr_co2(iCols,:)))
          call check_error_msg('rrtmgp_sw_main_set_vmr_ch4', &
               gas_concs%set_vmr(trim(active_gases_array(istr_ch4)),vmr_ch4(iCols,:)))
          call check_error_msg('rrtmgp_sw_main_set_vmr_n2o', &
               gas_concs%set_vmr(trim(active_gases_array(istr_n2o)),vmr_n2o(iCols,:)))
          call check_error_msg('rrtmgp_sw_main_set_vmr_h2o', &
               gas_concs%set_vmr(trim(active_gases_array(istr_h2o)),vmr_h2o(iCols,:)))
          call check_error_msg('rrtmgp_sw_main_set_vmr_o3',  &
               gas_concs%set_vmr(trim(active_gases_array(istr_o3)), vmr_o3(iCols,:)))

          ! ###################################################################################
          !
          ! Compute gas-optics
          !
          ! ###################################################################################

          call check_error_msg('rrtmgp_sw_main_gas_optics',sw_gas_props%gas_optics(&
               p_lay(iCols,:),          & ! IN  - Pressure @ layer-centers (Pa)
               p_lev(iCols,:),          & ! IN  - Pressure @ layer-interfaces (Pa)
               t_lay(iCols,:),          & ! IN  - Temperature @ layer-centers (K)
               gas_concs,               & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
               sw_optical_props_accum,  & ! OUT - RRTMGP DDT: Shortwave optical properties, by
                                          !                   spectral point (tau,ssa,g)
               toa_src_sw))               ! OUT - TOA incident shortwave radiation (spectral)
          ! Scale incident flux
          do iblck = 1, rrtmgp_phys_blksz
             toa_src_sw(iblck,:) = toa_src_sw(iblck,:)*solcon / sum(toa_src_sw(iblck,:))
          enddo

          ! ###################################################################################
          !
          ! Set surface albedo
          !
          ! Use near-IR albedo for bands with wavenumbers extending to 12850cm-1
          ! Use uv-vis albedo for bands with wavenumbers greater than 16000cm-1
          ! For overlapping band, average near-IR and us-vis albedos.
          !
          ! ###################################################################################
          do iblck = 1, rrtmgp_phys_blksz
             do iBand=1,sw_gas_props%get_nband()
                if (bandlimits(1,iBand) .lt. nIR_uvvis_bnd(1)) then
                   sfc_alb_dir(iBand,iblck) = sfc_alb_nir_dir(iCols(iblck))
                   sfc_alb_dif(iBand,iblck) = sfc_alb_nir_dif(iCols(iblck))
                endif
                if (bandlimits(1,iBand) .eq. nIR_uvvis_bnd(1)) then
                   sfc_alb_dir(iBand,iblck) = 0.5_kind_phys*(sfc_alb_nir_dir(iCols(iblck)) +    &
                                                             sfc_alb_uvvis_dir(iCols(iblck)))
                   sfc_alb_dif(iBand,iblck) = 0.5_kind_phys*(sfc_alb_nir_dif(iCols(iblck)) +    &
                                                             sfc_alb_uvvis_dif(iCols(iblck)))
                   ibd = iBand
                endif
                if (bandlimits(1,iBand) .ge. nIR_uvvis_bnd(2)) then
                   sfc_alb_dir(iBand,iblck) = sfc_alb_uvvis_dir(iCols(iblck))
                   sfc_alb_dif(iBand,iblck) = sfc_alb_uvvis_dif(iCols(iblck))
                endif
                if (bandlimits(1,iBand) .eq. uvb_bnd(1)) ibd_uv = iBand
             enddo
          enddo

          ! ###################################################################################
          !
          ! Compute optics for cloud(s) and precipitation, sample clouds...
          !
          ! ###################################################################################
          if (cloudy_column) then
             ! Gridmean/mp-clouds
             call check_error_msg('rrtmgp_sw_main_cloud_optics',sw_cloud_props%cloud_optics(&
                  cld_lwp(iCols,:),                     & ! IN  - Cloud liquid water path
                  cld_iwp(iCols,:),                     & ! IN  - Cloud ice water path
                  cld_reliq(iCols,:),                   & ! IN  - Cloud liquid effective radius
                  cld_reice(iCols,:),                   & ! IN  - Cloud ice effective radius
                  sw_optical_props_cloudsByBand))         ! OUT - RRTMGP DDT: Shortwave optical properties, 
                                                          !       in each band (tau,ssa,g)
             cldtausw(iCols,:) = sw_optical_props_cloudsByBand%tau(:,:,11)
          
             ! Include convective clouds?
             if (doGP_sgs_cnv) then
                ! Compute
                call check_error_msg('rrtmgp_sw_main_cnv_cloud_optics',sw_cloud_props%cloud_optics(&
                     cld_cnv_lwp(iCols,:),              & ! IN  - Convective cloud liquid water path (g/m2)
                     cld_cnv_iwp(iCols,:),              & ! IN  - Convective cloud ice water path (g/m2)
                     cld_cnv_reliq(iCols,:),            & ! IN  - Convective cloud liquid effective radius (microns)
                     cld_cnv_reice(iCols,:),            & ! IN  - Convective cloud ice effective radius (microns)
                     sw_optical_props_cnvcloudsByBand))   ! OUT - RRTMGP DDT containing convective cloud radiative properties
                                                          !       in each band
                ! Increment
                call check_error_msg('rrtmgp_sw_main_increment_cnvclouds_to_clouds',&
                     sw_optical_props_cnvcloudsByBand%increment(sw_optical_props_cloudsByBand))
             endif

             ! Include PBL clouds?
             if (doGP_sgs_pbl) then
                ! Compute
                call check_error_msg('rrtmgp_sw_main_pbl_cloud_optics',sw_cloud_props%cloud_optics(&
                     cld_pbl_lwp(iCols,:),              & ! IN  - PBL cloud liquid water path (g/m2)
                     cld_pbl_iwp(iCols,:),              & ! IN  - PBL cloud ice water path (g/m2)
                     cld_pbl_reliq(iCols,:),            & ! IN  - PBL cloud liquid effective radius (microns)
                     cld_pbl_reice(iCols,:),            & ! IN  - PBL cloud ice effective radius (microns)
                     sw_optical_props_pblcloudsByBand))   ! OUT - RRTMGP DDT containing PBL cloud radiative properties
                                                          !       in each band
                ! Increment
                call check_error_msg('rrtmgp_sw_main_increment_pblclouds_to_clouds',&
                     sw_optical_props_pblcloudsByBand%increment(sw_optical_props_cloudsByBand))
             endif
          
             ! Cloud precipitation optics: rain and snow(+groupel)
             do iblck = 1, rrtmgp_phys_blksz
                do iLay=1,nLay
                   if (cld_frac(iCols(iblck),iLay) .gt. ftiny) then
                      ! Rain/Snow optical depth (No band dependence)
                      tau_rain = cld_rwp(iCols(iblck),iLay)*a0r
                      if (cld_swp(iCols(iblck),iLay) .gt. 0. .and. cld_resnow(iCols(iblck),iLay) .gt. 10._kind_phys) then
                         tau_snow = cld_swp(iCols(iblck),iLay)*1.09087*(a0s + a1s/(1.0315*cld_resnow(iCols(iblck),iLay)))     ! fu's formula 
                      else
                         tau_snow = 0._kind_phys
                      endif
                      
                      ! Rain/Snow single-scattering albedo and asymmetry (Band dependent)
                      do iBand=1,sw_gas_props%get_nband()
                         ! By species
                         ssa_rain = tau_rain*(1.-b0r(iBand))
                         asy_rain = ssa_rain*c0r(iBand)
                         ssa_snow = tau_snow*(1.-(b0s(iBand)+b1s(iBand)*1.0315*cld_resnow(iCols(iblck),iLay)))
                         asy_snow = ssa_snow*c0s(iBand)
                         ! Combine
                         tau_prec = max(1.e-12_kind_phys, tau_rain + tau_snow)
                         ssa_prec = max(1.e-12_kind_phys, ssa_rain + ssa_snow)
                         asy_prec = max(1.e-12_kind_phys, asy_rain + asy_snow)
                         asyw     = asy_prec/max(1.e-12_kind_phys, ssa_prec)
                         ssaw     = min(1._kind_phys-0.000001, ssa_prec/tau_prec)
                         za1      = asyw * asyw
                         za2      = ssaw * za1                      
                         sw_optical_props_precipByBand%tau(iblck,iLay,iBand) = (1._kind_phys - za2) * tau_prec
                         sw_optical_props_precipByBand%ssa(iblck,iLay,iBand) = (ssaw - za2) / (1._kind_phys - za2)
                         sw_optical_props_precipByBand%g(iblck,iLay,iBand)   = asyw/(1+asyw)
                      enddo
                   endif
                enddo
             enddo
             ! Increment 
             call check_error_msg('rrtmgp_sw_main_increment_precip_to_clouds',&
                  sw_optical_props_precipByBand%increment(sw_optical_props_cloudsByBand))
          
             ! ###################################################################################
             !
             ! Cloud-sampling
             !
             ! ###################################################################################
             ! Change random number seed value for each radiation invocation (isubc_sw =1 or 2).
             if(isubc_sw == 1) then      ! advance prescribed permutation seed
                do iblck = 1, rrtmgp_phys_blksz
                   ipseed_sw(iblck) = sw_gas_props%get_ngpt() + iCols(iblck)
                enddo
             elseif (isubc_sw == 2) then ! use input array of permutaion seeds
                do iblck = 1, rrtmgp_phys_blksz
                   ipseed_sw(iblck) = icseed_sw(iCols(iblck))
                enddo
             endif

             ! Call RNG
             do iblck = 1, rrtmgp_phys_blksz
                call random_setseed(ipseed_sw(iblck),rng_stat)
                ! Use same rng for each layer
                if (iovr == iovr_max) then
                   call random_number(rng1D,rng_stat)
                   do iLay=1,nLay
                      rng3D(:,iLay,iblck) = rng1D
                   enddo
                else
                   do iLay=1,nLay
                      call random_number(rng1D,rng_stat)
                      rng3D(:,iLay,iblck) = rng1D
                   enddo
                endif
             enddo
             
             ! Cloud-overlap.
             ! Maximum-random, random or maximum.
             if (iovr == iovr_maxrand .or. iovr == iovr_rand .or. iovr == iovr_max) then
                call sampled_mask(real(rng3D, kind=kind_phys), cld_frac(iCols,:), maskMCICA)
             endif
             ! Exponential decorrelation length overlap
             if (iovr == iovr_dcorr) then
                do iblck = 1, rrtmgp_phys_blksz
                   ! Generate second RNG
                   call random_setseed(ipseed_sw(iblck),rng_stat)
                   call random_number(rng2D,rng_stat)
                   rng3D2(:,:,iblck) = reshape(source = rng2D,shape=[sw_gas_props%get_ngpt(),nLay])
                enddo
                !
                call sampled_mask(real(rng3D, kind=kind_phys), cld_frac(iCols,:), maskMCICA,                    &
                     overlap_param = cloud_overlap_param(iCols,1:nLay-1), randoms2 = real(rng3D2, kind=kind_phys))
             endif
             ! Exponential or Exponential-random
             if (iovr == iovr_exp .or. iovr == iovr_exprand) then
                call sampled_mask(real(rng3D, kind=kind_phys), cld_frac(iCols,:), maskMCICA,  &
                     overlap_param = cloud_overlap_param(iCols,1:nLay-1))
             endif
             ! Sampling. Map band optical depth to each g-point using McICA
             call check_error_msg('rrtmgp_sw_main_cloud_sampling',&
                  draw_samples(maskMCICA, .true., &
                  sw_optical_props_cloudsByBand, sw_optical_props_clouds))
          endif ! cloudy_column

          ! ###################################################################################
          !
          ! Compute clear-sky fluxes (gaseous+aerosol)
          !
          ! ###################################################################################
          ! Increment optics (always)
          sw_optical_props_aerosol_local%tau = aersw_tau(iCols,:,:)
          sw_optical_props_aerosol_local%ssa = aersw_ssa(iCols,:,:)
          sw_optical_props_aerosol_local%g   = aersw_g(iCols,:,:)
          call check_error_msg('rrtmgp_sw_main_increment_aerosol_to_clrsky', & 
               sw_optical_props_aerosol_local%increment(sw_optical_props_accum))

          ! Compute clear-sky fluxes (Yes for no-clouds. Optional for cloudy scenes)
          if (clear_column .or. doSWclrsky) then
             call check_error_msg('rrtmgp_sw_main_rte_sw_clrsky',rte_sw(     &
                  sw_optical_props_accum,    & ! IN  - optical-properties
                  top_at_1,                  & ! IN  - veritcal ordering flag
                  coszen(iCols),             & ! IN  - Cosine of solar zenith angle
                  toa_src_sw,                & ! IN  - incident solar flux at TOA
                  sfc_alb_dir,               & ! IN  - Shortwave surface albedo (direct)
                  sfc_alb_dif,               & ! IN  - Shortwave surface albedo (diffuse)
                  flux_clrsky))                ! OUT - Fluxes, clear-sky, 3D (1,nLay,nBand) 
             
             ! Store fluxes
             fluxswUP_clrsky(iCols,:)   = sum(flux_clrsky%bnd_flux_up, dim=3)
             fluxswDOWN_clrsky(iCols,:) = sum(flux_clrsky%bnd_flux_dn, dim=3)

             ! Compute surface downward beam/diffused flux components
             do iblck = 1, rrtmgp_phys_blksz
                do iBand=1,sw_gas_props%get_nband()
                   flux_dir = flux_clrsky%bnd_flux_dn(iblck,iSFC,iBand)
                   flux_dif = 0._kind_phys
                   ! Near-IR bands
                   if (iBand < ibd) then
                      scmpsw_clrsky(iblck)%nirbm = scmpsw_clrsky(iblck)%nirbm + flux_dir
                      scmpsw_clrsky(iblck)%nirdf = scmpsw_clrsky(iblck)%nirdf + flux_dif
                   endif
                   ! Transition band
                   if (iBand == ibd) then
                      scmpsw_clrsky(iblck)%nirbm = scmpsw_clrsky(iblck)%nirbm + flux_dir*0.5_kind_phys
                      scmpsw_clrsky(iblck)%nirdf = scmpsw_clrsky(iblck)%nirdf + flux_dif*0.5_kind_phys
                      scmpsw_clrsky(iblck)%visbm = scmpsw_clrsky(iblck)%visbm + flux_dir*0.5_kind_phys
                      scmpsw_clrsky(iblck)%visdf = scmpsw_clrsky(iblck)%visdf + flux_dif*0.5_kind_phys
                   endif
                   ! UV-VIS bands
                   if (iBand > ibd) then
                      scmpsw_clrsky(iblck)%visbm = scmpsw_clrsky(iblck)%visbm + flux_dir
                      scmpsw_clrsky(iblck)%visdf = scmpsw_clrsky(iblck)%visdf + flux_dif
                   endif
                   ! uv-b surface downward flux
                   scmpsw_clrsky(iblck)%uvbfc    = flux_clrsky%bnd_flux_dn(iblck,iSFC,ibd_uv)
                enddo
             enddo
          else
             fluxswUP_clrsky(iCols,:)   = 0._kind_phys
             fluxswDOWN_clrsky(iCols,:) = 0._kind_phys
             scmpsw                     = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
          endif

          ! ###################################################################################
          !
          ! All-sky fluxes (clear-sky + clouds + precipitation)
          !
          ! ###################################################################################
          if (cloudy_column) then
             ! Delta scale
             !call check_error_msg('rrtmgp_sw_main_delta_scale',sw_optical_props_clouds%delta_scale())

             ! Increment
             call check_error_msg('rrtmgp_sw_main_increment_clouds_to_clrsky', & 
                  sw_optical_props_clouds%increment(sw_optical_props_accum))

             ! Compute fluxes
             call check_error_msg('rrtmgp_sw_main_rte_sw_allsky',rte_sw(     &
                  sw_optical_props_accum,    & ! IN  - optical-properties
                  top_at_1,                  & ! IN  - veritcal ordering flag
                  coszen(iCols),             & ! IN  - Cosine of solar zenith angle
                  toa_src_sw,                & ! IN  - incident solar flux at TOA
                  sfc_alb_dir,               & ! IN  - Shortwave surface albedo (direct)
                  sfc_alb_dif,               & ! IN  - Shortwave surface albedo (diffuse)
                  flux_allsky))                ! OUT - Fluxes, clear-sky, 3D (1,nLay,nBand)
             
             ! Store fluxes
             fluxswUP_allsky(iCols,:)   = sum(flux_allsky%bnd_flux_up, dim=3)
             fluxswDOWN_allsky(iCols,:) = sum(flux_allsky%bnd_flux_dn, dim=3)
             
             ! Compute and store downward beam/diffused flux components
             do iblck = 1, rrtmgp_phys_blksz
                ! Loop over bands, sum fluxes...
                do iBand=1,sw_gas_props%get_nband()
                   flux_dir = flux_allsky%bnd_flux_dn_dir(iblck,iSFC,iBand) 
                   flux_dif = flux_allsky%bnd_flux_dn(iblck,iSFC,iBand) - flux_allsky%bnd_flux_dn_dir(iblck,iSFC,iBand)
                   ! Near-IR bands
                   if (iBand < ibd) then
                      scmpsw_allsky(iblck)%nirbm = scmpsw_allsky(iblck)%nirbm + flux_dir
                      scmpsw_allsky(iblck)%nirdf = scmpsw_allsky(iblck)%nirdf + flux_dif
                   endif
                   ! Transition band
                   if (iBand == ibd) then
                      scmpsw_allsky(iblck)%nirbm = scmpsw_allsky(iblck)%nirbm + flux_dir*0.5_kind_phys
                      scmpsw_allsky(iblck)%nirdf = scmpsw_allsky(iblck)%nirdf + flux_dif*0.5_kind_phys
                      scmpsw_allsky(iblck)%visbm = scmpsw_allsky(iblck)%visbm + flux_dir*0.5_kind_phys
                      scmpsw_allsky(iblck)%visdf = scmpsw_allsky(iblck)%visdf + flux_dif*0.5_kind_phys
                   endif
                   ! UV-VIS bands
                   if (iBand > ibd) then
                      scmpsw_allsky(iblck)%visbm = scmpsw_allsky(iblck)%visbm + flux_dir
                      scmpsw_allsky(iblck)%visdf = scmpsw_allsky(iblck)%visdf + flux_dif
                   endif
                   ! uv-b surface downward flux 
                   scmpsw_allsky(iblck)%uvbfc    = flux_allsky%bnd_flux_dn(iblck,iSFC,ibd_uv)
                enddo
                ! Store surface downward beam/diffused flux components
                if (zcf1(iblck) .gt. eps) then
                   scmpsw(iCols(iblck))%nirbm = scmpsw_allsky(iblck)%nirbm
                   scmpsw(iCols(iblck))%nirdf = scmpsw_allsky(iblck)%nirdf
                   scmpsw(iCols(iblck))%visbm = scmpsw_allsky(iblck)%visbm
                   scmpsw(iCols(iblck))%visdf = scmpsw_allsky(iblck)%visdf
                   scmpsw(iCols(iblck))%uvbfc = flux_allsky%bnd_flux_dn(iblck,iSFC,ibd_uv)
                else
                   scmpsw(iCols(iblck))%nirbm = scmpsw_clrsky(iblck)%nirbm
                   scmpsw(iCols(iblck))%nirdf = scmpsw_clrsky(iblck)%nirdf
                   scmpsw(iCols(iblck))%visbm = scmpsw_clrsky(iblck)%visbm
                   scmpsw(iCols(iblck))%visdf = scmpsw_clrsky(iblck)%visdf
                   scmpsw(iCols(iblck))%uvbfc = flux_clrsky%bnd_flux_dn(iblck,iSFC,ibd_uv)
                endif
                scmpsw(iCols(iblck))%uvbf0    = flux_clrsky%bnd_flux_dn(iblck,iSFC,ibd_uv)
             enddo
          else ! No clouds
             fluxswUP_allsky(iCols,:)   = sum(flux_clrsky%bnd_flux_up, dim=3)
             fluxswDOWN_allsky(iCols,:) = sum(flux_clrsky%bnd_flux_dn, dim=3)
             do iblck = 1, rrtmgp_phys_blksz
                scmpsw(iCols(iblck))%nirbm = scmpsw_clrsky(iblck)%nirbm
                scmpsw(iCols(iblck))%nirdf = scmpsw_clrsky(iblck)%nirdf
                scmpsw(iCols(iblck))%visbm = scmpsw_clrsky(iblck)%visbm
                scmpsw(iCols(iblck))%visdf = scmpsw_clrsky(iblck)%visdf
                scmpsw(iCols(iblck))%uvbfc = flux_clrsky%bnd_flux_dn(iblck,iSFC,ibd_uv)
                scmpsw(iCols(iblck))%uvbf0 = flux_clrsky%bnd_flux_dn(iblck,iSFC,ibd_uv)
             enddo
          endif
          !
       enddo ! nday
    else
       fluxswUP_allsky(:,:)   = 0._kind_phys
       fluxswDOWN_allsky(:,:) = 0._kind_phys
       fluxswUP_clrsky(:,:)   = 0._kind_phys
       fluxswDOWN_clrsky(:,:) = 0._kind_phys
       scmpsw                 = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
    endif
  end subroutine rrtmgp_sw_main_run
end module rrtmgp_sw_main
