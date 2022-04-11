! ###########################################################################################
! ###########################################################################################
module rrtmgp_lw_main
  use machine,                only: kind_phys
  use mo_optical_props,       only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_cloud_optics,        only: ty_cloud_optics
  use mo_rte_lw,              only: rte_lw
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,  only: ty_gas_concs
  use mo_fluxes_byband,       only: ty_fluxes_byband
  use mo_source_functions,    only: ty_source_func_lw
  use radiation_tools,        only: check_error_msg
  use rrtmgp_lw_gas_optics,   only: lw_gas_props,rrtmgp_lw_gas_optics_init
  use rrtmgp_lw_cloud_optics, only: lw_cloud_props, rrtmgp_lw_cloud_optics_init, abssnow0,  &
                                    abssnow1,absrain
  use module_radiation_gases, only: NF_VGAS, getgases, getozn
  use GFS_rrtmgp_pre,         only: iStr_h2o, iStr_co2, iStr_o3, iStr_n2o, iStr_ch4,        &
                                    iStr_o2, iStr_ccl4, iStr_cfc11, iStr_cfc12, iStr_cfc22 
  use mersenne_twister,       only: random_setseed, random_number, random_stat 
  use rrtmgp_sampling,        only: sampled_mask, draw_samples
  implicit none

  public rrtmgp_lw_main_init, rrtmgp_lw_main_run
contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_main_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_main_init
!! \htmlinclude rrtmgp_lw_main_int.html
!!
  subroutine rrtmgp_lw_main_init(rrtmgp_root_dir, rrtmgp_lw_file_gas, mpicomm, mpirank,     &
       mpiroot, minGPpres, maxGPpres, minGPtemp, maxGPtemp, active_gases_array, nrghice,    &
       doG_cldoptics, doGP_cldoptics_PADE, doGP_cldoptics_LUT,rrtmgp_lw_file_clouds, errmsg,&
       errflg)

    ! Inputs
    logical, intent(in) :: &
         doG_cldoptics,         & ! Use legacy RRTMG cloud-optics?
         doGP_cldoptics_PADE,   & ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT       ! Use RRTMGP cloud-optics: LUTs?
    integer, intent(inout) :: &
         nrghice                  ! Number of ice-roughness categories
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,       & ! RTE-RRTMGP root directory
         rrtmgp_lw_file_clouds, & ! RRTMGP file containing coefficients used to compute clouds optical properties
         rrtmgp_lw_file_gas       ! RRTMGP file containing coefficients used to compute gaseous optical properties
    integer,intent(in) :: &
         mpicomm,               & ! MPI communicator
         mpirank,               & ! Current MPI rank
         mpiroot                  ! Master MPI rank
    character(len=*), dimension(:), intent(in) :: &
         active_gases_array ! List of active gases from namelist as array)
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                   ! CCPP error message
    integer,          intent(out) :: &
         errflg                   ! CCPP error code
    real(kind_phys), intent(out) :: &
         minGPtemp,             & ! Minimum temperature allowed by RRTMGP.
         maxGPtemp,             & ! Maximum ...
         minGPpres,             & ! Minimum pressure allowed by RRTMGP. 
         maxGPpres                ! Maximum pressure allowed by RRTMGP. 

    ! Initialize CCPP error handling variables 
    errmsg = ''
    errflg = 0

    ! RRTMGP longwave gas-optics (k-distribution) initialization
    call rrtmgp_lw_gas_optics_init(rrtmgp_root_dir, rrtmgp_lw_file_gas, mpicomm, mpirank,   &
         mpiroot, minGPpres, maxGPpres, minGPtemp, maxGPtemp, active_gases_array, errmsg,   &
         errflg)

    ! RRTMGP longwave cloud-optics initialization
    call  rrtmgp_lw_cloud_optics_init(nrghice, mpicomm, mpirank, mpiroot, doG_cldoptics,    &
         doGP_cldoptics_PADE, doGP_cldoptics_LUT, rrtmgp_root_dir, rrtmgp_lw_file_clouds,   &
         errmsg, errflg)

  end subroutine rrtmgp_lw_main_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_main_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_main_run
!! \htmlinclude rrtmgp_lw_main_run.html
!!
  subroutine rrtmgp_lw_main_run(doLWrad, doLWclrsky, top_at_1, doGP_lwscat, use_LW_jacobian,&
       doGP_sgs_cnv, doGP_sgs_pbl, nCol, nLay, nGases, nGauss_angles, i_o3, icseed_lw, iovr,&
       iovr_convcld, iovr_max, iovr_maxrand, iovr_rand, iovr_dcorr, iovr_exp, iovr_exprand, &
       isubc_lw, semis, tsfg, p_lay, p_lev, t_lay, t_lev, vmr_o2, vmr_h2o, vmr_o3, vmr_ch4, &
       vmr_n2o, vmr_co2, cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp,         &
       cld_resnow, cld_rwp, cld_rerain, precip_frac, cld_cnv_lwp, cld_cnv_reliq,            &
       cld_cnv_iwp, cld_cnv_reice, cld_pbl_lwp, cld_pbl_reliq, cld_pbl_iwp, cld_pbl_reice,  &
       cloud_overlap_param, active_gases_array, lw_optical_props_aerosol,                   &
       fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky, fluxlwDOWN_clrsky, fluxlwUP_jac,&
       fluxlwUP_radtime, fluxlwDOWN_radtime, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doLWrad,             & ! Flag to calculate LW irradiances
         doLWclrsky,          & ! Flag to compute clear-sky fluxes (diagnostic)
         top_at_1,            & ! Vertical ordering flag
         use_LW_jacobian,     & ! Compute Jacobian of LW to update radiative fluxes between radiation calls?
         doGP_sgs_pbl,        & ! Flag for sgs MYNN-EDMF PBL cloud scheme
         doGP_sgs_cnv,        & ! Flag for sgs convective cloud scheme
         doGP_lwscat            ! Include scattering in LW cloud-optics?
    integer,intent(in) :: &
         nCol,                & ! Number of horizontal points
         nLay,                & ! Number of vertical grid points.
         nGases,              & ! Number of active gases in RRTMGP
         nGauss_angles,       & !
         i_o3,                & !
         iovr,                & ! Choice of cloud-overlap method
         iovr_convcld,        & ! Choice of convective cloud-overlap
         iovr_max,            & ! Flag for maximum cloud overlap method
         iovr_maxrand,        & ! Flag for maximum-random cloud overlap method
         iovr_rand,           & ! Flag for random cloud overlap method
         iovr_dcorr,          & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,            & ! Flag for exponential cloud overlap method
         iovr_exprand,        & ! Flag for exponential-random cloud overlap method
         isubc_lw               !
    integer,intent(in),dimension(:) :: &
         icseed_lw              ! Seed for random number generation for longwave radiation
    real(kind_phys), dimension(:), intent(in) :: &
         semis,              & ! Surface-emissivity
         tsfg                  !
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
         cld_pbl_lwp,         & ! Water path for       SGS PBL liquid cloud-particles
         cld_pbl_reliq,       & ! Effective radius for SGS PBL liquid cloud-particles
         cld_pbl_iwp,         & ! Water path for       SGS PBL ice    cloud-particles
         cld_pbl_reice,       & ! Effective radius for SGS PBL ice    cloud-particles
         cloud_overlap_param
    character(len=*), dimension(:), intent(in) :: &
         active_gases_array     ! List of active gases from namelist as array
    type(ty_optical_props_1scl),intent(inout) :: &
         lw_optical_props_aerosol ! RRTMGP DDT: Longwave aerosol optical properties (tau)

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
    type(ty_gas_concs) :: &
         gas_concentrations                   ! RRTMGP DDT: trace gas concentrations (vmr)
    type(ty_optical_props_1scl) :: &
         lw_optical_props_clrsky,           & ! RRTMGP DDT: longwave clear-sky radiative properties
         lw_optical_props_aerosol_local,    & ! RRTMGP DDT: longwave aerosol radiative properties
         lw_optical_props_cloudsByBand,     & ! RRTMGP DDT: Longwave optical properties in each band (clouds)
         lw_optical_props_cnvcloudsByBand,  & ! RRTMGP DDT: Longwave optical properties in each band (convective cloud)
         lw_optical_props_pblcloudsByBand,  & ! RRTMGP DDT: Longwave optical properties in each band (PBL cloud)
         lw_optical_props_precipByBand        ! RRTMGP DDT: Longwave optical properties in each band (precipitation)
    type(ty_optical_props_2str)  :: &
         lw_optical_props_clouds              ! RRTMGP DDT: Longwave optical properties in each band (sampled clouds)
    type(ty_source_func_lw)  :: &
         sources                              ! RRTMGP DDT: longwave source functions
    type(ty_fluxes_byband) :: &
         flux_allsky, flux_clrsky             ! RRTMGP DDT: Longwave flux profiles
    integer :: iCol, iLay, iGas, iBand, ipseed_lw
    type(random_stat) :: rng_stat
    real(kind_phys) :: tau_rain, tau_snow
    real(kind_phys), dimension(lw_gas_props%get_ngpt()) :: rng1D
    real(kind_phys), dimension(lw_gas_props%get_ngpt(),nLay,1) :: rng3D,rng3D2
    real(kind_phys), dimension(lw_gas_props%get_ngpt()*nLay) :: rng2D
    logical, dimension(1,nLay,lw_gas_props%get_ngpt()) :: maskMCICA
    real(kind_phys), dimension(1,nLay+1,lw_gas_props%get_nband()),target :: &
         fluxLW_up_allsky, fluxLW_up_clrsky, fluxLW_dn_allsky, fluxLW_dn_clrsky
    real(kind_phys), dimension(1,lw_gas_props%get_ngpt()) :: lw_Ds
    real(kind_phys), dimension(lw_gas_props%get_nband(),1) :: sfc_emiss_byband

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return

    ! ######################################################################################
    !
    ! Allocate/initialize RRTMGP DDT's
    !
    ! ######################################################################################
    !
    ! ty_gas_concs
    !
    gas_concentrations%ncol = 1
    gas_concentrations%nlay = nLay
    allocate(gas_concentrations%gas_name(nGases))
    allocate(gas_concentrations%concs(nGases))
    do iGas=1,nGases
       allocate(gas_concentrations%concs(iGas)%conc(1, nLay))
    enddo
    gas_concentrations%gas_name(:) = active_gases_array(:)
    !
    ! ty_optical_props
    !
    call check_error_msg('rrtmgp_lw_main_gas_optics_init',&
         lw_optical_props_clrsky%alloc_1scl(1, nLay, lw_gas_props))
    call check_error_msg('rrtmgp_lw_main_sources_init',&
         sources%alloc(1, nLay, lw_gas_props))
    call check_error_msg('rrtmgp_lw_main_cloud_optics_init',&
         lw_optical_props_cloudsByBand%alloc_1scl(1, nLay, lw_gas_props%get_band_lims_wavenumber()))
    call check_error_msg('rrtmgp_lw_main_precip_optics_init',&
         lw_optical_props_precipByBand%alloc_1scl(1, nLay, lw_gas_props%get_band_lims_wavenumber()))
    call check_error_msg('rrtmgp_lw_mian_cloud_sampling_init', & 
         lw_optical_props_clouds%alloc_2str(1, nLay, lw_gas_props))
    call check_error_msg('rrtmgp_lw_main_aerosol_optics_init',&
         lw_optical_props_aerosol_local%alloc_1scl(1, nLay, lw_gas_props%get_band_lims_wavenumber()))
    if (doGP_sgs_cnv) then
       call check_error_msg('rrtmgp_lw_main_cnv_cloud_optics_init',&
            lw_optical_props_cnvcloudsByBand%alloc_1scl(1, nLay, lw_gas_props%get_band_lims_wavenumber()))
    endif
    if (doGP_sgs_pbl) then
       call check_error_msg('rrtmgp_lw_main_pbl_cloud_optics_init',&
            lw_optical_props_pblcloudsByBand%alloc_1scl(1, nLay, lw_gas_props%get_band_lims_wavenumber()))
    endif
    !
    ! ty_fluxes_byband
    !
    flux_allsky%bnd_flux_up => fluxLW_up_allsky
    flux_allsky%bnd_flux_dn => fluxLW_dn_allsky
    flux_clrsky%bnd_flux_up => fluxLW_up_clrsky
    flux_clrsky%bnd_flux_dn => fluxLW_dn_clrsky

    ! Loop over all columns...
    do iCol=1,nCol
       ! Initialize/reset
       lw_optical_props_clrsky%tau       = 0._kind_phys
       lw_optical_props_precipByBand%tau = 0._kind_phys
       lw_optical_props_cloudsByBand%tau = 0._kind_phys
       lw_optical_props_clouds%tau       = 0._kind_phys
       lw_optical_props_clouds%ssa       = 1._kind_phys
       lw_optical_props_clouds%g         = 0._kind_phys
       if (doGP_sgs_cnv) lw_optical_props_cnvcloudsByBand%tau = 0._kind_phys
       if (doGP_sgs_pbl) lw_optical_props_pblcloudsByBand%tau = 0._kind_phys

       ! ###################################################################################
       !
       ! Set gas-concentrations
       !
       ! ###################################################################################
       gas_concentrations%concs(istr_o2)%conc(1,:)   = vmr_o2(iCol,:)
       gas_concentrations%concs(istr_co2)%conc(1,:)  = vmr_co2(iCol,:)
       gas_concentrations%concs(istr_ch4)%conc(1,:)  = vmr_ch4(iCol,:)
       gas_concentrations%concs(istr_n2o)%conc(1,:)  = vmr_n2o(iCol,:)
       gas_concentrations%concs(istr_h2o)%conc(1,:)  = vmr_h2o(iCol,:)
       gas_concentrations%concs(istr_o3)%conc(1,:)   = vmr_o3(iCol,:)

       ! ###################################################################################
       !
       ! Surface emissity in each band
       !
       ! ###################################################################################
       ! Assign same emissivity to all band
       if (semis(iCol) > 1e-6 .and. semis(iCol) <= 1.0) then
          do iBand=1,lw_gas_props%get_nband()
             sfc_emiss_byband(iBand,1) = semis(iCol)
          enddo
       else
          sfc_emiss_byband(1:lw_gas_props%get_nband(),1) = 1.0
       endif

       ! ###################################################################################
       !
       ! Gas-optics
       !
       ! ###################################################################################
       call check_error_msg('rrtmgp_lw_main_gas_optics',lw_gas_props%gas_optics(&
            p_lay(iCol:iCol,:),               & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(iCol:iCol,:),               & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(iCol:iCol,:),               & ! IN  - Temperature @ layer-centers (K)
            tsfg(iCol:iCol),                  & ! IN  - Skin-temperature (K)
            gas_concentrations,               & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
            lw_optical_props_clrsky,          & ! OUT - RRTMGP DDT: longwave optical properties
            sources,                          & ! OUT - RRTMGP DDT: source functions
            tlev=t_lev(iCol:iCol,:)))           ! IN  - Temperature @ layer-interfaces (K) (optional)

       ! ###################################################################################
       !
       ! Cloud-optics
       !
       ! ###################################################################################
       call check_error_msg('rrtmgp_lw_main_cloud_optics',lw_cloud_props%cloud_optics(&
            cld_lwp(iCol:iCol,:),              & ! IN  - Cloud liquid water path (g/m2)
            cld_iwp(iCol:iCol,:),              & ! IN  - Cloud ice water path (g/m2)
            cld_reliq(iCol:iCol,:),            & ! IN  - Cloud liquid effective radius (microns)
            cld_reice(iCol:iCol,:),            & ! IN  - Cloud ice effective radius (microns)
            lw_optical_props_cloudsByBand))      ! OUT - RRTMGP DDT containing cloud radiative properties
                                                 !       in each band
       
       ! Convective cloud-optics?
       if (doGP_sgs_cnv) then
          call check_error_msg('rrtmgp_lw_main_cnv_cloud_optics',lw_cloud_props%cloud_optics(&
               cld_cnv_lwp(iCol:iCol,:),          & ! IN  - Convective cloud liquid water path (g/m2)
               cld_cnv_iwp(iCol:iCol,:),          & ! IN  - Convective cloud ice water path (g/m2)
               cld_cnv_reliq(iCol:iCol,:),        & ! IN  - Convective cloud liquid effective radius (microns)
               cld_cnv_reice(iCol:iCol,:),        & ! IN  - Convective cloud ice effective radius (microns)
               lw_optical_props_cnvcloudsByBand))   ! OUT - RRTMGP DDT containing convective cloud radiative properties
                                                    !       in each band
          !call check_error_msg('rrtmgp_lw_main_increment_cnvclouds_to_clouds',&
          !     lw_optical_props_cnvcloudsByBand%increment(lw_optical_props_cloudsByBand))
       endif

       ! MYNN PBL cloud-optics?
       if (doGP_sgs_pbl) then
          call check_error_msg('rrtmgp_lw_main_pbl_cloud_optics',lw_cloud_props%cloud_optics(&
               cld_pbl_lwp(iCol:iCol,:),          & ! IN  - MYNN-EDMF PBL cloud liquid water path (g/m2)
               cld_pbl_iwp(iCol:iCol,:),          & ! IN  - MYNN-EDMF PBL cloud ice water path (g/m2)
               cld_pbl_reliq(iCol:iCol,:),        & ! IN  - MYNN-EDMF PBL cloud liquid effective radius (microns)
               cld_pbl_reice(iCol:iCol,:),        & ! IN  - MYNN-EDMF PBL cloud ice effective radius (microns)
               lw_optical_props_pblcloudsByBand))   ! OUT - RRTMGP DDT containing MYNN-EDMF PBL cloud radiative properties
                                                    !       in each band
          !call check_error_msg('rrtmgp_lw_main_increment_pblclouds_to_clouds',&
          !     lw_optical_props_pblcloudsByBand%increment(lw_optical_props_cloudsByBand))
       endif

       ! Cloud precipitation optics: rain and snow(+groupel)
       do iLay=1,nLay
          if (cld_frac(iCol,iLay) .gt. 0.) then
             ! Rain optical-depth (No band dependence)
             tau_rain = absrain*cld_rwp(iCol,iLay)

             ! Snow (+groupel) optical-depth (No band dependence)
             if (cld_swp(iCol,iLay) .gt. 0. .and. cld_resnow(iCol,iLay) .gt. 10._kind_phys) then
                tau_snow = abssnow0*1.05756*cld_swp(iCol,iLay)/cld_resnow(iCol,iLay)
             else
                tau_snow = 0.0
             endif
             do iBand=1,lw_gas_props%get_nband()
                lw_optical_props_precipByBand%tau(1,iLay,iBand) = tau_rain + tau_snow
             enddo
          endif
       enddo
       !call check_error_msg('rrtmgp_lw_main_increment_precip_to_clouds',&
       !     lw_optical_props_precipByBand%increment(lw_optical_props_cloudsByBand))

       ! ###################################################################################
       !
       ! Cloud-sampling
       !
       ! ###################################################################################
       ! Change random number seed value for each radiation invocation (isubc_lw =1 or 2).
       if(isubc_lw == 1) then      ! advance prescribed permutation seed
          ipseed_lw = lw_gas_props%get_ngpt() + iCol
       elseif (isubc_lw == 2) then ! use input array of permutaion seeds
          ipseed_lw = icseed_lw(iCol)
       endif
       ! Call RNG
       call random_setseed(ipseed_lw,rng_stat)
       ! Use same rng for each layer
       if (iovr == iovr_max) then
          call random_number(rng1D,rng_stat)
          do iLay=1,nLay
             rng3D(:,iLay,1) = rng1D
          enddo
       else
          do iLay=1,nLay
             call random_number(rng1D,rng_stat)
             rng3D(:,iLay,1) = rng1D
          enddo
       endif
       ! Cloud-overlap.
       ! Maximum-random, random or maximum.
       if (iovr == iovr_maxrand .or. iovr == iovr_rand .or. iovr == iovr_max) then
          call sampled_mask(rng3D, cld_frac(iCol:iCol,:), maskMCICA) 
       endif
       ! Exponential decorrelation length overlap
       if (iovr == iovr_dcorr) then
          ! Generate second RNG
          call random_setseed(ipseed_lw,rng_stat)
          call random_number(rng2D,rng_stat)
          rng3D2(:,:,1) = reshape(source = rng2D,shape=[lw_gas_props%get_ngpt(),nLay])
          !
          call sampled_mask(rng3D, cld_frac(iCol:iCol,:), maskMCICA,                    &
               overlap_param = cloud_overlap_param(iCol:iCol,1:nLay-1), randoms2 = rng3D2)
       endif
       ! Exponential or Exponential-random
       if (iovr == iovr_exp .or. iovr == iovr_exprand) then
          call sampled_mask(rng3D, cld_frac(iCol:iCol,:), maskMCICA,  &
               overlap_param = cloud_overlap_param(iCol:iCol,1:nLay-1))
       endif
       ! Sampling. Map band optical depth to each g-point using McICA
       call check_error_msg('rrtmgp_lw_main_cloud_sampling',&
            draw_samples(maskMCICA, .true., &
            lw_optical_props_cloudsByBand, lw_optical_props_clouds))

       ! ###################################################################################
       !
       ! Compute clear-sky fluxes (gaseous+aerosol) (optional)
       !
       ! ###################################################################################
       ! Add aerosol optics to gas optics
       lw_optical_props_aerosol_local%tau = lw_optical_props_aerosol%tau(iCol:iCol,:,:)
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
          fluxlwUP_clrsky(iCol,:)   = sum(flux_clrsky%bnd_flux_up(1,:,:),dim=2)
          fluxlwDOWN_clrsky(iCol,:) = sum(flux_clrsky%bnd_flux_dn(1,:,:),dim=2)
       else
          fluxlwUP_clrsky(iCol,:)   = 0.0
          fluxlwDOWN_clrsky(iCol,:) = 0.0   
       endif

       ! ###################################################################################
       !
       ! All-sky fluxes (clear-sky + clouds + precipitation)
       !
       ! ###################################################################################

       ! Include convective cloud?
       if (doGP_sgs_cnv) then
          call check_error_msg('rrtmgp_lw_main_increment_cnvclouds_to_clrsky',&
               lw_optical_props_cnvcloudsByBand%increment(lw_optical_props_clouds))
       endif
       
       ! Include MYNN-EDMF PBL clouds?
       if (doGP_sgs_pbl) then
          call check_error_msg('rrtmgp_lw_main_increment_pblclouds_to_clrsky',&
               lw_optical_props_pblcloudsByBand%increment(lw_optical_props_clouds))
       endif
       
       ! Add in precipitation
       call check_error_msg('rrtmgp_lw_main_increment_precip_to_clrsky',&
            lw_optical_props_precipByBand%increment(lw_optical_props_clouds))
       
       ! Include LW cloud-scattering?
       if (doGP_lwscat) then 
          ! Add clear-sky optics to cloud-optics (2-stream)
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
          ! Add cloud optics to clear-sky optics (scalar)
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
       fluxlwUP_allsky(iCol,:)   = sum(flux_allsky%bnd_flux_up(1,:,:),dim=2)
       fluxlwDOWN_allsky(iCol,:) = sum(flux_allsky%bnd_flux_dn(1,:,:),dim=2)
       
       ! Save fluxes for coupling
       fluxlwUP_radtime(iCol,:)   = fluxlwUP_allsky(iCol,:)
       fluxlwDOWN_radtime(iCol,:) = fluxlwDOWN_allsky(iCol,:)

    enddo

  end subroutine rrtmgp_lw_main_run

end module rrtmgp_lw_main
