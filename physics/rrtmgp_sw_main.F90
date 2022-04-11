! ###########################################################################################
! ###########################################################################################
module rrtmgp_sw_main
  use machine,                only: kind_phys
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
  use module_radiation_gases, only: NF_VGAS, getgases, getozn
  use GFS_rrtmgp_pre,         only: iStr_h2o, iStr_co2, iStr_o3, iStr_n2o, iStr_ch4,        &
                                    iStr_o2, iStr_ccl4, iStr_cfc11, iStr_cfc12, iStr_cfc22
  use mersenne_twister,       only: random_setseed, random_number, random_stat
  use rrtmgp_sampling,        only: sampled_mask, draw_samples
  implicit none

  public rrtmgp_sw_main_init, rrtmgp_sw_main_run
contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_main_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_main_init
!! \htmlinclude rrtmgp_sw_main_init.html
!!
  subroutine rrtmgp_sw_main_init(rrtmgp_root_dir, rrtmgp_sw_file_gas, mpicomm, mpirank,     &
       mpiroot, active_gases_array, nrghice, doG_cldoptics, doGP_cldoptics_PADE,            &
       doGP_cldoptics_LUT,rrtmgp_sw_file_clouds, errmsg, errflg)
    ! Inputs
    logical, intent(in) :: &
         doG_cldoptics,         & ! Use legacy RRTMG cloud-optics?
         doGP_cldoptics_PADE,   & ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT       ! Use RRTMGP cloud-optics: LUTs?
    integer, intent(inout) :: &
         nrghice                  ! Number of ice-roughness categories
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,       & ! RTE-RRTMGP root directory
         rrtmgp_sw_file_clouds, & ! RRTMGP file containing coefficients used to compute clouds optical properties
         rrtmgp_sw_file_gas       ! RRTMGP file containing coefficients used to compute gaseous optical properties
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

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! RRTMGP shortwave gas-optics (k-distribution) initialization
    call rrtmgp_sw_gas_optics_init(rrtmgp_root_dir, rrtmgp_sw_file_gas, mpicomm, mpirank,   &
         mpiroot, active_gases_array, errmsg, errflg)

    ! RRTMGP shortwave cloud-optics initialization 
    call rrtmgp_sw_cloud_optics_init(nrghice, mpicomm, mpirank, mpiroot, doG_cldoptics,     &
         doGP_cldoptics_PADE, doGP_cldoptics_LUT, rrtmgp_root_dir, rrtmgp_sw_file_clouds,   &
         errmsg, errflg)

  end subroutine rrtmgp_sw_main_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_main_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_main_run
!! \htmlinclude rrtmgp_sw_main_run.html
!!
  subroutine rrtmgp_sw_main_run(doSWrad, doSWclrsky, top_at_1, doGP_sgs_cnv, doGP_sgs_pbl,  &
       nCol, nDay, nLay, nGases, i_o3, idxday, icseed_sw, iovr, iovr_convcld, iovr_max,     &
       iovr_maxrand, iovr_rand, iovr_dcorr, iovr_exp, iovr_exprand, isubc_sw, iSFC,         &
       sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, coszen,      &
       p_lay, p_lev, t_lay, t_lev, vmr_o2, vmr_h2o, vmr_o3, vmr_ch4, vmr_n2o, vmr_co2,      &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,      &
       cld_rerain, precip_frac, cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp, cld_cnv_reice,     &
       cld_pbl_lwp, cld_pbl_reliq, cld_pbl_iwp, cld_pbl_reice, cloud_overlap_param,         &
       active_gases_array, sw_optical_props_aerosol, scmpsw, fluxswUP_allsky,               &
       fluxswDOWN_allsky, fluxswUP_clrsky, fluxswDOWN_clrsky, cldtausw, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doSWrad,             & ! Flag to calculate SW irradiances
         doSWclrsky,          & ! Flag to compute clear-sky fluxes (diagnostic)
         top_at_1,            & ! Vertical ordering flag
         doGP_sgs_pbl,        & ! Flag for sgs MYNN-EDMF PBL cloud scheme
         doGP_sgs_cnv           ! Flag for sgs convective cloud scheme
    integer,intent(in) :: &
         nCol,                & ! Number of horizontal points
         nDay,                & ! Number of daytime points
         nLay,                & ! Number of vertical grid points.
         nGases,              & ! Number of active gases in RRTMGP
         i_o3,                & !
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
         idxday,              & ! Index array for daytime points
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
         cld_pbl_lwp,         & ! Water path for       SGS PBL liquid cloud-particles
         cld_pbl_reliq,       & ! Effective radius for SGS PBL liquid cloud-particles
         cld_pbl_iwp,         & ! Water path for       SGS PBL ice    cloud-particles
         cld_pbl_reice,       & ! Effective radius for SGS PBL ice    cloud-particles
         cloud_overlap_param    !
    character(len=*), dimension(:), intent(in) :: &
         active_gases_array     ! List of active gases from namelist as array
    type(ty_optical_props_2str),intent(in) :: &
         sw_optical_props_aerosol ! RRTMGP DDT: Shortwave aerosol optical properties (tau,ssa,g)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                     ! CCPP error message
    integer, intent(out) :: &
         errflg                     ! CCPP error flag
    real(kind_phys), dimension(:,:), intent(out) :: &
         cldtausw                   ! Approx 10.mu band layer cloud optical depth  
    real(kind_phys), dimension(:,:), intent(inout) :: &
         fluxswUP_allsky,         & ! RRTMGP upward all-sky flux profiles (W/m2)
         fluxswDOWN_allsky,       & ! RRTMGP downward all-sky flux profiles (W/m2)
         fluxswUP_clrsky,         & ! RRTMGP upward clear-sky flux profiles (W/m2)
         fluxswDOWN_clrsky          ! RRTMGP downward clear-sky flux profiles (W/m2)
    type(cmpfsw_type), dimension(:), intent(inout) :: &
         scmpsw                     ! 2D surface fluxes, components:
                                    ! uvbfc - total sky downward uv-b flux (W/m2)
                                    ! uvbf0 - clear sky downward uv-b flux (W/m2)
                                    ! nirbm - downward nir direct beam flux (W/m2)
                                    ! nirdf - downward nir diffused flux (W/m2)
                                    ! visbm - downward uv+vis direct beam flux (W/m2)
                                    ! visdf - downward uv+vis diffused flux (W/m2)

    ! Local variables
    type(ty_gas_concs) :: &
         gas_concentrations                   ! RRTMGP DDT: trace gas concentrations (vmr)
    type(ty_optical_props_2str) :: &
         sw_optical_props_clrsky,           & ! RRTMGP DDT: Shortwave clear-sky radiative properties
         sw_optical_props_aerosol_local,    & ! RRTMGP DDT: Shortave aerosol radiative properties
         sw_optical_props_cloudsByBand,     & ! RRTMGP DDT: Shortwave optical properties in each band (clouds)
         sw_optical_props_cnvcloudsByBand,  & ! RRTMGP DDT: Shortwave optical properties in each band (convective cloud)
         sw_optical_props_pblcloudsByBand,  & ! RRTMGP DDT: Shortwave optical properties in each band (PBL cloud)
         sw_optical_props_precipByBand,     & ! RRTMGP DDT: Shortwave optical properties in each band (precipitation)
         sw_optical_props_clouds              ! RRTMGP DDT: Shortwave optical properties in each band (sampled clouds)
    type(ty_fluxes_byband) :: &
         flux_allsky,                       & ! RRTMGP DDT: All-sky flux (W/m2)
         flux_clrsky                          ! RRTMGP DDT: Clear-sky flux (W/m2)
    real(kind_phys) :: &
         tau_rain, tau_snow, ssa_rain, ssa_snow, asy_rain, asy_snow, &
         tau_prec, asy_prec, ssa_prec, asyw, ssaw, za1, za2
    real(kind_phys), dimension(sw_gas_props%get_ngpt()) :: rng1D
    real(kind_phys), dimension(sw_gas_props%get_ngpt(),nLay,1) :: rng3D,rng3D2
    real(kind_phys), dimension(sw_gas_props%get_ngpt()*nLay) :: rng2D
    logical, dimension(1,nLay,sw_gas_props%get_ngpt()) :: maskMCICA
    real(kind_phys), dimension(sw_gas_props%get_nband(),1) :: &
         sfc_alb_dir, sfc_alb_dif
    real(kind_phys), dimension(1,nLay+1,sw_gas_props%get_nband()),target :: &
         fluxSW_up_allsky, fluxSW_up_clrsky, fluxSW_dn_allsky, fluxSW_dn_clrsky, fluxSW_dn_dir_allsky
    integer :: iBand, ibd, iCol, iGas, iLay, ipseed_sw
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(2,sw_gas_props%get_nband()) :: bandlimits
    real(kind_phys), dimension(2), parameter :: nIR_uvvis_bnd = (/12850,16000/)
    real(kind_phys), dimension(1,sw_gas_props%get_ngpt()) :: toa_src_sw
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doSWrad) return
    if (nDay .le. 0) then
       fluxswUP_allsky(:,:)   = 0._kind_phys
       fluxswDOWN_allsky(:,:) = 0._kind_phys
       fluxswUP_clrsky(:,:)   = 0._kind_phys
       fluxswDOWN_clrsky(:,:) = 0._kind_phys
       scmpsw                 = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
       return
    endif

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
    call check_error_msg('rrtmgp_sw_main_gas_optics_init',&
         sw_optical_props_clrsky%alloc_2str(1, nLay, sw_gas_props))
    call check_error_msg('rrtmgp_sw_main_cloud_optics_init',&
         sw_optical_props_cloudsByBand%alloc_2str(1, nLay, sw_gas_props%get_band_lims_wavenumber()))
    call check_error_msg('rrtmgp_sw_main_precip_optics_init',&
         sw_optical_props_precipByBand%alloc_2str(1, nLay, sw_gas_props%get_band_lims_wavenumber()))
    call check_error_msg('rrtmgp_sw_mian_cloud_sampling_init', & 
         sw_optical_props_clouds%alloc_2str(1, nLay, sw_gas_props))
    call check_error_msg('rrtmgp_sw_main_aerosol_optics_init',&
         sw_optical_props_aerosol_local%alloc_2str(1, nLay, sw_gas_props%get_band_lims_wavenumber()))
    if (doGP_sgs_cnv) then
       call check_error_msg('rrtmgp_sw_main_cnv_cloud_optics_init',&
            sw_optical_props_cnvcloudsByBand%alloc_2str(1, nLay, sw_gas_props%get_band_lims_wavenumber()))
    endif
    if (doGP_sgs_pbl) then
       call check_error_msg('rrtmgp_sw_main_pbl_cloud_optics_init',&
            sw_optical_props_pblcloudsByBand%alloc_2str(1, nLay, sw_gas_props%get_band_lims_wavenumber()))
    endif
    !
    ! ty_fluxes_byband
    !
    flux_allsky%bnd_flux_up     => fluxSW_up_allsky
    flux_allsky%bnd_flux_dn     => fluxSW_dn_allsky
    flux_allsky%bnd_flux_dn_dir => fluxSW_dn_dir_allsky
    flux_clrsky%bnd_flux_up     => fluxSW_up_clrsky
    flux_clrsky%bnd_flux_dn     => fluxSW_dn_clrsky

    ! Loop over all (daylit)columns...
    do iCol=1,nDay
       ! Initialize/reset
       sw_optical_props_clouds%tau        = 0._kind_phys
       sw_optical_props_clouds%ssa        = 1._kind_phys
       sw_optical_props_clouds%g          = 0._kind_phys
       sw_optical_props_clrsky%tau        = 0._kind_phys
       sw_optical_props_clrsky%ssa        = 1._kind_phys
       sw_optical_props_clrsky%g          = 0._kind_phys
       sw_optical_props_cloudsByBand%tau  = 0._kind_phys
       sw_optical_props_cloudsByBand%ssa  = 1._kind_phys
       sw_optical_props_cloudsByBand%g    = 0._kind_phys
       sw_optical_props_precipByBand%tau  = 0._kind_phys
       sw_optical_props_precipByBand%ssa  = 1._kind_phys
       sw_optical_props_precipByBand%g    = 0._kind_phys
       sw_optical_props_aerosol_local%tau = 0._kind_phys
       sw_optical_props_aerosol_local%ssa = 1._kind_phys
       sw_optical_props_aerosol_local%g   = 0._kind_phys
       if (doGP_sgs_cnv) then
          sw_optical_props_cnvcloudsByBand%tau = 0._kind_phys
          sw_optical_props_cnvcloudsByBand%ssa = 1._kind_phys
          sw_optical_props_cnvcloudsByBand%g   = 0._kind_phys
       endif
       if (doGP_sgs_pbl) then
          sw_optical_props_pblcloudsByBand%tau = 0._kind_phys
          sw_optical_props_pblcloudsByBand%ssa = 1._kind_phys
          sw_optical_props_pblcloudsByBand%g   = 0._kind_phys
       endif

       ! ###################################################################################
       !
       ! Set gas-concentrations
       !
       ! ###################################################################################
       gas_concentrations%concs(istr_o2)%conc(1,:)   = vmr_o2(idxday(iCol),:)
       gas_concentrations%concs(istr_co2)%conc(1,:)  = vmr_co2(idxday(iCol),:)
       gas_concentrations%concs(istr_ch4)%conc(1,:)  = vmr_ch4(idxday(iCol),:)
       gas_concentrations%concs(istr_n2o)%conc(1,:)  = vmr_n2o(idxday(iCol),:)
       gas_concentrations%concs(istr_h2o)%conc(1,:)  = vmr_h2o(idxday(iCol),:)
       gas_concentrations%concs(istr_o3)%conc(1,:)   = vmr_o3(idxday(iCol),:)

       ! ###################################################################################
       !
       ! Set surface albedo
       !
       ! Use near-IR albedo for bands with wavenumbers extending to 12850cm-1
       ! Use uv-vis albedo for bands with wavenumbers greater than 16000cm-1
       ! For overlapping band, average near-IR and us-vis albedos.
       !
       ! ###################################################################################
       bandlimits = sw_gas_props%get_band_lims_wavenumber()
       do iBand=1,sw_gas_props%get_nband()
          if (bandlimits(1,iBand) .lt. nIR_uvvis_bnd(1)) then
             sfc_alb_dir(iBand,1) = sfc_alb_nir_dir(idxday(iCol))
             sfc_alb_dif(iBand,1) = sfc_alb_nir_dif(idxday(iCol))
          endif
          if (bandlimits(1,iBand) .eq. nIR_uvvis_bnd(1)) then
             sfc_alb_dir(iBand,1) = 0.5_kind_phys*(sfc_alb_nir_dir(idxday(iCol)) + sfc_alb_uvvis_dir(idxday(iCol)))
             sfc_alb_dif(iBand,1) = 0.5_kind_phys*(sfc_alb_nir_dif(idxday(iCol)) + sfc_alb_uvvis_dif(idxday(iCol)))
             ibd = iBand
          endif
          if (bandlimits(1,iBand) .ge. nIR_uvvis_bnd(2)) then
             sfc_alb_dir(iBand,1) = sfc_alb_uvvis_dir(idxday(iCol))
             sfc_alb_dif(iBand,1) = sfc_alb_uvvis_dif(idxday(iCol))
          endif
       enddo

       ! ###################################################################################
       !
       ! Gas-optics
       !
       ! ###################################################################################
       call check_error_msg('rrtmgp_sw_main_gas_optics',sw_gas_props%gas_optics(&
            p_lay(idxday(iCol:iCol),:),     & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(idxday(iCol:iCol),:),     & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(idxday(iCol:iCol),:),     & ! IN  - Temperature @ layer-centers (K)
            gas_concentrations,             & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
            sw_optical_props_clrsky,        & ! OUT - RRTMGP DDT: Shortwave optical properties, by
                                              !                   spectral point (tau,ssa,g)
            toa_src_sw))                      ! OUT - TOA incident shortwave radiation (spectral)

       ! ###################################################################################
       !
       ! Cloud-optics
       !
       ! ###################################################################################
       call check_error_msg('rrtmgp_sw_main_cloud_optics',sw_cloud_props%cloud_optics(&
            cld_lwp(idxday(iCol:iCol),:),   & ! IN  - Cloud liquid water path
            cld_iwp(idxday(iCol:iCol),:),   & ! IN  - Cloud ice water path
            cld_reliq(idxday(iCol:iCol),:), & ! IN  - Cloud liquid effective radius
            cld_reice(idxday(iCol:iCol),:), & ! IN  - Cloud ice effective radius
            sw_optical_props_cloudsByBand))   ! OUT - RRTMGP DDT: Shortwave optical properties, 
                                              !       in each band (tau,ssa,g)
       cldtausw(idxday(iCol),:) = sw_optical_props_cloudsByBand%tau(1,:,11)

       ! Convective cloud-optics?
       if (doGP_sgs_cnv) then
          call check_error_msg('rrtmgp_sw_main_cnv_cloud_optics',sw_cloud_props%cloud_optics(&
               cld_cnv_lwp(idxday(iCol:iCol),:),          & ! IN  - Convective cloud liquid water path (g/m2)
               cld_cnv_iwp(idxday(iCol:iCol),:),          & ! IN  - Convective cloud ice water path (g/m2)
               cld_cnv_reliq(idxday(iCol:iCol),:),        & ! IN  - Convective cloud liquid effective radius (microns)
               cld_cnv_reice(idxday(iCol:iCol),:),        & ! IN  - Convective cloud ice effective radius (microns)
               sw_optical_props_cnvcloudsByBand))           ! OUT - RRTMGP DDT containing convective cloud radiative properties
                                                            !       in each band
          !call check_error_msg('rrtmgp_sw_main_increment_cnvclouds_to_clouds',&
          !     sw_optical_props_cnvcloudsByBand%increment(sw_optical_props_cloudsByBand))
       endif

       ! MYNN PBL cloud-optics?
       if (doGP_sgs_pbl) then
          call check_error_msg('rrtmgp_sw_main_pbl_cloud_optics',sw_cloud_props%cloud_optics(&
               cld_pbl_lwp(idxday(iCol:iCol),:),          & ! IN  - MYNN-EDMF PBL cloud liquid water path (g/m2)
               cld_pbl_iwp(idxday(iCol:iCol),:),          & ! IN  - MYNN-EDMF PBL cloud ice water path (g/m2)
               cld_pbl_reliq(idxday(iCol:iCol),:),        & ! IN  - MYNN-EDMF PBL cloud liquid effective radius (microns)
               cld_pbl_reice(idxday(iCol:iCol),:),        & ! IN  - MYNN-EDMF PBL cloud ice effective radius (microns)
               sw_optical_props_pblcloudsByBand))           ! OUT - RRTMGP DDT containing MYNN-EDMF PBL cloud radiative properties
                                                            !       in each band
          !call check_error_msg('rrtmgp_sw_main_increment_pblclouds_to_clouds',&
          !     sw_optical_props_pblcloudsByBand%increment(sw_optical_props_cloudsByBand))
       endif

       ! Cloud precipitation optics: rain and snow(+groupel)
       do iLay=1,nLay
          if (cld_frac(idxday(iCol),iLay) .gt. 1.e-12_kind_phys) then
             ! Rain/Snow optical depth (No band dependence)
             tau_rain = cld_rwp(idxday(iCol),iLay)*a0r                   
             if (cld_swp(idxday(iCol),iLay) .gt. 0. .and. cld_resnow(idxday(iCol),iLay) .gt. 10._kind_phys) then
                tau_snow = cld_swp(idxday(iCol),iLay)*1.09087*(a0s + a1s/(1.0315*cld_resnow(idxday(iCol),iLay)))     ! fu's formula 
             else
                tau_snow = 0._kind_phys
             endif
             
             ! Rain/Snow single-scattering albedo and asymmetry (Band dependent)
             do iBand=1,sw_gas_props%get_nband()
                ! By species
                ssa_rain = tau_rain*(1.-b0r(iBand))
                asy_rain = ssa_rain*c0r(iBand)
                ssa_snow = tau_snow*(1.-(b0s(iBand)+b1s(iBand)*1.0315*cld_resnow(idxday(iCol),iLay)))
                asy_snow = ssa_snow*c0s(iBand)
                ! Combine
                tau_prec = max(1.e-12_kind_phys, tau_rain + tau_snow)
                ssa_prec = max(1.e-12_kind_phys, ssa_rain + ssa_snow)
                asy_prec = max(1.e-12_kind_phys, asy_rain + asy_snow)
                asyw     = asy_prec/max(1.e-12_kind_phys, ssa_prec)
                ssaw     = min(1._kind_phys-0.000001, ssa_prec/tau_prec)
                za1      = asyw * asyw
                za2      = ssaw * za1                      
                sw_optical_props_precipByBand%tau(1,iLay,iBand) = (1._kind_phys - za2) * tau_prec
                sw_optical_props_precipByBand%ssa(1,iLay,iBand) = (ssaw - za2) / (1._kind_phys - za2)
                sw_optical_props_precipByBand%g(1,iLay,iBand)   = asyw/(1+asyw)
             enddo
          endif
       enddo

       ! ###################################################################################
       !
       ! Cloud-sampling
       !
       ! ###################################################################################
       ! Change random number seed value for each radiation invocation (isubc_sw =1 or 2).
       if(isubc_sw == 1) then      ! advance prescribed permutation seed
          ipseed_sw = sw_gas_props%get_ngpt() + iCol
       elseif (isubc_sw == 2) then ! use input array of permutaion seeds
          ipseed_sw = icseed_sw(idxday(iCol))
       endif
       ! Call RNG
       call random_setseed(ipseed_sw,rng_stat)
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
          call sampled_mask(rng3D, cld_frac(idxday(iCol:iCol),:), maskMCICA) 
       endif
       ! Exponential decorrelation length overlap
       if (iovr == iovr_dcorr) then
          ! Generate second RNG
          call random_setseed(ipseed_sw,rng_stat)
          call random_number(rng2D,rng_stat)
          rng3D2(:,:,1) = reshape(source = rng2D,shape=[sw_gas_props%get_ngpt(),nLay])
          !
          call sampled_mask(rng3D, cld_frac(idxday(iCol:iCol),:), maskMCICA,                    &
               overlap_param = cloud_overlap_param(idxday(iCol:iCol),1:nLay-1), randoms2 = rng3D2)
       endif
       ! Exponential or Exponential-random
       if (iovr == iovr_exp .or. iovr == iovr_exprand) then
          call sampled_mask(rng3D, cld_frac(idxday(iCol:iCol),:), maskMCICA,  &
               overlap_param = cloud_overlap_param(idxday(iCol:iCol),1:nLay-1))
       endif
       ! Sampling. Map band optical depth to each g-point using McICA
       call check_error_msg('rrtmgp_sw_main_cloud_sampling',&
            draw_samples(maskMCICA, .true., &
            sw_optical_props_cloudsByBand, sw_optical_props_clouds))

       ! ###################################################################################
       !
       ! Compute clear-sky fluxes (gaseous+aerosol) (optional)
       !
       ! ###################################################################################
       ! Add aerosol optics to gas optics
       sw_optical_props_aerosol_local%tau = sw_optical_props_aerosol%tau(iCol:iCol,:,:)
       sw_optical_props_aerosol_local%ssa = sw_optical_props_aerosol%ssa(iCol:iCol,:,:)
       sw_optical_props_aerosol_local%g   = sw_optical_props_aerosol%g(iCol:iCol,:,:)
       call check_error_msg('rrtmgp_sw_main_increment_aerosol_to_clrsky',&
            sw_optical_props_aerosol_local%increment(sw_optical_props_clrsky))

       ! Delta-scale optical properties
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clrsky%delta_scale())
       if (doSWclrsky) then
          call check_error_msg('rrtmgp_sw_main_rte_sw_clrsky',rte_sw(     &
               sw_optical_props_clrsky,   & ! IN  - optical-properties
               top_at_1,                  & ! IN  - veritcal ordering flag
               coszen(idxday(iCol:iCol)), & ! IN  - Cosine of solar zenith angle
               toa_src_sw,                & ! IN  - incident solar flux at TOA
               sfc_alb_dir,               & ! IN  - Shortwave surface albedo (direct)
               sfc_alb_dif,               & ! IN  - Shortwave surface albedo (diffuse)
               flux_clrsky))                ! OUT - Fluxes, clear-sky, 3D (1,nLay,nBand) 
          ! Store fluxes
          fluxswUP_clrsky(idxday(iCol),:)   = sum(flux_clrsky%bnd_flux_up(1,:,:),dim=2)
          fluxswDOWN_clrsky(idxday(iCol),:) = sum(flux_clrsky%bnd_flux_dn(1,:,:),dim=2)
       else
          fluxswUP_clrsky(idxday(iCol),:)   = 0.0
          fluxswDOWN_clrsky(idxday(iCol),:) = 0.0
       endif

       ! ###################################################################################
       !
       ! All-sky fluxes (clear-sky + clouds + precipitation)
       !
       ! ###################################################################################

       ! Include convective cloud?
       if (doGP_sgs_cnv) then
          call check_error_msg('rrtmgp_sw_main_increment_cnvclouds_to_clrsky',&
               sw_optical_props_cnvcloudsByBand%increment(sw_optical_props_clouds))
       endif
       
       ! Include MYNN-EDMF PBL clouds?
       if (doGP_sgs_pbl) then
          call check_error_msg('rrtmgp_sw_main_increment_pblclouds_to_clrsky',&
               sw_optical_props_pblcloudsByBand%increment(sw_optical_props_clouds))
       endif
       
       ! Add in precipitation
       call check_error_msg('rrtmgp_sw_main_increment_precip_to_clrsky',&
            sw_optical_props_precipByBand%increment(sw_optical_props_clouds))

       ! Delta-scale optical properties
       call check_error_msg('rrtmgp_sw_main_delta_scale',sw_optical_props_clrsky%delta_scale())
       call check_error_msg('rrtmgp_sw_main_rte_sw_allsky',rte_sw(     &
            sw_optical_props_clouds,   & ! IN  - optical-properties
            top_at_1,                  & ! IN  - veritcal ordering flag
            coszen(idxday(iCol:iCol)), & ! IN  - Cosine of solar zenith angle
            toa_src_sw,                & ! IN  - incident solar flux at TOA
            sfc_alb_dir,               & ! IN  - Shortwave surface albedo (direct)
            sfc_alb_dif,               & ! IN  - Shortwave surface albedo (diffuse)
            flux_allsky))                ! OUT - Fluxes, clear-sky, 3D (1,nLay,nBand)
       
       ! Store fluxes
       fluxswUP_allsky(idxday(iCol),:)   = sum(flux_allsky%bnd_flux_up(1,:,:),dim=2)
       fluxswDOWN_allsky(idxday(iCol),:) = sum(flux_allsky%bnd_flux_dn(1,:,:),dim=2)
       ! Near IR
       scmpsw(idxday(iCol))%nirbm = sum(flux_allsky%bnd_flux_dn_dir(1,iSFC,1:ibd-1))  + &
                                        flux_allsky%bnd_flux_dn_dir(1,iSFC,ibd)/2.
       scmpsw(idxday(iCol))%nirdf = (sum(flux_allsky%bnd_flux_dn(1,iSFC,1:ibd-1))     + &
                                         flux_allsky%bnd_flux_dn(1,iSFC,ibd)/2.)      - &
                                    (sum(flux_allsky%bnd_flux_dn_dir(1,iSFC,1:ibd-1)) + &
                                         flux_allsky%bnd_flux_dn_dir(1,iSFC,ibd)/2.)
       ! UV-VIS
       scmpsw(idxday(iCol))%visbm = sum(flux_allsky%bnd_flux_dn_dir(1,iSFC,ibd+1:sw_gas_props%get_nband()))  + &
                                        flux_allsky%bnd_flux_dn_dir(1,iSFC,ibd)/2.
       scmpsw(idxday(iCol))%visdf = (sum(flux_allsky%bnd_flux_dn(1,iSFC,ibd+1:sw_gas_props%get_nband()))     + &
                                         flux_allsky%bnd_flux_dn(1,iSFC,ibd)/2. )                            - &
                                    (sum(flux_allsky%bnd_flux_dn_dir(1,iSFC,ibd+1:sw_gas_props%get_nband())) + &
                                         flux_allsky%bnd_flux_dn_dir(1,iSFC,ibd)/2.)
    enddo
  end subroutine rrtmgp_sw_main_run
end module rrtmgp_sw_main
