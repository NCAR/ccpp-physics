module rrtmgp_sw_gas_optics
  use machine,                only: kind_phys
  use module_radiation_gases, only: NF_VGAS
  use mo_rte_kind,            only: wl
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,  only: ty_gas_concs
  use rrtmgp_aux,             only: check_error_msg
  use mo_optical_props,       only: ty_optical_props_2str
  use mo_compute_bc,          only: compute_bc
  use GFS_rrtmgp_pre,         only: active_gases_array
  use netcdf

  implicit none

  ! RRTMGP k-distribution LUTs.
  type(ty_gas_optics_rrtmgp) :: sw_gas_props
  integer, dimension(:), allocatable :: &
       kminor_start_lowerSW,              & ! Starting index in the [1, nContributors] vector for a contributor
                                            ! given by \"minor_gases_lower\" (lower atmosphere)
       kminor_start_upperSW                 ! Starting index in the [1, nContributors] vector for a contributor
                                            ! given by \"minor_gases_upper\" (upper atmosphere)
  integer, dimension(:,:), allocatable :: &
       band2gptSW,                        & ! Beginning and ending gpoint for each band
       minor_limits_gpt_lowerSW,          & ! Beginning and ending gpoint for each minor interval in lower atmosphere
       minor_limits_gpt_upperSW             ! Beginning and ending gpoint for each minor interval in upper atmosphere
  integer, dimension(:,:,:), allocatable :: &
       key_speciesSW                        ! Key species pair for each band
  real(kind_phys) :: &
       press_ref_tropSW,                  & ! Reference pressure separating the lower and upper atmosphere [Pa]
       temp_ref_pSW,                      & ! Standard spectroscopic reference pressure [Pa]
       temp_ref_tSW,                      & ! Standard spectroscopic reference temperature [K]
       tsi_defaultSW,                     & ! 
       mg_defaultSW,                      & ! Mean value of Mg2 index over the average solar cycle from the NRLSSI2 model of solar variability
       sb_defaultSW                         ! Mean value of sunspot index over the average solar cycle from the NRLSSI2 model of solar variability
  real(kind_phys), dimension(:), allocatable :: &
       press_refSW,                       & ! Pressures for reference atmosphere; press_ref(# reference layers) [Pa]
       temp_refSW,                        & ! Temperatures for reference atmosphere; temp_ref(# reference layers) [K]
       solar_quietSW,                     & ! Spectrally-dependent quiet sun irradiance from the NRLSSI2 model of solar variability
       solar_facularSW,                   & ! Spectrally-dependent facular term from the NRLSSI2 model of solar variability
       solar_sunspotSW                      ! Spectrally-dependent sunspot term from the NRLSSI2 model of solar variability
  real(kind_phys), dimension(:,:), allocatable :: &
       band_limsSW                          ! Beginning and ending wavenumber [cm -1] for each band
  real(kind_phys), dimension(:,:,:), allocatable :: &
       vmr_refSW,                         & ! Volume mixing ratios for reference atmosphere
       kminor_lowerSW,                    & ! (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to
                                            ! [nTemp x nEta x nContributors] array)
       kminor_upperSW,                    & ! (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to
                                            ! [nTemp x nEta x nContributors] array)
       rayl_lowerSW,                      & ! Stored coefficients due to rayleigh scattering contribution
       rayl_upperSW                         ! Stored coefficients due to rayleigh scattering contribution
  real(kind_phys), dimension(:,:,:,:), allocatable :: &
       kmajorSW                             ! Stored absorption coefficients due to major absorbing gases
  character(len=32),  dimension(:), allocatable :: &
       gas_namesSW,                       & ! Names of absorbing gases
       gas_minorSW,                       & ! Name of absorbing minor gas
       identifier_minorSW,                & ! Unique string identifying minor gas
       minor_gases_lowerSW,               & ! Names of minor absorbing gases in lower atmosphere
       minor_gases_upperSW,               & ! Names of minor absorbing gases in upper atmosphere
       scaling_gas_lowerSW,               & ! Absorption also depends on the concentration of this gas
       scaling_gas_upperSW                  ! Absorption also depends on the concentration of this gas
  logical(wl), dimension(:), allocatable :: &
       minor_scales_with_density_lowerSW, & ! Density scaling is applied to minor absorption coefficients
       minor_scales_with_density_upperSW, & ! Density scaling is applied to minor absorption coefficients
       scale_by_complement_lowerSW,       & ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T)
       scale_by_complement_upperSW          ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T)
contains

  ! #########################################################################################
  ! SUBROUTINE sw_gas_optics_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_gas_optics_init
!! \htmlinclude rrtmgp_sw_gas_optics.html
!!
  subroutine rrtmgp_sw_gas_optics_init(nCol, nLev, nThreads, rrtmgp_root_dir,               &
       rrtmgp_sw_file_gas, gas_concentrations, mpicomm, mpirank, mpiroot, errmsg, errflg)

    ! Inputs
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,  & ! RTE-RRTMGP root directory
         rrtmgp_sw_file_gas  ! RRTMGP file containing coefficients used to compute gaseous optical properties
    integer,intent(in) :: &
         nCol,             & ! Number of horizontal gridpoints.
         nLev,             & ! Number of vertical levels.
         nThreads,         & ! Number of openMP threads
         mpicomm,          & ! MPI communicator
         mpirank,          & ! Current MPI rank
         mpiroot             ! Master MPI rank
    type(ty_gas_concs),intent(inout)  :: &
         gas_concentrations  ! RRTMGP DDT containing active trace gases.
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg              ! CCPP error message
    integer,          intent(out) :: &
         errflg              ! CCPP error code

    ! Local variables
    integer :: status, ncid, dimid, varID, iGas, ntemps, npress, ngptsSW, nabsorbers,    &
         nextrabsorbers, nminorabsorbers, nmixingfracs, nlayers, nbnds, npairs,          &
         nminor_absorber_intervals_lower, nminor_absorber_intervals_upper,               &
         ncontributors_lower, ncontributors_upper
    integer,dimension(:),allocatable :: temp1, temp2, temp3, temp4
    character(len=264) :: sw_gas_props_file

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the gphysics_nml
    sw_gas_props_file   = trim(rrtmgp_root_dir)//trim(rrtmgp_sw_file_gas)

    ! Read dimensions for k-distribution fields (only on master processor(0))
!    if (mpirank .eq. mpiroot) then
       ! Open file
       status = nf90_open(trim(sw_gas_props_file), NF90_NOWRITE, ncid)

       ! Read dimensions for k-distribution fields
       status = nf90_inq_dimid(ncid, 'temperature', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ntemps)
       status = nf90_inq_dimid(ncid, 'pressure', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=npress)
       status = nf90_inq_dimid(ncid, 'absorber', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nabsorbers)
       status = nf90_inq_dimid(ncid, 'minor_absorber',dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nminorabsorbers)
       status = nf90_inq_dimid(ncid, 'absorber_ext', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nextrabsorbers)
       status = nf90_inq_dimid(ncid, 'mixing_fraction', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nmixingfracs)
       status = nf90_inq_dimid(ncid, 'atmos_layer', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nlayers)
       status = nf90_inq_dimid(ncid, 'bnd', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nbnds)
       status = nf90_inq_dimid(ncid, 'gpt', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ngptsSW)
       status = nf90_inq_dimid(ncid, 'pair', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=npairs)
       status = nf90_inq_dimid(ncid, 'contributors_lower',dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ncontributors_lower)
       status = nf90_inq_dimid(ncid, 'contributors_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ncontributors_upper)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_lower', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nminor_absorber_intervals_lower)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nminor_absorber_intervals_upper)
 
       ! Allocate space for arrays
       if (.not. allocated(gas_namesSW))         &
            allocate(gas_namesSW(nabsorbers))
       if (.not. allocated(scaling_gas_lowerSW)) &
            allocate(scaling_gas_lowerSW(nminor_absorber_intervals_lower))
       if (.not. allocated(scaling_gas_upperSW)) &
            allocate(scaling_gas_upperSW(nminor_absorber_intervals_upper))
       if (.not. allocated(gas_minorSW))         &
            allocate(gas_minorSW(nminorabsorbers))
       if (.not. allocated(identifier_minorSW))  &
            allocate(identifier_minorSW(nminorabsorbers))
       if (.not. allocated(minor_gases_lowerSW)) &
            allocate(minor_gases_lowerSW(nminor_absorber_intervals_lower))
       if (.not. allocated(minor_gases_upperSW)) &
            allocate(minor_gases_upperSW(nminor_absorber_intervals_upper))
       if (.not. allocated(minor_limits_gpt_lowerSW)) &
            allocate(minor_limits_gpt_lowerSW(npairs,nminor_absorber_intervals_lower))
       if (.not. allocated(minor_limits_gpt_upperSW)) &
            allocate(minor_limits_gpt_upperSW(npairs,nminor_absorber_intervals_upper))
       if (.not. allocated(band2gptSW)) &
            allocate(band2gptSW(2,nbnds))
       if (.not. allocated(key_speciesSW)) &
            allocate(key_speciesSW(2,nlayers,nbnds))
       if (.not. allocated(band_limsSW)) &
            allocate(band_limsSW(2,nbnds))
       if (.not. allocated(press_refSW)) &
            allocate(press_refSW(npress))
       if (.not. allocated(temp_refSW)) &
            allocate(temp_refSW(ntemps))
       if (.not. allocated(vmr_refSW)) &
            allocate(vmr_refSW(nlayers, nextrabsorbers, ntemps))
       if (.not. allocated(kminor_lowerSW)) &
            allocate(kminor_lowerSW(ncontributors_lower, nmixingfracs, ntemps))
       if (.not. allocated(kmajorSW)) &
            allocate(kmajorSW(ngptsSW, nmixingfracs,  npress+1, ntemps))
       if (.not. allocated(kminor_start_lowerSW)) &
            allocate(kminor_start_lowerSW(nminor_absorber_intervals_lower))
       if (.not. allocated(kminor_upperSW)) &
            allocate(kminor_upperSW(ncontributors_upper, nmixingfracs, ntemps))
       if (.not. allocated(kminor_start_upperSW)) &
            allocate(kminor_start_upperSW(nminor_absorber_intervals_upper))
       if (.not. allocated(minor_scales_with_density_lowerSW)) &
            allocate(minor_scales_with_density_lowerSW(nminor_absorber_intervals_lower))
       if (.not. allocated(minor_scales_with_density_upperSW)) &
            allocate(minor_scales_with_density_upperSW(nminor_absorber_intervals_upper))
       if (.not. allocated(scale_by_complement_lowerSW)) &
            allocate(scale_by_complement_lowerSW(nminor_absorber_intervals_lower))
       if (.not. allocated(scale_by_complement_upperSW)) &
            allocate(scale_by_complement_upperSW(nminor_absorber_intervals_upper))
       if (.not. allocated(rayl_upperSW)) &
            allocate(rayl_upperSW(ngptsSW, nmixingfracs, ntemps))
       if (.not. allocated(rayl_lowerSW)) &
            allocate(rayl_lowerSW(ngptsSW, nmixingfracs, ntemps))
       if (.not. allocated(solar_quietSW)) &
            allocate(solar_quietSW(ngptsSW))
       if (.not. allocated(solar_facularSW)) &
            allocate(solar_facularSW(ngptsSW))	
       if (.not. allocated(solar_sunspotSW)) &
            allocate(solar_sunspotSW(ngptsSW))
       if (.not. allocated(temp1)) &
            allocate(temp1(nminor_absorber_intervals_lower))
       if (.not. allocated(temp2)) &
            allocate(temp2(nminor_absorber_intervals_upper))
       if (.not. allocated(temp3)) &
            allocate(temp3(nminor_absorber_intervals_lower))
       if (.not. allocated(temp4)) &
            allocate(temp4(nminor_absorber_intervals_upper))

       ! Read in fields from file
       if (mpirank==mpiroot) write (*,*) 'Reading RRTMGP shortwave k-distribution data ... '
       status = nf90_inq_varid(ncid, 'gas_names', varID)
       status = nf90_get_var(  ncid, varID, gas_namesSW)       
       status = nf90_inq_varid(ncid, 'scaling_gas_lower', varID)
       status = nf90_get_var(  ncid, varID, scaling_gas_lowerSW)       
       status = nf90_inq_varid(ncid, 'scaling_gas_upper', varID)
       status = nf90_get_var(  ncid, varID, scaling_gas_upperSW)       
       status = nf90_inq_varid(ncid, 'gas_minor', varID)
       status = nf90_get_var(  ncid, varID, gas_minorSW)       
       status = nf90_inq_varid(ncid, 'identifier_minor', varID)
       status = nf90_get_var(  ncid, varID, identifier_minorSW)       
       status = nf90_inq_varid(ncid, 'minor_gases_lower', varID)
       status = nf90_get_var(  ncid, varID, minor_gases_lowerSW)       
       status = nf90_inq_varid(ncid, 'minor_gases_upper', varID)
       status = nf90_get_var(  ncid, varID, minor_gases_upperSW)       
       status = nf90_inq_varid(ncid, 'minor_limits_gpt_lower', varID)
       status = nf90_get_var(  ncid, varID, minor_limits_gpt_lowerSW)       
       status = nf90_inq_varid(ncid, 'minor_limits_gpt_upper', varID)
       status = nf90_get_var(  ncid, varID, minor_limits_gpt_upperSW)       
       status = nf90_inq_varid(ncid, 'bnd_limits_gpt', varID)
       status = nf90_get_var(  ncid, varID, band2gptSW)       
       status = nf90_inq_varid(ncid, 'key_species', varID)
       status = nf90_get_var(  ncid, varID, key_speciesSW)       
       status = nf90_inq_varid(ncid,'bnd_limits_wavenumber', varID)
       status = nf90_get_var(  ncid, varID, band_limsSW)       
       status = nf90_inq_varid(ncid, 'press_ref', varID)
       status = nf90_get_var(  ncid, varID, press_refSW)       
       status = nf90_inq_varid(ncid, 'temp_ref', varID)
       status = nf90_get_var(  ncid, varID, temp_refSW)       
       status = nf90_inq_varid(ncid, 'absorption_coefficient_ref_P', varID)
       status = nf90_get_var(  ncid, varID, temp_ref_pSW)
       status = nf90_inq_varid(ncid, 'absorption_coefficient_ref_T', varID)
       status = nf90_get_var(  ncid, varID, temp_ref_tSW) 
       status = nf90_inq_varid(ncid, 'tsi_default', varID)
       status = nf90_get_var(  ncid, varID, tsi_defaultSW)
       status = nf90_inq_varid(ncid, 'mg_default', varID)
       status =	nf90_get_var(  ncid, varID, mg_defaultSW)
       status = nf90_inq_varid(ncid, 'sb_default', varID)
       status =	nf90_get_var(  ncid, varID, sb_defaultSW)
       status = nf90_inq_varid(ncid, 'press_ref_trop', varID)
       status = nf90_get_var(  ncid, varID, press_ref_tropSW)       
       status = nf90_inq_varid(ncid, 'kminor_lower', varID)
       status = nf90_get_var(  ncid, varID, kminor_lowerSW)       
       status = nf90_inq_varid(ncid, 'kminor_upper', varID)
       status = nf90_get_var(  ncid, varID, kminor_upperSW)       
       status = nf90_inq_varid(ncid, 'vmr_ref', varID)
       status = nf90_get_var(  ncid, varID, vmr_refSW)       
       status = nf90_inq_varid(ncid, 'kmajor', varID)
       status = nf90_get_var(  ncid, varID, kmajorSW)      
       status = nf90_inq_varid(ncid, 'kminor_start_lower', varID)
       status = nf90_get_var(  ncid, varID, kminor_start_lowerSW)       
       status = nf90_inq_varid(ncid, 'kminor_start_upper', varID)
       status = nf90_get_var(  ncid, varID, kminor_start_upperSW)       
       status = nf90_inq_varid(ncid, 'solar_source_quiet', varID)
       status = nf90_get_var(  ncid, varID, solar_quietSW)
       status = nf90_inq_varid(ncid, 'solar_source_facular', varID)
       status = nf90_get_var(  ncid, varID, solar_facularSW)
       status = nf90_inq_varid(ncid, 'solar_source_sunspot', varID)
       status = nf90_get_var(  ncid, varID, solar_sunspotSW)       
       status = nf90_inq_varid(ncid, 'rayl_lower', varID)
       status = nf90_get_var(  ncid, varID, rayl_lowerSW)
       status = nf90_inq_varid(ncid, 'rayl_upper', varID)
       status = nf90_get_var(  ncid, varID, rayl_upperSW)

       ! Logical fields are read in as integers and then converted to logicals.
       status = nf90_inq_varid(ncid,'minor_scales_with_density_lower', varID)
       status = nf90_get_var(  ncid, varID,temp1)
       minor_scales_with_density_lowerSW(:) = .false.
       where(temp1 .eq. 1) minor_scales_with_density_lowerSW(:) = .true.
       status = nf90_inq_varid(ncid,'minor_scales_with_density_upper', varID)
       status = nf90_get_var(  ncid, varID,temp2)
       minor_scales_with_density_upperSW(:) = .false.
       where(temp2 .eq. 1) minor_scales_with_density_upperSW(:) = .true.
       status = nf90_inq_varid(ncid,'scale_by_complement_lower', varID)
       status = nf90_get_var(  ncid, varID,temp3)
       scale_by_complement_lowerSW(:) = .false.
       where(temp3 .eq. 1) scale_by_complement_lowerSW(:) = .true.       
       status = nf90_inq_varid(ncid,'scale_by_complement_upper', varID)
       status = nf90_get_var(  ncid, varID,temp4)
       scale_by_complement_upperSW(:) = .false.
       where(temp4 .eq. 1) scale_by_complement_upperSW(:) = .true.
       
       ! Close
       status = nf90_close(ncid)
!    endif

    !   
    ! Initialize RRTMGP DDT's...
    !
    ! Shortwave k-distribution data
!$omp critical (load_sw_gas_optics)
    gas_concentrations%gas_name(:) = active_gases_array(:)
    call check_error_msg('sw_gas_optics_init',sw_gas_props%load(gas_concentrations,         &
         gas_namesSW, key_speciesSW, band2gptSW, band_limsSW, press_refSW, press_ref_tropSW,&
         temp_refSW, temp_ref_pSW, temp_ref_tSW, vmr_refSW, kmajorSW, kminor_lowerSW,       &
         kminor_upperSW, gas_minorSW, identifier_minorSW, minor_gases_lowerSW,              &
         minor_gases_upperSW, minor_limits_gpt_lowerSW, minor_limits_gpt_upperSW,           &
         minor_scales_with_density_lowerSW, minor_scales_with_density_upperSW,              &
         scaling_gas_lowerSW, scaling_gas_upperSW, scale_by_complement_lowerSW,             &
         scale_by_complement_upperSW, kminor_start_lowerSW, kminor_start_upperSW,           &
         solar_quietSW, solar_facularSW, solar_sunspotSW, tsi_defaultSW, mg_defaultSW,      &
         sb_defaultSW, rayl_lowerSW, rayl_upperSW))
!$omp end critical (load_sw_gas_optics)

  end subroutine rrtmgp_sw_gas_optics_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_gas_optics_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_gas_optics_run
!! \htmlinclude rrtmgp_sw_gas_optics.html
!!
  subroutine rrtmgp_sw_gas_optics_run(doSWrad, nCol, nLev, ngptsGPsw, nday, idxday,  p_lay, &
       p_lev, toa_src_sw, t_lay, t_lev, gas_concentrations, solcon, sw_optical_props_clrsky,&
       errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doSWrad                 ! Flag to calculate SW irradiances
    integer,intent(in) :: &
         ngptsGPsw,            & ! Number of spectral (g) points.
         nDay,                 & ! Number of daylit points.
         nCol,                 & ! Number of horizontal points
         nLev                    ! Number of vertical levels
    integer,intent(in),dimension(ncol) :: &
         idxday                  ! Indices for daylit points.
    real(kind_phys), dimension(ncol,nLev), intent(in) :: &
         p_lay,                & ! Pressure @ model layer-centers (hPa)
         t_lay                   ! Temperature (K)
    real(kind_phys), dimension(ncol,nLev+1), intent(in) :: &
         p_lev,                & ! Pressure @ model layer-interfaces (hPa)
         t_lev                   ! Temperature @ model levels
    type(ty_gas_concs),intent(in) :: &
         gas_concentrations      ! RRTMGP DDT: trace gas concentrations (vmr)
    real(kind_phys), intent(in) :: &
         solcon                  ! Solar constant

    ! Output
    character(len=*), intent(out) :: &
         errmsg                  ! CCPP error message
    integer,          intent(out) :: &
         errflg                  ! CCPP error code
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_clrsky ! RRTMGP DDT: clear-sky shortwave optical properties, spectral (tau,ssa,g) 
    real(kind_phys), dimension(nCol,ngptsGPsw), intent(out) :: &
         toa_src_sw              ! TOA incident spectral flux (W/m2)

    ! Local variables
    integer :: ij,iGas
    real(kind_phys), dimension(ncol,nLev) :: vmrTemp
    real(kind_phys), dimension(nday,ngptsGPsw) :: toa_src_sw_temp
    type(ty_gas_concs)  ::  gas_concentrations_daylit

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doSWrad) return

    toa_src_sw(:,:) = 0._kind_phys
    if (nDay .gt. 0) then
       !active_gases = gas_concentrations%get_gas_names()
       ! Allocate space
       call check_error_msg('rrtmgp_sw_gas_optics_run_alloc_2str',&
            sw_optical_props_clrsky%alloc_2str(nday, nLev, sw_gas_props))

       gas_concentrations_daylit%ncol = nDay
       gas_concentrations_daylit%nlay = nLev
       allocate(gas_concentrations_daylit%gas_name(gas_concentrations%get_num_gases()))
       allocate(gas_concentrations_daylit%concs(gas_concentrations%get_num_gases()))
       do iGas=1,gas_concentrations%get_num_gases()
          allocate(gas_concentrations_daylit%concs(iGas)%conc(nDay, nLev))
       enddo
       gas_concentrations_daylit%gas_name(:) = active_gases_array(:)

       ! Subset the gas concentrations.
       do iGas=1,gas_concentrations%get_num_gases()
          call check_error_msg('rrtmgp_sw_gas_optics_run_get_vmr',&
               gas_concentrations%get_vmr(trim(gas_concentrations_daylit%gas_name(iGas)),vmrTemp))
          call check_error_msg('rrtmgp_sw_gas_optics_run_set_vmr',&
               gas_concentrations_daylit%set_vmr(trim(gas_concentrations_daylit%gas_name(iGas)),vmrTemp(idxday(1:nday),:)))
       enddo

       ! Call SW gas-optics
       call check_error_msg('rrtmgp_sw_gas_optics_run',sw_gas_props%gas_optics(&
            p_lay(idxday(1:nday),:),          & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(idxday(1:nday),:),          & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(idxday(1:nday),:),          & ! IN  - Temperature @ layer-centers (K)
            gas_concentrations_daylit,        & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
            sw_optical_props_clrsky,          & ! OUT - RRTMGP DDT: Shortwave optical properties, by
                                                !                   spectral point (tau,ssa,g)
            toa_src_sw_temp))                   ! OUT - TOA incident shortwave radiation (spectral)
       toa_src_sw(idxday(1:nday),:) = toa_src_sw_temp

       ! Scale incident flux
       do ij=1,nday
          toa_src_sw(idxday(ij),:) = toa_src_sw(idxday(ij),:)*solcon/ &
                                     sum(toa_src_sw(idxday(ij),:))
       enddo
    endif

  end subroutine rrtmgp_sw_gas_optics_run

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_gas_optics_finalize
  ! #########################################################################################
  subroutine rrtmgp_sw_gas_optics_finalize()
  end subroutine rrtmgp_sw_gas_optics_finalize

end module rrtmgp_sw_gas_optics
 
