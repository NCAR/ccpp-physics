module rrtmgp_lw_gas_optics
  use machine,               only: kind_phys
  use mo_rte_kind,           only: wl
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs  
  use mo_source_functions,   only: ty_source_func_lw
  use mo_optical_props,      only: ty_optical_props_1scl
  use mo_compute_bc,         only: compute_bc
  use rrtmgp_aux,            only: check_error_msg
  use netcdf
  implicit none

  type(ty_gas_optics_rrtmgp) :: lw_gas_props
  integer, dimension(:), allocatable :: &
       kminor_start_lowerLW,              & ! Starting index in the [1, nContributors] vector for a contributor
                                            ! given by \"minor_gases_lower\" (lower atmosphere)
       kminor_start_upperLW                 ! Starting index in the [1, nContributors] vector for a contributor
                                            ! given by \"minor_gases_upper\" (upper atmosphere)
  integer, dimension(:,:), allocatable :: &
       band2gptLW,                        & ! Beginning and ending gpoint for each band
       minor_limits_gpt_lowerLW,          & ! Beginning and ending gpoint for each minor interval in lower atmosphere
       minor_limits_gpt_upperLW             ! Beginning and ending gpoint for each minor interval in upper atmosphere
  integer, dimension(:,:,:), allocatable :: &
       key_speciesLW                        ! Key species pair for each band
  real(kind_phys) :: &
       press_ref_tropLW,                  & ! Reference pressure separating the lower and upper atmosphere [Pa]
       temp_ref_pLW,                      & ! Standard spectroscopic reference pressure [Pa]
       temp_ref_tLW                         ! Standard spectroscopic reference temperature [K]
  real(kind_phys), dimension(:), allocatable :: &
       press_refLW,                       & ! Pressures for reference atmosphere; press_ref(# reference layers) [Pa]
       temp_refLW                           ! Temperatures for reference atmosphere; temp_ref(# reference layers) [K]
  real(kind_phys), dimension(:,:), allocatable :: &
       band_limsLW,                       & ! Beginning and ending wavenumber [cm -1] for each band
       totplnkLW,                         & ! Integrated Planck function by band
       optimal_angle_fitLW
  real(kind_phys), dimension(:,:,:), allocatable :: &
       vmr_refLW,                         & ! volume mixing ratios for reference atmospherer
       kminor_lowerLW,                    & ! (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to
                                            ! [nTemp x nEta x nContributors] array)
       kminor_upperLW,                    & ! (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to
                                            ! [nTemp x nEta x nContributors] array)
       rayl_lowerLW,                      & ! Not used in LW, rather allocated(rayl_lower) is used
       rayl_upperLW                         ! Not used in LW, rather allocated(rayl_upper) is used
  real(kind_phys), dimension(:,:,:,:), allocatable :: &
       kmajorLW,                          & ! Stored absorption coefficients due to major absorbing gases
       planck_fracLW                        ! Planck fractions   
  character(len=32),  dimension(:), allocatable :: &
       gas_namesLW,                       & ! Names of absorbing gases
       gas_minorLW,                       & ! Name of absorbing minor gas
       identifier_minorLW,                & ! Unique string identifying minor gas
       minor_gases_lowerLW,               & ! Names of minor absorbing gases in lower atmosphere
       minor_gases_upperLW,               & ! Names of minor absorbing gases in upper atmosphere
       scaling_gas_lowerLW,               & ! Absorption also depends on the concentration of this gas
       scaling_gas_upperLW                  ! Absorption also depends on the concentration of this gas
  logical(wl), dimension(:), allocatable :: &
       minor_scales_with_density_lowerLW, & ! Density scaling is applied to minor absorption coefficients
       minor_scales_with_density_upperLW, & ! Density scaling is applied to minor absorption coefficients
       scale_by_complement_lowerLW,       & ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T)
       scale_by_complement_upperLW          ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T)

contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_gas_optics_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_gas_optics_init
!! \htmlinclude rrtmgp_lw_gas_optics_init.html
!!
  subroutine rrtmgp_lw_gas_optics_init(rrtmgp_root_dir, rrtmgp_lw_file_gas, gas_concentrations,&
       nCol, nLev, mpicomm, mpirank, mpiroot, minGPpres, errmsg, errflg)

    ! Inputs
    type(ty_gas_concs), intent(in) :: &
         gas_concentrations  ! RRTMGP DDT: trace gas concentrations (vmr)
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,  & ! RTE-RRTMGP root directory
         rrtmgp_lw_file_gas  ! RRTMGP file containing coefficients used to compute gaseous optical properties
    integer,intent(in) :: &
         nCol,             & ! Number of horizontal points
         nLev,             & ! Number of vertical levels
         mpicomm,          & ! MPI communicator
         mpirank,          & ! Current MPI rank
         mpiroot             ! Master MPI rank
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                  ! CCPP error message
    integer,          intent(out) :: &
         errflg                  ! CCPP error code
    real(kind_phys), intent(out) :: &
         minGPpres               ! Minimum pressure allowed by RRTMGP. 
    ! Dimensions
    integer :: &
         ntemps, npress, ngpts, nabsorbers, nextrabsorbers, nminorabsorbers,&    
         nmixingfracs, nlayers, nbnds, npairs, ninternalSourcetemps,           &
         nminor_absorber_intervals_lower, nminor_absorber_intervals_upper,     &
         ncontributors_lower, ncontributors_upper,nfit_coeffs

    ! Local variables
    integer :: ncid, dimID, varID, status, iGas, ierr, ii
    integer,dimension(:),allocatable :: temp1, temp2, temp3, temp4, &
         temp_log_array1, temp_log_array2, temp_log_array3, temp_log_array4
    character(len=264) :: lw_gas_props_file

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the physics_nml
    lw_gas_props_file  = trim(rrtmgp_root_dir)//trim(rrtmgp_lw_file_gas)

    ! On master processor only...
!    if (mpirank .eq. mpiroot) then
       ! Open file
       status = nf90_open(trim(lw_gas_props_file), NF90_NOWRITE, ncid)

       ! Read dimensions for k-distribution fields
       status = nf90_inq_dimid(ncid, 'temperature', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ntemps)
       status = nf90_inq_dimid(ncid, 'pressure', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = npress)
       status = nf90_inq_dimid(ncid, 'absorber', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nabsorbers)
       status = nf90_inq_dimid(ncid, 'minor_absorber', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nminorabsorbers)
       status = nf90_inq_dimid(ncid, 'absorber_ext', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nextrabsorbers)
       status = nf90_inq_dimid(ncid, 'mixing_fraction', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nmixingfracs)
       status = nf90_inq_dimid(ncid, 'atmos_layer', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nlayers)
       status = nf90_inq_dimid(ncid, 'bnd', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nbnds)
       status = nf90_inq_dimid(ncid, 'gpt', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ngpts)
       status = nf90_inq_dimid(ncid, 'pair', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = npairs)
       status = nf90_inq_dimid(ncid, 'contributors_lower', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ncontributors_lower)
       status = nf90_inq_dimid(ncid, 'contributors_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ncontributors_upper)
       status = nf90_inq_dimid(ncid, 'fit_coeffs',  dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nfit_coeffs)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_lower', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nminor_absorber_intervals_lower)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nminor_absorber_intervals_upper)
       status = nf90_inq_dimid(ncid, 'temperature_Planck', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ninternalSourcetemps)

       ! Allocate space for arrays
       allocate(gas_namesLW(nabsorbers))
       allocate(scaling_gas_lowerLW(nminor_absorber_intervals_lower))
       allocate(scaling_gas_upperLW(nminor_absorber_intervals_upper))
       allocate(gas_minorLW(nminorabsorbers))
       allocate(identifier_minorLW(nminorabsorbers))
       allocate(minor_gases_lowerLW(nminor_absorber_intervals_lower))
       allocate(minor_gases_upperLW(nminor_absorber_intervals_upper))
       allocate(minor_limits_gpt_lowerLW(npairs,nminor_absorber_intervals_lower))
       allocate(minor_limits_gpt_upperLW(npairs,nminor_absorber_intervals_upper))
       allocate(band2gptLW(2,nbnds))
       allocate(key_speciesLW(2,nlayers,nbnds))
       allocate(band_limsLW(2,nbnds))
       allocate(press_refLW(npress))
       allocate(temp_refLW(ntemps))
       allocate(vmr_refLW(nlayers, nextrabsorbers, ntemps))
       allocate(kminor_lowerLW(ncontributors_lower, nmixingfracs, ntemps))
       allocate(kmajorLW(ngpts, nmixingfracs,  npress+1, ntemps))
       allocate(kminor_start_lowerLW(nminor_absorber_intervals_lower))
       allocate(kminor_upperLW(ncontributors_upper, nmixingfracs, ntemps))
       allocate(kminor_start_upperLW(nminor_absorber_intervals_upper))
       allocate(optimal_angle_fitLW(nfit_coeffs,nbnds))
       allocate(minor_scales_with_density_lowerLW(nminor_absorber_intervals_lower))
       allocate(minor_scales_with_density_upperLW(nminor_absorber_intervals_upper))
       allocate(scale_by_complement_lowerLW(nminor_absorber_intervals_lower))
       allocate(scale_by_complement_upperLW(nminor_absorber_intervals_upper))
       allocate(temp1(nminor_absorber_intervals_lower))
       allocate(temp2(nminor_absorber_intervals_upper))
       allocate(temp3(nminor_absorber_intervals_lower))
       allocate(temp4(nminor_absorber_intervals_upper))
       allocate(totplnkLW(ninternalSourcetemps, nbnds))
       allocate(planck_fracLW(ngpts, nmixingfracs, npress+1, ntemps))

       ! Read in fields from file
       if (mpirank==mpiroot) write (*,*) 'Reading RRTMGP longwave k-distribution data ... '
       status = nf90_inq_varid(ncid, 'gas_names', varID)
       status = nf90_get_var(  ncid, varID, gas_namesLW)
       status = nf90_inq_varid(ncid, 'scaling_gas_lower', varID)
       status = nf90_get_var(  ncid, varID, scaling_gas_lowerLW)
       status = nf90_inq_varid(ncid, 'scaling_gas_upper', varID)
       status = nf90_get_var(  ncid, varID, scaling_gas_upperLW)
       status = nf90_inq_varid(ncid, 'gas_minor', varID)
       status = nf90_get_var(  ncid, varID, gas_minorLW)
       status = nf90_inq_varid(ncid, 'identifier_minor', varID)
       status = nf90_get_var(  ncid, varID, identifier_minorLW)
       status = nf90_inq_varid(ncid, 'minor_gases_lower', varID)
       status = nf90_get_var(  ncid, varID, minor_gases_lowerLW)
       status = nf90_inq_varid(ncid, 'minor_gases_upper', varID)
       status = nf90_get_var(  ncid, varID, minor_gases_upperLW)
       status = nf90_inq_varid(ncid, 'minor_limits_gpt_lower', varID)
       status = nf90_get_var(  ncid, varID, minor_limits_gpt_lowerLW)
       status = nf90_inq_varid(ncid, 'minor_limits_gpt_upper', varID)
       status = nf90_get_var(  ncid, varID, minor_limits_gpt_upperLW)
       status = nf90_inq_varid(ncid, 'bnd_limits_gpt', varID)
       status = nf90_get_var(  ncid, varID, band2gptLW)
       status = nf90_inq_varid(ncid, 'key_species', varID)
       status = nf90_get_var(  ncid, varID, key_speciesLW)
       status = nf90_inq_varid(ncid, 'bnd_limits_wavenumber', varID)
       status = nf90_get_var(  ncid, varID, band_limsLW)
       status = nf90_inq_varid(ncid, 'press_ref', varID)
       status = nf90_get_var(  ncid, varID, press_refLW)
       status = nf90_inq_varid(ncid, 'temp_ref', varID)
       status = nf90_get_var(  ncid, varID, temp_refLW)
       status = nf90_inq_varid(ncid, 'absorption_coefficient_ref_P', varID)
       status = nf90_get_var(  ncid, varID, temp_ref_pLW)
       status = nf90_inq_varid(ncid, 'absorption_coefficient_ref_T', varID)
       status = nf90_get_var(  ncid, varID, temp_ref_tLW)
       status = nf90_inq_varid(ncid, 'press_ref_trop', varID)
       status = nf90_get_var(  ncid, varID, press_ref_tropLW)
       status = nf90_inq_varid(ncid, 'kminor_lower', varID)
       status = nf90_get_var(  ncid, varID, kminor_lowerLW)
       status = nf90_inq_varid(ncid, 'kminor_upper', varID)
       status = nf90_get_var(  ncid, varID, kminor_upperLW)
       status = nf90_inq_varid(ncid, 'vmr_ref', varID)
       status = nf90_get_var(  ncid, varID, vmr_refLW)
       status = nf90_inq_varid(ncid, 'optimal_angle_fit',varID)
       status = nf90_get_var(  ncid, varID, optimal_angle_fitLW)
       status = nf90_inq_varid(ncid, 'kmajor', varID)
       status = nf90_get_var(  ncid, varID, kmajorLW)
       status = nf90_inq_varid(ncid, 'kminor_start_lower', varID)
       status = nf90_get_var(  ncid, varID, kminor_start_lowerLW)
       status = nf90_inq_varid(ncid, 'kminor_start_upper', varID)
       status = nf90_get_var(  ncid, varID, kminor_start_upperLW)
       status = nf90_inq_varid(ncid, 'totplnk', varID)
       status = nf90_get_var(  ncid, varID, totplnkLW)
       status = nf90_inq_varid(ncid, 'plank_fraction', varID)
       status = nf90_get_var(  ncid, varID, planck_fracLW)

       ! Logical fields are read in as integers and then converted to logicals.
       status = nf90_inq_varid(ncid,'minor_scales_with_density_lower', varID)
       status = nf90_get_var(  ncid, varID,temp1)
       status = nf90_inq_varid(ncid,'minor_scales_with_density_upper', varID)
       status = nf90_get_var(  ncid, varID,temp2)
       status = nf90_inq_varid(ncid,'scale_by_complement_lower', varID)
       status = nf90_get_var(  ncid, varID,temp3)
       status = nf90_inq_varid(ncid,'scale_by_complement_upper', varID)
       status = nf90_get_var(  ncid, varID,temp4)
       status = nf90_close(ncid)

       do ii=1,nminor_absorber_intervals_lower
          if (temp1(ii) .eq. 0) minor_scales_with_density_lowerLW(ii) = .false.
          if (temp1(ii) .eq. 1) minor_scales_with_density_lowerLW(ii) = .true.
          if (temp3(ii) .eq. 0) scale_by_complement_lowerLW(ii)       = .false.
          if (temp3(ii) .eq. 1) scale_by_complement_lowerLW(ii)       = .true.
       enddo 
       do ii=1,nminor_absorber_intervals_upper
          if (temp2(ii) .eq. 0) minor_scales_with_density_upperLW(ii) = .false.
          if (temp2(ii) .eq. 1) minor_scales_with_density_upperLW(ii) = .true.
          if (temp4(ii) .eq. 0) scale_by_complement_upperLW(ii)       = .false.
          if (temp4(ii) .eq. 1) scale_by_complement_upperLW(ii)       = .true.
       enddo
!    endif

    !
    ! Initialize RRTMGP DDT's...
    !
!$omp critical (load_lw_gas_optics)
    ! Longwave k-distribution data.
    call check_error_msg('rrtmgp_lw_gas_optics_init',lw_gas_props%load(gas_concentrations,  &
         gas_namesLW, key_speciesLW, band2gptLW, band_limsLW, press_refLW, press_ref_tropLW,&
         temp_refLW,  temp_ref_pLW, temp_ref_tLW, vmr_refLW, kmajorLW, kminor_lowerLW,      &
         kminor_upperLW, gas_minorLW, identifier_minorLW, minor_gases_lowerLW,              &
         minor_gases_upperLW, minor_limits_gpt_lowerLW, minor_limits_gpt_upperLW,           &
         minor_scales_with_density_lowerLW, minor_scales_with_density_upperLW,              &
         scaling_gas_lowerLW, scaling_gas_upperLW, scale_by_complement_lowerLW,             &
         scale_by_complement_upperLW, kminor_start_lowerLW, kminor_start_upperLW, totplnkLW,&
         planck_fracLW, rayl_lowerLW, rayl_upperLW, optimal_angle_fitLW))
!$omp end critical (load_lw_gas_optics) 

    ! The minimum pressure allowed in GP RTE calculations. Used to bound uppermost layer
    ! temperature (GFS_rrtmgp_pre.F90)
    minGPpres = lw_gas_props%get_press_min()

  end subroutine rrtmgp_lw_gas_optics_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_gas_optics_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_gas_optics_run
!! \htmlinclude rrtmgp_lw_gas_optics_run.html
!!
  subroutine rrtmgp_lw_gas_optics_run(doLWrad, nCol, nLev, p_lay, p_lev, t_lay,&
       t_lev, tsfg, gas_concentrations, lw_optical_props_clrsky, sources,  errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doLWrad                 ! Flag to calculate LW irradiances
    integer,intent(in) :: &
         ncol,                &  ! Number of horizontal points
         nLev                    ! Number of vertical levels
    real(kind_phys), dimension(ncol,nLev), intent(in) :: &
         p_lay,                & ! Pressure @ model layer-centers (hPa)
         t_lay                   ! Temperature (K)
    real(kind_phys), dimension(ncol,nLev+1), intent(in) :: &
         p_lev,                & ! Pressure @ model layer-interfaces (hPa)
         t_lev                   ! Temperature @ model levels
    real(kind_phys), dimension(ncol), intent(in) :: &
         tsfg                    ! Surface ground temperature (K)
    type(ty_gas_concs),intent(in) :: &
         gas_concentrations      ! RRTMGP DDT: trace gas concentrations (vmr)

    ! Output
    character(len=*), intent(out) :: &
         errmsg                  ! CCPP error message
    integer,          intent(out) :: &
         errflg                  ! CCPP error code
    type(ty_optical_props_1scl),intent(out) :: &
         lw_optical_props_clrsky ! RRTMGP DDT: longwave clear-sky radiative properties
    type(ty_source_func_lw),intent(out) :: &
         sources                 ! RRTMGP DDT: longwave source functions

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return

    call check_error_msg('rrtmgp_lw_gas_optics_run',&
         lw_optical_props_clrsky%alloc_1scl(ncol, nLev, lw_gas_props))
    call check_error_msg('rrtmgp_lw_gas_optics_run',&
         sources%alloc(ncol, nLev, lw_gas_props))

    ! Gas-optics 
    call check_error_msg('rrtmgp_lw_gas_optics_run',lw_gas_props%gas_optics(&
         p_lay,                   & ! IN  - Pressure @ layer-centers (Pa)
         p_lev,                   & ! IN  - Pressure @ layer-interfaces (Pa)
         t_lay,                   & ! IN  - Temperature @ layer-centers (K)
         tsfg,                    & ! IN  - Skin-temperature (K)
         gas_concentrations,      & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
         lw_optical_props_clrsky, & ! OUT - RRTMGP DDT: longwave optical properties
         sources,                 & ! OUT - RRTMGP DDT: source functions
         tlev=t_lev))               ! IN  - Temperature @ layer-interfaces (K) (optional)

  end subroutine rrtmgp_lw_gas_optics_run

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_gas_optics_finalize
  ! #########################################################################################
  subroutine rrtmgp_lw_gas_optics_finalize()
  end subroutine rrtmgp_lw_gas_optics_finalize
  
end module rrtmgp_lw_gas_optics
