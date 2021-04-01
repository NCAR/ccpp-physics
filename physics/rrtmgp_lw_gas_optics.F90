module rrtmgp_lw_gas_optics
  use machine,               only: kind_phys
  use mo_rte_kind,           only: wl
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs  
  use mo_source_functions,   only: ty_source_func_lw
  use mo_optical_props,      only: ty_optical_props_1scl
  use mo_compute_bc,         only: compute_bc
  use rrtmgp_aux,            only: check_error_msg
  use GFS_rrtmgp_pre,        only: active_gases_array
  use netcdf
#ifdef MPI
  use mpi
#endif

  implicit none

  type(ty_gas_optics_rrtmgp) :: lw_gas_props
  integer :: &
       ntempsLW, npressLW, ngptsLW, nabsorbersLW, nextrabsorbersLW, nminorabsorbersLW,&
       nmixingfracsLW, nlayersLW, nbndsLW, npairsLW, ninternalSourcetempsLW,          &
       nminor_absorber_intervals_lowerLW, nminor_absorber_intervals_upperLW,          &
       ncontributors_lowerLW, ncontributors_upperLW, nfit_coeffsLW
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
  subroutine rrtmgp_lw_gas_optics_init(rrtmgp_root_dir, rrtmgp_lw_file_gas,                 &
       gas_concentrations, mpicomm, mpirank, mpiroot, minGPpres, minGPtemp, errmsg, errflg)

    ! Inputs
    type(ty_gas_concs), intent(inout) :: &
         gas_concentrations  ! RRTMGP DDT: trace gas concentrations (vmr)
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,  & ! RTE-RRTMGP root directory
         rrtmgp_lw_file_gas  ! RRTMGP file containing coefficients used to compute gaseous optical properties
    integer,intent(in) :: &
         mpicomm,          & ! MPI communicator
         mpirank,          & ! Current MPI rank
         mpiroot             ! Master MPI rank
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg              ! CCPP error message
    integer,          intent(out) :: &
         errflg              ! CCPP error code
    real(kind_phys), intent(out) :: &
         minGPtemp,        & ! Minimum temperature allowed by RRTMGP.
         minGPpres           ! Minimum pressure allowed by RRTMGP. 

    ! Local variables
    integer :: ncid, dimID, varID, status, iGas, ierr, ii, mpierr, iChar
    integer,dimension(:),allocatable :: temp1, temp2, temp3, temp4, &
         temp_log_array1, temp_log_array2, temp_log_array3, temp_log_array4
    character(len=264) :: lw_gas_props_file

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the physics_nml
    lw_gas_props_file  = trim(rrtmgp_root_dir)//trim(rrtmgp_lw_file_gas)

    ! #######################################################################################
    !
    ! Read dimensions for k-distribution fields...
    ! (ONLY master processor(0), if MPI enabled)
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif
       write (*,*) 'Reading RRTMGP longwave k-distribution metadata ... '

       ! Open file
       status = nf90_open(trim(lw_gas_props_file), NF90_NOWRITE, ncid)

       ! Read dimensions for k-distribution fields
       status = nf90_inq_dimid(ncid, 'temperature', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ntempsLW)
       status = nf90_inq_dimid(ncid, 'pressure', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = npressLW)
       status = nf90_inq_dimid(ncid, 'absorber', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nabsorbersLW)
       status = nf90_inq_dimid(ncid, 'minor_absorber', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nminorabsorbersLW)
       status = nf90_inq_dimid(ncid, 'absorber_ext', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nextrabsorbersLW)
       status = nf90_inq_dimid(ncid, 'mixing_fraction', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nmixingfracsLW)
       status = nf90_inq_dimid(ncid, 'atmos_layer', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nlayersLW)
       status = nf90_inq_dimid(ncid, 'bnd', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nbndsLW)
       status = nf90_inq_dimid(ncid, 'gpt', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ngptsLW)
       status = nf90_inq_dimid(ncid, 'pair', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = npairsLW)
       status = nf90_inq_dimid(ncid, 'contributors_lower', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ncontributors_lowerLW)
       status = nf90_inq_dimid(ncid, 'contributors_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ncontributors_upperLW)
       status = nf90_inq_dimid(ncid, 'fit_coeffs',  dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nfit_coeffsLW)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_lower', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nminor_absorber_intervals_lowerLW)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = nminor_absorber_intervals_upperLW)
       status = nf90_inq_dimid(ncid, 'temperature_Planck', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = ninternalSourcetempsLW)
#ifdef MPI
    endif ! On master processor

    ! Other processors waiting...
    call mpi_barrier(mpicomm, mpierr)

    ! #######################################################################################
    !
    ! Broadcast dimensions...
    ! (ALL processors)
    !
    ! #######################################################################################
    call mpi_bcast(ntempsLW,                          1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(npressLW,                          1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(ngptsLW,                           1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nabsorbersLW,                      1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nextrabsorbersLW,                  1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nminorabsorbersLW,                 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nmixingfracsLW,                    1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nlayersLW,                         1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nbndsLW,                           1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(npairsLW,                          1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(ninternalSourcetempsLW,            1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nminor_absorber_intervals_lowerLW, 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nminor_absorber_intervals_upperLW, 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(ncontributors_lowerLW,             1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(ncontributors_upperLW,             1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nfit_coeffsLW,                     1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
#endif

    ! Allocate space for arrays
    if (.not. allocated(gas_namesLW))                       &
         allocate(gas_namesLW(nabsorbersLW))
    if (.not. allocated(scaling_gas_lowerLW))               &
         allocate(scaling_gas_lowerLW(nminor_absorber_intervals_lowerLW))
    if (.not. allocated(scaling_gas_upperLW))               &
         allocate(scaling_gas_upperLW(nminor_absorber_intervals_upperLW))
    if (.not. allocated(gas_minorLW))                       &
         allocate(gas_minorLW(nminorabsorbersLW))
    if (.not. allocated(identifier_minorLW))                &
         allocate(identifier_minorLW(nminorabsorbersLW))
    if (.not. allocated(minor_gases_lowerLW))               &
         allocate(minor_gases_lowerLW(nminor_absorber_intervals_lowerLW))
    if (.not. allocated(minor_gases_upperLW))               &
         allocate(minor_gases_upperLW(nminor_absorber_intervals_upperLW))
    if (.not. allocated(minor_limits_gpt_lowerLW))          &
         allocate(minor_limits_gpt_lowerLW(npairsLW, nminor_absorber_intervals_lowerLW))
    if (.not. allocated(minor_limits_gpt_upperLW))          &
         allocate(minor_limits_gpt_upperLW(npairsLW, nminor_absorber_intervals_upperLW))
    if (.not. allocated(band2gptLW))                        &
         allocate(band2gptLW(2, nbndsLW))
    if (.not. allocated(key_speciesLW))                     &
         allocate(key_speciesLW(2, nlayersLW, nbndsLW))
    if (.not. allocated(band_limsLW))                       &
         allocate(band_limsLW(2, nbndsLW))
    if (.not. allocated(press_refLW))                       &
         allocate(press_refLW(npressLW))
    if (.not. allocated(temp_refLW))                        &
         allocate(temp_refLW(ntempsLW))
    if (.not. allocated(vmr_refLW))                         &
         allocate(vmr_refLW(nlayersLW, nextrabsorbersLW, ntempsLW))
    if (.not. allocated(kminor_lowerLW))                    &
         allocate(kminor_lowerLW(ncontributors_lowerLW, nmixingfracsLW, ntempsLW))
    if (.not. allocated(kmajorLW))                          &
         allocate(kmajorLW(ngptsLW, nmixingfracsLW,  npressLW+1, ntempsLW))
    if (.not. allocated(kminor_start_lowerLW))              &
         allocate(kminor_start_lowerLW(nminor_absorber_intervals_lowerLW))
    if (.not. allocated(kminor_upperLW))                    &
         allocate(kminor_upperLW(ncontributors_upperLW, nmixingfracsLW, ntempsLW))
    if (.not. allocated(kminor_start_upperLW))              &
         allocate(kminor_start_upperLW(nminor_absorber_intervals_upperLW))
    if (.not. allocated(optimal_angle_fitLW))               &
         allocate(optimal_angle_fitLW(nfit_coeffsLW, nbndsLW))
    if (.not. allocated(minor_scales_with_density_lowerLW)) &
         allocate(minor_scales_with_density_lowerLW(nminor_absorber_intervals_lowerLW))
    if (.not. allocated(minor_scales_with_density_upperLW)) &
         allocate(minor_scales_with_density_upperLW(nminor_absorber_intervals_upperLW))
    if (.not. allocated(scale_by_complement_lowerLW))       &
         allocate(scale_by_complement_lowerLW(nminor_absorber_intervals_lowerLW))
    if (.not. allocated(scale_by_complement_upperLW))       &
         allocate(scale_by_complement_upperLW(nminor_absorber_intervals_upperLW))
    if (.not. allocated(temp1))                             &
         allocate(temp1(nminor_absorber_intervals_lowerLW))
    if (.not. allocated(temp2))                             &
         allocate(temp2(nminor_absorber_intervals_upperLW))
    if (.not. allocated(temp3))                             &
         allocate(temp3(nminor_absorber_intervals_lowerLW))
    if (.not. allocated(temp4))                             &
         allocate(temp4(nminor_absorber_intervals_upperLW))
    if (.not. allocated(totplnkLW))                         &
         allocate(totplnkLW(ninternalSourcetempsLW, nbndsLW))
    if (.not. allocated(planck_fracLW))                     &
         allocate(planck_fracLW(ngptsLW, nmixingfracsLW, npressLW+1, ntempsLW))

    ! #######################################################################################
    !
    ! Read in data ...
    ! (ONLY master processor(0), if MPI enabled) 
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif
       write (*,*) 'Reading RRTMGP longwave k-distribution data ... '
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

       do ii=1,nminor_absorber_intervals_lowerLW
          if (temp1(ii) .eq. 0) minor_scales_with_density_lowerLW(ii) = .false.
          if (temp1(ii) .eq. 1) minor_scales_with_density_lowerLW(ii) = .true.
          if (temp3(ii) .eq. 0) scale_by_complement_lowerLW(ii)       = .false.
          if (temp3(ii) .eq. 1) scale_by_complement_lowerLW(ii)       = .true.
       enddo 
       do ii=1,nminor_absorber_intervals_upperLW
          if (temp2(ii) .eq. 0) minor_scales_with_density_upperLW(ii) = .false.
          if (temp2(ii) .eq. 1) minor_scales_with_density_upperLW(ii) = .true.
          if (temp4(ii) .eq. 0) scale_by_complement_upperLW(ii)       = .false.
          if (temp4(ii) .eq. 1) scale_by_complement_upperLW(ii)       = .true.
       enddo
#ifdef MPI
    endif ! Master process

    ! Other processors waiting...
    call mpi_barrier(mpicomm, mpierr)

    ! #######################################################################################
    !
    ! Broadcast data... 
    ! (ALL processors)
    !
    ! #######################################################################################

    ! Real scalars
    call mpi_bcast(press_ref_tropLW, 1,      MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(temp_ref_pLW,     1,      MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(temp_ref_tLW,     1,      MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)

    ! Integer arrays
    call mpi_bcast(kminor_start_lowerLW,               &
         size(kminor_start_lowerLW),         MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(kminor_start_upperLW,               &
         size(kminor_start_upperLW),         MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(band2gptLW,                         &
         size(band2gptLW),                   MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(minor_limits_gpt_lowerLW,           &
         size(minor_limits_gpt_lowerLW),     MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(minor_limits_gpt_upperLW,           &
         size(minor_limits_gpt_upperLW),     MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(key_speciesLW,                      &
         size(key_speciesLW),                MPI_INTEGER,          mpiroot, mpicomm, mpierr)

    ! Real arrays
    call mpi_bcast(press_refLW,                        &
         size(press_refLW),                  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(temp_refLW,                         &
         size(temp_refLW),                   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(band_limsLW,                        &
         size(band_limsLW),                  MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(totplnkLW,                          &
         size(totplnkLW),                    MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(optimal_angle_fitLW,                &
         size(optimal_angle_fitLW),          MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(vmr_refLW,                          &
         size(vmr_refLW),                    MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(kminor_lowerLW,                     &
         size(kminor_lowerLW),               MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(kminor_upperLW,                     &
         size(kminor_upperLW),               MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(kmajorLW,                           &
         size(kmajorLW),                     MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(planck_fracLW,                      &
         size(planck_fracLW),                MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)


    ! Characters
    do iChar=1,nabsorbersLW
       call mpi_bcast(gas_namesLW(iChar),              &
         len(gas_namesLW(iChar)),                 MPI_CHARACTER,   mpiroot, mpicomm, mpierr)
    enddo
    do iChar=1,nminorabsorbersLW
       call mpi_bcast(gas_minorLW(iChar),              &
         len(gas_minorLW(iChar)),                 MPI_CHARACTER,   mpiroot, mpicomm, mpierr)
       call mpi_bcast(identifier_minorLW(iChar),       &
         len(identifier_minorLW(iChar)),          MPI_CHARACTER,   mpiroot, mpicomm, mpierr)
    enddo
    do iChar=1,nminor_absorber_intervals_lowerLW
       call mpi_bcast(minor_gases_lowerLW(iChar),      &
         len(minor_gases_lowerLW(iChar)),          MPI_CHARACTER,  mpiroot, mpicomm, mpierr)
       call mpi_bcast(scaling_gas_lowerLW(iChar),      &
         len(scaling_gas_lowerLW(iChar)),          MPI_CHARACTER,  mpiroot, mpicomm, mpierr)
    enddo
    do iChar=1,nminor_absorber_intervals_upperLW
       call mpi_bcast(minor_gases_upperLW(iChar),      &
         len(minor_gases_upperLW(iChar)),          MPI_CHARACTER,  mpiroot, mpicomm, mpierr)
       call mpi_bcast(scaling_gas_upperLW(iChar),      &
         len(scaling_gas_upperLW(iChar)),          MPI_CHARACTER,  mpiroot, mpicomm, mpierr)
    enddo

    ! Logicals
    call mpi_bcast(minor_scales_with_density_lowerLW,  &
         size(minor_scales_with_density_lowerLW),  MPI_LOGICAL,    mpiroot, mpicomm, mpierr)
    call mpi_bcast(minor_scales_with_density_upperLW,  &
         size(minor_scales_with_density_upperLW),  MPI_LOGICAL,    mpiroot, mpicomm, mpierr)
    call mpi_bcast(scale_by_complement_lowerLW,        &
         size(scale_by_complement_lowerLW),        MPI_LOGICAL,    mpiroot, mpicomm, mpierr)
    call mpi_bcast(scale_by_complement_upperLW,        &
         size(scale_by_complement_upperLW),        MPI_LOGICAL,    mpiroot, mpicomm, mpierr)

    call mpi_barrier(mpicomm, mpierr)
#endif

    ! #######################################################################################
    !   
    ! Initialize RRTMGP DDT's...
    !
    ! #######################################################################################
    gas_concentrations%gas_name(:) = active_gases_array(:)
    call check_error_msg('rrtmgp_lw_gas_optics_init',lw_gas_props%load(gas_concentrations,  &
         gas_namesLW, key_speciesLW, band2gptLW, band_limsLW, press_refLW, press_ref_tropLW,&
         temp_refLW,  temp_ref_pLW, temp_ref_tLW, vmr_refLW, kmajorLW, kminor_lowerLW,      &
         kminor_upperLW, gas_minorLW, identifier_minorLW, minor_gases_lowerLW,              &
         minor_gases_upperLW, minor_limits_gpt_lowerLW, minor_limits_gpt_upperLW,           &
         minor_scales_with_density_lowerLW, minor_scales_with_density_upperLW,              &
         scaling_gas_lowerLW, scaling_gas_upperLW, scale_by_complement_lowerLW,             &
         scale_by_complement_upperLW, kminor_start_lowerLW, kminor_start_upperLW, totplnkLW,&
         planck_fracLW, rayl_lowerLW, rayl_upperLW, optimal_angle_fitLW))

    ! The minimum pressure allowed in GP RTE calculations. Used to bound uppermost layer
    ! temperature (GFS_rrtmgp_pre.F90)
    minGPpres = lw_gas_props%get_press_min()
    minGPtemp = lw_gas_props%get_temp_min() 

  end subroutine rrtmgp_lw_gas_optics_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_gas_optics_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_gas_optics_run
!! \htmlinclude rrtmgp_lw_gas_optics_run.html
!!
  subroutine rrtmgp_lw_gas_optics_run(doLWrad, nCol, nLev, p_lay, p_lev, t_lay, t_lev, tsfg,&
       gas_concentrations, lw_optical_props_clrsky, sources,  errmsg, errflg)

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
