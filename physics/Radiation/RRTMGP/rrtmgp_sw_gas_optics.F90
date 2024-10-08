!> \file rrtmgp_sw_gas_optics.F90
!!

!> This module contains a routine to initialize the k-distribution data used
!! by the RRTMGP shortwave radiation scheme.
module rrtmgp_sw_gas_optics
  use machine,                only: kind_phys
  use mo_rte_kind,            only: wl
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,  only: ty_gas_concs
  use radiation_tools,        only: check_error_msg
  use netcdf
#ifdef MPI
  use mpi_f08
#endif

  implicit none
  real(kind_phys),parameter :: &
       tsi_default = 1360.85767381726, &
       mg_default  = 0.1567652,        &
       sb_default  = 902.7126

  ! RRTMGP k-distribution LUTs.
  type(ty_gas_optics_rrtmgp) :: sw_gas_props
  integer :: &
       ntempsSW, npressSW, ngptsSW, nabsorbersSW, nextrabsorbersSW, nminorabsorbersSW, &
       nmixingfracsSW, nlayersSW, nbndsSW, npairsSW, nminor_absorber_intervals_lowerSW,&
       nminor_absorber_intervals_upperSW, ncontributors_lowerSW, ncontributors_upperSW
  integer, dimension(:), allocatable :: &
       kminor_start_lowerSW,              & !< Starting index in the [1, nContributors] vector for a contributor
                                            !< given by \"minor_gases_lower\" (lower atmosphere)
       kminor_start_upperSW                 !< Starting index in the [1, nContributors] vector for a contributor
                                            !< given by \"minor_gases_upper\" (upper atmosphere)
  integer, dimension(:,:), allocatable :: &
       band2gptSW,                        & !< Beginning and ending gpoint for each band
       minor_limits_gpt_lowerSW,          & !< Beginning and ending gpoint for each minor interval in lower atmosphere
       minor_limits_gpt_upperSW             !< Beginning and ending gpoint for each minor interval in upper atmosphere
  integer, dimension(:,:,:), allocatable :: &
       key_speciesSW                        !< Key species pair for each band
  real(kind_phys) :: &
       press_ref_tropSW,                  & !< Reference pressure separating the lower and upper atmosphere [Pa]
       temp_ref_pSW,                      & !< Standard spectroscopic reference pressure [Pa]
       temp_ref_tSW,                      & !< Standard spectroscopic reference temperature [K]
       tsi_defaultSW,                     & !< 
       mg_defaultSW,                      & !< Mean value of Mg2 index over the average solar cycle from the NRLSSI2 model of solar variability
       sb_defaultSW                         !< Mean value of sunspot index over the average solar cycle from the NRLSSI2 model of solar variability
  real(kind_phys), dimension(:), allocatable :: &
       press_refSW,                       & !< Pressures for reference atmosphere; press_ref(# reference layers) [Pa]
       temp_refSW,                        & !< Temperatures for reference atmosphere; temp_ref(# reference layers) [K]
       solar_quietSW,                     & !< Spectrally-dependent quiet sun irradiance from the NRLSSI2 model of solar variability
       solar_facularSW,                   & !< Spectrally-dependent facular term from the NRLSSI2 model of solar variability
       solar_sunspotSW                      !< Spectrally-dependent sunspot term from the NRLSSI2 model of solar variability
  real(kind_phys), dimension(:,:), allocatable :: &
       band_limsSW                          !< Beginning and ending wavenumber [cm -1] for each band
  real(kind_phys), dimension(:,:,:), allocatable :: &
       vmr_refSW,                         & !< Volume mixing ratios for reference atmosphere
       kminor_lowerSW,                    & !< (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to
                                            !< [nTemp x nEta x nContributors] array)
       kminor_upperSW,                    & !< (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to
                                            !< [nTemp x nEta x nContributors] array)
       rayl_lowerSW,                      & !< Stored coefficients due to rayleigh scattering contribution
       rayl_upperSW                         !< Stored coefficients due to rayleigh scattering contribution
  real(kind_phys), dimension(:,:,:,:), allocatable :: &
       kmajorSW                             !< Stored absorption coefficients due to major absorbing gases
  character(len=32),  dimension(:), allocatable :: &
       gas_namesSW,                       & !< Names of absorbing gases
       gas_minorSW,                       & !< Name of absorbing minor gas
       identifier_minorSW,                & !< Unique string identifying minor gas
       minor_gases_lowerSW,               & !< Names of minor absorbing gases in lower atmosphere
       minor_gases_upperSW,               & !< Names of minor absorbing gases in upper atmosphere
       scaling_gas_lowerSW,               & !< Absorption also depends on the concentration of this gas
       scaling_gas_upperSW                  !< Absorption also depends on the concentration of this gas
  logical(wl), dimension(:), allocatable :: &
       minor_scales_with_density_lowerSW, & !< Density scaling is applied to minor absorption coefficients
       minor_scales_with_density_upperSW, & !< Density scaling is applied to minor absorption coefficients
       scale_by_complement_lowerSW,       & !< Absorption is scaled by concentration of scaling_gas (F) or its complement (T)
       scale_by_complement_upperSW          !< Absorption is scaled by concentration of scaling_gas (F) or its complement (T)
contains

!> \section arg_table_rrtmgp_sw_gas_optics_init Argument Table
!! \htmlinclude rrtmgp_sw_gas_optics.html
!!
!! RRTMGP relies heavility on derived-data-types, which contain type-bound procedures 
!! that are referenced frequently throughout the RRTMGP shortwave scheme. The data needed
!! for the correlated k-distribution is also contained within this type. Within this module,
!! the full k-distribution data is read in, reduced by the "active gases" provided, and
!! loaded into the RRTMGP DDT, ty_gas_optics_rrtmgp.
  subroutine rrtmgp_sw_gas_optics_init(rrtmgp_root_dir, rrtmgp_sw_file_gas,              &
       active_gases_array, mpicomm, mpirank, mpiroot, errmsg, errflg)

    ! Inputs
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,  & !< RTE-RRTMGP root directory
         rrtmgp_sw_file_gas  !< RRTMGP file containing K-distribution data
    character(len=*), dimension(:), intent(in) :: &
         active_gases_array  !< List of active gases from namelist as array
    type(MPI_Comm),intent(in) :: &
         mpicomm             !< MPI communicator
    integer,intent(in) :: &
         mpirank,          & !< Current MPI rank
         mpiroot             !< Master MPI rank

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg              !< CCPP error message
    integer,          intent(out) :: &
         errflg              !< CCPP error code

    ! Local variables
    integer :: status, ncid, dimid, varID, mpierr, iChar
    integer,dimension(:),allocatable :: temp1, temp2, temp3, temp4
    character(len=264) :: sw_gas_props_file
    type(ty_gas_concs) :: gas_concs  ! RRTMGP DDT containing active trace gases

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the gfphysics_nml
    sw_gas_props_file   = trim(rrtmgp_root_dir)//trim(rrtmgp_sw_file_gas)

    ! #######################################################################################
    !
    ! Read dimensions for k-distribution fields...
    ! (ONLY master processor(0), if MPI enabled)
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif
       write (*,*) 'Reading RRTMGP shortwave k-distribution metadata ... '

       ! Open file
       status = nf90_open(trim(sw_gas_props_file), NF90_NOWRITE, ncid)

       ! Read dimensions for k-distribution fields
       status = nf90_inq_dimid(ncid, 'temperature', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ntempsSW)
       status = nf90_inq_dimid(ncid, 'pressure', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=npressSW)
       status = nf90_inq_dimid(ncid, 'absorber', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nabsorbersSW)
       status = nf90_inq_dimid(ncid, 'minor_absorber',dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nminorabsorbersSW)
       status = nf90_inq_dimid(ncid, 'absorber_ext', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nextrabsorbersSW)
       status = nf90_inq_dimid(ncid, 'mixing_fraction', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nmixingfracsSW)
       status = nf90_inq_dimid(ncid, 'atmos_layer', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nlayersSW)
       status = nf90_inq_dimid(ncid, 'bnd', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nbndsSW)
       status = nf90_inq_dimid(ncid, 'gpt', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ngptsSW)
       status = nf90_inq_dimid(ncid, 'pair', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=npairsSW)
       status = nf90_inq_dimid(ncid, 'contributors_lower',dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ncontributors_lowerSW)
       status = nf90_inq_dimid(ncid, 'contributors_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ncontributors_upperSW)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_lower', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nminor_absorber_intervals_lowerSW)
       status = nf90_inq_dimid(ncid, 'minor_absorber_intervals_upper', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nminor_absorber_intervals_upperSW)

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
    call mpi_bcast(nbndsSW,                           1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(ngptsSW,                           1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nmixingfracsSW,                    1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(ntempsSW,                          1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(npressSW,                          1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nabsorbersSW,                      1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nextrabsorbersSW,                  1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nminorabsorbersSW,                 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nlayersSW,                         1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(npairsSW,                          1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(ncontributors_upperSW,             1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(ncontributors_lowerSW,             1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nminor_absorber_intervals_upperSW, 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nminor_absorber_intervals_lowerSW, 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
#endif

    ! #######################################################################################
    !
    ! Allocate space for arrays...
    ! (ALL processors)
    !
    ! #######################################################################################
    if (.not. allocated(gas_namesSW))         &
         allocate(gas_namesSW(nabsorbersSW))
    if (.not. allocated(scaling_gas_lowerSW)) &
         allocate(scaling_gas_lowerSW(nminor_absorber_intervals_lowerSW))
    if (.not. allocated(scaling_gas_upperSW)) &
         allocate(scaling_gas_upperSW(nminor_absorber_intervals_upperSW))
    if (.not. allocated(gas_minorSW))         &
         allocate(gas_minorSW(nminorabsorbersSW))
    if (.not. allocated(identifier_minorSW))  &
         allocate(identifier_minorSW(nminorabsorbersSW))
    if (.not. allocated(minor_gases_lowerSW)) &
         allocate(minor_gases_lowerSW(nminor_absorber_intervals_lowerSW))
    if (.not. allocated(minor_gases_upperSW)) &
         allocate(minor_gases_upperSW(nminor_absorber_intervals_upperSW))
    if (.not. allocated(minor_limits_gpt_lowerSW)) &
         allocate(minor_limits_gpt_lowerSW(npairsSW,nminor_absorber_intervals_lowerSW))
    if (.not. allocated(minor_limits_gpt_upperSW)) &
         allocate(minor_limits_gpt_upperSW(npairsSW,nminor_absorber_intervals_upperSW))
    if (.not. allocated(band2gptSW)) &
         allocate(band2gptSW(2,nbndsSW))
    if (.not. allocated(key_speciesSW)) &
         allocate(key_speciesSW(2,nlayersSW,nbndsSW))
    if (.not. allocated(band_limsSW)) &
         allocate(band_limsSW(2,nbndsSW))
    if (.not. allocated(press_refSW)) &
         allocate(press_refSW(npressSW))
    if (.not. allocated(temp_refSW)) &
         allocate(temp_refSW(ntempsSW))
    if (.not. allocated(vmr_refSW)) &
         allocate(vmr_refSW(nlayersSW, nextrabsorbersSW, ntempsSW))
    if (.not. allocated(kminor_lowerSW)) &
         allocate(kminor_lowerSW(ncontributors_lowerSW, nmixingfracsSW, ntempsSW))
    if (.not. allocated(kmajorSW)) &
         allocate(kmajorSW(ngptsSW, nmixingfracsSW,  npressSW+1, ntempsSW))
    if (.not. allocated(kminor_start_lowerSW)) &
         allocate(kminor_start_lowerSW(nminor_absorber_intervals_lowerSW))
    if (.not. allocated(kminor_upperSW)) &
         allocate(kminor_upperSW(ncontributors_upperSW, nmixingfracsSW, ntempsSW))
    if (.not. allocated(kminor_start_upperSW)) &
         allocate(kminor_start_upperSW(nminor_absorber_intervals_upperSW))
    if (.not. allocated(minor_scales_with_density_lowerSW)) &
         allocate(minor_scales_with_density_lowerSW(nminor_absorber_intervals_lowerSW))
    if (.not. allocated(minor_scales_with_density_upperSW)) &
         allocate(minor_scales_with_density_upperSW(nminor_absorber_intervals_upperSW))
    if (.not. allocated(scale_by_complement_lowerSW)) &
         allocate(scale_by_complement_lowerSW(nminor_absorber_intervals_lowerSW))
    if (.not. allocated(scale_by_complement_upperSW)) &
         allocate(scale_by_complement_upperSW(nminor_absorber_intervals_upperSW))
    if (.not. allocated(rayl_upperSW)) &
         allocate(rayl_upperSW(ngptsSW, nmixingfracsSW, ntempsSW))
    if (.not. allocated(rayl_lowerSW)) &
         allocate(rayl_lowerSW(ngptsSW, nmixingfracsSW, ntempsSW))
    if (.not. allocated(solar_quietSW)) &
         allocate(solar_quietSW(ngptsSW))
    if (.not. allocated(solar_facularSW)) &
         allocate(solar_facularSW(ngptsSW))	
    if (.not. allocated(solar_sunspotSW)) &
         allocate(solar_sunspotSW(ngptsSW))
    if (.not. allocated(temp1)) &
         allocate(temp1(nminor_absorber_intervals_lowerSW))
    if (.not. allocated(temp2)) &
         allocate(temp2(nminor_absorber_intervals_upperSW))
    if (.not. allocated(temp3)) &
         allocate(temp3(nminor_absorber_intervals_lowerSW))
    if (.not. allocated(temp4)) &
         allocate(temp4(nminor_absorber_intervals_upperSW))

    ! #######################################################################################
    !
    ! Read in data ...
    ! (ONLY master processor(0), if MPI enabled) 
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif
       write (*,*) 'Reading RRTMGP shortwave k-distribution data ... '
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
       if (status .eq. 0) then
          status = nf90_get_var(  ncid, varID, tsi_defaultSW)
       else
          tsi_defaultSW = tsi_default
       endif
       status = nf90_inq_varid(ncid, 'mg_default', varID)
       if (status .eq. 0) then
          status =	nf90_get_var(  ncid, varID, mg_defaultSW)
       else
          mg_defaultSW = mg_default
       endif
       status = nf90_inq_varid(ncid, 'sb_default', varID)
       if (status .eq. 0) then
          status =	nf90_get_var(  ncid, varID, sb_defaultSW)
       else
          sb_defaultSW = sb_default
       endif
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
    call mpi_bcast(press_ref_tropSW, 1,           MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(temp_ref_pSW,     1,           MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(temp_ref_tSW,     1,           MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(tsi_defaultSW,    1,           MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(mg_defaultSW,     1,           MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(sb_defaultSW,     1,           MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    
    ! Integer arrays
    call mpi_bcast(kminor_start_lowerSW,               &
         size(kminor_start_lowerSW),              MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(kminor_start_upperSW,               &
         size(kminor_start_upperSW),              MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(band2gptSW,                         &
         size(band2gptSW),                        MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(minor_limits_gpt_lowerSW,           &
         size(minor_limits_gpt_lowerSW),          MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(minor_limits_gpt_upperSW,           &
         size(minor_limits_gpt_upperSW),          MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(key_speciesSW,                      &
         size(key_speciesSW),                     MPI_INTEGER,          mpiroot, mpicomm, mpierr)
    
    ! Real arrays
    call mpi_bcast(press_refSW,                        &
         size(press_refSW),                       MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(temp_refSW,                         &
         size(temp_refSW),                        MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(solar_quietSW,                      &
         size(solar_quietSW),                     MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(solar_facularSW,                    &
         size(solar_facularSW),                   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(solar_sunspotSW,                    &
         size(solar_sunspotSW),                   MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(band_limsSW,                        &
         size(band_limsSW),                       MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(vmr_refSW,                          &
         size(vmr_refSW),                         MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(kminor_lowerSW,                     &
         size(kminor_lowerSW),                    MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(kminor_upperSW,                     &
         size(kminor_upperSW),                    MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(rayl_lowerSW,                       &
         size(rayl_lowerSW),                      MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(rayl_upperSW,                       &
         size(rayl_upperSW),                      MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(kmajorSW,                           &
         size(kmajorSW),                          MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    
    ! Characters
    do iChar=1,nabsorbersSW
       call mpi_bcast(gas_namesSW(iChar),              &
         len(gas_namesSW(iChar)),                 MPI_CHARACTER,        mpiroot, mpicomm, mpierr)
    enddo
    do iChar=1,nminorabsorbersSW
       call mpi_bcast(gas_minorSW(iChar),              &
            len(gas_minorSW(iChar)),              MPI_CHARACTER,        mpiroot, mpicomm, mpierr)
       call mpi_bcast(identifier_minorSW(iChar),       &
            len(identifier_minorSW(iChar)),       MPI_CHARACTER,        mpiroot, mpicomm, mpierr)
    enddo
    do iChar=1,nminor_absorber_intervals_lowerSW
       call mpi_bcast(minor_gases_lowerSW(iChar),      &
            len(minor_gases_lowerSW(iChar)),      MPI_CHARACTER,        mpiroot, mpicomm, mpierr)
       call mpi_bcast(scaling_gas_lowerSW(iChar),      &
            len(scaling_gas_lowerSW(iChar)),      MPI_CHARACTER,        mpiroot, mpicomm, mpierr)
    enddo

    do iChar=1,nminor_absorber_intervals_upperSW
       call mpi_bcast(minor_gases_upperSW(iChar),      &
            len(minor_gases_upperSW(iChar)),      MPI_CHARACTER,        mpiroot, mpicomm, mpierr)
       call mpi_bcast(scaling_gas_upperSW(iChar),      &
            len(scaling_gas_upperSW(iChar)),      MPI_CHARACTER,        mpiroot, mpicomm, mpierr)
    enddo
    
    ! Logicals
    call mpi_bcast(minor_scales_with_density_lowerSW,  &
         size(minor_scales_with_density_lowerSW), MPI_LOGICAL,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(minor_scales_with_density_upperSW,  &
         size(minor_scales_with_density_upperSW), MPI_LOGICAL,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(scale_by_complement_lowerSW,        &
         size(scale_by_complement_lowerSW),       MPI_LOGICAL,          mpiroot, mpicomm, mpierr)
    call mpi_bcast(scale_by_complement_upperSW,        &
         size(scale_by_complement_upperSW),       MPI_LOGICAL,          mpiroot, mpicomm, mpierr)

    call mpi_barrier(mpicomm, mpierr)
#endif

    ! #######################################################################################
    !   
    ! Initialize RRTMGP DDT's...
    !
    ! #######################################################################################
    call check_error_msg('rrtmgp_sw_gas_optics_init_gas_concs',gas_concs%init(active_gases_array))
    call check_error_msg('rrtmgp_sw_gas_optics_init_load',sw_gas_props%load(gas_concs,      &
         gas_namesSW, key_speciesSW, band2gptSW, band_limsSW, press_refSW, press_ref_tropSW,&
         temp_refSW, temp_ref_pSW, temp_ref_tSW, vmr_refSW, kmajorSW, kminor_lowerSW,       &
         kminor_upperSW, gas_minorSW, identifier_minorSW, minor_gases_lowerSW,              &
         minor_gases_upperSW, minor_limits_gpt_lowerSW, minor_limits_gpt_upperSW,           &
         minor_scales_with_density_lowerSW, minor_scales_with_density_upperSW,              &
         scaling_gas_lowerSW, scaling_gas_upperSW, scale_by_complement_lowerSW,             &
         scale_by_complement_upperSW, kminor_start_lowerSW, kminor_start_upperSW,           &
         solar_quietSW, solar_facularSW, solar_sunspotSW, tsi_defaultSW, mg_defaultSW,      &
         sb_defaultSW, rayl_lowerSW, rayl_upperSW))

  end subroutine rrtmgp_sw_gas_optics_init
end module rrtmgp_sw_gas_optics
 
