! ###########################################################################################
! ###########################################################################################
module rrtmgp_sw
  use GFS_typedefs,              only: GFS_control_type
  use physparam,                 only: iovrsw, icldflg, iswcliq, isubcsw
  use machine,                   only: kind_phys
  use mo_rte_kind,               only: wl
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,     only: ty_gas_concs
  use mo_fluxes,                 only: ty_fluxes_broadband
  use mo_optical_props,          only: ty_optical_props_1scl,ty_optical_props_2str  
  use mo_rte_sw,                 only: rte_sw
  use mo_heating_rates,          only: compute_heating_rate
  use mo_rrtmgp_constants,       only: grav, avogad
  use module_radsw_parameters,   only: topfsw_type, sfcfsw_type, profsw_type, cmpfsw_type


  implicit none

  ! Parameters
  integer,parameter :: nGases = 6
  real(kind_phys),parameter :: epsilon=1.0e-6
  character(len=3),parameter, dimension(nGases) :: &
       active_gases = (/ 'h2o', 'co2', 'o3 ', 'n2o', 'ch4', 'o2 '/)

  ! Molecular weight ratios (for converting mmr to vmr)
  real(kind_phys), parameter :: &
       amd   = 28.9644_kind_phys,  & ! Molecular weight of dry-air     (g/mol)
       amw   = 18.0154_kind_phys,  & ! Molecular weight of water vapor (g/mol)
       amo3  = 47.9982_kind_phys,  & ! Modelular weight of ozone       (g/mol)
       amdw  = amd/amw,            & ! Molecular weight of dry air / water vapor
       amdo3 = amd/amo3              ! Molecular weight of dry air / ozone
  !
  real (kind_phys), parameter :: &
       s0 = 1368.22                  ! Solar constant (W/m2)

  ! Logical flags for optional output fields in rrtmgp_sw_run(), default=.false.
  logical :: &
       l_AllSky_HR_byband  = .false., & ! 2D [ncol,nlay] all-sky heating rates, in each band [ncol,nlay,nBandsSW]?
       l_ClrSky_HR         = .false., & ! 2D [ncol,nlay] clear-sky heating rate?
       l_fluxes2D          = .false., & ! 2D [ncol,nlay] radiative fluxes *Note* fluxes is a DDT w/ 4 fields.
       l_sfcFluxes1D       = .false.    ! 1D [ncol] surface fluxes  *Note* fluxes is a DDT w/ 6 fields.

  ! Module parameters (set during rrtmgp_sw_init())
  integer :: &
       nGptsSW,            & ! Number of SW spectral g-points
       nBandsSW,           & ! Number of SW bands
       ipsdsw0               ! Initial seed for McICA

  ! Classes used by rte+rrtmgp
  type(ty_gas_optics_rrtmgp) :: &
       kdist_sw_clr
  type(ty_gas_concs)  :: &
       gas_concs_sw

  public rrtmgp_sw_init, rrtmgp_sw_run, rrtmgp_sw_finalize
contains
  ! #########################################################################################
  ! rrtmgp_sw_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_init Argument Table
!! | local_name      | standard_name             | long_name                                               | units | rank | type             |    kind   | intent | optional |
!! |-----------------|---------------------------|---------------------------------------------------------|-------|------|------------------|-----------|--------|----------|
!! | Model           | GFS_control_type_instance | Fortran DDT containing FV3-GFS model control parameters | DDT   |    0 | GFS_control_type |           | in     | F        |
!! | mpirank         | mpi_rank                  | current MPI rank                                        | index |    0 | integer          |           | in     | F        |
!! | mpiroot         | mpi_root                  | master MPI rank                                         | index |    0 | integer          |           | in     | F        |
!! | mpicomm         | mpi_comm                  | MPI communicator                                        | index |    0 | integer          |           | in     | F        |
!! | errmsg          | ccpp_error_message        | error message for error handling in CCPP                | none  |    0 | character        | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag           | error flag for error handling in CCPP                   | flag  |    0 | integer          |           | out    | F        |
!!
  ! #########################################################################################
  subroutine rrtmgp_sw_init(Model,mpicomm, mpirank, mpiroot, errmsg, errflg)
    use netcdf
#ifdef MPI
    use mpi
#endif
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model      ! DDT containing model control parameters
    integer,intent(in) :: &
         mpicomm, & ! MPI communicator
         mpirank, & ! Current MPI rank
         mpiroot    ! Master MPI rank
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer,          intent(out) :: &
         errflg     ! Error code

    ! Fields from the K-distribution files
    ! Variables that will be passed to gas_optics%load()
    integer, dimension(:), allocatable :: &
         kminor_start_lower_sw,              & ! used by RRTGMP gas optics 
         kminor_start_upper_sw                 ! used by RRTGMP gas optics 
    integer, dimension(:,:), allocatable :: &
         band2gpt_sw,                        & ! used by RRTGMP gas optics 
         minor_limits_gpt_lower_sw,          & ! used by RRTGMP gas optics 
         minor_limits_gpt_upper_sw             ! used by RRTGMP gas optics 
    integer, dimension(:,:,:), allocatable :: &
         key_species_sw                        ! used by RRTGMP gas optics 
    real(kind_phys) :: &
         press_ref_trop_sw,                  & ! used by RRTGMP gas optics 
         temp_ref_p_sw,                      & ! used by RRTGMP gas optics 
         temp_ref_t_sw                         ! used by RRTGMP gas optics 
    real(kind_phys), dimension(:), allocatable :: &
         press_ref_sw,                       & ! used by RRTGMP gas optics 
         temp_ref_sw,                        & ! used by RRTGMP gas optics 
         solar_source_sw                       ! used by RRTGMP gas optics 
    real(kind_phys), dimension(:,:), allocatable :: &
         band_lims_sw                          ! used by RRTGMP gas optics 

    real(kind_phys), dimension(:,:,:), allocatable :: &
         vmr_ref_sw,                         & ! used by RRTGMP gas optics 
         kminor_lower_sw,                    & ! used by RRTGMP gas optics 
         kminor_upper_sw,                    & ! used by RRTGMP gas optics 
         rayl_lower_sw,                      & ! used by RRTGMP gas optics 
         rayl_upper_sw                         ! used by RRTGMP gas optics 
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         kmajor_sw                             ! used by RRTGMP gas optics 
    character(len=32),  dimension(:), allocatable :: &
         gas_names_sw,                       & ! used by RRTGMP gas optics 
         gas_minor_sw,                       & ! used by RRTGMP gas optics 
         identifier_minor_sw,                & ! used by RRTGMP gas optics 
         minor_gases_lower_sw,               & ! used by RRTGMP gas optics 
         minor_gases_upper_sw,               & ! used by RRTGMP gas optics 
         scaling_gas_lower_sw,               & ! used by RRTGMP gas optics 
         scaling_gas_upper_sw                  ! used by RRTGMP gas optics 
    logical(wl), dimension(:), allocatable :: &
         minor_scales_with_density_lower_sw, & ! used by RRTGMP gas optics 
         minor_scales_with_density_upper_sw, & ! used by RRTGMP gas optics 
         scale_by_complement_lower_sw,       & ! used by RRTGMP gas optics 
         scale_by_complement_upper_sw          ! used by RRTGMP gas optics 
    ! Dimensions (to be broadcast across all processors)
    integer :: &
         ntemps_sw,                          & ! used by RRTGMP gas optics 
         npress_sw,                          & ! used by RRTGMP gas optics 
         nabsorbers_sw,                      & ! used by RRTGMP gas optics 
         nextrabsorbers_sw,                  & ! used by RRTGMP gas optics 
         nminorabsorbers_sw,                 & ! used by RRTGMP gas optics 
         nmixingfracs_sw,                    & ! used by RRTGMP gas optics 
         nlayers_sw,                         & ! used by RRTGMP gas optics 
         nbnds_sw,                           & ! used by RRTGMP gas optics 
         ngpts_sw,                           & ! used by RRTGMP gas optics 
         npairs_sw,                          & ! used by RRTGMP gas optics 
         nminor_absorber_intervals_lower_sw, & ! used by RRTGMP gas optics 
         nminor_absorber_intervals_upper_sw, & ! used by RRTGMP gas optics 
         ncontributors_lower_sw,             & ! used by RRTGMP gas optics 
         ncontributors_upper_sw                ! used by RRTGMP gas optics 

    ! Local variables
    integer :: status,ncid_sw,dimid,varID,ij,iGas
    character(len=264) :: kdist_file
    integer,dimension(:),allocatable :: temp1,temp2,temp3,temp4,temp_log_array1,&
         temp_log_array2, temp_log_array3, temp_log_array4

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Ensure that requested cloud overlap is reasonable.
    if ( iovrsw<0 .or. iovrsw>3 ) then
       print *,'  *** Error in specification of cloud overlap flag',   &
            ' IOVRSW=',iovrsw,' in RSWINIT !!'
       stop
    endif

    ! Check cloud flags for consistency.
    if ((icldflg == 0 .and. iswcliq /= 0) .or.                        &
         (icldflg == 1 .and. iswcliq == 0)) then
       print *,'  *** Model cloud scheme inconsistent with SW',        &
            ' radiation cloud radiative property setup !!'
       stop
    endif
    if ( isubcsw==0 .and. iovrsw>2 ) then
       print *,'  *** IOVRSW=',iovrsw,' is not available for ISUBCSW=0 setting!!'
       print *,'      The program will use maximum/random overlap instead.'
       iovrsw = 1
    endif

    ! Filenames are set in the gfs_physics_nml (scm/src/GFS_typedefs.F90)
    kdist_file      = trim(Model%rrtmgp_root)//trim(Model%kdist_sw_file_gas)

    ! Read dimensions for k-distribution fields (only on master processor(0))
    if (mpirank .eq. mpiroot) then
       if(nf90_open(trim(kdist_file), NF90_WRITE, ncid_sw) .eq. NF90_NOERR) then
          status = nf90_inq_dimid(ncid_sw, 'temperature', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=ntemps_sw)
          status = nf90_inq_dimid(ncid_sw, 'pressure', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=npress_sw)
          status = nf90_inq_dimid(ncid_sw, 'absorber', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=nabsorbers_sw)
          status = nf90_inq_dimid(ncid_sw, 'minor_absorber', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=nminorabsorbers_sw)
          status = nf90_inq_dimid(ncid_sw, 'absorber_ext', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=nextrabsorbers_sw)
          status = nf90_inq_dimid(ncid_sw, 'mixing_fraction', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=nmixingfracs_sw)
          status = nf90_inq_dimid(ncid_sw, 'atmos_layer', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=nlayers_sw)
          status = nf90_inq_dimid(ncid_sw, 'bnd', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=nbnds_sw)
          status = nf90_inq_dimid(ncid_sw, 'gpt', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=ngpts_sw)
          status = nf90_inq_dimid(ncid_sw, 'pair', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=npairs_sw)
          status = nf90_inq_dimid(ncid_sw, 'contributors_lower', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=ncontributors_lower_sw)
          status = nf90_inq_dimid(ncid_sw, 'contributors_upper', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=ncontributors_upper_sw)
          status = nf90_inq_dimid(ncid_sw, 'minor_absorber_intervals_lower', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=nminor_absorber_intervals_lower_sw)
          status = nf90_inq_dimid(ncid_sw, 'minor_absorber_intervals_upper', dimid)
          status = nf90_inquire_dimension(ncid_sw, dimid, len=nminor_absorber_intervals_upper_sw)
          status = nf90_close(ncid_sw)
       endif
    endif

    ! Broadcast dimensions to all processors
#ifdef MPI
    call MPI_BCAST(ntemps_sw,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(npress_sw,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nabsorbers_sw,                      1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nminorabsorbers_sw,                 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nextraabsorbers_sw,                 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nmixingfracs_sw,                    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nlayers_sw,                         1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nbnds_sw,                           1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(ngpts_sw,                           1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(npairs_sw,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(ncontributors_lower_sw,             1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(ncontributors_upper_sw,             1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nminor_absorber_intervals_lower_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nminor_absorber_intervals_upper_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
#endif

    ! On master processor, allocate space, read in fields, broadcast to all processors
    if (mpirank .eq. mpiroot) then
       ! Allocate space for arrays
       allocate(gas_names_sw(nabsorbers_sw))
       allocate(scaling_gas_lower_sw(nminor_absorber_intervals_lower_sw))
       allocate(scaling_gas_upper_sw(nminor_absorber_intervals_upper_sw))
       allocate(gas_minor_sw(nminorabsorbers_sw))
       allocate(identifier_minor_sw(nminorabsorbers_sw))
       allocate(minor_gases_lower_sw(nminor_absorber_intervals_lower_sw))
       allocate(minor_gases_upper_sw(nminor_absorber_intervals_upper_sw))
       allocate(minor_limits_gpt_lower_sw(npairs_sw,nminor_absorber_intervals_lower_sw))
       allocate(minor_limits_gpt_upper_sw(npairs_sw,nminor_absorber_intervals_upper_sw))
       allocate(band2gpt_sw(2,nbnds_sw))
       allocate(key_species_sw(2,nlayers_sw,nbnds_sw))
       allocate(band_lims_sw(2,nbnds_sw))
       allocate(press_ref_sw(npress_sw))
       allocate(temp_ref_sw(ntemps_sw))
       allocate(vmr_ref_sw(nlayers_sw, nextrabsorbers_sw, ntemps_sw))
       allocate(kminor_lower_sw(ncontributors_lower_sw, nmixingfracs_sw, ntemps_sw))
       allocate(kmajor_sw(ngpts_sw, nmixingfracs_sw,  npress_sw+1, ntemps_sw))
       allocate(kminor_start_lower_sw(nminor_absorber_intervals_lower_sw))
       allocate(kminor_upper_sw(ncontributors_upper_sw, nmixingfracs_sw, ntemps_sw))
       allocate(kminor_start_upper_sw(nminor_absorber_intervals_upper_sw))
       allocate(minor_scales_with_density_lower_sw(nminor_absorber_intervals_lower_sw))
       allocate(minor_scales_with_density_upper_sw(nminor_absorber_intervals_upper_sw))
       allocate(scale_by_complement_lower_sw(nminor_absorber_intervals_lower_sw))
       allocate(scale_by_complement_upper_sw(nminor_absorber_intervals_upper_sw))
       allocate(rayl_upper_sw(ngpts_sw, nmixingfracs_sw, ntemps_sw))
       allocate(rayl_lower_sw(ngpts_sw, nmixingfracs_sw, ntemps_sw))
       allocate(solar_source_sw(ngpts_sw))
       allocate(temp1(nminor_absorber_intervals_lower_sw))
       allocate(temp2(nminor_absorber_intervals_upper_sw))
       allocate(temp3(nminor_absorber_intervals_lower_sw))
       allocate(temp4(nminor_absorber_intervals_upper_sw))

       ! Read in fields from file
       if(nf90_open(trim(kdist_file), NF90_WRITE, ncid_sw) .eq. NF90_NOERR) then
          status = nf90_inq_varid(ncid_sw,'gas_names',varID)
          status = nf90_get_var(ncid_sw,varID,gas_names_sw)
          !
          status = nf90_inq_varid(ncid_sw,'scaling_gas_lower',varID)
          status = nf90_get_var(ncid_sw,varID,scaling_gas_lower_sw)
          !
          status = nf90_inq_varid(ncid_sw,'scaling_gas_upper',varID)
          status = nf90_get_var(ncid_sw,varID,scaling_gas_upper_sw)
          !
          status = nf90_inq_varid(ncid_sw,'gas_minor',varID)
          status = nf90_get_var(ncid_sw,varID,gas_minor_sw)
          !
          status = nf90_inq_varid(ncid_sw,'identifier_minor',varID)
          status = nf90_get_var(ncid_sw,varID,identifier_minor_sw)
          !
          status = nf90_inq_varid(ncid_sw,'minor_gases_lower',varID)
          status = nf90_get_var(ncid_sw,varID,minor_gases_lower_sw)
          !
          status = nf90_inq_varid(ncid_sw,'minor_gases_upper',varID)
          status = nf90_get_var(ncid_sw,varID,minor_gases_upper_sw)
          !
          status = nf90_inq_varid(ncid_sw,'minor_limits_gpt_lower',varID)
          status = nf90_get_var(ncid_sw,varID,minor_limits_gpt_lower_sw)
          !
          status = nf90_inq_varid(ncid_sw,'minor_limits_gpt_upper',varID)
          status = nf90_get_var(ncid_sw,varID,minor_limits_gpt_upper_sw)
          !
          status = nf90_inq_varid(ncid_sw,'bnd_limits_gpt',varID)
          status = nf90_get_var(ncid_sw,varID,band2gpt_sw)
          !
          status = nf90_inq_varid(ncid_sw,'key_species',varID)
          status = nf90_get_var(ncid_sw,varID,key_species_sw)
          !
          status = nf90_inq_varid(ncid_sw,'bnd_limits_wavenumber',varID)
          status = nf90_get_var(ncid_sw,varID,band_lims_sw)
          !
          status = nf90_inq_varid(ncid_sw,'press_ref',varID)
          status = nf90_get_var(ncid_sw,varID,press_ref_sw)
          !
          status = nf90_inq_varid(ncid_sw,'temp_ref',varID)
          status = nf90_get_var(ncid_sw,varID,temp_ref_sw)
          !
          status = nf90_inq_varid(ncid_sw,'absorption_coefficient_ref_P',varID)
          status = nf90_get_var(ncid_sw,varID,temp_ref_p_sw)
          !
          status = nf90_inq_varid(ncid_sw,'absorption_coefficient_ref_T',varID)
          status = nf90_get_var(ncid_sw,varID,temp_ref_t_sw)
          !
          status = nf90_inq_varid(ncid_sw,'press_ref_trop',varID)
          status = nf90_get_var(ncid_sw,varID,press_ref_trop_sw)
          !
          status = nf90_inq_varid(ncid_sw,'kminor_lower',varID)
          status = nf90_get_var(ncid_sw,varID,kminor_lower_sw)
          !
          status = nf90_inq_varid(ncid_sw,'kminor_upper',varID)
          status = nf90_get_var(ncid_sw,varID,kminor_upper_sw)
          !
          status = nf90_inq_varid(ncid_sw,'vmr_ref',varID)
          status = nf90_get_var(ncid_sw,varID,vmr_ref_sw)
          !
          status = nf90_inq_varid(ncid_sw,'kmajor',varID)
          status = nf90_get_var(ncid_sw,varID,kmajor_sw)
          !
          status = nf90_inq_varid(ncid_sw,'kminor_start_lower',varID)
          status = nf90_get_var(ncid_sw,varID,kminor_start_lower_sw)
          !
          status = nf90_inq_varid(ncid_sw,'kminor_start_upper',varID)
          status = nf90_get_var(ncid_sw,varID,kminor_start_upper_sw)
          !
          status = nf90_inq_varid(ncid_sw,'solar_source',varID)
          status = nf90_get_var(ncid_sw,varID,solar_source_sw)
          !
          status = nf90_inq_varid(ncid_sw,'rayl_lower',varID)
          status = nf90_get_var(ncid_sw,varID,rayl_lower_sw)

          status = nf90_inq_varid(ncid_sw,'rayl_upper',varID)
          status = nf90_get_var(ncid_sw,varID,rayl_upper_sw)

          ! Logical fields are read in as integers and then converted to logicals.
          status = nf90_inq_varid(ncid_sw,'minor_scales_with_density_lower',varID)
          status = nf90_get_var(ncid_sw,varID,temp1)
          minor_scales_with_density_lower_sw(:) = .false.
          where(temp1 .eq. 1) minor_scales_with_density_lower_sw(:) = .true.
          !
          status = nf90_inq_varid(ncid_sw,'minor_scales_with_density_upper',varID)
          status = nf90_get_var(ncid_sw,varID,temp2)
          minor_scales_with_density_upper_sw(:) = .false.
          where(temp2 .eq. 1) minor_scales_with_density_upper_sw(:) = .true.
          !
          status = nf90_inq_varid(ncid_sw,'scale_by_complement_lower',varID)
          status = nf90_get_var(ncid_sw,varID,temp3)
          scale_by_complement_lower_sw(:) = .false.
          where(temp3 .eq. 1) scale_by_complement_lower_sw(:) = .true.
          !
          status = nf90_inq_varid(ncid_sw,'scale_by_complement_upper',varID)
          status = nf90_get_var(ncid_sw,varID,temp4)
          scale_by_complement_upper_sw(:) = .false.
          where(temp4 .eq. 1) scale_by_complement_upper_sw(:) = .true.
          
          ! Close
          status = nf90_close(ncid_sw)
       endif
    endif

    ! Broadcast arrays to all processors
#ifdef MPI
    call MPI_BCAST(minor_limits_gpt_upper_sw, size(minor_limits_gpt_upper_sw), MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(minor_limits_gpt_lower_sw, size(minor_limits_gpt_lower_sw), MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(kminor_start_upper_sw,     size(kminor_start_upper_sw),     MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(kminor_start_lower_sw,     size(kminor_start_lower_sw),     MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(key_species_sw,            size(key_species_sw),            MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(band2gpt_sw,               size(band2gpt_sw),               MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(band_lims_sw,              size(band_lims_sw),              kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(press_ref_sw,              size(press_ref_sw),              kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(temp_ref_sw,               size(temp_ref_sw),               kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(kminor_lower_sw,           size(kminor_lower_sw),           kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(kminor_upper_sw,           size(kminor_upper_sw),           kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(scaling_gas_lower_sw,      size(scaling_gas_lower_sw),      kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(scaling_gas_upper_sw,      size(scaling_gas_upper_sw),      kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(vmr_ref_sw,                size(vmr_ref_sw),                kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(kmajor_sw,                 size(kmajor_sw),                 kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(temp_ref_p_sw,             1,                               kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(temp_ref_t_sw,             1,                               kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(press_ref_trop_sw,         1,                               kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(solar_source_sw,           size(solar_source_sw),           kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(rayl_lower_sw,             size(rayl_lower_sw),             kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(rayl_upper_sw,             size(rayl_upper_sw),             kind_phys,   mpiroot, mpicomm, ierr)
    ! Character arrays
    do ij=1,nabsorbers_sw
       call MPI_BCAST(gas_names_sw(ij),         32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
    enddo
    do ij=1,nminorabsorbers_sw
       call MPI_BCAST(gas_minor_sw(ij),         32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(identifier_minor_sw(ij),  32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
    enddo
    do ij=1,nminor_absorber_intervals_lower_sw
       call MPI_BCAST(minor_gases_lower_sw(ij), 32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
    enddo
    do ij=1,nminor_absorber_intervals_upper_sw
       call MPI_BCAST(minor_gases_upper_sw(ij), 32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
    enddo
    ! Logical arrays (First convert to integer-array, then broadcast)
    !
    allocate(temp_log_array1(nminor_absorber_intervals_lower_sw))
    where(minor_scales_with_density_lower_sw)
       temp_log_array1 = 1
    elsewhere
       temp_log_array1 = 0
    end where
    call MPI_BCAST(temp_log_array1, size(temp_log_array1), MPI_INTEGER,  mpiroot, mpicomm, ierr)
    !
    allocate(temp_log_array2(nminor_absorber_intervals_lower_sw))
    where(scale_by_complement_lower_sw)
       temp_log_array2 = 1
    elsewhere
       temp_log_array2 = 0
    end where
    call MPI_BCAST(temp_log_array2, size(temp_log_array2), MPI_INTEGER,  mpiroot, mpicomm, ierr)
    !
    allocate(temp_log_array3(nminor_absorber_intervals_upper_sw))
    where(minor_scales_with_density_upper_sw)
       temp_log_array3 = 1
    elsewhere
       temp_log_array3 = 0
    end where
    call MPI_BCAST(temp_log_array3, size(temp_log_array3), MPI_INTEGER,  mpiroot, mpicomm, ierr)
    !
    allocate(temp_log_array4(nminor_absorber_intervals_upper_sw))
    where(scale_by_complement_upper_sw)
       temp_log_array4 = 1
    elsewhere
       temp_log_array4 = 0
    end where
    call MPI_BCAST(temp_log_array4, size(temp_log_array4), MPI_INTEGER,  mpiroot, mpicomm, ierr)
#endif

    ! Initialize gas concentrations and gas optics class with data
    do iGas=1,nGases
       call check_error_msg(gas_concs_sw%set_vmr(active_gases(iGas), 0._kind_phys))
    enddo    
    call check_error_msg(kdist_sw_clr%load(gas_concs_sw, gas_names_sw, key_species_sw, band2gpt_sw, &
         band_lims_sw, press_ref_sw, press_ref_trop_sw, temp_ref_sw,  temp_ref_p_sw, temp_ref_t_sw, &
         vmr_ref_sw, kmajor_sw, kminor_lower_sw, kminor_upper_sw, gas_minor_sw,identifier_minor_sw, &
         minor_gases_lower_sw, minor_gases_upper_sw, minor_limits_gpt_lower_sw,                     &
         minor_limits_gpt_upper_sw, minor_scales_with_density_lower_sw,                             &
         minor_scales_with_density_upper_sw, scaling_gas_lower_sw,                                  &
         scaling_gas_upper_sw, scale_by_complement_lower_sw,                                        &
         scale_by_complement_upper_sw, kminor_start_lower_sw, kminor_start_upper_sw,                &
         solar_source_sw, rayl_lower_sw, rayl_upper_sw))

    ! Set band index by g-point array 
    nBandsSW = kdist_sw_clr%get_nband()
    nGptsSW  = kdist_sw_clr%get_ngpt()

    ! Set initial permutation seed for McICA, initially set to number of G-points
    ipsdsw0 = kdist_sw_clr%get_ngpt()

  end subroutine rrtmgp_sw_init
  ! #########################################################################################
  ! RRTMGP_SW_RUN
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_run Argument Table
!! | local_name      | standard_name                                                                                  | long_name                                                                | units   | rank | type        |    kind   | intent | optional |
!! |-----------------|------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|---------|------|-------------|-----------|--------|----------|
!! | p_lay           | air_pressure_at_layer_for_radiation_in_hPa                                                     | air pressure layer                                                       | hPa     |    2 | real        | kind_phys | in     | F        |
!! | p_lev           | air_pressure_at_interface_for_radiation_in_hPa                                                 | air pressure level                                                       | hPa     |    2 | real        | kind_phys | in     | F        |
!! | t_lay           | air_temperature_at_layer_for_radiation                                                         | air temperature layer                                                    | K       |    2 | real        | kind_phys | in     | F        |
!! | t_lev           | air_temperature_at_interface_for_radiation                                                     | air temperature level                                                    | K       |    2 | real        | kind_phys | in     | F        |
!! | q_lay           | water_vapor_specific_humidity_at_layer_for_radiation                                           | specific humidity layer                                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | o3_lay          | ozone_concentration_at_layer_for_radiation                                                     | ozone concentration layer                                                | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_co2         | volume_mixing_ratio_co2                                                                        | volume mixing ratio co2                                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_n2o         | volume_mixing_ratio_n2o                                                                        | volume mixing ratio no2                                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_ch4         | volume_mixing_ratio_ch4                                                                        | volume mixing ratio ch4                                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_o2          | volume_mixing_ratio_o2                                                                         | volume mixing ratio o2                                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_co          | volume_mixing_ratio_co                                                                         | volume mixing ratio co                                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_cfc11       | volume_mixing_ratio_cfc11                                                                      | volume mixing ratio cfc11                                                | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_cfc12       | volume_mixing_ratio_cfc12                                                                      | volume mixing ratio cfc12                                                | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_cfc22       | volume_mixing_ratio_cfc22                                                                      | volume mixing ratio cfc22                                                | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_ccl4        | volume_mixing_ratio_ccl4                                                                       | volume mixing ratio ccl4                                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | icseed          | seed_random_numbers_sw                                                                         | seed for random number generation for shortwave radiation                | none    |    1 | integer     |           | in     | F        |
!! | tau_aer         | aerosol_optical_depth_for_longwave_bands_01-16                                                 | aerosol optical depth for shortwave bands 01-16                          | none    |    3 | real        | kind_phys | in     | F        |
!! | ssa_aer         | aerosol_single_scattering_albedo_for_longwave_bands_01-16                                      | aerosol single scattering albedo for shortwave bands 01-16               | frac    |    3 | real        | kind_phys | in     | F        |
!! | asy_aer         | aerosol_asymmetry_parameter_for_shortwave_bands_01-16                                          | aerosol asymmetry paramter for shortwave bands 01-16                     | none    |    3 | real        | kind_phys | in     | F        |
!! | sfcalb_nir_dir  | surface_albedo_due_to_near_IR_direct                                                           | surface albedo due to near IR direct beam                                | frac    |    1 | real        | kind_phys | in     | F        |
!! | sfcalb_nir_dif  | surface_albedo_due_to_near_IR_diffused                                                         | surface albedo due to near IR diffused beam                              | frac    |    1 | real        | kind_phys | in     | F        |
!! | sfcalb_uvis_dir | surface_albedo_due_to_UV_and_VIS_direct                                                        | surface albedo due to UV+VIS direct beam                                 | frac    |    1 | real        | kind_phys | in     | F        |
!! | sfcalb_uvis_dif | surface_albedo_due_to_UV_and_VIS_diffused                                                      | surface albedo due to UV+VIS diffused beam                               | frac    |    1 | real        | kind_phys | in     | F        |
!! | dzlyr           | layer_thickness_for_radiation                                                                  | layer thickness                                                          | km      |    2 | real        | kind_phys | in     | F        |
!! | delpin          | layer_pressure_thickness_for_radiation                                                         | layer pressure thickness                                                 | hPa     |    2 | real        | kind_phys | in     | F        | 
!! | de_lgth         | cloud_decorrelation_length                                                                     | cloud decorrelation length                                               | km      |    1 | real        | kind_phys | in     | F        | 
!! | cossza          | cosine_of_zenith_angle                                                                         | cosine of the solar zenit angle                                          | none    |    1 | real        | kind_phys | in     | F        |
!! | solcon          | solar_constant                                                                                 | solar constant                                                           | W m-2   |    0 | real        | kind_phys | in     | F        |
!! | nday            | daytime_points_dimension                                                                       | daytime points dimension                                                 | count   |    0 | integer     |           | in     | F        |
!! | idxday          | daytime_points                                                                                 | daytime points                                                           | index   |    1 | integer     |           | in     | F        |
!! | ncol            | horizontal_loop_extent                                                                         | horizontal dimension                                                     | count   |    0 | integer     |           | in     | F        |
!! | nlay            | adjusted_vertical_layer_dimension_for_radiation                                                | number of vertical layers for radiation                                  | count   |    0 | integer     |           | in     | F        |
!! | lprint          | flag_print                                                                                     | flag to print                                                            | flag    |    0 | logical     |           | in     | F        |
!! | cldfrac         | total_cloud_fraction                                                                           | total cloud fraction                                                     | frac    |    2 | real        | kind_phys | in     | F        |
!! | lsswr           | flag_to_calc_sw                                                                                | flag to calculate SW irradiances                                         | flag    |    0 | logical     |           | in     | F        |
!! | hswc            | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step                    | shortwave total sky heating rate                                         | K s-1   |    2 | real        | kind_phys | inout  | F        |
!! | topflx          | sw_fluxes_top_atmosphere                                                                       | shortwave total sky fluxes at the top of the atm                         | W m-2   |    1 | topfsw_type |           | inout  | F        |
!! | sfcflx          | sw_fluxes_sfc                                                                                  | shortwave total sky fluxes at the Earth surface                          | W m-2   |    1 | sfcfsw_type |           | inout  | F        |
!! | cldtau          | cloud_optical_depth_layers_at_0.55mu_band                                                      | approx .55mu band layer cloud optical depth                              | none    |    2 | real        | kind_phys | inout  | F        |
!! | hsw0            | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | shortwave clear sky heating rate                                         | K s-1   |    2 | real        | kind_phys | inout  | T        |
!! | hswb            | sw_heating_rate_spectral                                                                       | shortwave total sky heating rate (spectral)                              | K s-1   |    3 | real        | kind_phys | inout  | T        |
!! | flxprf          | sw_fluxes                                                                                      | sw fluxes total sky / csk and up / down at levels                        | W m-2   |    2 | profsw_type |           | inout  | T        |
!! | fdncmp          | components_of_surface_downward_shortwave_fluxes                                                | derived type for special components of surface downward shortwave fluxes | W m-2   |    1 | cmpfsw_type |           | inout  | T        |
!! | cld_lwp         | cloud_liquid_water_path                                                                        | cloud liquid water path                                                  | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_liq     | mean_effective_radius_for_liquid_cloud                                                         | mean effective radius for liquid cloud                                   | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_iwp         | cloud_ice_water_path                                                                           | cloud ice water path                                                     | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_ice     | mean_effective_radius_for_ice_cloud                                                            | mean effective radius for ice cloud                                      | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_rwp         | cloud_rain_water_path                                                                          | cloud rain water path                                                    | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_rain    | mean_effective_radius_for_rain_drop                                                            | mean effective radius for rain drop                                      | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_swp         | cloud_snow_water_path                                                                          | cloud snow water path                                                    | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_snow    | mean_effective_radius_for_snow_flake                                                           | mean effective radius for snow flake                                     | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_od          | cloud_optical_depth                                                                            | cloud optical depth                                                      | none    |    2 | real        | kind_phys | in     | T        |
!! | cld_ssa         | cloud_single_scattering_albedo                                                                 | cloud single scattering albedo                                           | frac    |    2 | real        | kind_phys | in     | T        |
!! | cld_asy         | cloud_asymmetry_parameter                                                                      | cloud asymmetry parameter                                                | none    |    2 | real        | kind_phys | in     | T        |
!! | errmsg          | ccpp_error_message                                                                             | error message for error handling in CCPP                                 | none    |    0 | character   | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                                                                | error flag for error handling in CCPP                                    | flag    |    0 | integer     |           | out    | F        |
!!
  subroutine rrtmgp_sw_run(p_lay, p_lev, t_lay, t_lev, q_lay, o3_lay, vmr_co2, vmr_n2o,     & ! IN
       vmr_ch4, vmr_o2, vmr_co, vmr_cfc11, vmr_cfc12, vmr_cfc22, vmr_ccl4, icseed, tau_aer, & ! IN
       ssa_aer, asy_aer, sfcalb_nir_dir, sfcalb_nir_dif,  sfcalb_uvis_dir, sfcalb_uvis_dif, & ! IN
       dzlyr, delpin, de_lgth, cossza, solcon, nday, idxday, ncol, nlay, lprint, cldfrac,   & ! IN
       lsswr,                                                                               & ! IN
       hswc, topflx, sfcflx, cldtau,                                                        & ! OUT
       hsw0, hswB, flxprf, fdncmp, cld_lwp, cld_ref_liq, cld_iwp, cld_ref_ice, cld_rwp,     & ! OUT(optional)
       cld_ref_rain, cld_swp, cld_ref_snow, cld_od, cld_ssa, cld_asy,                       & ! OUT(optional)
       errmsg, errflg)

    ! Inputs
    integer, intent(in) :: &
         ncol,            & ! Number of horizontal grid-points
         nlay,            & ! Number of vertical layers
         nday               ! Number of daytime points
    integer, intent(in), dimension(ncol) :: &
         icseed             ! auxiliary special cloud related array when module 
                            ! variable isubcsw=2, it provides permutation seed 
                            ! for each column profile that are used for generating 
                            ! random numbers. when isubcsw /=2, it will not be used.
    integer, intent(in), dimension(nday) :: &
         idxday             ! Index array for daytime points
    logical, intent(in) :: &
         lprint,          & ! Control flag for diagnostics
         lsswr              ! Flag to calculate RRTMGP SW? 
    real(kind_phys), intent(in) :: &
         solcon             ! Solar constant                          (W/m2)
    real(kind_phys), dimension(ncol), intent(in) :: &
         sfcalb_nir_dir,  & ! Surface albedo direct (near-IR)         (1)
         sfcalb_nir_dif,  & ! Surface albedo diffuse (near-IR)        (1)
         sfcalb_uvis_dir, & ! Surface albedo direct (UV and Visible)  (1)
         sfcalb_uvis_dif, & ! Surface albedo diffuse (UV and Visible) (1)
         de_lgth,         & ! Cloud decorrelation length              (km)
         cossza             ! Cosine of solar zenith angle            (1)
    real(kind_phys), dimension(ncol,nlay), intent(in) :: &
         dzlyr,           & ! layer thinkness                         (km)
         delpin,          & ! layer thickness                         (mb)
         cldfrac,         & ! Cloud-fraction                          (1)
         p_lay,           & ! Pressure @ model layer-centers          (mb)
         t_lay,           & ! Temperature                             (K)
         q_lay,           & ! Specific humidity                       (kg/kg)
         o3_lay,          & ! O3 mass mixing-ratio                    (kg/kg)
         vmr_co2,         & ! Co2 volume-mixing ratio                 (kg/kg)
         vmr_n2o,         & ! N2o volume-mixing ratio                 (kg/kg)
         vmr_ch4,         & ! Ch4 volume-mixing ratio                 (kg/kg)
         vmr_o2,          & ! O2 volume-mixing ratio                  (kg/kg)
         vmr_co,          & ! Co volume-mixing ratio                  (kg/kg)
         vmr_cfc11,       & ! CFC11 volume-mixing ratio               (kg/kg)
         vmr_cfc12,       & ! CFC12 volume-mixing ratio               (kg/kg)
         vmr_cfc22,       & ! CFC22 volume-mixing ratio               (kg/kg)
         vmr_ccl4           ! CCl4 volume-mixing ratio                (kg/kg)
    real(kind_phys), dimension(ncol,nlay+1), intent(in) :: &
         p_lev,           & ! Pressure @ model layer-interfaces       (mb)
         t_lev              ! Temperature                             (K)
    real(kind_phys), dimension(ncol,nlay,nbandsSW), intent(in) :: &
         tau_aer,         & ! Aerosol optical depth                   (1) 
         ssa_aer,         & ! Aerosol single-scattering albedo        (1)
         asy_aer            ! Aerosol asymmetry parameter             (1)
    ! Inputs (optional)
    real(kind_phys), dimension(ncol,nlay), intent(in), optional:: &
         cld_lwp,      & ! Cloud liquid water path                (g/m2)
         cld_ref_liq,  & ! Effective radius (liquid)              (micron)
         cld_iwp,      & ! Cloud ice water path                   (g/m2)
         cld_ref_ice,  & ! Effective radius (ice)                 (micron)
         cld_rwp,      & ! Cloud rain water path                  (g/m2)
         cld_ref_rain, & ! Effective radius (rain-drop)           (micron)
         cld_swp,      & ! Cloud snow-water path                  (g/m2)
         cld_ref_snow, & ! Effective radius (snow-flake)          (micron) 
         cld_od,       & ! Cloud optical-depth                    (1)
         cld_ssa,      & ! Cloud single-scattering albedo         (1)
         cld_asy         ! Cloud asymmetry parameter              (1)

    ! Outputs (mandatory)
    character(len=*), intent(out) :: &
         errmsg             ! Error message
    integer, intent(out) :: &
         errflg             ! Error code
    real(kind_phys),dimension(ncol,nlay), intent(inout) :: &
         hswc,            & ! All-sky heating-rate                    (K/sec)
         cldtau             ! ~0.55mu band layer tau                  (1)
    type(topfsw_type), dimension(ncol), intent(inout) :: &
         topflx             ! radiation fluxes at top, components:
                            ! upfxc - total sky upward flux at top    (w/m2)
                            ! upfx0 - clear sky upward flux at top    (w/m2)
    type(sfcfsw_type), dimension(ncol), intent(inout) :: &
         sfcflx             ! radiation fluxes at sfc, components:
                            ! upfxc - total sky upward flux at sfc    (w/m2)  
                            ! upfx0 - clear sky upward flux at sfc    (w/m2)
                            ! dnfxc - total sky downward flux at sfc  (w/m2)
                            ! dnfx0 - clear sky downward flux at sfc  (w/m2)
    ! Outputs (optional)
    real(kind_phys), dimension(ncol,nlay), intent(inout), optional :: &
         hsw0               ! Clear-sky heating-rate                  (K/sec)
    real(kind_phys), dimension(ncol,nlay,nBandsSW), intent(inout), optional :: &
         hswb               ! All-sky heating rate, in each band      (K/sec)
    type(profsw_type), dimension(ncol,nlay+1), intent(inout), optional :: &
         flxprf             ! 2D radiative fluxes, components:
                            ! upfxc - total sky upward flux           (W/m2)
                            ! dnfxc - total sky dnward flux           (W/m2)
                            ! upfx0 - clear sky upward flux           (W/m2)
                            ! dnfx0 - clear sky dnward flux           (W/m2)
    type(cmpfsw_type), dimension(ncol), intent(inout), optional :: &
         fdncmp             ! 2D surface fluxes, components:
                            ! uvbfc - total sky downward uv-b flux at (W/m2)
                            ! uvbf0 - clear sky downward uv-b flux at (W/m2)
                            ! nirbm - downward nir direct beam flux   (W/m2)
                            ! nirdf - downward nir diffused flux      (W/m2)
                            ! visbm - downward uv+vis direct beam flux(W/m2)
                            ! visdf - downward uv+vis diffused flux   (W/m2)

    ! RTE+RRTMGP classes
    type(ty_optical_props_2str) :: &
         optical_props_clr  ! Optical properties for gaseous atmosphere
    type(ty_optical_props_2str) :: &
         optical_props_aer  ! Optical properties for aerosols
    type(ty_fluxes_broadband) :: &
         fluxAllSky,         & ! All-sky flux                      (W/m2)
         fluxClrSky            ! Clear-sky flux                    (W/m2)

    ! Local variables
    integer :: iCol, iBand, iGpt, iDay
    integer,dimension(ncol) :: ipseed
    real(kind_phys) :: solAdjFac
    logical :: top_at_1=.false.
    real(kind_phys), dimension(ncol,nlay) :: vmr_o3, vmr_h2o, thetaTendClrSky, &
         thetaTendAllSky,coldry,tem0
    real(kind_phys), dimension(nday,nlay+1),target :: &
         flux_up_allSky, flux_up_clrSky, flux_dn_allSky, flux_dn_clrSky, p_lev2
    real(kind_phys), dimension(nday,nGptsSW) :: toa_flux

    ! Initialize
    errmsg = ''
    errflg = 0
    if (.not. lsswr) return
    if (nday <= 0) return

    ! Are any optional outputs requested?
    l_ClrSky_HR        = present(hsw0)
    l_AllSky_HR_byband = present(hswb)
    l_fluxes2D         = present(flxprf)
    l_sfcFluxes1D      = present(fdncmp)

    ! Check for optional input arguments, this depends on cloud method
    if (iswcliq > 0) then    ! use prognostic cloud method
       if ( .not.present(cld_lwp) .or. .not.present(cld_ref_liq) .or.  &
            .not.present(cld_iwp) .or. .not.present(cld_ref_ice) .or.  &
            .not.present(cld_rwp) .or. .not.present(cld_ref_rain) .or. &
            .not.present(cld_swp) .or. .not.present(cld_ref_snow) )then
          write(errmsg,'(*(a))')                                        &
               'Logic error: iswcliq>0 requires the following',   &
               ' optional arguments to be present:',              &
               ' cld_lwp, cld_ref_liq, cld_iwp, cld_ref_ice,',    &
               ' cld_rwp, cld_ref_rain, cld_swp, cld_ref_snow'
          errflg = 1
          return
       end if
    else                     ! use diagnostic cloud method
       if ( .not.present(cld_od) .or. .not.present(cld_ssa) .or.       &
            .not.present(cld_asy)) then
          write(errmsg,'(*(a))')                                        &
               'Logic error: iswcliq<=0 requires the following',  &
               ' optional arguments to be present:',              &
               ' cld_od, cld_ssa, cld_asy'
          errflg = 1
          return
       end if
    endif

    ! What is vertical ordering?
    top_at_1 = (p_lay(1,1) .lt. p_lay(1,nlay))

    ! Change random number seed value for each radiation invocation (isubcsw =1 or 2).
    if(isubcsw == 1) then      ! advance prescribed permutation seed
       do iCol = 1, ncol
          ipseed(iCol) = ipsdsw0 + iCol
       enddo
    elseif (isubcsw == 2) then ! use input array of permutaion seeds
       do iCol = 1, ncol
          ipseed(iCol) = icseed(iCol)
       enddo
    endif

    ! Compute volume mixing-ratios for ozone (mmr) and specific-humidity.
    vmr_h2o = merge((q_lay/(1-q_lay))*amdw, 0., q_lay  .ne. 1.)
    vmr_o3  = merge(o3_lay*amdo3,           0., o3_lay .gt. 0.)

    ! Compute dry air column amount
    tem0   = (1. - vmr_h2o)*amd + vmr_h2o*amw
    coldry = ( 1.0e-20 * 1.0e3 *avogad)*delpin / (100.*grav*tem0*(1. + vmr_h2o))

    ! Input model-level pressure @ the top-of-model is set to 1Pa, whereas RRTMGP minimum 
    ! pressure needs to be slightly greater than that, ~1.00518Pa
    p_lev2=p_lev
    p_lev2(:,nlay+1) = kdist_sw_clr%get_press_min()/100.

    ! Compute solar constant adjustment factor..
    solAdjFac = solcon / s0

    ! Initialize outputs
    hswc(:,:)   = 0.
    cldtau(:,:) = 0.
    topflx      = topfsw_type ( 0., 0., 0. )
    sfcflx      = sfcfsw_type ( 0., 0., 0., 0. )
    if (l_ClrSky_HR) then
       hsw0(:,:) = 0.
    endif
    if(l_AllSky_HR_byband) then
       hswb(:,:,:) = 0.
    endif
    if (l_fluxes2D) then
       flxprf = profsw_type ( 0., 0., 0., 0. )
    endif
    if (l_sfcFluxes1D) then
       fdncmp = cmpfsw_type (0.,0.,0.,0.,0.,0.)
    endif

    ! #######################################################################################
    ! CALL RRTMGP (Only for daylit (idxday) points)
    ! #######################################################################################
    if (nDay .gt. 0) then

       ! Allocate space for source functions and gas optical properties
       call check_error_msg(optical_props_clr%alloc_2str(nday, nlay, kdist_sw_clr))
       call check_error_msg(optical_props_aer%alloc_2str(nday, nlay, kdist_sw_clr))

       ! Initialize RRTMGP files
       fluxAllSky%flux_up => flux_up_allSky
       fluxAllsky%flux_dn => flux_dn_allSky
       fluxClrSky%flux_up => flux_up_clrSky
       fluxClrsky%flux_dn => flux_dn_clrSky
       
       ! #######################################################################################
       ! 1) Clear-sky fluxes (gaseous-atmosphere + aerosols)
       ! #######################################################################################
       ! 1a) Set gas concentrations
       print*,'Clear-Sky(SW): Set Gas Concentrations'
       call gas_concs_sw%reset()
       call check_error_msg(gas_concs_sw%set_vmr('o2',  vmr_o2(idxday,1:nlay)))
       call check_error_msg(gas_concs_sw%set_vmr('co2', vmr_co2(idxday,1:nlay)))
       call check_error_msg(gas_concs_sw%set_vmr('ch4', vmr_ch4(idxday,1:nlay)))
       call check_error_msg(gas_concs_sw%set_vmr('n2o', vmr_n2o(idxday,1:nlay)))
       call check_error_msg(gas_concs_sw%set_vmr('h2o', vmr_h2o(idxday,1:nlay)))
       call check_error_msg(gas_concs_sw%set_vmr('o3',  vmr_o3(idxday,1:nlay)))
       
       ! 1b) Compute the optical properties of the atmosphere and the Planck source functions
       !    from pressures, temperatures, and gas concentrations...
       print*,'Clear-Sky(SW): Optics'
       call check_error_msg(kdist_sw_clr%gas_optics(      &
            p_lay(idxday,1:nlay)*100.,      &
            p_lev2(idxday,1:nlay+1)*100.,   &
            t_lay(idxday,1:nlay),           &
            gas_concs_sw,                   &
            optical_props_clr,              &
            toa_flux))

       ! 1c) Add contribution from aerosols. Also, apply diffusivity angle
       print*,'Clear-Sky(SW): Increment Aerosol'
       do iDay=1,nDay
          do iGpt=1,nGptsSW
             iBand = kdist_sw_clr%convert_gpt2band(iGpt)
             optical_props_aer%tau(iDay,1:nlay,iGpt) = tau_aer(idxday(iDay),1:nlay,iBand)
             optical_props_aer%ssa(iDay,1:nlay,iGpt) = ssa_aer(idxday(iDay),1:nlay,iBand)
             optical_props_aer%g(iDay,1:nlay,iGpt)   = asy_aer(idxday(iDay),1:nlay,iBand)
          enddo
       enddo
       call check_error_msg(optical_props_aer%increment(optical_props_clr))

       ! 1c) Compute the clear-sky broadband fluxes
       print*,'Clear-Sky(SW): Fluxes'
       call check_error_msg(rte_sw(optical_props_clr, top_at_1, cossza(idxday), toa_flux,&
            spread(sfcalb_nir_dir(idxday),1, ncopies = nBandsSW), & 
            spread(sfcalb_nir_dif(idxday),1, ncopies = nBandsSW), &
            fluxClrSky))

       ! 1d) Compute heating rates
       print*,'Clear-Sky(SW): Heating-rates'
       call check_error_msg(compute_heating_rate(   &
            fluxClrSky%flux_up,                     &
            fluxClrSky%flux_dn,                     &
            p_lev2(idxday,1:nlay+1)*100., &
            thetaTendClrSky))
       thetaTendAllSky = thetaTendClrSky
   
    end if ! Daylit days

    ! #######################################################################################
    ! Copy fluxes from RRTGMP types into model radiation types.
    ! #######################################################################################
    ! Mandatory outputs
    topflx(idxday)%upfxc = fluxAllSky%flux_up(:,nlay+1)
    topflx(idxday)%upfx0 = fluxClrSky%flux_up(:,nlay+1)
    sfcflx(idxday)%upfxc = fluxAllSky%flux_up(:,1)
    sfcflx(idxday)%upfx0 = fluxClrSky%flux_up(:,1)
    sfcflx(idxday)%dnfxc = fluxAllSky%flux_dn(:,1)
    sfcflx(idxday)%dnfx0 = fluxClrSky%flux_dn(:,1)
    cldtau(idxday,:)     = 0.
    hswc(idxday,:)       = thetaTendAllSky

    ! Optional output
    if(l_fluxes2D) then
       flxprf(idxday,:)%upfxc = fluxAllSky%flux_up
       flxprf(idxday,:)%dnfxc = fluxAllSky%flux_dn
       flxprf(idxday,:)%upfx0 = fluxClrSky%flux_up
       flxprf(idxday,:)%dnfx0 = fluxClrSky%flux_dn
    endif
    if (l_AllSky_HR_byband) then
       do iBand=1,nBandsSW
          hswb(idxday,:,iBand) = thetaTendAllSky
       end do
    endif
    if (l_ClrSky_HR) then
       hsw0(idxday,:) = thetaTendClrSky
    endif

  end subroutine rrtmgp_sw_run
  ! #########################################################################################
  ! #########################################################################################
  subroutine rrtmgp_sw_finalize()
  end subroutine rrtmgp_sw_finalize

  ! #########################################################################################
  ! Ancillary functions
  ! #########################################################################################
  subroutine check_error_msg(error_msg)
    character(len=*), intent(in) :: error_msg
    
    if(error_msg /= "") then
       print*,"ERROR(rrtmgp_sw_main.F90): "
       print*,trim(error_msg)
       return
    end if
  end subroutine check_error_msg  
  ! #########################################################################################
  ! #########################################################################################
end module rrtmgp_sw
