module rrtmgp_sw_gas_optics
  use machine,               only: kind_phys
  use GFS_typedefs,          only: GFS_control_type, GFS_radtend_type
  use mo_rte_kind,           only: wl
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use netcdf

  ! Parameters
  integer :: ipsdsw0

contains
!! \section arg_table_rrtmgp_sw_gas_optics_init Argument Table
!! | local_name   | standard_name                                | long_name                                                          | units | rank | type                 | kind  | intent | optional |
!! |--------------|----------------------------------------------|--------------------------------------------------------------------|-------|------|----------------------|-------|--------|----------|
!! | Model        | GFS_control_type_instance                    | Fortran DDT containing FV3-GFS model control parameters            | DDT   |    0 | GFS_control_type     |       | in     | F        |
!! | Radtend      | GFS_radtend_type_instance                    | Fortran DDT containing FV3-GFS radiation tendencies                | DDT   |    0 | GFS_radtend_type     |       | in     | F        |
!! | mpirank      | mpi_rank                                     | current MPI rank                                                   | index |    0 | integer              |       | in     | F        |
!! | mpiroot      | mpi_root                                     | master MPI rank                                                    | index |    0 | integer              |       | in     | F        |
!! | mpicomm      | mpi_comm                                     | MPI communicator                                                   | index |    0 | integer              |       | in     | F        |
!! | errmsg       | ccpp_error_message                           | error message for error handling in CCPP                           | none  |    0 | character            | len=* | out    | F        |
!! | errflg       | ccpp_error_flag                              | error flag for error handling in CCPP                              | flag  |    0 | integer              |       | out    | F        |
!! | sw_gas_props | coefficients_for_sw_gas_optics               | DDT containing spectral information for RRTMGP SW radiation scheme | DDT   |    0 | ty_gas_optics_rrtmgp |       | out    | F        |
!! | ngpts_sw     | number_of_spectral_points_for_SW_calculation | Number of spectral points for SW RRTMGP calculation                | none  |    0 | integer              |       | out    | F        |
!!
  ! #########################################################################################
  ! SUBROUTINE sw_gas_optics_init
  ! #########################################################################################
  subroutine rrtmgp_sw_gas_optics_init(Model, Radtend, mpicomm, mpirank, mpiroot, sw_gas_props, &
       ngpts_sw, errmsg, errflg)
    use netcdf
#ifdef MPI
    use mpi
#endif

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model      ! DDT containing model control parameters
    type(GFS_radtend_type), intent(in) :: &
         Radtend
    integer,intent(in) :: &
         mpicomm, & ! MPI communicator
         mpirank, & ! Current MPI rank
         mpiroot    ! Master MPI rank

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer,          intent(out) :: &
         errflg     ! Error code
    integer, intent(out) :: &
         ngpts_sw
    type(ty_gas_optics_rrtmgp),intent(out) :: &
         sw_gas_props

    ! Fields from the K-distribution files
    ! Variables that will be passed to gas_optics%load()
    type(ty_gas_concs)  :: &
         gas_concentrations
    integer, dimension(:), allocatable :: &
         kminor_start_lower_sw,              & ! used by RRTMGP gas optics 
         kminor_start_upper_sw                 ! used by RRTMGP gas optics 
    integer, dimension(:,:), allocatable :: &
         band2gpt_sw,                        & ! used by RRTMGP gas optics 
         minor_limits_gpt_lower_sw,          & ! used by RRTMGP gas optics 
         minor_limits_gpt_upper_sw             ! used by RRTMGP gas optics 
    integer, dimension(:,:,:), allocatable :: &
         key_species_sw                        ! used by RRTMGP gas optics 
    real(kind_phys) :: &
         press_ref_trop_sw,                  & ! used by RRTMGP gas optics 
         temp_ref_p_sw,                      & ! used by RRTMGP gas optics 
         temp_ref_t_sw                         ! used by RRTMGP gas optics 
    real(kind_phys), dimension(:), allocatable :: &
         press_ref_sw,                       & ! used by RRTMGP gas optics 
         temp_ref_sw,                        & ! used by RRTMGP gas optics 
         solar_source_sw                       ! used by RRTMGP gas optics 
    real(kind_phys), dimension(:,:), allocatable :: &
         band_lims_sw                          ! used by RRTMGP gas optics                         

    real(kind_phys), dimension(:,:,:), allocatable :: &
         vmr_ref_sw,                         & ! used by RRTMGP gas optics 
         kminor_lower_sw,                    & ! used by RRTMGP gas optics 
         kminor_upper_sw,                    & ! used by RRTMGP gas optics 
         rayl_lower_sw,                      & ! used by RRTMGP gas optics 
         rayl_upper_sw                         ! used by RRTMGP gas optics 
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         kmajor_sw                             ! used by RRTMGP gas optics 
    character(len=32),  dimension(:), allocatable :: &
         gas_names_sw,                       & ! used by RRTMGP gas optics 
         gas_minor_sw,                       & ! used by RRTMGP gas optics 
         identifier_minor_sw,                & ! used by RRTMGP gas optics 
         minor_gases_lower_sw,               & ! used by RRTMGP gas optics 
         minor_gases_upper_sw,               & ! used by RRTMGP gas optics 
         scaling_gas_lower_sw,               & ! used by RRTMGP gas optics 
         scaling_gas_upper_sw                  ! used by RRTMGP gas optics 
    logical(wl), dimension(:), allocatable :: &
         minor_scales_with_density_lower_sw, & ! used by RRTMGP gas optics 
         minor_scales_with_density_upper_sw, & ! used by RRTMGP gas optics 
         scale_by_complement_lower_sw,       & ! used by RRTMGP gas optics 
         scale_by_complement_upper_sw          ! used by RRTMGP gas optics 
    ! Dimensions (to be broadcast across all processors)
    integer :: &
         ntemps_sw,                          & ! used by RRTMGP gas optics 
         npress_sw,                          & ! used by RRTMGP gas optics 
         nabsorbers_sw,                      & ! used by RRTMGP gas optics 
         nextrabsorbers_sw,                  & ! used by RRTMGP gas optics 
         nminorabsorbers_sw,                 & ! used by RRTMGP gas optics 
         nmixingfracs_sw,                    & ! used by RRTMGP gas optics 
         nlayers_sw,                         & ! used by RRTMGP gas optics 
         nbnds_sw,                           & ! used by RRTMGP gas optics 
         npairs_sw,                          & ! used by RRTMGP gas optics 
         nminor_absorber_intervals_lower_sw, & ! used by RRTMGP gas optics 
         nminor_absorber_intervals_upper_sw, & ! used by RRTMGP gas optics 
         ncontributors_lower_sw,             & ! used by RRTMGP gas optics 
         ncontributors_upper_sw                ! used by RRTMGP gas optics 

    ! Local variables
    integer :: status,ncid_sw,ncid_sw_clds,dimid,varID,ij,iGas
    integer,dimension(:),allocatable :: temp1,temp2,temp3,temp4,temp_log_array1,&
         temp_log_array2, temp_log_array3, temp_log_array4
    character(len=264) :: sw_gas_props_file

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the gfs_physics_nml (scm/src/GFS_typedefs.F90)
    sw_gas_props_file   = trim(Model%rrtmgp_root)//trim(Model%sw_file_gas)

    ! Read dimensions for k-distribution fields (only on master processor(0))
    if (mpirank .eq. mpiroot) then
       if(nf90_open(trim(sw_gas_props_file), NF90_WRITE, ncid_sw) .eq. NF90_NOERR) then
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

    ! On master processor,  read in fields, broadcast to all processors
    if (mpirank .eq. mpiroot) then
       ! Read in fields from file
       if(nf90_open(trim(sw_gas_props_file), NF90_WRITE, ncid_sw) .eq. NF90_NOERR) then
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
    do iGas=1,Model%nGases
       call check_error_msg('sw_gas_optics_init',gas_concentrations%set_vmr(Radtend%active_gases(iGas), 0._kind_phys))
    enddo
    call check_error_msg('sw_gas_optics_init',sw_gas_props%load(gas_concentrations, gas_names_sw,   &
         key_species_sw, band2gpt_sw, band_lims_sw, press_ref_sw, press_ref_trop_sw, temp_ref_sw,   &
         temp_ref_p_sw, temp_ref_t_sw, vmr_ref_sw, kmajor_sw, kminor_lower_sw, kminor_upper_sw,     &
         gas_minor_sw,identifier_minor_sw, minor_gases_lower_sw, minor_gases_upper_sw,              &
         minor_limits_gpt_lower_sw,minor_limits_gpt_upper_sw, minor_scales_with_density_lower_sw,   &
         minor_scales_with_density_upper_sw, scaling_gas_lower_sw,                                  &
         scaling_gas_upper_sw, scale_by_complement_lower_sw,                                        &
         scale_by_complement_upper_sw, kminor_start_lower_sw, kminor_start_upper_sw,                &
         solar_source_sw, rayl_lower_sw, rayl_upper_sw))

    ! Set initial permutation seed for McICA, initially set to number of G-points
    ipsdsw0 = sw_gas_props%get_ngpt()
  end subroutine rrtmgp_sw_gas_optics_init

  subroutine rrtmgp_sw_gas_optics_run()
  end subroutine rrtmgp_sw_gas_optics_run
  subroutine rrtmgp_sw_gas_optics_finalize()
  end subroutine rrtmgp_sw_gas_optics_finalize
  ! #########################################################################################
  ! SUBROUTINE check_error_msg
  ! #########################################################################################
  subroutine check_error_msg(routine_name, error_msg)
    character(len=*), intent(in) :: &
         error_msg, routine_name
    
    if(error_msg /= "") then
       print*,"ERROR("//trim(routine_name)//"): "
       print*,trim(error_msg)
       return
    end if
  end subroutine check_error_msg  
end module rrtmgp_sw_gas_optics
 
