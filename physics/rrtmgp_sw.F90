! ###########################################################################################
! ###########################################################################################
module rrtmgp_sw
  use machine,                 only: kind_phys
  use GFS_typedefs,            only: GFS_control_type
  use mo_rte_kind,             only: wl
  use mo_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,         only: ty_cloud_optics
  use mo_optical_props,        only: ty_optical_props_2str
  use mo_rrtmgp_clr_all_sky,   only: rte_sw
  use mo_gas_concentrations,   only: ty_gas_concs
  use mo_fluxes_byband,        only: ty_fluxes_byband
  use module_radsw_parameters, only: cmpfsw_type

  ! Parameters
  integer,parameter :: nGases = 6
  real(kind_phys),parameter :: epsilon=1.0e-6
  character(len=3),parameter, dimension(nGases) :: &
       active_gases = (/ 'h2o', 'co2', 'o3 ', 'n2o', 'ch4', 'o2 '/)
  integer :: nrghice, ipsdsw0

  public rrtmgp_sw_init, rrtmgp_sw_run, rrtmgp_sw_finalize
contains


!! \section arg_table_rrtmgp_sw_init Argument Table
!! | local_name         | standard_name                                   | long_name                                                                 | units | rank | type                 |    kind   | intent | optional |
!! |--------------------|-------------------------------------------------|---------------------------------------------------------------------------|-------|------|----------------------|-----------|--------|----------|
!! | Model              | GFS_control_type_instance                       | Fortran DDT containing FV3-GFS model control parameters                   | DDT   |    0 | GFS_control_type     |           | in     | F        |
!! | mpirank            | mpi_rank                                        | current MPI rank                                                          | index |    0 | integer              |           | in     | F        |
!! | mpiroot            | mpi_root                                        | master MPI rank                                                           | index |    0 | integer              |           | in     | F        |
!! | mpicomm            | mpi_comm                                        | MPI communicator                                                          | index |    0 | integer              |           | in     | F        |
!! | errmsg             | ccpp_error_message                              | error message for error handling in CCPP                                  | none  |    0 | character            | len=*     | out    | F        |
!! | errflg             | ccpp_error_flag                                 | error flag for error handling in CCPP                                     | flag  |    0 | integer              |           | out    | F        |
!! | kdist_sw           | K_distribution_file_for_RRTMGP_SW_scheme        | DDT containing spectral information for RRTMGP SW radiation scheme        | DDT   |    0 | ty_gas_optics_rrtmgp |           | inout  | F        |
!! | kdist_cldy_sw      | K_distribution_file_for_cloudy_RRTMGP_SW_scheme | DDT containing spectral information for cloudy RRTMGP SW radiation scheme | DDT   |    0 | ty_cloud_optics      |           | inout  | F        |
!!
  ! #########################################################################################
  subroutine rrtmgp_sw_init(Model,mpicomm, mpirank, mpiroot, kdist_sw, kdist_cldy_sw,   &
        errmsg, errflg)
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
    type(ty_gas_optics_rrtmgp),intent(inout) :: &
         kdist_sw
    type(ty_cloud_optics),intent(inout) :: &
         kdist_cldy_sw

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer,          intent(out) :: &
         errflg     ! Error code

    ! Fields from the K-distribution files
    ! Variables that will be passed to gas_optics%load()
    type(ty_gas_concs)  :: &
         gas_concentrations
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
         temp_ref_t_sw,                      & ! used by RRTGMP gas optics 
         radliq_lwr_sw,                      & ! used by RRTGMP cloud optics 
         radliq_upr_sw,                      & ! used by RRTGMP cloud optics 
         radliq_fac_sw,                      & ! used by RRTGMP cloud optics 
         radice_lwr_sw,                      & ! used by RRTGMP cloud optics 
         radice_upr_sw,                      & ! used by RRTGMP cloud optics 
         radice_fac_sw                         ! used by RRTGMP cloud optics 

    real(kind_phys), dimension(:), allocatable :: &
         press_ref_sw,                       & ! used by RRTGMP gas optics 
         temp_ref_sw,                        & ! used by RRTGMP gas optics 
         solar_source_sw,                    & ! used by RRTGMP gas optics 
         pade_sizereg_extliq_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_ssaliq_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_asyliq_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_extice_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_ssaice_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_asyice_sw                ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:,:), allocatable :: &
         band_lims_sw,                       &  ! used by RRTGMP gas optics 
         lut_extliq_sw,                      & ! used by RRTGMP cloud optics 
         lut_ssaliq_sw,                      & ! used by RRTGMP cloud optics 
         lut_asyliq_sw,                      & ! used by RRTGMP cloud optics 
         band_lims_cldy_sw                     ! used by RRTGMP cloud optics                          

    real(kind_phys), dimension(:,:,:), allocatable :: &
         vmr_ref_sw,                         & ! used by RRTGMP gas optics 
         kminor_lower_sw,                    & ! used by RRTGMP gas optics 
         kminor_upper_sw,                    & ! used by RRTGMP gas optics 
         rayl_lower_sw,                      & ! used by RRTGMP gas optics 
         rayl_upper_sw,                      & ! used by RRTGMP gas optics 
         lut_extice_sw,                      & ! used by RRTGMP cloud optics 
         lut_ssaice_sw,                      & ! used by RRTGMP cloud optics 
         lut_asyice_sw,                      & ! used by RRTGMP cloud optics 
         pade_extliq_sw,                     & ! used by RRTGMP cloud optics 
         pade_ssaliq_sw,                     & ! used by RRTGMP cloud optics 
         pade_asyliq_sw                        ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         kmajor_sw,                          & ! used by RRTGMP gas optics 
         pade_extice_sw,                     & ! used by RRTGMP cloud optics 
         pade_ssaice_sw,                     & ! used by RRTGMP cloud optics 
         pade_asyice_sw                        ! used by RRTGMP cloud optics 
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
         ncontributors_upper_sw,             & ! used by RRTGMP gas optics 
         nbandSWcldy_sw,                     & ! used by RRTGMP cloud optics 
         nsize_liq_sw,                       & ! used by RRTGMP cloud optics 
         nsize_ice_sw,                       & ! used by RRTGMP cloud optics 
         nsizereg_sw,                        & ! used by RRTGMP cloud optics 
         ncoeff_ext_sw,                      & ! used by RRTGMP cloud optics 
         ncoeff_ssa_g_sw,                    & ! used by RRTGMP cloud optics 
         nbound_sw,                          & ! used by RRTGMP cloud optics  
         npairsSWcldy_sw                       ! used by RRTGMP cloud optics 

    ! Local variables
    integer :: status,ncid_sw,ncid_sw_clds,dimid,varID,ij,iGas
    character(len=264) :: kdist_file, kdist_cldy_file
    integer,dimension(:),allocatable :: temp1,temp2,temp3,temp4,temp_log_array1,&
         temp_log_array2, temp_log_array3, temp_log_array4

    ! Initialize
    errmsg = ''
    errflg = 0

    ! How are we handling cloud-optics?
    rrtmgp_sw_cld_phys = Model%rrtmgp_cld_phys

    ! Filenames are set in the gfs_physics_nml (scm/src/GFS_typedefs.F90)
    kdist_file      = trim(Model%rrtmgp_root)//trim(Model%kdist_sw_file_gas)
    kdist_cldy_file = trim(Model%rrtmgp_root)//trim(Model%kdist_sw_file_clouds)

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
       call check_error_msg('rrtmgp_sw_init',gas_concentrations%set_vmr(active_gases(iGas), 0._kind_phys))
    enddo
    call check_error_msg('rrtmgp_sw_init',kdist_sw%load(gas_concentrations, gas_names_sw, key_species_sw, band2gpt_sw, &
         band_lims_sw, press_ref_sw, press_ref_trop_sw, temp_ref_sw,  temp_ref_p_sw, temp_ref_t_sw, &
         vmr_ref_sw, kmajor_sw, kminor_lower_sw, kminor_upper_sw, gas_minor_sw,identifier_minor_sw, &
         minor_gases_lower_sw, minor_gases_upper_sw, minor_limits_gpt_lower_sw,                     &
         minor_limits_gpt_upper_sw, minor_scales_with_density_lower_sw,                             &
         minor_scales_with_density_upper_sw, scaling_gas_lower_sw,                                  &
         scaling_gas_upper_sw, scale_by_complement_lower_sw,                                        &
         scale_by_complement_upper_sw, kminor_start_lower_sw, kminor_start_upper_sw,                &
         solar_source_sw, rayl_lower_sw, rayl_upper_sw))

    ! Set band index by g-point array 
    nBandsSW = kdist_sw%get_nband()
    nGptsSW  = kdist_sw%get_ngpt()

    ! Set initial permutation seed for McICA, initially set to number of G-points
    ipsdsw0 = kdist_sw%get_ngpt()

    ! #######################################################################################
    ! If RRTMGP cloud-optics are requested, read tables and broadcast.
    ! #######################################################################################
    ! Read dimensions for k-distribution fields (only on master processor(0))
    if (mpirank .eq. mpiroot) then
       if(nf90_open(trim(kdist_cldy_file), NF90_WRITE, ncid_sw_clds) == NF90_NOERR) then
          status = nf90_inq_dimid(ncid_sw_clds, 'nband', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=nbandSWcldy_sw)
          status = nf90_inq_dimid(ncid_sw_clds, 'nrghice', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=nrghice)
          status = nf90_inq_dimid(ncid_sw_clds, 'nsize_liq', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=nsize_liq_sw)
          status = nf90_inq_dimid(ncid_sw_clds, 'nsize_ice', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=nsize_ice_sw)
          status = nf90_inq_dimid(ncid_sw_clds, 'nsizereg', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=nsizereg_sw)
          status = nf90_inq_dimid(ncid_sw_clds, 'ncoeff_ext', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=ncoeff_ext_sw)
          status = nf90_inq_dimid(ncid_sw_clds, 'ncoeff_ssa_g', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=ncoeff_ssa_g_sw)
          status = nf90_inq_dimid(ncid_sw_clds, 'nbound', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=nbound_sw)
          status = nf90_inq_dimid(ncid_sw_clds, 'pair', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=npairsSWcldy_sw)
          status = nf90_close(ncid_sw_clds)
       endif
    endif

    ! Broadcast dimensions to all processors
#ifdef MPI
    if (rrtmgp_sw_cld_phys .eq. 1 .or. rrtmgp_sw_cld_phys .eq. 2) then
       call MPI_BCAST(nbandSWcldy_sw,  1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nrghice,         1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_liq_sw,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_ice_sw,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsizereg_sw,     1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ext_sw,   1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ssa_g_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nbound_sw,       1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(npairsSWcldy_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    endif
#endif

    if (rrtmgp_sw_cld_phys .eq. 1) then
       allocate(lut_extliq_sw(nsize_liq_sw, nBandSWcldy_sw))
       allocate(lut_ssaliq_sw(nsize_liq_sw, nBandSWcldy_sw))
       allocate(lut_asyliq_sw(nsize_liq_sw, nBandSWcldy_sw))
       allocate(lut_extice_sw(nsize_ice_sw, nBandSWcldy_sw, nrghice))
       allocate(lut_ssaice_sw(nsize_ice_sw, nBandSWcldy_sw, nrghice))
       allocate(lut_asyice_sw(nsize_ice_sw, nBandSWcldy_sw, nrghice))
       allocate(band_lims_cldy_sw(2, nBandSWcldy_sw))
    endif
    if (rrtmgp_sw_cld_phys .eq. 2) then
       allocate(pade_extliq_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ext_sw ))
       allocate(pade_ssaliq_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ssa_g_sw))
       allocate(pade_asyliq_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ssa_g_sw))
       allocate(pade_extice_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ext_sw,   nrghice))
       allocate(pade_ssaice_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ssa_g_sw, nrghice))
       allocate(pade_asyice_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ssa_g_sw, nrghice))
       allocate(pade_sizereg_extliq_sw(nbound_sw))
       allocate(pade_sizereg_ssaliq_sw(nbound_sw))
       allocate(pade_sizereg_asyliq_sw(nbound_sw))
       allocate(pade_sizereg_extice_sw(nbound_sw))
       allocate(pade_sizereg_ssaice_sw(nbound_sw))
       allocate(pade_sizereg_asyice_sw(nbound_sw))
       allocate(band_lims_cldy_sw(2,nbandSWcldy_sw))
    endif

    ! On master processor, allocate space, read in fields, broadcast to all processors
    if (mpirank .eq. mpiroot) then
       ! 
       if (rrtmgp_sw_cld_phys .eq. 1) then
          !
          if(nf90_open(trim(kdist_cldy_file), NF90_WRITE, ncid_sw_clds) == NF90_NOERR) then
             status = nf90_inq_varid(ncid_sw_clds,'radliq_lwr',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radliq_lwr_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radliq_upr',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radliq_upr_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radliq_fac',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radliq_fac_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radice_lwr',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radice_lwr_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radice_upr',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radice_upr_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radice_fac',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radice_fac_sw)
             status = nf90_inq_varid(ncid_sw_clds,'lut_extliq',varID)
             status = nf90_get_var(ncid_sw_clds,varID,lut_extliq_sw)
             status = nf90_inq_varid(ncid_sw_clds,'lut_ssaliq',varID)
             status = nf90_get_var(ncid_sw_clds,varID,lut_ssaliq_sw)
             status = nf90_inq_varid(ncid_sw_clds,'lut_asyliq',varID)
             status = nf90_get_var(ncid_sw_clds,varID,lut_asyliq_sw)
             status = nf90_inq_varid(ncid_sw_clds,'lut_extice',varID)
             status = nf90_get_var(ncid_sw_clds,varID,lut_extice_sw)
             status = nf90_inq_varid(ncid_sw_clds,'lut_ssaice',varID)
             status = nf90_get_var(ncid_sw_clds,varID,lut_ssaice_sw)
             status = nf90_inq_varid(ncid_sw_clds,'lut_asyice',varID)
             status = nf90_get_var(ncid_sw_clds,varID,lut_asyice_sw)
             status = nf90_inq_varid(ncid_sw_clds,'bnd_limits_wavenumber',varID)
             status = nf90_get_var(ncid_sw_clds,varID,band_lims_cldy_sw)
             status = nf90_close(ncid_sw_clds)
          endif
       endif
       !
       if (rrtmgp_sw_cld_phys .eq. 2) then
          !
          if(nf90_open(trim(kdist_cldy_file), NF90_WRITE, ncid_sw_clds) == NF90_NOERR) then
             status = nf90_inq_varid(ncid_sw_clds,'radliq_lwr',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radliq_lwr_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radliq_upr',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radliq_upr_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radliq_fac',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radliq_fac_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radice_lwr',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radice_lwr_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radice_upr',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radice_upr_sw)
             status = nf90_inq_varid(ncid_sw_clds,'radice_fac',varID)
             status = nf90_get_var(ncid_sw_clds,varID,radice_fac_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_extliq',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_extliq_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_ssaliq',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_ssaliq_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_asyliq',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_asyliq_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_extice',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_extice_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_ssaice',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_ssaice_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_asyice',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_asyice_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_sizreg_extliq',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_sizereg_extliq_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_sizreg_ssaliq',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_sizereg_ssaliq_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_sizreg_asyliq',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_sizereg_asyliq_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_sizreg_extice',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_sizereg_extice_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_sizreg_ssaice',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_sizereg_ssaice_sw)
             status = nf90_inq_varid(ncid_sw_clds,'pade_sizreg_asyice',varID)
             status = nf90_get_var(ncid_sw_clds,varID,pade_sizereg_asyice_sw)
             status = nf90_inq_varid(ncid_sw_clds,'bnd_limits_wavenumber',varID)
             status = nf90_get_var(ncid_sw_clds,varID,band_lims_cldy_sw)
             status = nf90_close(ncid_sw_clds)
          endif
       endif
    endif

    ! Broadcast arrays to all processors
#ifdef MPI
    if (rrtmgp_sw_cld_phys .eq. 1) then
       call MPI_BCAST(radliq_lwr_sw,           size(radliq_lwr_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radliq_upr_sw,           size(radliq_upr_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radliq_fac_sw,           size(radliq_fac_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_lwr_sw,           size(radice_lwr_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_upr_sw,           size(radice_upr_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_fac_sw,           size(radice_fac_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_extliq_sw,           size(lut_extliq_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_ssaliq_sw,           size(lut_ssaliq_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_asyliq_sw,           size(lut_asyliq_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_extice_sw,           size(lut_extice_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_ssaice_sw,           size(lut_ssaice_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_asyice_sw,           size(lut_asyice_sw),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(band_lims_cldy_sw),      size(band_lims_cldy_sw),      kind_phys,   mpiroot, mpicomm, ierr)    
    endif
    if (rrtmgp_sw_cld_phys .eq. 2) then
       call MPI_BCAST(pade_extliq_sw,          size(pade_extliq_sw),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_ssaliq_sw,          size(pade_ssaliq_sw),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_asyliq_sw,          size(pade_asyliq_sw),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_extice_sw,          size(pade_extice_sw),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_ssaice_sw,          size(pade_ssaice_sw),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_asyice_sw,          size(pade_asyice_sw),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_extliq_sw), size(pade_sizereg_extliq_sw), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_ssaliq_sw), size(pade_sizereg_ssaliq_sw), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_asyliq_sw), size(pade_sizereg_asyliq_sw), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_extice_sw), size(pade_sizereg_extice_sw), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_ssaice_sw), size(pade_sizereg_ssaice_sw), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_asyice_sw), size(pade_sizereg_asyice_sw), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(band_lims_cldy_sw),      size(band_lims_cldy_sw),      kind_phys,   mpiroot, mpicomm, ierr)    
    endif
#endif

    ! Load tables data for RRTGMP cloud-optics  
    if (rrtmgp_sw_cld_phys .eq. 1) then
       call check_error_msg('rrtmgp_sw_init',kdist_cldy_sw%set_ice_roughness(nrghice))
       call check_error_msg('rrtmgp_sw_init',kdist_cldy_sw%load(band_lims_cldy_sw, radliq_lwr_sw,            &
            radliq_upr_sw, radliq_fac_sw, radice_lwr_sw, radice_upr_sw, radice_fac_sw,      &
            lut_extliq_sw, lut_ssaliq_sw, lut_asyliq_sw, lut_extice_sw, lut_ssaice_sw,      &
            lut_asyice_sw))
    endif
    if (rrtmgp_sw_cld_phys .eq. 2) then
       call check_error_msg('rrtmgp_sw_init',kdist_cldy_sw%set_ice_roughness(nrghice))
       call check_error_msg('rrtmgp_sw_init', kdist_cldy_sw%load(band_lims_cldy_sw, pade_extliq_sw,           &
            pade_ssaliq_sw, pade_asyliq_sw, pade_extice_sw, pade_ssaice_sw, pade_asyice_sw, &
            pade_sizereg_extliq_sw, pade_sizereg_ssaliq_sw, pade_sizereg_asyliq_sw,         &
            pade_sizereg_extice_sw, pade_sizereg_ssaice_sw, pade_sizereg_asyice_sw))
    endif

  end subroutine rrtmgp_sw_init

  ! #########################################################################################
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_run Argument Table
!! | local_name              | standard_name                                                                                  | long_name                                                                | units | rank | type                  |    kind   | intent | optional |
!! |-------------------------|------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|-------|------|-----------------------|-----------|--------|----------|
!! | ncol                    | horizontal_loop_extent                                                                         | horizontal dimension                                                     | count |    0 | integer               |           | in     | F        |
!! | nlay                    | adjusted_vertical_layer_dimension_for_radiation                                                | number of vertical layers for radiation                                  | count |    0 | integer               |           | in     | F        |
!! | p_lay                   | air_pressure_at_layer_for_radiation_in_hPa                                                     | air pressure layer                                                       | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | p_lev                   | air_pressure_at_interface_for_radiation_in_hPa                                                 | air pressure level                                                       | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | t_lay                   | air_temperature_at_layer_for_radiation                                                         | air temperature layer                                                    | K     |    2 | real                  | kind_phys | in     | F        |
!! | kdist_sw                | K_distribution_file_for_RRTMGP_SW_scheme                                                       | DDT containing spectral information for RRTMGP SW radiation scheme       | DDT   |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | optical_propsSW_clds    | shortwave_optical_properties_for_cloudy_atmosphere                                             | Fortran DDT containing RRTMGP optical properties                         | DDT   |    0 | ty_optical_props_2str |           | in     | F        |
!! | optical_propsSW_aerosol | shortwave_optical_properties_for_aerosols                                                      | Fortran DDT containing RRTMGP optical properties                         | DDT   |    0 | ty_optical_props_2str |           | in     | F        |
!! | gas_concentrations      | Gas_concentrations_for_RRTMGP_suite_sw                                                         | DDT containing gas concentrations for RRTMGP radiation scheme            | DDT   |    0 | ty_gas_concs          |           | in     | F        |
!! | lsswr                   | flag_to_calc_sw                                                                                | flag to calculate SW irradiances                                         | flag  |    0 | logical               |           | in     | F        |
!! | sfcalb_nir_dir          | surface_shortwave_albedo_near_infrared_direct_in_each_band                                     | surface sw near-infrared direct albedo in each SW band                   | frac  |    2 | real                  | kind_phys | out    | F        |
!! | sfcalb_nir_dif          | surface_shortwave_albedo_near_infrared_diffuse_in_each_band                                    | surface sw near-infrared diffuse albedo in each SW band                  | frac  |    2 | real                  | kind_phys | out    | F        |
!! | cossza                  | cosine_of_zenith_angle                                                                         | cosine of the solar zenit angle                                          | none  |    1 | real                  | kind_phys | in     | F        |
!! | nday                    | daytime_points_dimension                                                                       | daytime points dimension                                                 | count |    0 | integer               |           | in     | F        |
!! | idxday                  | daytime_points                                                                                 | daytime points                                                           | index |    1 | integer               |           | in     | F        |
!! | fluxSW_allsky           | sw_flux_profiles_byband_allsky                                                                 | Fortran DDT containing RRTMGP 3D fluxes                                  | DDT   |    0 | ty_fluxes_byband      |           | out    | F        |
!! | fluxSW_clrsky           | sw_flux_profiles_byband_clrsky                                                                 | Fortran DDT containing RRTMGP 3D fluxes                                  | DDT   |    0 | ty_fluxes_byband      |           | out    | F        |
!! | hsw0                    | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | shortwave clear sky heating rate                                         | K s-1 |    2 | real                  | kind_phys | in     | T        |
!! | hswb                    | sw_heating_rate_spectral                                                                       | shortwave total sky heating rate (spectral)                              | K s-1 |    3 | real                  | kind_phys | in     | T        |
!! | scmpsw                  | components_of_surface_downward_shortwave_fluxes                                                | derived type for special components of surface downward shortwave fluxes | W m-2 |    1 | cmpfsw_type           |           | inout  | F        |
!! | errmsg                  | ccpp_error_message                                                                             | error message for error handling in CCPP                                 | none  |    0 | character             | len=*     | out    | F        |
!! | errflg                  | ccpp_error_flag                                                                                | error flag for error handling in CCPP                                    | flag  |    0 | integer               |           | out    | F        |
!!
  subroutine rrtmgp_sw_run(ncol, nlay, kdist_sw, p_lay, t_lay, p_lev, gas_concentrations, &
       optical_propsSW_clds, optical_propsSW_aerosol,&
       lsswr, sfcalb_nir_dir, sfcalb_nir_dif, cossza,  nday, idxday, fluxSW_allsky, fluxSW_clrsky, hsw0, hswb, scmpsw, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: &
         ncol,                 & ! Number of horizontal gridpoints
         nlay,                 & ! Number of vertical layers
         nday                    ! Number of daytime points
    integer, intent(in), dimension(nday) :: &
         idxday                  ! Index array for daytime points
    real(kind_phys), dimension(ncol,nlay), intent(in) :: &
         p_lay,                & ! Pressure @ model layer-centers         (hPa)
         t_lay                   ! Temperature                            (K)
    real(kind_phys), dimension(ncol,nlay+1), intent(in) :: &
         p_lev                   ! Pressure @ model layer-interfaces      (hPa)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         kdist_sw                ! DDT containing SW spectral information
    real(kind_phys), dimension(kdist_sw%get_nband(),ncol), intent(in) :: &
         sfcalb_nir_dir,  & ! Surface albedo direct (near-IR)         (1)
         sfcalb_nir_dif     ! Surface albedo diffuse (near-IR)        (1)
    real(kind_phys), dimension(ncol), intent(in) :: &
         cossza             ! Cosine of solar zenith angle            (1)
    type(ty_optical_props_2str),intent(in) :: &
         optical_propsSW_clds, & ! RRTMGP DDT: longwave cloud radiative properties 
         optical_propsSW_aerosol ! RRTMGP DDT: longwave aerosol radiative properties

    type(ty_gas_concs),intent(in) :: &
         gas_concentrations      ! RRTMGP DDT: trace gas concentrations   (vmr)
    logical, intent(in) :: &
         lsswr                   ! Flag to calculate SW irradiances

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg
    type(ty_fluxes_byband),intent(out) :: &
         fluxSW_allsky, & ! All-sky flux                      (W/m2)
         fluxSW_clrsky    ! Clear-sky flux                    (W/m2)

    ! Inputs (optional) (NOTE. We only need the optional arguments to know what fluxes to output, HR's are computed later)
    real(kind_phys), dimension(ncol,nlay), optional, intent(in) :: &
         hsw0             ! Clear-sky heating rate            (K/sec)
    real(kind_phys), dimension(ncol,nlay,kdist_sw%get_nband()), intent(in), optional :: &
         hswb             ! All-sky heating rate, by band     (K/sec)
    ! Outputs (optional)
    type(cmpfsw_type), dimension(ncol), intent(out),optional :: &
         scmpsw           ! 2D surface fluxes, components:
                          ! uvbfc - total sky downward uv-b flux at  (W/m2)
                          ! uvbf0 - clear sky downward uv-b flux at  (W/m2)
                          ! nirbm - downward nir direct beam flux    (W/m2)
                          ! nirdf - downward nir diffused flux       (W/m2)
                          ! visbm - downward uv+vis direct beam flux (W/m2)
                          ! visdf - downward uv+vis diffused flux    (W/m2)


    ! Local variables
    real(kind_phys), dimension(nday,nlay+1),target :: &
         fluxSW_up_allsky, fluxSW_up_clrsky, fluxSW_dn_allsky, fluxSW_dn_clrsky
    real(kind_phys), dimension(nday,nlay+1,kdist_sw%get_nband()),target :: &
         fluxSWBB_up_allsky, fluxSWBB_dn_allsky
    logical :: l_ClrSky_HR=.false., l_AllSky_HR_byband=.false., l_scmpsw=.false.

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    if (.not. lsswr) return

    ! Are any optional outputs requested? Need to know now to compute correct fluxes.
    l_ClrSky_HR        = present(hsw0)
    l_AllSky_HR_byband = present(hswb)
    l_scmpsw           = present(scmpsw)
    if ( l_scmpsw ) then
       scmpsw = cmpfsw_type (0., 0., 0., 0., 0., 0.)
    endif

    if (nDay .gt. 0) then
       ! Initialize RRTMGP DDT containing 2D(3D) fluxes
       fluxSW_allsky%flux_up => fluxSW_up_allsky
       fluxSW_allsky%flux_dn => fluxSW_dn_allsky
       fluxSW_clrsky%flux_up => fluxSW_up_clrsky
       fluxSW_clrsky%flux_dn => fluxSW_dn_clrsky
       ! Only calculate fluxes by-band, only when heating-rate profiles by band are requested.
       if (l_AllSky_HR_byband) then
          fluxSW_allsky%bnd_flux_up => fluxSWBB_up_allsky
          fluxSW_allsky%bnd_flux_dn => fluxSWBB_dn_allsky
       endif
       
       ! Call RRTMGP SW scheme
       call check_error_msg('rrtmgp_sw_run',rte_sw(               &
            kdist_sw,                             & ! IN  - spectral information 
            gas_concentrations,                   & ! IN  - gas concentrations (vmr)
            p_lay(idxday,1:nlay),                 & ! IN  - pressure at layer interfaces (Pa)  
            t_lay(idxday,1:nlay),                 & ! IN  - temperature at layer interfaes (K)
            p_lev(idxday,1:nlay+1),               & ! IN  - pressure at layer centers (Pa)
            cossza(idxday),                       & ! IN  - Cosine of solar zenith angle
            sfcalb_nir_dir(:,idxday),             & ! IN  - Shortwave surface albedo (direct)
            sfcalb_nir_dif(:,idxday),             & ! IN  - Shortwave surface albedo (diffuse)
            optical_propsSW_clds,                 & ! IN  - DDT containing cloud optical information 
            fluxSW_allsky,                        & ! OUT - Fluxes, all-sky, 3D (nCol,nLay,nBand) 
            fluxSW_clrsky,                        & ! OUT - Fluxes, clear-sky, 3D (nCol,nLay,nBand) 
            aer_props = optical_propsSW_aerosol))   ! IN(optional) - DDT containing aerosol optical information
    endif

  end subroutine rrtmgp_sw_run
  
  subroutine rrtmgp_sw_finalize()
  end subroutine rrtmgp_sw_finalize
  subroutine check_error_msg(routine_name, error_msg)
    character(len=*), intent(in) :: &
         error_msg, routine_name
    
    if(error_msg /= "") then
       print*,"ERROR("//trim(routine_name)//"): "
       print*,trim(error_msg)
       return
    end if
  end subroutine check_error_msg    


end module rrtmgp_sw
