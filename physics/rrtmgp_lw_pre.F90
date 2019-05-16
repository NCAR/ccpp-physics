!>\file rrtmgp_lw_pre.f90
!! This file contains a call to module_radiation_surface::setemis() to
!! setup surface emissivity for LW radiation.
module rrtmgp_lw_pre
  use machine,                   only: kind_phys
  use GFS_typedefs,              only: GFS_control_type
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,     only: ty_gas_concs
  use mo_cloud_optics,           only: ty_cloud_optics
  use mo_rte_kind,               only: wl
  use rrtmgp_lw,                 only: check_error_msg
        
  ! Parameters
  integer,parameter :: nGases = 6
  real(kind_phys),parameter :: epsilon=1.0e-6
  character(len=3),parameter, dimension(nGases) :: &
       active_gases = (/ 'h2o', 'co2', 'o3 ', 'n2o', 'ch4', 'o2 '/)
  integer :: nrghice, ipsdlw0
  
contains

!! \section arg_table_rrtmgp_lw_pre_init Argument Table
!! | local_name         | standard_name                                   | long_name                                                                 | units | rank | type                 |    kind   | intent | optional |
!! |--------------------|-------------------------------------------------|---------------------------------------------------------------------------|-------|------|----------------------|-----------|--------|----------|
!! | Model              | GFS_control_type_instance                       | Fortran DDT containing FV3-GFS model control parameters                   | DDT   |    0 | GFS_control_type     |           | in     | F        |
!! | mpirank            | mpi_rank                                        | current MPI rank                                                          | index |    0 | integer              |           | in     | F        |
!! | mpiroot            | mpi_root                                        | master MPI rank                                                           | index |    0 | integer              |           | in     | F        |
!! | mpicomm            | mpi_comm                                        | MPI communicator                                                          | index |    0 | integer              |           | in     | F        |
!! | errmsg             | ccpp_error_message                              | error message for error handling in CCPP                                  | none  |    0 | character            | len=*     | out    | F        |
!! | errflg             | ccpp_error_flag                                 | error flag for error handling in CCPP                                     | flag  |    0 | integer              |           | out    | F        |
!! | kdist_lw           | K_distribution_file_for_RRTMGP_LW_scheme        | DDT containing spectral information for RRTMGP LW radiation scheme        | DDT   |    0 | ty_gas_optics_rrtmgp |           | out    | F        |
!! | kdist_cldy_lw      | K_distribution_file_for_cloudy_RRTMGP_LW_scheme | DDT containing spectral information for cloudy RRTMGP LW radiation scheme | DDT   |    0 | ty_cloud_optics      |           | out    | F        |
!! | gas_concentrations | Gas_concentrations_for_RRTMGP_suite             | DDT containing gas concentrations for RRTMGP radiation scheme             | DDT   |    0 | ty_gas_concs         |           | out    | F        |
!!
  ! #########################################################################################
  subroutine rrtmgp_lw_pre_init(Model,mpicomm, mpirank, mpiroot, kdist_lw, kdist_cldy_lw,   &
       gas_concentrations, errmsg, errflg)
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
    type(ty_gas_optics_rrtmgp),intent(out) :: &
         kdist_lw
    type(ty_cloud_optics),intent(out) :: &
         kdist_cldy_lw
    type(ty_gas_concs),intent(out) :: &
         gas_concentrations
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer,          intent(out) :: &
         errflg     ! Error code

    ! Variables that will be passed to gas_optics%load()
    integer, dimension(:), allocatable :: &
         kminor_start_lower,              & ! used by RRTGMP gas optics 
         kminor_start_upper                 ! used by RRTGMP gas optics 
    integer, dimension(:,:), allocatable :: &
         band2gpt,                        & ! used by RRTGMP gas optics 
         minor_limits_gpt_lower,          & ! used by RRTGMP gas optics 
         minor_limits_gpt_upper             ! used by RRTGMP gas optics 
    integer, dimension(:,:,:), allocatable :: &
         key_species                        ! used by RRTGMP gas optics 
    real(kind_phys) :: &
         press_ref_trop,                  & ! used by RRTGMP gas optics 
         temp_ref_p,                      & ! used by RRTGMP gas optics 
         temp_ref_t,                      & ! used by RRTGMP gas optics 
         radliq_lwr,                      & ! used by RRTGMP cloud optics 
         radliq_upr,                      & ! used by RRTGMP cloud optics 
         radliq_fac,                      & ! used by RRTGMP cloud optics 
         radice_lwr,                      & ! used by RRTGMP cloud optics 
         radice_upr,                      & ! used by RRTGMP cloud optics 
         radice_fac                         ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:), allocatable :: &
         press_ref,                       & ! used by RRTGMP gas optics 
         temp_ref,                        & ! used by RRTGMP gas optics 
         pade_sizereg_extliq,             & ! used by RRTGMP cloud optics 
         pade_sizereg_ssaliq,             & ! used by RRTGMP cloud optics 
         pade_sizereg_asyliq,             & ! used by RRTGMP cloud optics 
         pade_sizereg_extice,             & ! used by RRTGMP cloud optics 
         pade_sizereg_ssaice,             & ! used by RRTGMP cloud optics 
         pade_sizereg_asyice                ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:,:), allocatable :: &
         band_lims,                       & ! used by RRTGMP gas optics 
         totplnk,                         & ! used by RRTGMP gas optics 
         lut_extliq,                      & ! used by RRTGMP cloud optics 
         lut_ssaliq,                      & ! used by RRTGMP cloud optics 
         lut_asyliq,                      & ! used by RRTGMP cloud optics 
         band_lims_cldy                     ! used by RRTGMP cloud optics 

    real(kind_phys), dimension(:,:,:), allocatable :: &
         vmr_ref,                         & ! used by RRTGMP gas optics 
         kminor_lower,                    & ! used by RRTGMP gas optics 
         kminor_upper,                    & ! used by RRTGMP gas optics 
         rayl_lower,                      & ! used by RRTGMP gas optics 
         rayl_upper,                      & ! used by RRTGMP gas optics 
         lut_extice,                      & ! used by RRTGMP cloud optics 
         lut_ssaice,                      & ! used by RRTGMP cloud optics 
         lut_asyice,                      & ! used by RRTGMP cloud optics 
         pade_extliq,                     & ! used by RRTGMP cloud optics 
         pade_ssaliq,                     & ! used by RRTGMP cloud optics 
         pade_asyliq                        ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         kmajor,                          & ! used by RRTGMP gas optics 
         planck_frac,                     & ! used by RRTGMP gas optics 
         pade_extice,                     & ! used by RRTGMP cloud optics 
         pade_ssaice,                     & ! used by RRTGMP cloud optics 
         pade_asyice                        ! used by RRTGMP cloud optics 
    character(len=32),  dimension(:), allocatable :: &
         gas_names,                       & ! used by RRTGMP gas optics 
         gas_minor,                       & ! used by RRTGMP gas optics 
         identifier_minor,                & ! used by RRTGMP gas optics 
         minor_gases_lower,               & ! used by RRTGMP gas optics 
         minor_gases_upper,               & ! used by RRTGMP gas optics 
         scaling_gas_lower,               & ! used by RRTGMP gas optics 
         scaling_gas_upper                  ! used by RRTGMP gas optics 
    logical(wl), dimension(:), allocatable :: &
         minor_scales_with_density_lower, & ! used by RRTGMP gas optics 
         minor_scales_with_density_upper, & ! used by RRTGMP gas optics 
         scale_by_complement_lower,       & ! used by RRTGMP gas optics 
         scale_by_complement_upper          ! used by RRTGMP gas optics 

    ! Dimensions (to be broadcast across all processors)
    integer :: &
         ntemps,                          & ! used by RRTGMP gas optics 
         npress,                          & ! used by RRTGMP gas optics 
         nabsorbers,                      & ! used by RRTGMP gas optics 
         nextrabsorbers,                  & ! used by RRTGMP gas optics 
         nminorabsorbers,                 & ! used by RRTGMP gas optics 
         nmixingfracs,                    & ! used by RRTGMP gas optics 
         nlayers,                         & ! used by RRTGMP gas optics 
         nbnds,                           & ! used by RRTGMP gas optics 
         ngpts,                           & ! used by RRTGMP gas optics 
         npairs,                          & ! used by RRTGMP gas optics 
         ninternalSourcetemps,            & ! used by RRTGMP gas optics 
         nminor_absorber_intervals_lower, & ! used by RRTGMP gas optics 
         nminor_absorber_intervals_upper, & ! used by RRTGMP gas optics 
         ncontributors_lower,             & ! used by RRTGMP gas optics 
         ncontributors_upper,             & ! used by RRTGMP gas optics 
         nbandLWcldy,                     & ! used by RRTGMP cloud optics 
         nsize_liq,                       & ! used by RRTGMP cloud optics 
         nsize_ice,                       & ! used by RRTGMP cloud optics 
         nsizereg,                        & ! used by RRTGMP cloud optics 
         ncoeff_ext,                      & ! used by RRTGMP cloud optics 
         ncoeff_ssa_g,                    & ! used by RRTGMP cloud optics 
         nbound,                          & ! used by RRTGMP cloud optics  
         npairsLWcldy                       ! used by RRTGMP cloud optics 

    ! Local variables
    integer :: ncid_lw,dimID,varID,status,igpt,iGas,ij,ierr,ncid_lw_clds
    integer,dimension(:),allocatable :: temp1,temp2,temp3,temp4,temp_log_array1,&
    	temp_log_array2, temp_log_array3, temp_log_array4
    character(len=264) :: kdist_file,kdist_cldy_file
    integer,parameter :: max_strlen=256

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Ensure that requested cloud overlap is reasonable.
    if ( iovrlw .lt. 0 .or. iovrlw .gt. 3 ) then
       print *,'  *** Error in specification of cloud overlap flag',   &
            ' IOVRLW=',iovrlw,' in RLWINIT !!'
       stop
    elseif ( iovrlw .ge. 2 .and. isubclw .eq. 0 ) then 
       print *,'  *** IOVRLW=',iovrlw,' is not available for',         &
            ' ISUBCLW=0 setting!!'
       print *,'      The program uses maximum/random overlap',        &
            ' instead.'
       iovrlw = 1
    endif

    ! Check cloud flags for consistency.
    if ((icldflg .eq. 0 .and. ilwcliq .ne. 0) .or. &
        (icldflg .eq. 1 .and. ilwcliq .eq. 0)) then
       print *,'  *** Model cloud scheme inconsistent with LW',        &
            ' radiation cloud radiative property setup !!'
       stop
    endif

    ! How are we handling cloud-optics?
    rrtmgp_lw_cld_phys = Model%rrtmgp_cld_phys

    ! Filenames are set in the gfs_physics_nml (scm/src/GFS_typedefs.F90)
    kdist_file      = trim(Model%rrtmgp_root)//trim(Model%kdist_lw_file_gas)
    kdist_cldy_file = trim(Model%rrtmgp_root)//trim(Model%kdist_lw_file_clouds)

    ! Read dimensions for k-distribution fields (only on master processor(0))
    if (mpirank .eq. mpiroot) then
       if(nf90_open(trim(kdist_file), NF90_WRITE, ncid_lw) .eq. NF90_NOERR) then
          status = nf90_inq_dimid(ncid_lw, 'temperature', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=ntemps)
          status = nf90_inq_dimid(ncid_lw, 'pressure', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=npress)
          status = nf90_inq_dimid(ncid_lw, 'absorber', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=nabsorbers)
          status = nf90_inq_dimid(ncid_lw, 'minor_absorber', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=nminorabsorbers)
          status = nf90_inq_dimid(ncid_lw, 'absorber_ext', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=nextrabsorbers)
          status = nf90_inq_dimid(ncid_lw, 'mixing_fraction', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=nmixingfracs)
          status = nf90_inq_dimid(ncid_lw, 'atmos_layer', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=nlayers)
          status = nf90_inq_dimid(ncid_lw, 'bnd', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=nbnds)
          status = nf90_inq_dimid(ncid_lw, 'gpt', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=ngpts)
          status = nf90_inq_dimid(ncid_lw, 'pair', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=npairs)
          status = nf90_inq_dimid(ncid_lw, 'contributors_lower', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=ncontributors_lower)
          status = nf90_inq_dimid(ncid_lw, 'contributors_upper', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=ncontributors_upper)
          status = nf90_inq_dimid(ncid_lw, 'minor_absorber_intervals_lower', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=nminor_absorber_intervals_lower)
          status = nf90_inq_dimid(ncid_lw, 'minor_absorber_intervals_upper', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=nminor_absorber_intervals_upper)
          status = nf90_inq_dimid(ncid_lw, 'temperature_Planck', dimid)
          status = nf90_inquire_dimension(ncid_lw, dimid, len=ninternalSourcetemps)
          status = nf90_close(ncid_lw)
       endif
    endif
    
    ! Broadcast dimensions to all processors
#ifdef MPI
    call MPI_BCAST(ntemps,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(npress,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nabsorbers,                      1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nminorabsorbers,                 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nextraabsorbers,                 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nmixingfracs,                    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nlayers,                         1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nbnds,                           1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(ngpts,                           1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(npairs,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(ncontributors_lower,             1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(ncontributors_upper,             1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nminor_absorber_intervals_lower, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(nminor_absorber_intervals_upper, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(ninternalSourcetemps,            1, MPI_INTEGER, mpiroot, mpicomm, ierr)
#endif
    
    !if (mpirank .eq. mpiroot) then
    ! Allocate space for arrays
    allocate(gas_names(nabsorbers))
    allocate(scaling_gas_lower(nminor_absorber_intervals_lower))
    allocate(scaling_gas_upper(nminor_absorber_intervals_upper))
    allocate(gas_minor(nminorabsorbers))
    allocate(identifier_minor(nminorabsorbers))
    allocate(minor_gases_lower(nminor_absorber_intervals_lower))
    allocate(minor_gases_upper(nminor_absorber_intervals_upper))
    allocate(minor_limits_gpt_lower(npairs,nminor_absorber_intervals_lower))
    allocate(minor_limits_gpt_upper(npairs,nminor_absorber_intervals_upper))
    allocate(band2gpt(2,nbnds))
    allocate(key_species(2,nlayers,nbnds))
    allocate(band_lims(2,nbnds))
    allocate(press_ref(npress))
    allocate(temp_ref(ntemps))
    allocate(vmr_ref(nlayers, nextrabsorbers, ntemps))
    allocate(kminor_lower(ncontributors_lower, nmixingfracs, ntemps))
    allocate(kmajor(ngpts, nmixingfracs,  npress+1, ntemps))
    allocate(kminor_start_lower(nminor_absorber_intervals_lower))
    allocate(kminor_upper(ncontributors_upper, nmixingfracs, ntemps))
    allocate(kminor_start_upper(nminor_absorber_intervals_upper))
    allocate(minor_scales_with_density_lower(nminor_absorber_intervals_lower))
    allocate(minor_scales_with_density_upper(nminor_absorber_intervals_upper))
    allocate(scale_by_complement_lower(nminor_absorber_intervals_lower))
    allocate(scale_by_complement_upper(nminor_absorber_intervals_upper))
    allocate(temp1(nminor_absorber_intervals_lower))
    allocate(temp2(nminor_absorber_intervals_upper))
    allocate(temp3(nminor_absorber_intervals_lower))
    allocate(temp4(nminor_absorber_intervals_upper))
    allocate(totplnk(ninternalSourcetemps, nbnds))
    allocate(planck_frac(ngpts, nmixingfracs, npress+1, ntemps))

    if (mpirank .eq. mpiroot) then
       ! Read in fields from file
       if(nf90_open(trim(kdist_file), NF90_WRITE, ncid_lw) .eq. NF90_NOERR) then
          status = nf90_inq_varid(ncid_lw,'gas_names',varID)
          status = nf90_get_var(ncid_lw,varID,gas_names)
          !
          status = nf90_inq_varid(ncid_lw,'scaling_gas_lower',varID)
          status = nf90_get_var(ncid_lw,varID,scaling_gas_lower)
          !
          status = nf90_inq_varid(ncid_lw,'scaling_gas_upper',varID)
          status = nf90_get_var(ncid_lw,varID,scaling_gas_upper)
          !
          status = nf90_inq_varid(ncid_lw,'gas_minor',varID)
          status = nf90_get_var(ncid_lw,varID,gas_minor)
          !
          status = nf90_inq_varid(ncid_lw,'identifier_minor',varID)
          status = nf90_get_var(ncid_lw,varID,identifier_minor)
          !
          status = nf90_inq_varid(ncid_lw,'minor_gases_lower',varID)
          status = nf90_get_var(ncid_lw,varID,minor_gases_lower)
          !
          status = nf90_inq_varid(ncid_lw,'minor_gases_upper',varID)
          status = nf90_get_var(ncid_lw,varID,minor_gases_upper)
          !
          status = nf90_inq_varid(ncid_lw,'minor_limits_gpt_lower',varID)
          status = nf90_get_var(ncid_lw,varID,minor_limits_gpt_lower)
          !
          status = nf90_inq_varid(ncid_lw,'minor_limits_gpt_upper',varID)
          status = nf90_get_var(ncid_lw,varID,minor_limits_gpt_upper)
          !
          status = nf90_inq_varid(ncid_lw,'bnd_limits_gpt',varID)
          status = nf90_get_var(ncid_lw,varID,band2gpt)
          !
          status = nf90_inq_varid(ncid_lw,'key_species',varID)
          status = nf90_get_var(ncid_lw,varID,key_species)
          !
          status = nf90_inq_varid(ncid_lw,'bnd_limits_wavenumber',varID)
          status = nf90_get_var(ncid_lw,varID,band_lims)
          !
          status = nf90_inq_varid(ncid_lw,'press_ref',varID)
          status = nf90_get_var(ncid_lw,varID,press_ref)
          !
          status = nf90_inq_varid(ncid_lw,'temp_ref',varID)
          status = nf90_get_var(ncid_lw,varID,temp_ref)
          !
          status = nf90_inq_varid(ncid_lw,'absorption_coefficient_ref_P',varID)
          status = nf90_get_var(ncid_lw,varID,temp_ref_p)
          !
          status = nf90_inq_varid(ncid_lw,'absorption_coefficient_ref_T',varID)
          status = nf90_get_var(ncid_lw,varID,temp_ref_t)
          !
          status = nf90_inq_varid(ncid_lw,'press_ref_trop',varID)
          status = nf90_get_var(ncid_lw,varID,press_ref_trop)
          !
          status = nf90_inq_varid(ncid_lw,'kminor_lower',varID)
          status = nf90_get_var(ncid_lw,varID,kminor_lower)
          !
          status = nf90_inq_varid(ncid_lw,'kminor_upper',varID)
          status = nf90_get_var(ncid_lw,varID,kminor_upper)
          !
          status = nf90_inq_varid(ncid_lw,'vmr_ref',varID)
          status = nf90_get_var(ncid_lw,varID,vmr_ref)
          !
          status = nf90_inq_varid(ncid_lw,'kmajor',varID)
          status = nf90_get_var(ncid_lw,varID,kmajor)
          !
          status = nf90_inq_varid(ncid_lw,'kminor_start_lower',varID)
          status = nf90_get_var(ncid_lw,varID,kminor_start_lower)
          !
          status = nf90_inq_varid(ncid_lw,'kminor_start_upper',varID)
          status = nf90_get_var(ncid_lw,varID,kminor_start_upper)
          !
          status = nf90_inq_varid(ncid_lw,'totplnk',varID)
          status = nf90_get_var(ncid_lw,varID,totplnk)
          !
          status = nf90_inq_varid(ncid_lw,'plank_fraction',varID)
          status = nf90_get_var(ncid_lw,varID,planck_frac)
          
          ! Logical fields are read in as integers and then converted to logicals.
          status = nf90_inq_varid(ncid_lw,'minor_scales_with_density_lower',varID)
          status = nf90_get_var(ncid_lw,varID,temp1)
          minor_scales_with_density_lower(:) = .false.
          where(temp1 .eq. 1) minor_scales_with_density_lower(:) = .true.
          !
          status = nf90_inq_varid(ncid_lw,'minor_scales_with_density_upper',varID)
          status = nf90_get_var(ncid_lw,varID,temp2)
          minor_scales_with_density_upper(:) = .false.
          where(temp2 .eq. 1) minor_scales_with_density_upper(:) = .true.
          !
          status = nf90_inq_varid(ncid_lw,'scale_by_complement_lower',varID)
          status = nf90_get_var(ncid_lw,varID,temp3)
          scale_by_complement_lower(:) = .false.
          where(temp3 .eq. 1) scale_by_complement_lower(:) = .true.
          !
          status = nf90_inq_varid(ncid_lw,'scale_by_complement_upper',varID)
          status = nf90_get_var(ncid_lw,varID,temp4)
          scale_by_complement_upper(:) = .false.
          where(temp4 .eq. 1) scale_by_complement_upper(:) = .true.
          
          ! Close
          status = nf90_close(ncid_lw)
       endif
    endif

    ! Broadcast arrays to all processors
#ifdef MPI
    call MPI_BCAST(minor_limits_gpt_upper, size(minor_limits_gpt_upper), MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(minor_limits_gpt_lower, size(minor_limits_gpt_lower), MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(kminor_start_upper,     size(kminor_start_upper),     MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(kminor_start_lower,     size(kminor_start_lower),     MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(key_species,            size(key_species),            MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(band2gpt,               size(band2gpt),               MPI_INTEGER, mpiroot, mpicomm, ierr)
    call MPI_BCAST(band_lims,              size(band_lims),              kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(press_ref,              size(press_ref),              kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(temp_ref,               size(temp_ref),               kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(kminor_lower,           size(kminor_lower),           kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(kminor_upper,           size(kminor_upper),           kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(scaling_gas_lower,      size(scaling_gas_lower),      kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(scaling_gas_upper,      size(scaling_gas_upper),      kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(vmr_ref,                size(vmr_ref),                kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(kmajor,                 size(kmajor),                 kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(temp_ref_p,             1,                            kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(temp_ref_t,             1,                            kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(press_ref_trop,         1,                            kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(totplnk,                size(totplnk),                kind_phys,   mpiroot, mpicomm, ierr)
    call MPI_BCAST(planck_frac,            size(planck_frac),            kind_phys,   mpiroot, mpicomm, ierr)
    ! Character arrays
    do ij=1,nabsorbers
       call MPI_BCAST(gas_names(ij),         32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
    enddo
    do ij=1,nminorabsorbers
       call MPI_BCAST(gas_minor(ij),         32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(identifier_minor(ij),  32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
    enddo
    do ij=1,nminor_absorber_intervals_lower
       call MPI_BCAST(minor_gases_lower(ij), 32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
    enddo
    do ij=1,nminor_absorber_intervals_upper
       call MPI_BCAST(minor_gases_upper(ij), 32,  MPI_CHAR,   mpiroot, mpicomm, ierr)
    enddo
    ! Logical arrays (First convert to integer-array, then broadcast)
    !
    allocate(temp_log_array1(nminor_absorber_intervals_lower))
    where(minor_scales_with_density_lower)
       temp_log_array1 = 1
    elsewhere
       temp_log_array1 = 0
    end where
    call MPI_BCAST(temp_log_array1, size(temp_log_array1), MPI_INTEGER,  mpiroot, mpicomm, ierr)
    !
    allocate(temp_log_array2(nminor_absorber_intervals_lower))
    where(scale_by_complement_lower)
       temp_log_array2 = 1
    elsewhere
       temp_log_array2 = 0
    end where
    call MPI_BCAST(temp_log_array2, size(temp_log_array2), MPI_INTEGER,  mpiroot, mpicomm, ierr)
    !
    allocate(temp_log_array3(nminor_absorber_intervals_upper))
    where(minor_scales_with_density_upper)
       temp_log_array3 = 1
    elsewhere
       temp_log_array3 = 0
    end where
    call MPI_BCAST(temp_log_array3, size(temp_log_array3), MPI_INTEGER,  mpiroot, mpicomm, ierr)
    !
    allocate(temp_log_array4(nminor_absorber_intervals_upper))
    where(scale_by_complement_upper)
       temp_log_array4 = 1
    elsewhere
       temp_log_array4 = 0
    end where
    call MPI_BCAST(temp_log_array4, size(temp_log_array4), MPI_INTEGER,  mpiroot, mpicomm, ierr)
#endif

    ! Initialize gas concentrations and gas optics class with data
    do iGas=1,nGases
       call check_error_msg(gas_concentrations%set_vmr(active_gases(iGas), 0._kind_phys))
    enddo    
    call check_error_msg(kdist_lw%load(gas_concentrations, gas_names, key_species, band2gpt, &
         band_lims, press_ref, press_ref_trop, temp_ref,  temp_ref_p, temp_ref_t,   &
         vmr_ref, kmajor, kminor_lower, kminor_upper, gas_minor,identifier_minor,   &
         minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower,              &
         minor_limits_gpt_upper, minor_scales_with_density_lower,                   &
         minor_scales_with_density_upper, scaling_gas_lower,                        &
         scaling_gas_upper, scale_by_complement_lower,                              &
         scale_by_complement_upper, kminor_start_lower, kminor_start_upper,         &
         totplnk, planck_frac, rayl_lower, rayl_upper))

    ! Set initial permutation seed for McICA, initially set to number of G-points
    ipsdlw0 = kdist_lw%get_ngpt()

    ! #######################################################################################
    ! If RRTMGP cloud-optics are requested, read tables and broadcast.
    ! #######################################################################################
    ! Read dimensions for k-distribution fields (only on master processor(0))
    if (mpirank .eq. mpiroot) then
       if(nf90_open(trim(kdist_cldy_file), NF90_WRITE, ncid_lw_clds) == NF90_NOERR) then
          status = nf90_inq_dimid(ncid_lw_clds, 'nband', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=nbandLWcldy)
          status = nf90_inq_dimid(ncid_lw_clds, 'nrghice', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=nrghice)
          status = nf90_inq_dimid(ncid_lw_clds, 'nsize_liq', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=nsize_liq)
          status = nf90_inq_dimid(ncid_lw_clds, 'nsize_ice', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=nsize_ice)
          status = nf90_inq_dimid(ncid_lw_clds, 'nsizereg', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=nsizereg)
          status = nf90_inq_dimid(ncid_lw_clds, 'ncoeff_ext', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=ncoeff_ext)
          status = nf90_inq_dimid(ncid_lw_clds, 'ncoeff_ssa_g', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=ncoeff_ssa_g)
          status = nf90_inq_dimid(ncid_lw_clds, 'nbound', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=nbound)
          status = nf90_inq_dimid(ncid_lw_clds, 'pair', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=npairsLWcldy)
          status = nf90_close(ncid_lw_clds)
       endif
    endif

    ! Broadcast dimensions to all processors
#ifdef MPI
    if (rrtmgp_lw_cld_phys .eq. 1 .or. rrtmgp_lw_cld_phys .eq. 2) then
       call MPI_BCAST(nbandLWcldy,  1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nrghice,      1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_liq,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_ice,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsizereg,     1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ext,   1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ssa_g, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nbound,       1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(npairsLWcldy, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    endif
#endif

    if (rrtmgp_lw_cld_phys .eq. 1) then
       allocate(lut_extliq(nsize_liq, nBandLWcldy))
       allocate(lut_ssaliq(nsize_liq, nBandLWcldy))
       allocate(lut_asyliq(nsize_liq, nBandLWcldy))
       allocate(lut_extice(nsize_ice, nBandLWcldy, nrghice))
       allocate(lut_ssaice(nsize_ice, nBandLWcldy, nrghice))
       allocate(lut_asyice(nsize_ice, nBandLWcldy, nrghice))
       allocate(band_lims_cldy(2, nBandLWcldy))
    endif
    if (rrtmgp_lw_cld_phys .eq. 2) then
       allocate(pade_extliq(nbandLWcldy, nsizereg,  ncoeff_ext ))
       allocate(pade_ssaliq(nbandLWcldy, nsizereg,  ncoeff_ssa_g))
       allocate(pade_asyliq(nbandLWcldy, nsizereg,  ncoeff_ssa_g))
       allocate(pade_extice(nbandLWcldy, nsizereg,  ncoeff_ext,   nrghice))
       allocate(pade_ssaice(nbandLWcldy, nsizereg,  ncoeff_ssa_g, nrghice))
       allocate(pade_asyice(nbandLWcldy, nsizereg,  ncoeff_ssa_g, nrghice))
       allocate(pade_sizereg_extliq(nbound))
       allocate(pade_sizereg_ssaliq(nbound))
       allocate(pade_sizereg_asyliq(nbound))
       allocate(pade_sizereg_extice(nbound))
       allocate(pade_sizereg_ssaice(nbound))
       allocate(pade_sizereg_asyice(nbound))
       allocate(band_lims_cldy(2,nbandLWcldy))
    endif

    ! On master processor, allocate space, read in fields, broadcast to all processors
    if (mpirank .eq. mpiroot) then
       ! 
       if (rrtmgp_lw_cld_phys .eq. 1) then
          !
          if(nf90_open(trim(kdist_cldy_file), NF90_WRITE, ncid_lw_clds) == NF90_NOERR) then
             status = nf90_inq_varid(ncid_lw_clds,'radliq_lwr',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radliq_lwr)
             status = nf90_inq_varid(ncid_lw_clds,'radliq_upr',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radliq_upr)
             status = nf90_inq_varid(ncid_lw_clds,'radliq_fac',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radliq_fac)
             status = nf90_inq_varid(ncid_lw_clds,'radice_lwr',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radice_lwr)
             status = nf90_inq_varid(ncid_lw_clds,'radice_upr',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radice_upr)
             status = nf90_inq_varid(ncid_lw_clds,'radice_fac',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radice_fac)
             status = nf90_inq_varid(ncid_lw_clds,'lut_extliq',varID)
             status = nf90_get_var(ncid_lw_clds,varID,lut_extliq)
             status = nf90_inq_varid(ncid_lw_clds,'lut_ssaliq',varID)
             status = nf90_get_var(ncid_lw_clds,varID,lut_ssaliq)
             status = nf90_inq_varid(ncid_lw_clds,'lut_asyliq',varID)
             status = nf90_get_var(ncid_lw_clds,varID,lut_asyliq)
             status = nf90_inq_varid(ncid_lw_clds,'lut_extice',varID)
             status = nf90_get_var(ncid_lw_clds,varID,lut_extice)
             status = nf90_inq_varid(ncid_lw_clds,'lut_ssaice',varID)
             status = nf90_get_var(ncid_lw_clds,varID,lut_ssaice)
             status = nf90_inq_varid(ncid_lw_clds,'lut_asyice',varID)
             status = nf90_get_var(ncid_lw_clds,varID,lut_asyice)
             status = nf90_inq_varid(ncid_lw_clds,'bnd_limits_wavenumber',varID)
             status = nf90_get_var(ncid_lw_clds,varID,band_lims_cldy)
             status = nf90_close(ncid_lw_clds)
          endif
       endif
       !
       if (rrtmgp_lw_cld_phys .eq. 2) then
          !
          if(nf90_open(trim(kdist_cldy_file), NF90_WRITE, ncid_lw_clds) == NF90_NOERR) then
             status = nf90_inq_varid(ncid_lw_clds,'radliq_lwr',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radliq_lwr)
             status = nf90_inq_varid(ncid_lw_clds,'radliq_upr',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radliq_upr)
             status = nf90_inq_varid(ncid_lw_clds,'radliq_fac',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radliq_fac)
             status = nf90_inq_varid(ncid_lw_clds,'radice_lwr',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radice_lwr)
             status = nf90_inq_varid(ncid_lw_clds,'radice_upr',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radice_upr)
             status = nf90_inq_varid(ncid_lw_clds,'radice_fac',varID)
             status = nf90_get_var(ncid_lw_clds,varID,radice_fac)
             status = nf90_inq_varid(ncid_lw_clds,'pade_extliq',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_extliq)
             status = nf90_inq_varid(ncid_lw_clds,'pade_ssaliq',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_ssaliq)
             status = nf90_inq_varid(ncid_lw_clds,'pade_asyliq',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_asyliq)
             status = nf90_inq_varid(ncid_lw_clds,'pade_extice',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_extice)
             status = nf90_inq_varid(ncid_lw_clds,'pade_ssaice',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_ssaice)
             status = nf90_inq_varid(ncid_lw_clds,'pade_asyice',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_asyice)
             status = nf90_inq_varid(ncid_lw_clds,'pade_sizreg_extliq',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_sizereg_extliq)
             status = nf90_inq_varid(ncid_lw_clds,'pade_sizreg_ssaliq',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_sizereg_ssaliq)
             status = nf90_inq_varid(ncid_lw_clds,'pade_sizreg_asyliq',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_sizereg_asyliq)
             status = nf90_inq_varid(ncid_lw_clds,'pade_sizreg_extice',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_sizereg_extice)
             status = nf90_inq_varid(ncid_lw_clds,'pade_sizreg_ssaice',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_sizereg_ssaice)
             status = nf90_inq_varid(ncid_lw_clds,'pade_sizreg_asyice',varID)
             status = nf90_get_var(ncid_lw_clds,varID,pade_sizereg_asyice)
             status = nf90_inq_varid(ncid_lw_clds,'bnd_limits_wavenumber',varID)
             status = nf90_get_var(ncid_lw_clds,varID,band_lims_cldy)
             status = nf90_close(ncid_lw_clds)
          endif
       endif
    endif

    ! Broadcast arrays to all processors
#ifdef MPI
    if (rrtmgp_lw_cld_phys .eq. 1) then
       call MPI_BCAST(radliq_lwr,           size(radliq_lwr),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radliq_upr,           size(radliq_upr),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radliq_fac,           size(radliq_fac),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_lwr,           size(radice_lwr),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_upr,           size(radice_upr),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_fac,           size(radice_fac),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_extliq,           size(lut_extliq),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_ssaliq,           size(lut_ssaliq),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_asyliq,           size(lut_asyliq),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_extice,           size(lut_extice),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_ssaice,           size(lut_ssaice),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_asyice,           size(lut_asyice),          kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(band_lims_cldy),      size(band_lims_cldy),      kind_phys,   mpiroot, mpicomm, ierr)    
    endif
    if (rrtmgp_lw_cld_phys .eq. 2) then
       call MPI_BCAST(pade_extliq,          size(pade_extliq),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_ssaliq,          size(pade_ssaliq),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_asyliq,          size(pade_asyliq),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_extice,          size(pade_extice),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_ssaice,          size(pade_ssaice),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_asyice,          size(pade_asyice),         kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_extliq), size(pade_sizereg_extliq), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_ssaliq), size(pade_sizereg_ssaliq), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_asyliq), size(pade_sizereg_asyliq), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_extice), size(pade_sizereg_extice), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_ssaice), size(pade_sizereg_ssaice), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_asyice), size(pade_sizereg_asyice), kind_phys,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(band_lims_cldy),      size(band_lims_cldy),      kind_phys,   mpiroot, mpicomm, ierr)    
    endif
#endif

    ! Load tables data for RRTGMP cloud-optics  
    if (rrtmgp_lw_cld_phys .eq. 1) then
       call check_error_msg(kdist_cldy_lw%set_ice_roughness(nrghice))
       call check_error_msg(kdist_cldy_lw%load(band_lims_cldy, radliq_lwr, radliq_upr, &
            radliq_fac, radice_lwr, radice_upr, radice_fac, lut_extliq, lut_ssaliq,    &
            lut_asyliq, lut_extice, lut_ssaice, lut_asyice))
    endif
    if (rrtmgp_lw_cld_phys .eq. 2) then
       call check_error_msg(kdist_cldy_lw%set_ice_roughness(nrghice))
       call check_error_msg(kdist_cldy_lw%load(band_lims_cldy, pade_extliq,            &
            pade_ssaliq, pade_asyliq, pade_extice, pade_ssaice, pade_asyice,           &
            pade_sizereg_extliq, pade_sizereg_ssaliq, pade_sizereg_asyliq,             &
            pade_sizereg_extice, pade_sizereg_ssaice, pade_sizereg_asyice))
    endif

  end subroutine rrtmgp_lw_pre_init

!> \section arg_table_rrtmgp_lw_pre_run Argument Table
!! | local_name       | standard_name                             | long_name                                                          | units    | rank |  type                 |   kind    | intent | optional |
!! |------------------|-------------------------------------------|--------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model            | GFS_control_type_instance                 | Fortran DDT containing FV3-GFS model control parameters            | DDT      |    0 | GFS_control_type      |           | in     | F        |
!! | Grid             | GFS_grid_type_instance                    | Fortran DDT containing FV3-GFS grid and interpolation related data | DDT      |    0 | GFS_grid_type         |           | in     | F        |
!! | Sfcprop          | GFS_sfcprop_type_instance                 | Fortran DDT containing FV3-GFS surface fields                      | DDT      |    0 | GFS_sfcprop_type      |           | in     | F        |
!! | Radtend          | GFS_radtend_type_instance                 | Fortran DDT containing FV3-GFS radiation tendencies                | DDT      |    0 | GFS_radtend_type      |           | inout  | F        |
!! | im               | horizontal_loop_extent                    | horizontal loop extent                                             | count    |    0 | integer               |           | in     | F        |
!! | tsfg             | surface_ground_temperature_for_radiation  | surface ground temperature for radiation                           | K        |    1 | real                  | kind_phys | in     | F        |
!! | tsfa             | surface_air_temperature_for_radiation     | lowest model layer air temperature for radiation                   | K        |    1 | real                  | kind_phys | in     | F        |
!! | errmsg           | ccpp_error_message                        | error message for error handling in CCPP                           | none     |    0 | character             | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                           | error flag for error handling in CCPP                              | flag     |    0 | integer               |           | out    | F        |
!! | kdist_lw         | K_distribution_file_for_RRTMGP_LW_scheme  | DDT containing spectral information for RRTMGP LW radiation scheme | DDT      |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | sfc_emiss_byband | surface_longwave_emissivity_in_each_band  | surface lw emissivity in fraction in each LW band                  | frac     |    2 | real                  | kind_phys | inout  | F        |
!!
      subroutine rrtmgp_lw_pre_run (Model, Grid, Sfcprop, Radtend, im, tsfg, tsfa, kdist_lw, sfc_emiss_byband, errmsg, errflg)
    
      use machine,                   only: kind_phys

      use GFS_typedefs,              only: GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type,           &
                                           GFS_sfcprop_type         
      use module_radiation_surface,  only: setemis
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp

      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_sfcprop_type),         intent(in)    :: Sfcprop
      type(GFS_grid_type),            intent(in)    :: Grid
      integer, intent(in)                           :: im
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) ::  tsfa, tsfg
      type(ty_gas_optics_rrtmgp),intent(in) :: &
           kdist_lw        ! DDT containing LW spectral information
      real(kind_phys),dimension(kdist_lw%get_nband(),im),intent(inout) :: sfc_emiss_byband
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      integer :: ij

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      print*,'In RRTMGP_lW_PRE_RUN(): top',shape(sfc_emiss_byband),im
      print*,'In RRTMGP_lW_PRE_RUN(): top',shape(Radtend%semis)
      if (Model%lslwr) then
!>  - Call module_radiation_surface::setemis(),to setup surface
!! emissivity for LW radiation.
        call setemis (Grid%xlon, Grid%xlat, Sfcprop%slmsk,        &        !  ---  inputs
                     Sfcprop%snowd, Sfcprop%sncovr, Sfcprop%zorl, &
                     tsfg, tsfa, Sfcprop%hprim, IM,               &
                      Radtend%semis)                              !  ---  outputs
        do ij=1,kdist_lw%get_nband()
           print*,ij
           sfc_emiss_byband(ij,1:im) = Radtend%semis(1:im)
        enddo
      endif
      print*,'In RRTMGP_lW_PRE_RUN(): bottom'

      end subroutine rrtmgp_lw_pre_run

!> \section arg_table_rrtmgp_lw_pre_finalize Argument Table
!!
       subroutine rrtmgp_lw_pre_finalize ()
       end subroutine rrtmgp_lw_pre_finalize
!! @}
       end module rrtmgp_lw_pre
