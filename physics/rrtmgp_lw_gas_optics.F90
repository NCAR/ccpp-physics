module rrtmgp_lw_gas_optics
  use machine,               only: kind_phys
  use GFS_typedefs,          only: GFS_control_type, GFS_radtend_type
  use mo_rte_kind,           only: wl
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use netcdf

  ! Parameters
  integer :: ipsdlw0

contains
!! \section arg_table_rrtmgp_lw_gas_optics_init Argument Table
!! | local_name   | standard_name                                | long_name                                                          | units | rank | type                 | kind  | intent | optional |
!! |--------------|----------------------------------------------|--------------------------------------------------------------------|-------|------|----------------------|-------|--------|----------|
!! | Model        | GFS_control_type_instance                    | Fortran DDT containing FV3-GFS model control parameters            | DDT   |    0 | GFS_control_type     |       | in     | F        |
!! | Radtend      | GFS_radtend_type_instance                    | Fortran DDT containing FV3-GFS radiation tendencies                | DDT   |    0 | GFS_radtend_type     |       | in     | F        |
!! | mpirank      | mpi_rank                                     | current MPI rank                                                   | index |    0 | integer              |       | in     | F        |
!! | mpiroot      | mpi_root                                     | master MPI rank                                                    | index |    0 | integer              |       | in     | F        |
!! | mpicomm      | mpi_comm                                     | MPI communicator                                                   | index |    0 | integer              |       | in     | F        |
!! | errmsg       | ccpp_error_message                           | error message for error handling in CCPP                           | none  |    0 | character            | len=* | out    | F        |
!! | errflg       | ccpp_error_flag                              | error flag for error handling in CCPP                              | flag  |    0 | integer              |       | out    | F        |
!! | lw_gas_props | coefficients_for_lw_gas_optics               | DDT containing spectral information for RRTMGP LW radiation scheme | DDT   |    0 | ty_gas_optics_rrtmgp |       | out    | F        |
!! | ngpts_lw     | number_of_spectral_points_for_LW_calculation | Number of spectral points for LW RRTMGP calculation                | none  |    0 | integer              |       | out    | F        |
!!
  ! #########################################################################################
  ! #########################################################################################
  subroutine rrtmgp_lw_gas_optics_init(Model, Radtend, mpicomm, mpirank, mpiroot, lw_gas_props,      &
       ngpts_lw, errmsg, errflg)
    use netcdf
    
#ifdef MPI
    use mpi
#endif

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model        ! DDT containing model control parameters
    type(GFS_radtend_type), intent(in) :: &
         Radtend      ! DDT containing FV3-GFS radiation tendencies
    integer,intent(in) :: &
         mpicomm,   & ! MPI communicator
         mpirank,   & ! Current MPI rank
         mpiroot      ! Master MPI rank
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg       ! Error message
    integer,          intent(out) :: &
         errflg       ! Error code
    integer, intent(out) :: &
         ngpts_lw     ! Number of g-points
    type(ty_gas_optics_rrtmgp),intent(out) :: &
         lw_gas_props ! DDT containing spectral information for RRTMGP LW radiation scheme

    ! Variables that will be passed to gas_optics%load()
    type(ty_gas_concs) :: &
         gas_concentrations
    integer, dimension(:), allocatable :: &
         kminor_start_lower,              & ! used by RRTMGP gas optics 
         kminor_start_upper                 ! used by RRTMGP gas optics 
    integer, dimension(:,:), allocatable :: &
         band2gpt,                        & ! used by RRTMGP gas optics 
         minor_limits_gpt_lower,          & ! used by RRTMGP gas optics 
         minor_limits_gpt_upper             ! used by RRTMGP gas optics 
    integer, dimension(:,:,:), allocatable :: &
         key_species                        ! used by RRTMGP gas optics 
    real(kind_phys) :: &
         press_ref_trop,                  & ! used by RRTMGP gas optics 
         temp_ref_p,                      & ! used by RRTMGP gas optics 
         temp_ref_t                         ! used by RRTMGP gas optics 
    real(kind_phys), dimension(:), allocatable :: &
         press_ref,                       & ! used by RRTMGP gas optics 
         temp_ref                           ! used by RRTMGP gas optics 
    real(kind_phys), dimension(:,:), allocatable :: &
         band_lims,                       & ! used by RRTMGP gas optics 
         totplnk                            ! used by RRTMGP gas optics 
    real(kind_phys), dimension(:,:,:), allocatable :: &
         vmr_ref,                         & ! used by RRTMGP gas optics 
         kminor_lower,                    & ! used by RRTMGP gas optics 
         kminor_upper,                    & ! used by RRTMGP gas optics 
         rayl_lower,                      & ! used by RRTMGP gas optics 
         rayl_upper                         ! used by RRTMGP gas optics 
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         kmajor,                          & ! used by RRTMGP gas optics 
         planck_frac                        ! used by RRTMGP gas optics  
    character(len=32),  dimension(:), allocatable :: &
         gas_names,                       & ! used by RRTMGP gas optics 
         gas_minor,                       & ! used by RRTMGP gas optics 
         identifier_minor,                & ! used by RRTMGP gas optics 
         minor_gases_lower,               & ! used by RRTMGP gas optics 
         minor_gases_upper,               & ! used by RRTMGP gas optics 
         scaling_gas_lower,               & ! used by RRTMGP gas optics 
         scaling_gas_upper                  ! used by RRTMGP gas optics 
    logical(wl), dimension(:), allocatable :: &
         minor_scales_with_density_lower, & ! used by RRTMGP gas optics 
         minor_scales_with_density_upper, & ! used by RRTMGP gas optics 
         scale_by_complement_lower,       & ! used by RRTMGP gas optics 
         scale_by_complement_upper          ! used by RRTMGP gas optics 

    ! Dimensions (to be broadcast across all processors)
    integer :: &
         ntemps,                          & ! used by RRTMGP gas optics 
         npress,                          & ! used by RRTMGP gas optics 
         nabsorbers,                      & ! used by RRTMGP gas optics 
         nextrabsorbers,                  & ! used by RRTMGP gas optics 
         nminorabsorbers,                 & ! used by RRTMGP gas optics 
         nmixingfracs,                    & ! used by RRTMGP gas optics 
         nlayers,                         & ! used by RRTMGP gas optics 
         nbnds,                           & ! used by RRTMGP gas optics 
         npairs,                          & ! used by RRTMGP gas optics 
         ninternalSourcetemps,            & ! used by RRTMGP gas optics 
         nminor_absorber_intervals_lower, & ! used by RRTMGP gas optics 
         nminor_absorber_intervals_upper, & ! used by RRTMGP gas optics 
         ncontributors_lower,             & ! used by RRTMGP gas optics 
         ncontributors_upper                ! used by RRTMGP gas optics 

    ! Local variables
    integer :: ncid_lw,dimID,varID,status,igpt,iGas,ij,ierr
    integer,dimension(:),allocatable :: temp1,temp2,temp3,temp4,temp_log_array1,&
    	temp_log_array2, temp_log_array3, temp_log_array4
    character(len=264) :: lw_gas_props_file
    integer,parameter :: max_strlen=256

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the gfs_physics_nml (scm/src/GFS_typedefs.F90)
    lw_gas_props_file  = trim(Model%rrtmgp_root)//trim(Model%lw_file_gas)

    ! Read dimensions for k-distribution fields (only on master processor(0))
    if (mpirank .eq. mpiroot) then
       if(nf90_open(trim(lw_gas_props_file), NF90_WRITE, ncid_lw) .eq. NF90_NOERR) then
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
          status = nf90_inquire_dimension(ncid_lw, dimid, len=ngpts_lw)
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
    call MPI_BCAST(ngpts_lw,                        1, MPI_INTEGER, mpiroot, mpicomm, ierr)
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
    allocate(kmajor(ngpts_lw, nmixingfracs,  npress+1, ntemps))
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
    allocate(planck_frac(ngpts_lw, nmixingfracs, npress+1, ntemps))

    if (mpirank .eq. mpiroot) then
       ! Read in fields from file
       if(nf90_open(trim(lw_gas_props_file), NF90_WRITE, ncid_lw) .eq. NF90_NOERR) then
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
    do iGas=1,Model%nGases
       call check_error_msg('lw_gas_optics_init',gas_concentrations%set_vmr(Radtend%active_gases(iGas), 0._kind_phys))
    enddo    
    call check_error_msg('lw_gas_optics_init',lw_gas_props%load(gas_concentrations, gas_names, &
         key_species, band2gpt, band_lims, press_ref, press_ref_trop, temp_ref,  temp_ref_p,   &
         temp_ref_t,  vmr_ref, kmajor, kminor_lower, kminor_upper, gas_minor,identifier_minor, &
         minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower, minor_limits_gpt_upper, &
         minor_scales_with_density_lower,  minor_scales_with_density_upper, scaling_gas_lower, &
         scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper,              &
         kminor_start_lower, kminor_start_upper, totplnk, planck_frac, rayl_lower, rayl_upper))

    ! Set initial permutation seed for McICA, initially set to number of G-points
    ipsdlw0 = lw_gas_props%get_ngpt()
  end subroutine rrtmgp_lw_gas_optics_init

  subroutine rrtmgp_lw_gas_optics_run()
  end subroutine rrtmgp_lw_gas_optics_run
  subroutine rrtmgp_lw_gas_optics_finalize()
  end subroutine rrtmgp_lw_gas_optics_finalize

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
end module rrtmgp_lw_gas_optics
