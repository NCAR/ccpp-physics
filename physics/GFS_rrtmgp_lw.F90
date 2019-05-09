! ###########################################################################################
! ###########################################################################################
module GFS_rrtmgp_lw
  use mo_gas_concentrations,     only: ty_gas_concs
  use mo_fluxes,                 only: ty_fluxes_broadband
  use mo_fluxes_byband,          only: ty_fluxes_byband
  use mo_optical_props,          only: ty_optical_props_1scl,ty_optical_props_2str
  use mo_source_functions,       only: ty_source_func_lw
  use mo_rte_kind,               only: wl
  use mo_heating_rates,          only: compute_heating_rate
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp_type
  use mo_cloud_optics,           only: ty_cloud_optics_type
  use mo_cloud_sampling,         only: sampled_mask_max_ran, sampled_mask_exp_ran, draw_samples
  use machine,                   only: kind_phys
  use module_radlw_parameters,   only: topflw_type, sfcflw_type, proflw_type
  use physparam,                 only: ilwcliq,isubclw,iovrlw,ilwrgas,icldflg,ilwrate
  use GFS_typedefs,              only: GFS_control_type
  use mo_rrtmgp_constants,       only: grav, avogad
  use mo_rrtmgp_lw_cloud_optics, only: rrtmgp_lw_cloud_optics
  use mersenne_twister,          only: random_setseed, random_number, random_stat
  use mo_rrtmgp_clr_all_sky,     only: rte_lw

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

  ! Logical flags for optional output fields in rrtmgp_lw_run(), default=.false.
  logical :: &
       l_AllSky_HR_byband  = .false., & ! 2D [ncol,nlay] all-sky heating rates, in each band [ncol,nlay,nBandsLW]?
       l_ClrSky_HR         = .false., & ! 2D [ncol,nlay] clear-sky heating rate?
       l_fluxes2D          = .false.    ! 2D [ncol,nlay] radiative fluxes? *Note* fluxes is a DDT w/ 4 fields.

  ! Module parameters (set during rrtmgp_lw_init())
  integer :: &
       rrtmgp_lw_cld_phys, & ! RRTMGP cloud-physics (0-RRTMG, 1-RRTGMP(LUT), 2-RRTMGP(Pade))
       nGptsLW,            & ! Number of LW spectral g-points
       nBandsLW,           & ! Number of LW bands
       nrghice,            & ! Number of ice roughness categories
       ipsdlw0               ! Initial see for McICA

  integer,allocatable,dimension(:) :: &
       ngb_LW   ! Band index for each g-points

  ! Classes used by rte+rrtmgp
  !type(ty_gas_optics_rrtmgp) :: &
  !     kdist_lw
  !type(ty_cloud_optics) :: &
  !     kdist_lw_cldy
  type(ty_gas_concs)  :: &
       gas_concentrations

  public GFS_rrtmgp_lw_init, GFS_rrtmgp_lw_run, GFS_rrtmgp_lw_finalize
contains
  ! #########################################################################################
  ! GFS_rrtmgp_lw_init
  ! #########################################################################################
!! \section arg_table_GFS_rrtmgp_lw_init Argument Table
!! | local_name         | standard_name                                   | long_name                                                                 | units | rank | type                      |    kind   | intent | optional |
!! |--------------------|-------------------------------------------------|---------------------------------------------------------------------------|-------|------|---------------------------|-----------|--------|----------|
!! | Model              | GFS_control_type_instance                       | Fortran DDT containing FV3-GFS model control parameters                   | DDT   |    0 | GFS_control_type          |           | in     | F        |
!! | mpirank            | mpi_rank                                        | current MPI rank                                                          | index |    0 | integer                   |           | in     | F        |
!! | mpiroot            | mpi_root                                        | master MPI rank                                                           | index |    0 | integer                   |           | in     | F        |
!! | mpicomm            | mpi_comm                                        | MPI communicator                                                          | index |    0 | integer                   |           | in     | F        |
!! | errmsg             | ccpp_error_message                              | error message for error handling in CCPP                                  | none  |    0 | character                 | len=*     | out    | F        |
!! | errflg             | ccpp_error_flag                                 | error flag for error handling in CCPP                                     | flag  |    0 | integer                   |           | out    | F        |
!! | kdist_lw           | K_distribution_file_for_RRTMGP_LW_scheme        | DDT containing spectral information for RRTMGP LW radiation scheme        | DDT   |    0 | ty_gas_optics_rrtmgp_type |           | inout  | F        |
!! | kdist_lw_cldy      | K_distribution_file_for_cloudy_RRTMGP_LW_scheme | DDT containing spectral information for cloudy RRTMGP LW radiation scheme | DDT   |    0 | ty_cloud_optics_type      |           | inout  | F        |
!!
  ! #########################################################################################
  subroutine GFS_rrtmgp_lw_init(Model,mpicomm, mpirank, mpiroot, kdist_lw, kdist_lw_cldy,   &
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
    type(ty_gas_optics_rrtmgp_type),intent(inout) :: &
         kdist_lw
    type(ty_cloud_optics_type),intent(inout) :: &
         kdist_lw_cldy
!    type(ty_gas_concs_type),intent(inout)  :: &
!         gas_concentrations

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

    ! Set band index by g-point array 
    nBandsLW = kdist_lw%get_nband()
    nGptsLW  = kdist_lw%get_ngpt()
    ngb_LW   = kdist_lw%get_gpoint_bands()

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
       call check_error_msg(kdist_lw_cldy%set_ice_roughness(nrghice))
       call check_error_msg(kdist_lw_cldy%load(band_lims_cldy, radliq_lwr, radliq_upr, &
            radliq_fac, radice_lwr, radice_upr, radice_fac, lut_extliq, lut_ssaliq,    &
            lut_asyliq, lut_extice, lut_ssaice, lut_asyice))
    endif
    if (rrtmgp_lw_cld_phys .eq. 2) then
       call check_error_msg(kdist_lw_cldy%set_ice_roughness(nrghice))
       call check_error_msg(kdist_lw_cldy%load(band_lims_cldy, pade_extliq,            &
            pade_ssaliq, pade_asyliq, pade_extice, pade_ssaice, pade_asyice,           &
            pade_sizereg_extliq, pade_sizereg_ssaliq, pade_sizereg_asyliq,             &
            pade_sizereg_extice, pade_sizereg_ssaice, pade_sizereg_asyice))
    endif

  end subroutine GFS_rrtmgp_lw_init

  ! #########################################################################################
  ! rrtmg_lw_run
  ! #########################################################################################
!! \section arg_table_GFS_rrtmgp_lw_run Argument Table
!! | local_name      | standard_name                                                                                 | long_name                                                 | units   | rank | type        |    kind   | intent | optional |
!! |-----------------|-----------------------------------------------------------------------------------------------|-----------------------------------------------------------|---------|------|-------------|-----------|--------|----------|
!! | p_lay           | air_pressure_at_layer_for_radiation_in_hPa                                                    | air pressure layer                                        | hPa     |    2 | real        | kind_phys | in     | F        |
!! | p_lev           | air_pressure_at_interface_for_radiation_in_hPa                                                | air pressure level                                        | hPa     |    2 | real        | kind_phys | in     | F        |
!! | t_lay           | air_temperature_at_layer_for_radiation                                                        | air temperature layer                                     | K       |    2 | real        | kind_phys | in     | F        |
!! | t_lev           | air_temperature_at_interface_for_radiation                                                    | air temperature level                                     | K       |    2 | real        | kind_phys | in     | F        |
!! | q_lay           | water_vapor_specific_humidity_at_layer_for_radiation                                          | specific humidity layer                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | o3_lay          | ozone_concentration_at_layer_for_radiation                                                    | ozone concentration layer                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_co2         | volume_mixing_ratio_co2                                                                       | volume mixing ratio co2                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_n2o         | volume_mixing_ratio_n2o                                                                       | volume mixing ratio no2                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_ch4         | volume_mixing_ratio_ch4                                                                       | volume mixing ratio ch4                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_o2          | volume_mixing_ratio_o2                                                                        | volume mixing ratio o2                                    | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_co          | volume_mixing_ratio_co                                                                        | volume mixing ratio co                                    | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_cfc11       | volume_mixing_ratio_cfc11                                                                     | volume mixing ratio cfc11                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_cfc12       | volume_mixing_ratio_cfc12                                                                     | volume mixing ratio cfc12                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_cfc22       | volume_mixing_ratio_cfc22                                                                     | volume mixing ratio cfc22                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | vmr_ccl4        | volume_mixing_ratio_ccl4                                                                      | volume mixing ratio ccl4                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | icseed          | seed_random_numbers_lw                                                                        | seed for random number generation for longwave radiation  | none    |    1 | integer     |           | in     | F        |
!! | tau_aer         | aerosol_optical_depth_for_longwave_bands_01-16                                                | aerosol optical depth for longwave bands 01-16            | none    |    3 | real        | kind_phys | in     | F        |
!! | ssa_aer         | aerosol_single_scattering_albedo_for_longwave_bands_01-16                                     | aerosol single scattering albedo for longwave bands 01-16 | frac    |    3 | real        | kind_phys | in     | F        |
!! | sfc_emiss       | surface_longwave_emissivity                                                                   | surface emissivity                                        | frac    |    1 | real        | kind_phys | in     | F        |
!! | skt             | surface_ground_temperature_for_radiation                                                      | surface ground temperature for radiation                  | K       |    1 | real        | kind_phys | in     | F        |
!! | dzlyr           | layer_thickness_for_radiation                                                                 | layer thickness                                           | km      |    2 | real        | kind_phys | in     | F        |
!! | delpin          | layer_pressure_thickness_for_radiation                                                        | layer pressure thickness                                  | hPa     |    2 | real        | kind_phys | in     | F        | 
!! | de_lgth         | cloud_decorrelation_length                                                                    | cloud decorrelation length                                | km      |    1 | real        | kind_phys | in     | F        | 
!! | ncol            | horizontal_loop_extent                                                                        | horizontal dimension                                      | count   |    0 | integer     |           | in     | F        |
!! | nlay            | adjusted_vertical_layer_dimension_for_radiation                                               | number of vertical layers for radiation                   | count   |    0 | integer     |           | in     | F        |
!! | lprint          | flag_print                                                                                    | flag to print                                             | flag    |    0 | logical     |           | in     | F        |
!! | cldfrac         | total_cloud_fraction                                                                          | total cloud fraction                                      | frac    |    2 | real        | kind_phys | in     | F        |
!! | lslwr           | flag_to_calc_lw                                                                               | flag to calculate LW irradiances                          | flag    |    0 | logical     |           | in     | F        |
!! | hlwc            | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step                    | longwave total sky heating rate                           | K s-1   |    2 | real        | kind_phys | out    | F        |
!! | topflx          | lw_fluxes_top_atmosphere                                                                      | longwave total sky fluxes at the top of the atm           | W m-2   |    1 | topflw_type |           | inout  | F        |
!! | sfcflx          | lw_fluxes_sfc                                                                                 | longwave total sky fluxes at the Earth surface            | W m-2   |    1 | sfcflw_type |           | inout  | F        |
!! | cldtau          | cloud_optical_depth_layers_at_10mu_band                                                       | approx 10mu band layer cloud optical depth                | none    |    2 | real        | kind_phys | inout  | F        |
!! | hlw0            | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step | longwave clear sky heating rate                           | K s-1   |    2 | real        | kind_phys | inout  | T        |
!! | hlwb            | lw_heating_rate_spectral                                                                      | longwave total sky heating rate (spectral)                | K s-1   |    3 | real        | kind_phys | inout  | T        |
!! | flxprf          | lw_fluxes                                                                                     | lw fluxes total sky / csk and up / down at levels         | W m-2   |    2 | proflw_type |           | inout  | T        |
!! | cld_lwp         | cloud_liquid_water_path                                                                       | cloud liquid water path                                   | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_liq     | mean_effective_radius_for_liquid_cloud                                                        | mean effective radius for liquid cloud                    | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_iwp         | cloud_ice_water_path                                                                          | cloud ice water path                                      | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_ice     | mean_effective_radius_for_ice_cloud                                                           | mean effective radius for ice cloud                       | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_rwp         | cloud_rain_water_path                                                                         | cloud ice water path                                      | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_rain    | mean_effective_radius_for_rain_drop                                                           | mean effective radius for rain drop                       | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_swp         | cloud_snow_water_path                                                                         | cloud snow water path                                     | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_snow    | mean_effective_radius_for_snow_flake                                                          | mean effective radius for snow flake                      | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_od          | cloud_optical_depth                                                                           | cloud optical depth                                       | none    |    2 | real        | kind_phys | in     | T        |
!! | errmsg          | ccpp_error_message                                                                            | error message for error handling in CCPP                  | none    |    0 | character   | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                                                               | error flag for error handling in CCPP                     | flag    |    0 | integer     |           | out    | F        |
!! | kdist_lw        | K_distribution_file_for_RRTMGP_LW_scheme                                                      | DDT containing spectral information for RRTMGP LW radiation scheme | DDT   |    0 | ty_gas_optics_rrtmgp_type |           | in     | F        |
!! | kdist_lw_cldy   | K_distribution_file_for_cloudy_RRTMGP_LW_scheme                                               | DDT containing spectral information for cloudy RRTMGP LW radiation scheme | DDT   |    0 | ty_cloud_optics_type |           | inout  | F        |
!!
  ! #########################################################################################
  subroutine GFS_rrtmgp_lw_run(p_lay, p_lev, t_lay, t_lev, q_lay, o3_lay, vmr_co2, vmr_n2o,     & ! IN
       vmr_ch4, vmr_o2, vmr_co, vmr_cfc11, vmr_cfc12, vmr_cfc22, vmr_ccl4, icseed, tau_aer, & ! IN
       ssa_aer, sfc_emiss, skt, dzlyr, delpin, de_lgth, ncol, nlay, lprint, cldfrac, lslwr, & ! IN
       kdist_lw, kdist_lw_cldy, &
       hlwc, topflx, sfcflx, cldtau,                                                        & ! OUT
       hlw0, hlwb, flxprf,                                                                  & ! OPT(out)
       cld_lwp, cld_ref_liq, cld_iwp, cld_ref_ice, cld_rwp, cld_ref_rain, cld_swp,          & ! OPT(in)
       cld_ref_snow, cld_od, errmsg, errflg)                                                  ! OPT(in)

    ! Inputs
    integer,intent(in) :: &
         ncol,         & ! Number of horizontal grid-points
         nlay            ! Number of vertical layers
    integer,intent(in),dimension(ncol) :: &
         icseed          ! auxiliary special cloud related array when module 
                         ! variable isubclw=2, it provides permutation seed 
                         ! for each column profile that are used for generating 
                         ! random numbers. when isubclw /=2, it will not be used.
    logical,intent(in) :: &
         lprint,       & ! Control flag for diagnostics
         lslwr           ! Flag to calculate RRTMGP LW?           (1)
    type(ty_gas_optics_rrtmgp_type),intent(in) :: &
         kdist_lw        ! DDT containing LW spectral information
    type(ty_cloud_optics_type),intent(in) :: &
         kdist_lw_cldy
!    type(ty_gas_concs_type),intent(inout) :: &
!         gas_concentrations
    real(kind_phys), dimension(ncol), intent(in) :: &
         sfc_emiss,    & ! Surface emissivity                     (1)
         skt,          & ! Surface(skin) temperature              (K)
         de_lgth         ! Cloud decorrelation length             (km)
    real(kind_phys), dimension(ncol,nlay), intent(in) :: &
         dzlyr,        & ! layer thinkness                        (km)
         delpin,       & ! layer thickness                        (mb)
         cldfrac,      & ! Cloud-fraction                         (1)
         p_lay,        & ! Pressure @ model layer-centers         (mb)
         t_lay,        & ! Temperature                            (K)
         q_lay,        & ! Specific humidity                      (kg/kg)
         o3_lay,       & ! O3 mass mixing-ratio                   (kg/kg)
         vmr_co2,      & ! Co2 volume-mixing ratio                (kg/kg)
         vmr_n2o,      & ! N2o volume-mixing ratio                (kg/kg)
         vmr_ch4,      & ! Ch4 volume-mixing ratio                (kg/kg)
         vmr_o2,       & ! O2 volume-mixing ratio                 (kg/kg)
         vmr_co,       & ! Co volume-mixing ratio                 (kg/kg)
         vmr_cfc11,    & ! CFC11 volume-mixing ratio              (kg/kg)
         vmr_cfc12,    & ! CFC12 volume-mixing ratio              (kg/kg)
         vmr_cfc22,    & ! CFC22 volume-mixing ratio              (kg/kg)
         vmr_ccl4        ! CCl4 volume-mixing ratio               (kg/kg)
    real(kind_phys), dimension(ncol,nlay+1), intent(in) :: &
         p_lev,        & ! Pressure @ model layer-interfaces      (mb)
         t_lev           ! Temperature                            (K)
    real(kind_phys), dimension(ncol,nlay,nbandsLW),intent(in) :: &
         tau_aer,      & ! Aerosol optical depth                  (1) 
         ssa_aer         ! Aerosol single-scattering albedo       (1)
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
         cld_od          ! Cloud optical-depth                    (1)

    ! Outputs (mandatory)
    character(len=*), intent(out) :: &
         errmsg          ! Error message
    integer,          intent(out) :: &
         errflg          ! Error code
    real(kind_phys),dimension(ncol,nlay),intent(inout) :: &
         hlwc,         & ! All-sky heating-rate                   (K/sec)
         cldtau          ! ~10mu band layer tau                   (1)
    type(topflw_type), dimension(ncol), intent(inout) :: &
         topflx          ! radiation fluxes at top, components:
                         ! upfxc - total sky upward flux at top   (w/m2)
                         ! upfx0 - clear sky upward flux at top   (w/m2)
    type(sfcflw_type), dimension(ncol), intent(inout) :: &
         sfcflx          ! radiation fluxes at sfc, components:
                         ! upfxc - total sky upward flux at sfc   (w/m2)  
                         ! upfx0 - clear sky upward flux at sfc   (w/m2)
                         ! dnfxc - total sky downward flux at sfc (w/m2)
                         ! dnfx0 - clear sky downward flux at sfc (w/m2)
    ! Outputs (optional)
    real(kind_phys), dimension(ncol,nlay,nbandsLW), optional, intent(inout) :: &
         hlwb            ! All-sky heating rate, in each band     (K/sec)
    real(kind_phys), dimension(ncol,nlay), optional, intent(inout) :: &
         hlw0            ! Clear-sky heating rate                 (K/sec)
    type(proflw_type), dimension(ncol,nlay+1), optional, intent(inout) :: &
         flxprf          ! 2D radiative fluxes, components:
                         ! upfxc - total sky upward flux          (W/m2)
                         ! dnfxc - total sky dnward flux          (W/m2)
                         ! upfx0 - clear sky upward flux          (W/m2)
                         ! dnfx0 - clear sky dnward flux          (W/m2)

    ! Local variables
    integer :: iGpt,iCol,iLay,iBand,iTOA,iSFC
    integer,dimension(ncol) :: ipseed
    real(kind_phys), dimension(nBandsLW,ncol) :: &
         semiss        
    real(kind_phys), dimension(ncol,nlay+1),target :: &
         flux_up_allSky, flux_up_clrSky, flux_dn_allSky, flux_dn_clrSky
    real(kind_phys), dimension(ncol,nlay+1,nBandsLW),target :: &
         fluxBB_up_allSky, fluxBB_dn_allSky
    real(kind_phys), dimension(ncol,nlay) :: &
         vmr_o3, vmr_h2o, thetaTendClrSky,thetaTendAllSky, cld_ref_liq2, &
         cld_ref_ice2,tau_snow,tau_rain
    real(kind_phys), dimension(ncol,nlay,nBandsLW) :: &
         tau_cld,thetaTendByBandAllSky
    real(kind_phys), dimension(nGptsLW,nlay,ncol) :: &
         rng3D
    real(kind_phys), dimension(nGptsLW*nLay) :: &
         rng1D
    logical,dimension(ncol,nlay) :: &
         liqmask,icemask
    logical, dimension(ncol,nlay,nGptsLW) :: &
         cldfracMCICA
    logical :: &
         top_at_1=.false.

    ! Types used by Random Number Generator
    type(random_stat) :: rng_stat

    ! RTE+RRTMGP classes
    type(ty_optical_props_1scl) :: &
         optical_props_clr,  & ! Optical properties for gaseous atmosphere
         optical_props_cldy, & ! Optical properties for clouds (by band)
         optical_props_mcica,& ! Optical properties for clouds (sampled)
         optical_props_aer     ! Optical properties for aerosols

    type(ty_source_func_lw) :: &
         sources               ! source function
    type(ty_fluxes_byband) :: &
         fluxAllSky,         & ! All-sky flux                      (W/m2)
         fluxClrSky            ! Clear-sky flux                    (W/m2)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    if (.not. lslwr) return

    ! Some consistency checks...

    ! Are any optional outputs requested?
    l_ClrSky_HR        = present(hlw0)
    l_AllSky_HR_byband = present(hlwb)
    l_fluxes2D         = present(flxprf)

    ! Check for optional input arguments, this depends on cloud method
    if (ilwcliq > 0) then    ! use prognostic cloud method
       if (.not. present(cld_lwp) .or. .not. present(cld_ref_liq)  .or. &
           .not. present(cld_iwp) .or. .not. present(cld_ref_ice)  .or. &
           .not. present(cld_rwp) .or. .not. present(cld_ref_rain) .or. &
           .not. present(cld_swp) .or. .not. present(cld_ref_snow)) then
          write(errmsg,'(*(a))')                                        &
               'Logic error: ilwcliq>0 requires the following',         &
               ' optional arguments to be present:',                    &
               ' cld_lwp, cld_ref_liq, cld_iwp, cld_ref_ice,',          &
               ' cld_rwp, cld_ref_rain, cld_swp, cld_ref_snow'         
          errflg = 1
          return
       end if
    else                     ! use diagnostic cloud method
       if (.not. present(cld_od) ) then
          write(errmsg,'(*(a))')                                        &
               'Logic error: ilwcliq<=0 requires the following',        &
               ' optional argument to be present: cld_od'
          errflg = 1
          return
       end if
    end if

    ! What is vertical ordering?
    top_at_1 = (p_lay(1,1) .lt. p_lay(1,nlay))
    if (top_at_1) then 
       iSFC = nlay+1
       iTOA = 1
    else
       iSFC = 1
       iTOA = nlay+1
    endif

    ! Change random number seed value for each radiation invocation (isubclw =1 or 2).
    if(isubclw == 1) then      ! advance prescribed permutation seed
       do iCol = 1, ncol
          ipseed(iCol) = ipsdlw0 + iCol
       enddo
    elseif (isubclw == 2) then ! use input array of permutaion seeds
       do iCol = 1, ncol
          ipseed(iCol) = icseed(iCol)
       enddo
    endif

    ! Surface emissivity
    semiss(:,:) = 1.
    do iBand=1,nBandsLW
       where(sfc_emiss .gt. epsilon .and. sfc_emiss .le. 1) semiss(iBand,:) = sfc_emiss
    enddo

    ! Compute volume mixing-ratios for ozone (mmr) and specific-humidity.
    vmr_h2o = merge((q_lay/(1-q_lay))*amdw, 0., q_lay  .ne. 1.)
    vmr_o3  = merge(o3_lay*amdo3,           0., o3_lay .gt. 0.)

    ! Compute ice/liquid cloud masks, needed by rrtmgp_cloud_optics
    liqmask = (cldfrac .gt. 0 .and. cld_lwp .gt. 0)
    icemask = (cldfrac .gt. 0 .and. cld_iwp .gt. 0)

    ! RRTMGP cloud_optics expects particle size to be in a certain range. bound here
    if (rrtmgp_lw_cld_phys .gt. 0) then
       cld_ref_ice2 = cld_ref_ice
       where(cld_ref_ice2 .gt. kdist_lw_cldy%get_max_radius_ice()) cld_ref_ice2=kdist_lw_cldy%get_max_radius_ice()
       where(cld_ref_ice2 .lt. kdist_lw_cldy%get_min_radius_ice()) cld_ref_ice2=kdist_lw_cldy%get_min_radius_ice()
       cld_ref_liq2 = cld_ref_liq
       where(cld_ref_liq2 .gt. kdist_lw_cldy%get_max_radius_liq()) cld_ref_liq2=kdist_lw_cldy%get_max_radius_liq()
       where(cld_ref_liq2 .lt. kdist_lw_cldy%get_min_radius_liq()) cld_ref_liq2=kdist_lw_cldy%get_min_radius_liq()
    endif

    ! #######################################################################################
    ! Call RRTMGP
    ! #######################################################################################
    ! Allocate space for source functions and gas optical properties [ncol,nlay,ngpts]
    call check_error_msg(sources%alloc(                 nCol, nLay, kdist_lw))
    call check_error_msg(optical_props_clr%alloc_1scl(  nCol, nLay, kdist_lw))
    call check_error_msg(optical_props_mcica%alloc_1scl(nCol, nLay, kdist_lw))
    ! Cloud optics [nCol,nLay,nBands]
    call check_error_msg(optical_props_cldy%init(optical_props_clr%get_band_lims_wavenumber()))
    call check_error_msg(optical_props_cldy%alloc_1scl(ncol,nlay))
    ! Aerosol optics [Ccol,nLay,nBands]
    call check_error_msg(optical_props_aer%init(optical_props_clr%get_band_lims_wavenumber()))
    call check_error_msg(optical_props_aer%alloc_1scl(ncol,nlay))

    ! Initialize RRTMGP files 
    fluxAllSky%flux_up   => flux_up_allSky
    fluxAllsky%flux_dn   => flux_dn_allSky
    fluxClrSky%flux_up   => flux_up_clrSky
    fluxClrsky%flux_dn   => flux_dn_clrSky
    ! Only calculate fluxes by-band, only when heating-rate profiles by band are requested.
    if (l_AllSky_HR_byband) then
       fluxAllSky%bnd_flux_up => fluxBB_up_allSky
       fluxAllsky%bnd_flux_dn => fluxBB_dn_allSky
    endif

    ! #######################################################################################
    ! Set gas concentrations
    ! #######################################################################################
    call gas_concentrations%reset()
    call check_error_msg(gas_concentrations%set_vmr('o2',  vmr_o2))
    call check_error_msg(gas_concentrations%set_vmr('co2', vmr_co2))
    call check_error_msg(gas_concentrations%set_vmr('ch4', vmr_ch4))
    call check_error_msg(gas_concentrations%set_vmr('n2o', vmr_n2o))
    call check_error_msg(gas_concentrations%set_vmr('h2o', vmr_h2o))
    call check_error_msg(gas_concentrations%set_vmr('o3',  vmr_o3))

    ! #######################################################################################
    ! Copy aerosol to RRTMGP DDT
    ! #######################################################################################
    optical_props_aer%tau = tau_aer * (1. - ssa_aer)

    ! #######################################################################################
    ! Compute cloud-optics for RTE.
    ! #######################################################################################

    ! Compute in-cloud radiative properties 
    if (any(cldfrac .gt. 0)) then
       ! i) RRTMG cloud optics.
       ! If using RRTMG cloud-physics. Model can provide either cloud-optics (cld_od) or 
       ! cloud-properties by type (cloud LWP,snow effective radius, etc...)
       if (rrtmgp_lw_cld_phys .eq. 0) then
          ! Cloud-optical properties by type provided.
          if (.not. present(cld_od)) then
             call rrtmgp_lw_cloud_optics(ncol, nlay, nBandsLW, cld_lwp, cld_ref_liq, cld_iwp,  &
                  cld_ref_ice, cld_rwp, cld_ref_rain, cld_swp, cld_ref_snow, cldfrac, tau_cld)
             optical_props_cldy%tau = tau_cld
          else
          ! Cloud-optical depth provided.
             do iCol=1,ncol
                do iLay=1,nlay
                   if (cldfrac(iCol,iLay) .gt. 1e-20_kind_phys) then
                      optical_props_cldy%tau(iCol,iLay,:) = cld_od(iCol,iLay)
                   else
                      optical_props_cldy%tau(iCol,iLay,:) = 0._kind_phys
                   endif
                end do
             end do
          endif
       endif

       ! ii) Use RRTMGP cloud-optics.
       if (rrtmgp_lw_cld_phys .gt. 0) then
          call check_error_msg(kdist_lw_cldy%cloud_optics(ncol, nlay, nBandsLW, nrghice,     &
               liqmask, icemask, cld_lwp, cld_iwp, cld_ref_liq2, cld_ref_ice2, optical_props_cldy))        
       end if
    endif

    ! #######################################################################################
    ! Call McICA to generate subcolumns.
    ! #######################################################################################
    if (isubclw .gt. 0) then

       ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
       ! and layers. ([nGpts,nLayer,nColumn]-> [nGpts*nLayer]*nColumn)
       do iCol=1,nCol
          call random_setseed(ipseed(icol),rng_stat)
          call random_number(rng1D,rng_stat)
          rng3D(:,:,iCol) = reshape(source = rng1D,shape=[nGptsLW,nLay])
       enddo

       ! Call McICA
       select case ( iovrlw )
       ! Maximumn-random 
       case(1)
          call check_error_msg(sampled_mask_max_ran(rng3D,cldfrac,cldfracMCICA))       
       end select
          
       ! Map band optical depth to each g-point using McICA
       call check_error_msg(draw_samples(cldfracMCICA,optical_props_cldy,optical_props_mcica))
    endif

    ! #######################################################################################
    ! Compute fluxes
    ! #######################################################################################
    call check_error_msg(rte_lw(  &
         kdist_lw,                &
         gas_concentrations,      &
         p_lay(1:ncol,1:nlay),    &
         t_lay(1:ncol,1:nlay),    &
         p_lev(1:ncol,1:nlay+1),  &
         skt(1:ncol),             &
         semiss,                  &
         optical_props_mcica,     &
         fluxAllSky,              &
         fluxClrSky,              &
         aer_props = optical_props_aer))

    ! #######################################################################################
    ! Compute heating rates
    ! #######################################################################################
    if (l_ClrSky_HR) then
       call check_error_msg(compute_heating_rate(     &
            fluxClrSky%flux_up,                       &
            fluxClrSky%flux_dn,                       &
            p_lev(1:ncol,1:nlay+1),                   &
            thetaTendClrSky))
    endif
    if (l_AllSky_HR_byband) then
       do iBand=1,nBandsLW
          call check_error_msg(compute_heating_rate(  &
               fluxAllSky%bnd_flux_up(:,:,iBand),     &
               fluxAllSky%bnd_flux_dn(:,:,iBand),     &
               p_lev(1:ncol,1:nlay+1),                &
               thetaTendByBandAllSky(:,:,iBand)))
       enddo
    else
       call check_error_msg(compute_heating_rate(     &
            fluxAllSky%flux_up,                       &
            fluxAllSky%flux_dn,                       &
            p_lev(1:ncol,1:nlay+1),                   &
            thetaTendAllSky))
    endif

    ! #######################################################################################
    ! Copy fluxes from RRTGMP types into model radiation types.
    ! #######################################################################################
    ! Mandatory outputs
    topflx%upfxc = fluxAllSky%flux_up(:,iTOA)
    topflx%upfx0 = fluxClrSky%flux_up(:,iTOA)
    sfcflx%upfxc = fluxAllSky%flux_up(:,iSFC)
    sfcflx%upfx0 = fluxClrSky%flux_up(:,iSFC)
    sfcflx%dnfxc = fluxAllSky%flux_dn(:,iSFC)
    sfcflx%dnfx0 = fluxClrSky%flux_dn(:,iSFC)
    cldtau       = optical_props_cldy%tau(:,:,7)
    hlwc         = thetaTendAllSky

    ! Optional output
    if(l_fluxes2D) then
       flxprf%upfxc = fluxAllSky%flux_up
       flxprf%dnfxc = fluxAllSky%flux_dn
       flxprf%upfx0 = fluxClrSky%flux_up
       flxprf%dnfx0 = fluxClrSky%flux_dn
    endif
    if (l_AllSky_HR_byband) then
       hlwb = thetaTendByBandAllSky
    endif
    if (l_ClrSky_HR) then
       hlw0 = thetaTendClrSky
    endif

  end subroutine GFS_rrtmgp_lw_run
  !
  subroutine GFS_rrtmgp_lw_finalize()
  end subroutine GFS_rrtmgp_lw_finalize

  ! #########################################################################################
  ! Ancillary functions
  ! #########################################################################################
  subroutine check_error_msg(error_msg)
    character(len=*), intent(in) :: error_msg
    
    if(error_msg /= "") then
       print*,"ERROR(GFS_rrtmgp_lw_main.F90): "
       print*,trim(error_msg)
       return
    end if
  end subroutine check_error_msg  
  
  
end module GFS_rrtmgp_lw
