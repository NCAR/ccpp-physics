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

contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_gas_optics_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_gas_optics_init
!! \htmlinclude rrtmgp_lw_gas_optics.html
!!
  subroutine rrtmgp_lw_gas_optics_init(rrtmgp_root_dir, rrtmgp_lw_file_gas, rrtmgp_nGases,   &
       active_gases_array, mpicomm, mpirank, mpiroot, lw_gas_props, errmsg, errflg)
    use netcdf

#ifdef MPI
    use mpi
#endif

    ! Inputs
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,  & ! RTE-RRTMGP root directory
         rrtmgp_lw_file_gas  ! RRTMGP file containing coefficients used to compute gaseous optical properties
    integer, intent(in) :: &
         rrtmgp_nGases       ! Number of trace gases active in RRTMGP
    character(len=128),dimension(rrtmgp_nGases), intent(in) :: &
         active_gases_array  ! Character array containing trace gases to include in RRTMGP
    integer,intent(in) :: &
         mpicomm,          & ! MPI communicator
         mpirank,          & ! Current MPI rank
         mpiroot             ! Master MPI rank
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg              ! CCPP error message
    integer,          intent(out) :: &
         errflg              ! CCPP error code
    type(ty_gas_optics_rrtmgp),intent(out) :: &
         lw_gas_props        ! RRTMGP DDT: longwave spectral information

    ! Variables that will be passed to gas_optics%load()
    type(ty_gas_concs) :: &
         gas_concentrations                 ! RRTMGP DDT: trace gas concentrations (vmr)
    integer, dimension(:), allocatable :: &
         kminor_start_lower,              & ! Starting index in the [1, nContributors] vector for a contributor
                                            ! given by \"minor_gases_lower\" (lower atmosphere)  
         kminor_start_upper                 ! Starting index in the [1, nContributors] vector for a contributor 
                                            ! given by \"minor_gases_upper\" (upper atmosphere)  
    integer, dimension(:,:), allocatable :: &
         band2gpt,                        & ! Beginning and ending gpoint for each band   
         minor_limits_gpt_lower,          & ! Beginning and ending gpoint for each minor interval in lower atmosphere  
         minor_limits_gpt_upper             ! Beginning and ending gpoint for each minor interval in upper atmosphere
    integer, dimension(:,:,:), allocatable :: &
         key_species                        ! Key species pair for each band  
    real(kind_phys) :: &
         press_ref_trop,                  & ! Reference pressure separating the lower and upper atmosphere [Pa]   
         temp_ref_p,                      & ! Standard spectroscopic reference pressure [Pa]
         temp_ref_t                         ! Standard spectroscopic reference temperature [K]
    real(kind_phys), dimension(:), allocatable :: &
         press_ref,                       & ! Pressures for reference atmosphere; press_ref(# reference layers) [Pa]   
         temp_ref                           ! Temperatures for reference atmosphere; temp_ref(# reference layers) [K]  
    real(kind_phys), dimension(:,:), allocatable :: &
         band_lims,                       & ! Beginning and ending wavenumber [cm -1] for each band  
         totplnk                            ! Integrated Planck function by band  
    real(kind_phys), dimension(:,:,:), allocatable :: &
         vmr_ref,                         & ! volume mixing ratios for reference atmosphere   
         kminor_lower,                    & ! (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to 
                                            ! [nTemp x nEta x nContributors] array)  
         kminor_upper,                    & ! (transformed from [nTemp x nEta x nGpt x nAbsorbers] array to 
                                            ! [nTemp x nEta x nContributors] array)  
         rayl_lower,                      & ! Not used in LW, rather allocated(rayl_lower) is used  
         rayl_upper                         ! Not used in LW, rather allocated(rayl_upper) is used    
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         kmajor,                          & ! Stored absorption coefficients due to major absorbing gases  
         planck_frac                        ! Planck fractions   
    character(len=32),  dimension(:), allocatable :: &
         gas_names,                       & ! Names of absorbing gases  
         gas_minor,                       & ! Name of absorbing minor gas  
         identifier_minor,                & ! Unique string identifying minor gas  
         minor_gases_lower,               & ! Names of minor absorbing gases in lower atmosphere   
         minor_gases_upper,               & ! Names of minor absorbing gases in upper atmosphere   
         scaling_gas_lower,               & ! Absorption also depends on the concentration of this gas  
         scaling_gas_upper                  ! Absorption also depends on the concentration of this gas  
    logical(wl), dimension(:), allocatable :: &
         minor_scales_with_density_lower, & ! Density scaling is applied to minor absorption coefficients   
         minor_scales_with_density_upper, & ! Density scaling is applied to minor absorption coefficients 
         scale_by_complement_lower,       & ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T) 
         scale_by_complement_upper          ! Absorption is scaled by concentration of scaling_gas (F) or its complement (T)

    ! Dimensions
    integer :: &
         ntemps, npress, ngpts_lw, nabsorbers, nextrabsorbers, nminorabsorbers,&    
         nmixingfracs, nlayers, nbnds, npairs, ninternalSourcetemps,           &
         nminor_absorber_intervals_lower, nminor_absorber_intervals_upper,     &
         ncontributors_lower, ncontributors_upper

    ! Local variables
    integer :: ncid_lw, dimID, varID, status, iGas
    integer,dimension(:),allocatable :: temp1, temp2, temp3, temp4, &
         temp_log_array1, temp_log_array2, temp_log_array3, temp_log_array4
    character(len=264) :: lw_gas_props_file
#ifdef MPI
    integer :: ierr
#endif

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the physics_nml
    lw_gas_props_file  = trim(rrtmgp_root_dir)//trim(rrtmgp_lw_file_gas)

    ! Read dimensions for k-distribution fields (only on master processor(0))
!    if (mpirank .eq. mpiroot) then
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
!    endif

    ! Broadcast dimensions to all processors
!#ifdef MPI
!    call MPI_BARRIER(mpicomm, ierr)
!    call MPI_BCAST(ntemps,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(npress,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nabsorbers,                      1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nminorabsorbers,                 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nextrabsorbers,                  1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nmixingfracs,                    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nlayers,                         1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nbnds,                           1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(ngpts_lw,                        1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(npairs,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(ncontributors_lower,             1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(ncontributors_upper,             1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nminor_absorber_intervals_lower, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nminor_absorber_intervals_upper, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(ninternalSourcetemps,            1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BARRIER(mpicomm, ierr)
!#endif
    
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

!    if (mpirank .eq. mpiroot) then
       write (*,*) 'Reading RRTMGP longwave k-distribution data ... '
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
!    endif


    ! Broadcast arrays to all processors
!#ifdef MPI
!    call MPI_BARRIER(mpicomm, ierr)
!    if (mpirank==mpiroot) write (*,*) 'Broadcasting RRTMGP longwave k-distribution data ... '
!    call MPI_BCAST(minor_limits_gpt_upper,   size(minor_limits_gpt_upper), MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(minor_limits_gpt_lower,   size(minor_limits_gpt_lower), MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_start_upper,       size(kminor_start_upper),     MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_start_lower,       size(kminor_start_lower),     MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(key_species,              size(key_species),            MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(band2gpt,                 size(band2gpt),               MPI_INTEGER,            mpiroot, mpicomm, ierr)
!#ifndef SINGLE_PREC
!    call MPI_BCAST(band_lims,                size(band_lims),              MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(press_ref,                size(press_ref),              MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref,                 size(temp_ref),               MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_lower,             size(kminor_lower),           MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_upper,             size(kminor_upper),           MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scaling_gas_lower,        size(scaling_gas_lower),      MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scaling_gas_upper,        size(scaling_gas_upper),      MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(vmr_ref,                  size(vmr_ref),                MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kmajor,                   size(kmajor),                 MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_p,               1,                            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_t,               1,                            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(press_ref_trop,           1,                            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(totplnk,                  size(totplnk),                MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(planck_frac,              size(planck_frac),            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!#else
!    call MPI_BCAST(band_lims,                size(band_lims),              MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(press_ref,                size(press_ref),              MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref,                 size(temp_ref),               MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_lower,             size(kminor_lower),           MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_upper,             size(kminor_upper),           MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scaling_gas_lower,        size(scaling_gas_lower),      MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scaling_gas_upper,        size(scaling_gas_upper),      MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(vmr_ref,                  size(vmr_ref),                MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kmajor,                   size(kmajor),                 MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_p,               1,                            MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_t,               1,                            MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(press_ref_trop,           1,                            MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(totplnk,                  size(totplnk),                MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(planck_frac,              size(planck_frac),            MPI_REAL,               mpiroot, mpicomm, ierr)
!#endif
!    ! Character arrays
!    do ij=1,nabsorbers
!       call MPI_BCAST(gas_names(ij),         len(gas_names(ij)),           MPI_CHAR,               mpiroot, mpicomm, ierr)
!    enddo
!    do ij=1,nminorabsorbers
!       call MPI_BCAST(gas_minor(ij),         len(gas_minor(ij)),           MPI_CHAR,               mpiroot, mpicomm, ierr)
!       call MPI_BCAST(identifier_minor(ij),  len(identifier_minor(ij)),    MPI_CHAR,               mpiroot, mpicomm, ierr)
!    enddo
!    do ij=1,nminor_absorber_intervals_lower
!       call MPI_BCAST(minor_gases_lower(ij), len(minor_gases_lower(ij)),   MPI_CHAR,               mpiroot, mpicomm, ierr)
!    enddo
!    do ij=1,nminor_absorber_intervals_upper
!       call MPI_BCAST(minor_gases_upper(ij), len(minor_gases_upper(ij)),   MPI_CHAR,               mpiroot, mpicomm, ierr)
!    enddo
!    ! Logical arrays
!    !
!    call MPI_BCAST(minor_scales_with_density_lower, nminor_absorber_intervals_lower, MPI_LOGICAL,  mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scale_by_complement_lower,       nminor_absorber_intervals_lower, MPI_LOGICAL,  mpiroot, mpicomm, ierr)
!    call MPI_BCAST(minor_scales_with_density_upper, nminor_absorber_intervals_upper, MPI_LOGICAL,  mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scale_by_complement_upper,       nminor_absorber_intervals_upper, MPI_LOGICAL,  mpiroot, mpicomm, ierr)
!    call MPI_BARRIER(mpicomm, ierr)
!#endif

    ! Initialize gas concentrations and gas optics class with data
    do iGas=1,rrtmgp_nGases
       call check_error_msg('lw_gas_optics_init',gas_concentrations%set_vmr(active_gases_array(iGas), 0._kind_phys))
    enddo    
    call check_error_msg('lw_gas_optics_init',lw_gas_props%load(gas_concentrations, gas_names, &
         key_species, band2gpt, band_lims, press_ref, press_ref_trop, temp_ref,  temp_ref_p,   &
         temp_ref_t,  vmr_ref, kmajor, kminor_lower, kminor_upper, gas_minor,identifier_minor, &
         minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower, minor_limits_gpt_upper, &
         minor_scales_with_density_lower,  minor_scales_with_density_upper, scaling_gas_lower, &
         scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper,              &
         kminor_start_lower, kminor_start_upper, totplnk, planck_frac, rayl_lower, rayl_upper))

  end subroutine rrtmgp_lw_gas_optics_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_gas_optics_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_gas_optics_run
!! \htmlinclude rrtmgp_lw_gas_optics.html
!!
  subroutine rrtmgp_lw_gas_optics_run(doLWrad, nCol, nLev, lw_gas_props, p_lay, p_lev, t_lay,&
       t_lev, skt, gas_concentrations, lw_optical_props_clrsky, sources,  errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doLWrad                 ! Flag to calculate LW irradiances
    integer,intent(in) :: &
         ncol,                &  ! Number of horizontal points
         nLev                    ! Number of vertical levels
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props            ! RRTMGP DDT:
    real(kind_phys), dimension(ncol,nLev), intent(in) :: &
         p_lay,                & ! Pressure @ model layer-centers (hPa)
         t_lay                   ! Temperature (K)
    real(kind_phys), dimension(ncol,nLev+1), intent(in) :: &
         p_lev,                & ! Pressure @ model layer-interfaces (hPa)
         t_lev                   ! Temperature @ model levels
    real(kind_phys), dimension(ncol), intent(in) :: &
         skt                     ! Surface(skin) temperature (K)
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

    ! Allocate space
    call check_error_msg('rrtmgp_lw_gas_optics_run',lw_optical_props_clrsky%alloc_1scl(ncol, nLev, lw_gas_props))
    call check_error_msg('rrtmgp_lw_gas_optics_run',sources%init(lw_gas_props))
    call check_error_msg('rrtmgp_lw_gas_optics_run',sources%alloc(ncol, nLev))

    ! Gas-optics 
    call check_error_msg('rrtmgp_lw_gas_optics_run',lw_gas_props%gas_optics(&
         p_lay,                   & ! IN  - Pressure @ layer-centers (Pa)
         p_lev,                   & ! IN  - Pressure @ layer-interfaces (Pa)
         t_lay,                   & ! IN  - Temperature @ layer-centers (K)
         skt,                     & ! IN  - Skin-temperature (K)
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
