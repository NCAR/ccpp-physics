module rrtmgp_sw_gas_optics
  use machine,                only: kind_phys
  use module_radiation_gases, only: NF_VGAS
  use mo_rte_kind,            only: wl
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations,  only: ty_gas_concs
  use rrtmgp_aux,             only: check_error_msg
  use mo_optical_props,       only: ty_optical_props_2str
  use mo_compute_bc,          only: compute_bc
  use netcdf

contains

  ! #########################################################################################
  ! SUBROUTINE sw_gas_optics_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_gas_optics_init
!! \htmlinclude rrtmgp_sw_gas_optics.html
!!
  subroutine rrtmgp_sw_gas_optics_init(rrtmgp_root_dir, rrtmgp_sw_file_gas, rrtmgp_nGases,   &
       active_gases_array, mpicomm, mpirank, mpiroot, sw_gas_props, ipsdsw0, errmsg, errflg)
    use netcdf
!#ifdef MPI
!    use mpi
!#endif

    ! Inputs
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,  & ! RTE-RRTMGP root directory
         rrtmgp_sw_file_gas  ! RRTMGP file containing coefficients used to compute gaseous optical properties
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
         errmsg              ! Error message
    integer,          intent(out) :: &
         errflg,           & ! Error code
         ipsdsw0             !
    type(ty_gas_optics_rrtmgp),intent(out) :: &
         sw_gas_props        ! RRTMGP DDT: 

    ! Fields from the K-distribution files
    ! Variables that will be passed to gas_optics%load()
    type(ty_gas_concs)  :: &
         gas_concentrations
    integer, dimension(:), allocatable :: &
         kminor_start_lower_sw,              & !  
         kminor_start_upper_sw                 !  
    integer, dimension(:,:), allocatable :: &
         band2gpt_sw,                        & !  
         minor_limits_gpt_lower_sw,          & !  
         minor_limits_gpt_upper_sw             !  
    integer, dimension(:,:,:), allocatable :: &
         key_species_sw                        !  
    real(kind_phys) :: &
         press_ref_trop_sw,                  & !  
         temp_ref_p_sw,                      & !  
         temp_ref_t_sw                         !  
    real(kind_phys), dimension(:), allocatable :: &
         press_ref_sw,                       & !  
         temp_ref_sw,                        & !  
         solar_source_sw                       !  
    real(kind_phys), dimension(:,:), allocatable :: &
         band_lims_sw                          !                          

    real(kind_phys), dimension(:,:,:), allocatable :: &
         vmr_ref_sw,                         & !  
         kminor_lower_sw,                    & !  
         kminor_upper_sw,                    & !  
         rayl_lower_sw,                      & !  
         rayl_upper_sw                         !  
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         kmajor_sw                             !  
    character(len=32),  dimension(:), allocatable :: &
         gas_names_sw,                       & !  
         gas_minor_sw,                       & !  
         identifier_minor_sw,                & !  
         minor_gases_lower_sw,               & !  
         minor_gases_upper_sw,               & !  
         scaling_gas_lower_sw,               & !  
         scaling_gas_upper_sw                  !  
    logical(wl), dimension(:), allocatable :: &
         minor_scales_with_density_lower_sw, & !  
         minor_scales_with_density_upper_sw, & !  
         scale_by_complement_lower_sw,       & !  
         scale_by_complement_upper_sw          !  
    ! Dimensions (to be broadcast across all processors)
    integer :: &
         ntemps_sw,                          & !  
         npress_sw,                          & !  
         ngpts_sw,                           & ! 
         nabsorbers_sw,                      & !  
         nextrabsorbers_sw,                  & !  
         nminorabsorbers_sw,                 & !  
         nmixingfracs_sw,                    & !  
         nlayers_sw,                         & !  
         nbnds_sw,                           & !  
         npairs_sw,                          & !  
         nminor_absorber_intervals_lower_sw, & !  
         nminor_absorber_intervals_upper_sw, & !  
         ncontributors_lower_sw,             & !  
         ncontributors_upper_sw                !  

    ! Local variables
    integer :: status,ncid_sw,dimid,varID,iGas
    integer,dimension(:),allocatable :: temp1,temp2,temp3,temp4
    character(len=264) :: sw_gas_props_file
!#ifdef MPI
!    integer :: ierr
!#endif

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the gphysics_nml
    sw_gas_props_file   = trim(rrtmgp_root_dir)//trim(rrtmgp_sw_file_gas)

    ! Read dimensions for k-distribution fields (only on master processor(0))
!    if (mpirank .eq. mpiroot) then
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
!    endif

    ! Broadcast dimensions to all processors
!#ifdef MPI
!    call MPI_BARRIER(mpicomm, ierr)
!    call MPI_BCAST(ntemps_sw,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(npress_sw,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nabsorbers_sw,                      1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nminorabsorbers_sw,                 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nextrabsorbers_sw,                  1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nmixingfracs_sw,                    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nlayers_sw,                         1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nbnds_sw,                           1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(ngpts_sw,                           1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(npairs_sw,                          1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(ncontributors_lower_sw,             1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(ncontributors_upper_sw,             1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nminor_absorber_intervals_lower_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    call MPI_BCAST(nminor_absorber_intervals_upper_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!#endif

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
!    if (mpirank .eq. mpiroot) then
       write (*,*) 'Reading RRTMGP shortwave k-distribution data ... '
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
!    endif

    ! Broadcast arrays to all processors
!#ifdef MPI
!    call MPI_BARRIER(mpicomm, ierr)
!    write (*,*) 'Broadcasting RRTMGP shortwave k-distribution data ... '
!    call MPI_BCAST(minor_limits_gpt_upper_sw,   size(minor_limits_gpt_upper_sw), MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(minor_limits_gpt_lower_sw,   size(minor_limits_gpt_lower_sw), MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_start_upper_sw,       size(kminor_start_upper_sw),     MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_start_lower_sw,       size(kminor_start_lower_sw),     MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(key_species_sw,              size(key_species_sw),            MPI_INTEGER,            mpiroot, mpicomm, ierr)
!    call MPI_BCAST(band2gpt_sw,                 size(band2gpt_sw),               MPI_INTEGER,            mpiroot, mpicomm, ierr)
!#ifndef SINGLE_PREC
!    call MPI_BCAST(band_lims_sw,                size(band_lims_sw),              MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(press_ref_sw,                size(press_ref_sw),              MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_sw,                 size(temp_ref_sw),               MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_lower_sw,             size(kminor_lower_sw),           MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_upper_sw,             size(kminor_upper_sw),           MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scaling_gas_lower_sw,        size(scaling_gas_lower_sw),      MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scaling_gas_upper_sw,        size(scaling_gas_upper_sw),      MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(vmr_ref_sw,                  size(vmr_ref_sw),                MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kmajor_sw,                   size(kmajor_sw),                 MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_p_sw,               1,                               MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_t_sw,               1,                               MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(press_ref_trop_sw,           1,                               MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(solar_source_sw,             size(solar_source_sw),           MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(rayl_lower_sw,               size(rayl_lower_sw),             MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!    call MPI_BCAST(rayl_upper_sw,               size(rayl_upper_sw),             MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!#else
!    call MPI_BCAST(band_lims_sw,                size(band_lims_sw),              MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(press_ref_sw,                size(press_ref_sw),              MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_sw,                 size(temp_ref_sw),               MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_lower_sw,             size(kminor_lower_sw),           MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kminor_upper_sw,             size(kminor_upper_sw),           MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scaling_gas_lower_sw,        size(scaling_gas_lower_sw),      MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scaling_gas_upper_sw,        size(scaling_gas_upper_sw),      MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(vmr_ref_sw,                  size(vmr_ref_sw),                MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(kmajor_sw,                   size(kmajor_sw),                 MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_p_sw,               1,                               MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(temp_ref_t_sw,               1,                               MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(press_ref_trop_sw,           1,                               MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(solar_source_sw,             size(solar_source_sw),           MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(rayl_lower_sw,               size(rayl_lower_sw),             MPI_REAL,               mpiroot, mpicomm, ierr)
!    call MPI_BCAST(rayl_upper_sw,               size(rayl_upper_sw),             MPI_REAL,               mpiroot, mpicomm, ierr)
!#endif
!    ! Character arrays
!    do ij=1,nabsorbers_sw
!       call MPI_BCAST(gas_names_sw(ij),         len(gas_names_sw(ij)),           MPI_CHAR,               mpiroot, mpicomm, ierr)
!    enddo
!    do ij=1,nminorabsorbers_sw
!       call MPI_BCAST(gas_minor_sw(ij),         len(gas_minor_sw(ij)),           MPI_CHAR,               mpiroot, mpicomm, ierr)
!       call MPI_BCAST(identifier_minor_sw(ij),  len(identifier_minor_sw(ij)),    MPI_CHAR,               mpiroot, mpicomm, ierr)
!    enddo
!    do ij=1,nminor_absorber_intervals_lower_sw
!       call MPI_BCAST(minor_gases_lower_sw(ij), len(minor_gases_lower_sw(ij)),   MPI_CHAR,               mpiroot, mpicomm, ierr)
!    enddo
!    do ij=1,nminor_absorber_intervals_upper_sw
!       call MPI_BCAST(minor_gases_upper_sw(ij), len(minor_gases_upper_sw(ij)),   MPI_CHAR,               mpiroot, mpicomm, ierr)
!    enddo
!
!    ! Logical arrays
!    call MPI_BCAST(minor_scales_with_density_lower_sw, nminor_absorber_intervals_lower_sw, MPI_LOGICAL,  mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scale_by_complement_lower_sw,       nminor_absorber_intervals_lower_sw, MPI_LOGICAL,  mpiroot, mpicomm, ierr)
!    call MPI_BCAST(minor_scales_with_density_upper_sw, nminor_absorber_intervals_upper_sw, MPI_LOGICAL,  mpiroot, mpicomm, ierr)
!    call MPI_BCAST(scale_by_complement_upper_sw,       nminor_absorber_intervals_upper_sw, MPI_LOGICAL,  mpiroot, mpicomm, ierr)
!#endif

    ! Initialize gas concentrations and gas optics class with data
    do iGas=1,rrtmgp_nGases
       call check_error_msg('sw_gas_optics_init',gas_concentrations%set_vmr(active_gases_array(iGas), 0._kind_phys))
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

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_gas_optics_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_gas_optics_run
!! \htmlinclude rrtmgp_sw_gas_optics.html
!!
  subroutine rrtmgp_sw_gas_optics_run(doSWrad, nCol, nLev, nday, idxday, sw_gas_props, p_lay,&
       p_lev, toa_src_sw, t_lay, t_lev, gas_concentrations, solcon, rrtmgp_nGases,           &
       active_gases_array, sw_optical_props_clrsky, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doSWrad             ! Flag to calculate SW irradiances
    integer,intent(in) :: &
         nDay,                 & ! Number of daylit points.
         nCol,                 & ! Number of horizontal points
         nLev                    ! Number of vertical levels
    integer,intent(in),dimension(ncol) :: &
         idxday                  ! Indices for daylit points.
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props            ! RRTMGP DDT: spectral information for RRTMGP SW radiation scheme
    real(kind_phys), dimension(ncol,nLev), intent(in) :: &
         p_lay,                & ! Pressure @ model layer-centers         (hPa)
         t_lay                   ! Temperature                            (K)
    real(kind_phys), dimension(ncol,nLev+1), intent(in) :: &
         p_lev,                & ! Pressure @ model layer-interfaces      (hPa)
         t_lev                   ! Temperature @ model levels
    type(ty_gas_concs),intent(in) :: &
         gas_concentrations      ! RRTMGP DDT: trace gas concentrations   (vmr)
    real(kind_phys), intent(in) :: &
         solcon                  ! Solar constant
    integer, intent(in) :: &
         rrtmgp_nGases       ! Number of trace gases active in RRTMGP
    character(len=128),dimension(rrtmgp_nGases), intent(in) :: &
         active_gases_array  ! Character array containing trace gases to include in RRTMGP

    ! Output
    character(len=*), intent(out) :: &
         errmsg                  ! Error message
    integer,          intent(out) :: &
         errflg                  ! Error code
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_clrsky ! RRTMGP DDT: clear-sky shortwave optical properties, spectral (tau,ssa,g) 
    real(kind_phys), dimension(ncol,sw_gas_props%get_ngpt()), intent(out) :: &
         toa_src_sw              ! TOA incident spectral flux (W/m2)

    ! Local variables
    integer :: ij,iGas
    real(kind_phys), dimension(ncol,nLev) :: vmrTemp
    real(kind_phys), dimension(nday,sw_gas_props%get_ngpt()) :: toa_src_sw_temp
    type(ty_gas_concs) :: &
         gas_concentrations_daylit    ! RRTMGP DDT: trace gas concentrations   (vmr)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doSWrad) return

    if (nDay .gt. 0) then
       ! Allocate space
       call check_error_msg('rrtmgp_sw_gas_optics_run',sw_optical_props_clrsky%alloc_2str(nday, nLev, sw_gas_props))

       ! Subset the gas concentrations, only need daylit points.
       do iGas=1,rrtmgp_nGases
          call check_error_msg('rrtmgp_sw_rte_run',&
               gas_concentrations%get_vmr(trim(active_gases_array(iGas)),vmrTemp))
          call check_error_msg('rrtmgp_sw_rte_run',&
               gas_concentrations_daylit%set_vmr(trim(active_gases_array(iGas)),vmrTemp(idxday(1:nday),:)))
       enddo

       ! Gas-optics
       call check_error_msg('rrtmgp_sw_gas_optics_run',sw_gas_props%gas_optics(&
            p_lay(idxday(1:nday),:),   & ! IN  - Pressure @ layer-centers (Pa)
            p_lev(idxday(1:nday),:),   & ! IN  - Pressure @ layer-interfaces (Pa)
            t_lay(idxday(1:nday),:),   & ! IN  - Temperature @ layer-centers (K)
            gas_concentrations_daylit, & ! IN  - RRTMGP DDT: trace gas volumne mixing-ratios
            sw_optical_props_clrsky,   & ! OUT - RRTMGP DDT: Shortwave optical properties, by
                                         !                   spectral point (tau,ssa,g)
            toa_src_sw_temp))            ! OUT - TOA incident shortwave radiation (spectral)
       toa_src_sw(idxday(1:nday),:) = toa_src_sw_temp
       ! Scale incident flux
       do ij=1,nday
          toa_src_sw(idxday(ij),:) = toa_src_sw(idxday(ij),:)*solcon/ &
                                     sum(toa_src_sw(idxday(ij),:))
       enddo
    endif

  end subroutine rrtmgp_sw_gas_optics_run

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_gas_optics_finalize
  ! #########################################################################################
  subroutine rrtmgp_sw_gas_optics_finalize()
  end subroutine rrtmgp_sw_gas_optics_finalize

end module rrtmgp_sw_gas_optics
 
