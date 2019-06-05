module rrtmgp_lw_cloud_optics
  use machine,                  only: kind_phys
  use GFS_typedefs,             only: GFS_control_type
  use mo_rte_kind,              only: wl
  use mo_cloud_optics,          only: ty_cloud_optics
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use physparam,                only: isubclw, iovrlw
  use mo_optical_props,         only: ty_optical_props_1scl
  use mo_cloud_sampling,        only: sampled_mask_max_ran, sampled_mask_exp_ran, draw_samples
  use mersenne_twister,         only: random_setseed, random_number, random_stat
  use mo_rrtmg_lw_cloud_optics, only: rrtmg_lw_cloud_optics   
  use rrtmgp_lw_gas_optics,     only: ipsdlw0
  use netcdf

  integer :: nrghice_lw
  public rrtmgp_lw_cloud_optics_init, rrtmgp_lw_cloud_optics_run, rrtmgp_lw_cloud_optics_finalize
contains

!! \section arg_table_rrtmgp_lw_cloud_optics_init Argument Table
!! | local_name     | standard_name                    | long_name                                                          | units | rank | type             | kind  | intent | optional |
!! |----------------|----------------------------------|--------------------------------------------------------------------|-------|------|------------------|-------|--------|----------|
!! | Model          | GFS_control_type_instance        | Fortran DDT containing FV3-GFS model control parameters            | DDT   |    0 | GFS_control_type |       | in     | F        |
!! | mpirank        | mpi_rank                         | current MPI rank                                                   | index |    0 | integer          |       | in     | F        |
!! | mpiroot        | mpi_root                         | master MPI rank                                                    | index |    0 | integer          |       | in     | F        |
!! | mpicomm        | mpi_comm                         | MPI communicator                                                   | index |    0 | integer          |       | in     | F        |
!! | errmsg         | ccpp_error_message               | error message for error handling in CCPP                           | none  |    0 | character        | len=* | out    | F        |
!! | errflg         | ccpp_error_flag                  | error flag for error handling in CCPP                              | flag  |    0 | integer          |       | out    | F        |
!! | lw_cloud_props | coefficients_for_lw_cloud_optics | DDT containing spectral information for RRTMGP LW radiation scheme | DDT   |    0 | ty_cloud_optics  |       | out    | F        |
!!
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_init()
  ! #########################################################################################
  subroutine rrtmgp_lw_cloud_optics_init(Model, mpicomm, mpirank, mpiroot, lw_cloud_props,  &
       errmsg, errflg)
#ifdef MPI
    use mpi
#endif

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model          ! DDT containing model control parameters
    integer,intent(in) :: &
         mpicomm,     & ! MPI communicator
         mpirank,     & ! Current MPI rank
         mpiroot        ! Master MPI rank
 
    ! Outputs
    type(ty_cloud_optics),intent(out) :: &
         lw_cloud_props ! DDT containing spectral information for RRTMGP LW radiation scheme
    character(len=*), intent(out) :: &
         errmsg         ! Error message
    integer,          intent(out) :: &
         errflg         ! Error code

    ! Variables that will be passed to cloud_optics%load()
    real(kind_phys) :: &
         radliq_lwr,                      & ! used by RRTMGP cloud optics 
         radliq_upr,                      & ! used by RRTMGP cloud optics 
         radliq_fac,                      & ! used by RRTMGP cloud optics 
         radice_lwr,                      & ! used by RRTMGP cloud optics 
         radice_upr,                      & ! used by RRTMGP cloud optics 
         radice_fac                         ! used by RRTMGP cloud optics 
    real(kind_phys), dimension(:), allocatable :: &
         pade_sizereg_extliq,             & ! used by RRTMGP cloud optics 
         pade_sizereg_ssaliq,             & ! used by RRTMGP cloud optics 
         pade_sizereg_asyliq,             & ! used by RRTMGP cloud optics 
         pade_sizereg_extice,             & ! used by RRTMGP cloud optics 
         pade_sizereg_ssaice,             & ! used by RRTMGP cloud optics 
         pade_sizereg_asyice                ! used by RRTMGP cloud optics 
    real(kind_phys), dimension(:,:), allocatable :: &
         lut_extliq,                      & ! used by RRTMGP cloud optics 
         lut_ssaliq,                      & ! used by RRTMGP cloud optics 
         lut_asyliq,                      & ! used by RRTMGP cloud optics 
         band_lims_cldy                     ! used by RRTMGP cloud optics 

    real(kind_phys), dimension(:,:,:), allocatable :: &
         lut_extice,                      & ! used by RRTMGP cloud optics 
         lut_ssaice,                      & ! used by RRTMGP cloud optics 
         lut_asyice,                      & ! used by RRTMGP cloud optics 
         pade_extliq,                     & ! used by RRTMGP cloud optics 
         pade_ssaliq,                     & ! used by RRTMGP cloud optics 
         pade_asyliq                        ! used by RRTMGP cloud optics 
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         pade_extice,                     & ! used by RRTMGP cloud optics 
         pade_ssaice,                     & ! used by RRTMGP cloud optics 
         pade_asyice                        ! used by RRTMGP cloud optics 
    ! Dimensions
    integer :: &
         nbandLWcldy,                     & ! used by RRTMGP cloud optics 
         nsize_liq,                       & ! used by RRTMGP cloud optics 
         nsize_ice,                       & ! used by RRTMGP cloud optics 
         nsizereg,                        & ! used by RRTMGP cloud optics 
         ncoeff_ext,                      & ! used by RRTMGP cloud optics 
         ncoeff_ssa_g,                    & ! used by RRTMGP cloud optics 
         nbound,                          & ! used by RRTMGP cloud optics  
         npairsLWcldy                       ! used by RRTMGP cloud optics 

    ! Local variables
    integer :: dimID,varID,status,igpt,iGas,ij,ierr,ncid_lw_clds
    integer,dimension(:),allocatable :: temp1,temp2,temp3,temp4,temp_log_array1,&
    	temp_log_array2, temp_log_array3, temp_log_array4
    character(len=264) :: lw_cloud_props_file
    integer,parameter :: max_strlen=256

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the gfs_physics_nml (scm/src/GFS_typedefs.F90)
    lw_cloud_props_file = trim(Model%rrtmgp_root)//trim(Model%lw_file_clouds)
    ! Read dimensions for k-distribution fields (only on master processor(0))
    if (mpirank .eq. mpiroot) then
       if(nf90_open(trim(lw_cloud_props_file), NF90_WRITE, ncid_lw_clds) == NF90_NOERR) then
          status = nf90_inq_dimid(ncid_lw_clds, 'nband', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=nbandLWcldy)
          status = nf90_inq_dimid(ncid_lw_clds, 'nrghice', dimid)
          status = nf90_inquire_dimension(ncid_lw_clds, dimid, len=nrghice_lw)
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
    if (Model%rrtmgp_cld_optics .eq. 1 .or. Model%rrtmgp_cld_optics .eq. 2) then
       call MPI_BCAST(nbandLWcldy,  1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nrghice_lw,   1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_liq,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_ice,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsizereg,     1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ext,   1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ssa_g, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nbound,       1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(npairsLWcldy, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    endif
#endif

    if (Model%rrtmgp_cld_optics .eq. 1) then
       allocate(lut_extliq(nsize_liq, nBandLWcldy))
       allocate(lut_ssaliq(nsize_liq, nBandLWcldy))
       allocate(lut_asyliq(nsize_liq, nBandLWcldy))
       allocate(lut_extice(nsize_ice, nBandLWcldy, nrghice_lw))
       allocate(lut_ssaice(nsize_ice, nBandLWcldy, nrghice_lw))
       allocate(lut_asyice(nsize_ice, nBandLWcldy, nrghice_lw))
       allocate(band_lims_cldy(2, nBandLWcldy))
    endif
    if (Model%rrtmgp_cld_optics .eq. 2) then
       allocate(pade_extliq(nbandLWcldy, nsizereg,  ncoeff_ext ))
       allocate(pade_ssaliq(nbandLWcldy, nsizereg,  ncoeff_ssa_g))
       allocate(pade_asyliq(nbandLWcldy, nsizereg,  ncoeff_ssa_g))
       allocate(pade_extice(nbandLWcldy, nsizereg,  ncoeff_ext,   nrghice_lw))
       allocate(pade_ssaice(nbandLWcldy, nsizereg,  ncoeff_ssa_g, nrghice_lw))
       allocate(pade_asyice(nbandLWcldy, nsizereg,  ncoeff_ssa_g, nrghice_lw))
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
       if (Model%rrtmgp_cld_optics .eq. 1) then
          !
          if(nf90_open(trim(lw_cloud_props_file), NF90_WRITE, ncid_lw_clds) == NF90_NOERR) then
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
       if (Model%rrtmgp_cld_optics .eq. 2) then
          !
          if(nf90_open(trim(lw_cloud_props_file), NF90_WRITE, ncid_lw_clds) == NF90_NOERR) then
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
    if (Model%rrtmgp_cld_optics .eq. 1) then
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
    if (Model%rrtmgp_cld_optics .eq. 2) then
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

    ! Load tables data for RRTMGP cloud-optics  
    if (Model%rrtmgp_cld_optics .eq. 1) then
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%set_ice_roughness(nrghice_lw))
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%load(band_lims_cldy, &
            radliq_lwr, radliq_upr, radliq_fac, radice_lwr, radice_upr, radice_fac,    &
            lut_extliq, lut_ssaliq, lut_asyliq, lut_extice, lut_ssaice, lut_asyice))
    endif
    if (Model%rrtmgp_cld_optics .eq. 2) then
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%set_ice_roughness(nrghice_lw))
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%load(band_lims_cldy, &
            pade_extliq, pade_ssaliq, pade_asyliq, pade_extice, pade_ssaice,           &
            pade_asyice, pade_sizereg_extliq, pade_sizereg_ssaliq, pade_sizereg_asyliq,&
            pade_sizereg_extice, pade_sizereg_ssaice, pade_sizereg_asyice))
    endif
  end subroutine rrtmgp_lw_cloud_optics_init


!! \section arg_table_rrtmgp_lw_cloud_optics_run Argument Table
!! | local_name            | standard_name                                       | long_name                                                                    | units   | rank | type                  |    kind   | intent | optional |
!! |-----------------------|-----------------------------------------------------|------------------------------------------------------------------------------|---------|------|-----------------------|-----------|--------|----------|
!! | Model                 | GFS_control_type_instance                           | Fortran DDT containing FV3-GFS model control parameters                      | DDT     |    0 | GFS_control_type      |           | in     | F        |
!! | ncol                  | horizontal_loop_extent                              | horizontal dimension                                                         | count   |    0 | integer               |           | in     | F        |
!! | ngpts_lw              | number_of_spectral_points_for_LW_calculation        | Number of spectral points for LW RRTMGP calculation                          | none    |    0 | integer               |           | in     | F        |
!! | p_lay                 | air_pressure_at_layer_for_RRTMGP_in_hPa             | air pressure layer                                                           | hPa     |    2 | real                  | kind_phys | in     | F        |
!! | t_lay                 | air_temperature_at_layer_for_RRTMGP                 | air temperature layer                                                        | K       |    2 | real                  | kind_phys | in     | F        |
!! | p_lev                 | air_pressure_at_interface_for_RRTMGP_in_hPa         | air pressure level                                                           | hPa     |    2 | real                  | kind_phys | in     | F        |
!! | cld_frac              | total_cloud_fraction                                | layer total cloud fraction                                                   | frac    |    2 | real                  | kind_phys | in     | F        |
!! | cld_lwp               | cloud_liquid_water_path                             | layer cloud liquid water path                                                | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_reliq             | mean_effective_radius_for_liquid_cloud              | mean effective radius for liquid cloud                                       | micron  |    2 | real                  | kind_phys | in     | F        |
!! | cld_iwp               | cloud_ice_water_path                                | layer cloud ice water path                                                   | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_reice             | mean_effective_radius_for_ice_cloud                 | mean effective radius for ice cloud                                          | micron  |    2 | real                  | kind_phys | in     | F        |
!! | cld_swp               | cloud_snow_water_path                               | layer cloud snow water path                                                  | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_resnow            | mean_effective_radius_for_snow_flake                | mean effective radius for snow cloud                                         | micron  |    2 | real                  | kind_phys | in     | F        |
!! | cld_rwp               | cloud_rain_water_path                               | layer cloud rain water path                                                  | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_rerain            | mean_effective_radius_for_rain_drop                 | mean effective radius for rain cloud                                         | micron  |    2 | real                  | kind_phys | in     | F        |
!! | icseed_lw             | seed_random_numbers_sw                              | seed for random number generation for shortwave radiation                    | none    |    1 | integer               |           | in     | F        |
!! | aerosols              | aerosol_optical_properties_for_longwave_bands_01-16 | aerosol optical properties for longwave bands 01-16                          | various |    4 | real                  | kind_phys | in     | F        |
!! | lw_cloud_props        | coefficients_for_lw_cloud_optics                    | DDT containing spectral information for cloudy RRTMGP LW radiation scheme    | DDT     |    0 | ty_cloud_optics       |           | in     | F        |
!! | lw_gas_props          | coefficients_for_lw_gas_optics                      | DDT containing spectral information for RRTMGP LW radiation scheme           | DDT     |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | optical_props_clouds  | longwave_optical_properties_for_cloudy_atmosphere   | Fortran DDT containing RRTMGP optical properties                             | DDT     |    0 | ty_optical_props_1scl |           | out    | F        |
!! | optical_props_aerosol | longwave_optical_properties_for_aerosols            | Fortran DDT containing RRTMGP optical properties                             | DDT     |    0 | ty_optical_props_1scl |           | out    | F        |
!! | cldtaulw              | cloud_optical_depth_layers_at_10mu_band             | approx 10mu band layer cloud optical depth                                   | none    |    2 | real                  | kind_phys | out    | F        |
!! | errmsg                | ccpp_error_message                                  | error message for error handling in CCPP                                     | none    |    0 | character             | len=*     | out    | F        |
!! | errflg                | ccpp_error_flag                                     | error flag for error handling in CCPP                                        | flag    |    0 | integer               |           | out    | F        |
!!
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_run()
  ! #########################################################################################
  subroutine rrtmgp_lw_cloud_optics_run(Model, ncol, ngpts_lw, icseed_lw, p_lay, t_lay, p_lev, cld_frac,       &
       cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp, cld_rerain,    &
       aerosols, lw_cloud_props, lw_gas_props,                              &
       optical_props_clouds, optical_props_aerosol, cldtaulw, errmsg, errflg)
    
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model               ! DDT containing FV3-GFS model control parameters 
    integer, intent(in) :: &
         ncol,             & ! Number of horizontal gridpoints
         ngpts_lw            ! Number of spectral points
    integer,intent(in),dimension(ncol) :: &
         icseed_lw           ! auxiliary special cloud related array when module 
                             ! variable isubclw=2, it provides permutation seed 
                             ! for each column profile that are used for generating 
                             ! random numbers. when isubclw /=2, it will not be used.
    real(kind_phys), dimension(ncol,model%levs), intent(in) :: &
         p_lay,            & ! Pressure @ model layer-centers         (hPa)
         t_lay               ! Temperature                            (K)
    real(kind_phys), dimension(ncol,model%levs+1), intent(in) :: &
         p_lev               ! Pressure @ model layer-interfaces      (hPa)
    real(kind_phys), dimension(ncol,model%levs),intent(in) :: &
         cld_frac,         & ! Total cloud fraction by layer
         cld_lwp,          & ! Cloud liquid water path
         cld_reliq,        & ! Cloud liquid effective radius
         cld_iwp,          & ! Cloud ice water path
         cld_reice,        & ! Cloud ice effective radius
         cld_swp,          & ! Cloud snow water path       (used only fro RRTMG scheme)
         cld_resnow,       & ! Cloud snow effective radius (used only fro RRTMG scheme)
         cld_rwp,          & ! Cloud rain water path       (used only fro RRTMG scheme)
         cld_rerain          ! Cloud rain effective radius (used only fro RRTMG scheme)
    type(ty_cloud_optics),intent(in) :: &
         lw_cloud_props      !
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props
    real(kind_phys), intent(in),dimension(ncol, model%levs, lw_cloud_props%get_nband(),3) :: &
         aerosols            !
    real(kind_phys), dimension(ncol,Model%levs), intent(out) :: &
         cldtaulw            ! approx 10.mu band layer cloud optical depth  

    ! Outputs
    type(ty_optical_props_1scl),intent(out) :: &
         optical_props_clouds, & !
         optical_props_aerosol   !
    integer, intent(out) :: &
         errflg                  !
    character(len=*), intent(out) :: &
         errmsg                  !

    ! Local variables
    integer :: iCol
    integer,dimension(ncol) :: ipseed_lw
    logical,dimension(ncol,model%levs) :: liqmask, icemask
    type(ty_optical_props_1scl) :: optical_props_cloudsByBand
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(ngpts_lw,model%levs,ncol) :: rng3D
    real(kind_phys), dimension(ngpts_lw*model%levs) :: rng1D
    logical, dimension(ncol,model%levs,ngpts_lw) :: cldfracMCICA
    real(kind_phys), dimension(ncol,model%levs,lw_cloud_props%get_nband()) :: &
         tau_cld

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. Model%lslwr) return

    ! #######################################################################################
    ! Change random number seed value for each radiation invocation (isubclw =1 or 2).
    ! #######################################################################################
    if(isubclw == 1) then      ! advance prescribed permutation seed
       do iCol = 1, nCol
          ipseed_lw(iCol) = ipsdlw0 + iCol
       enddo
    elseif (isubclw == 2) then ! use input array of permutaion seeds
       do iCol = 1, nCol
          ipseed_lw(iCol) = icseed_lw(iCol)
       enddo
    endif
    
    ! #######################################################################################
    ! Compute ice/liquid cloud masks, needed by rrtmgp_cloud_optics
    ! #######################################################################################
    liqmask = (cld_frac .gt. 0 .and. cld_lwp .gt. 0)
    icemask = (cld_frac .gt. 0 .and. cld_iwp .gt. 0)

    ! #######################################################################################
    ! Allocate space for RRTMGP DDTs containing cloud and aerosol radiative properties
    ! #######################################################################################
    ! Cloud optics [nCol,model%levs,nBands]
    call check_error_msg('rrtmgp_lw_cloud_optics_run',optical_props_cloudsByBand%alloc_1scl(&
         ncol, model%levs, lw_cloud_props%get_band_lims_wavenumber()))
    ! Aerosol optics [nCol,model%levs,nBands]
    call check_error_msg('rrtmgp_lw_cloud_optics_run',optical_props_aerosol%alloc_1scl(     &
         ncol, model%levs, lw_cloud_props%get_band_lims_wavenumber()))
    ! Cloud optics [nCol,model%levs,nGpts]
    call check_error_msg('rrtmgp_lw_cloud_optics_run',optical_props_clouds%alloc_1scl(      &
         ncol, model%levs, lw_gas_props))

    ! #######################################################################################
    ! Copy aerosol optical information to RRTMGP DDT
    ! #######################################################################################
    optical_props_aerosol%tau = aerosols(:,:,:,1) * (1. - aerosols(:,:,:,2))

    ! #######################################################################################
    ! Compute cloud-optics for RTE.
    ! #######################################################################################
    if (Model%rrtmgp_cld_optics .gt. 0) then
       ! i) RRTMGP cloud-optics.
       call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_cloud_props%cloud_optics(&
            ncol,                       & ! IN  - Number of horizontal gridpoints 
            model%levs,                 & ! IN  - Number of vertical layers
            lw_cloud_props%get_nband(), & ! IN  - Number of LW bands
            nrghice_lw,                 & ! IN  - Number of ice-roughness categories
            liqmask,                    & ! IN  - Liquid-cloud mask
            icemask,                    & ! IN  - Ice-cloud mask
            cld_lwp,                    & ! IN  - Cloud liquid water path
            cld_iwp,                    & ! IN  - Cloud ice water path
            cld_reliq,                  & ! IN  - Cloud liquid effective radius
            cld_reice,                  & ! IN  - Cloud ice effective radius
            optical_props_cloudsByBand))  ! OUT - RRTMGP DDT containing cloud radiative properties
                                          !       in each band
    else
       ! ii) RRTMG cloud-optics.
       if (any(cld_frac .gt. 0)) then
          call rrtmg_lw_cloud_optics(ncol, model%levs, lw_cloud_props%get_nband(), cld_lwp,     &
               cld_reliq, cld_iwp, cld_reice, cld_rwp, cld_rerain, cld_swp, cld_resnow,    &
               cld_frac, tau_cld)
          optical_props_cloudsByBand%tau = tau_cld
       endif
    endif

    ! #######################################################################################
    ! Call McICA to generate subcolumns.
    ! #######################################################################################
    ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
    ! and layers. ([nGpts,model%levs,nColumn]-> [nGpts*model%levs]*nColumn)
    do iCol=1,ncol
       call random_setseed(ipseed_lw(icol),rng_stat)
       call random_number(rng1D,rng_stat)
       rng3D(:,:,iCol) = reshape(source = rng1D,shape=[ngpts_lw,model%levs])
    enddo
    
    ! Call McICA
    select case ( iovrlw )
       ! Maximumn-random 
    case(1)
       call check_error_msg('rrtmgp_lw_cloud_optics_run',sampled_mask_max_ran(rng3D,cld_frac,cldfracMCICA))       
    end select
    
    ! Map band optical depth to each g-point using McICA
    call check_error_msg('rrtmgp_lw_cloud_optics_run',draw_samples(cldfracMCICA,optical_props_cloudsByBand,optical_props_clouds))
    
    ! GFS_RRTMGP_POST_RUN() requires the LW optical depth ~10microns
    cldtaulw = optical_props_cloudsByBand%tau(:,:,7)

  end subroutine rrtmgp_lw_cloud_optics_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_finalize()
  ! #########################################################################################
  subroutine rrtmgp_lw_cloud_optics_finalize()
  end subroutine rrtmgp_lw_cloud_optics_finalize

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

end module rrtmgp_lw_cloud_optics
