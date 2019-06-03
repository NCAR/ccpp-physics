module rrtmgp_aux
  use machine,               only: kind_phys
  use GFS_typedefs,          only: GFS_control_type
  use mo_rte_kind,           only: wl
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,       only: ty_cloud_optics
  use mo_gas_concentrations, only: ty_gas_concs
  use netcdf

  ! Parameters
  integer,parameter :: nGases = 6
  character(len=3),parameter, dimension(nGases) :: &
       active_gases = (/ 'h2o', 'co2', 'o3 ', 'n2o', 'ch4', 'o2 '/)
  integer :: nrghice_lw, nrghice_sw, ipsdlw0, ipsdsw0
  
contains
  
  subroutine rrtmgp_aux_init()
  end subroutine rrtmgp_aux_init
  subroutine rrtmgp_aux_run()
  end subroutine rrtmgp_aux_run
  subroutine rrtmgp_aux_finalize()
  end subroutine rrtmgp_aux_finalize

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_init()
  ! #########################################################################################
  subroutine lw_cloud_optics_init(Model, mpicomm, mpirank, mpiroot, lw_cloud_props,  &
       errmsg, errflg)
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
    type(ty_cloud_optics),intent(inout) :: &
         lw_cloud_props
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer,          intent(out) :: &
         errflg     ! Error code

    ! Variables that will be passed to cloud_optics%load()
    real(kind_phys) :: &
         radliq_lwr,                      & ! used by RRTGMP cloud optics 
         radliq_upr,                      & ! used by RRTGMP cloud optics 
         radliq_fac,                      & ! used by RRTGMP cloud optics 
         radice_lwr,                      & ! used by RRTGMP cloud optics 
         radice_upr,                      & ! used by RRTGMP cloud optics 
         radice_fac                         ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:), allocatable :: &
         pade_sizereg_extliq,             & ! used by RRTGMP cloud optics 
         pade_sizereg_ssaliq,             & ! used by RRTGMP cloud optics 
         pade_sizereg_asyliq,             & ! used by RRTGMP cloud optics 
         pade_sizereg_extice,             & ! used by RRTGMP cloud optics 
         pade_sizereg_ssaice,             & ! used by RRTGMP cloud optics 
         pade_sizereg_asyice                ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:,:), allocatable :: &
         lut_extliq,                      & ! used by RRTGMP cloud optics 
         lut_ssaliq,                      & ! used by RRTGMP cloud optics 
         lut_asyliq,                      & ! used by RRTGMP cloud optics 
         band_lims_cldy                     ! used by RRTGMP cloud optics 

    real(kind_phys), dimension(:,:,:), allocatable :: &
         lut_extice,                      & ! used by RRTGMP cloud optics 
         lut_ssaice,                      & ! used by RRTGMP cloud optics 
         lut_asyice,                      & ! used by RRTGMP cloud optics 
         pade_extliq,                     & ! used by RRTGMP cloud optics 
         pade_ssaliq,                     & ! used by RRTGMP cloud optics 
         pade_asyliq                        ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         pade_extice,                     & ! used by RRTGMP cloud optics 
         pade_ssaice,                     & ! used by RRTGMP cloud optics 
         pade_asyice                        ! used by RRTGMP cloud optics 
    ! Dimensions
    integer :: &
         nbandLWcldy,                     & ! used by RRTGMP cloud optics 
         nsize_liq,                       & ! used by RRTGMP cloud optics 
         nsize_ice,                       & ! used by RRTGMP cloud optics 
         nsizereg,                        & ! used by RRTGMP cloud optics 
         ncoeff_ext,                      & ! used by RRTGMP cloud optics 
         ncoeff_ssa_g,                    & ! used by RRTGMP cloud optics 
         nbound,                          & ! used by RRTGMP cloud optics  
         npairsLWcldy                       ! used by RRTGMP cloud optics 

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

    ! Load tables data for RRTGMP cloud-optics  
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
  end subroutine lw_cloud_optics_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_gas_optics_init()
  ! #########################################################################################
  subroutine lw_gas_optics_init(Model, mpicomm, mpirank, mpiroot, lw_gas_props,      &
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
         lw_gas_props
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer,          intent(out) :: &
         errflg     ! Error code

    ! Variables that will be passed to gas_optics%load()
    type(ty_gas_concs) :: &
         gas_concentrations
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
         temp_ref_t                         ! used by RRTGMP gas optics 
    real(kind_phys), dimension(:), allocatable :: &
         press_ref,                       & ! used by RRTGMP gas optics 
         temp_ref                           ! used by RRTGMP gas optics 
    real(kind_phys), dimension(:,:), allocatable :: &
         band_lims,                       & ! used by RRTGMP gas optics 
         totplnk                            ! used by RRTGMP gas optics 
    real(kind_phys), dimension(:,:,:), allocatable :: &
         vmr_ref,                         & ! used by RRTGMP gas optics 
         kminor_lower,                    & ! used by RRTGMP gas optics 
         kminor_upper,                    & ! used by RRTGMP gas optics 
         rayl_lower,                      & ! used by RRTGMP gas optics 
         rayl_upper                         ! used by RRTGMP gas optics 
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         kmajor,                          & ! used by RRTGMP gas optics 
         planck_frac                        ! used by RRTGMP gas optics  
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
         ncontributors_upper                ! used by RRTGMP gas optics 

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
    do iGas=1,nGases
       call check_error_msg('lw_gas_optics_init',gas_concentrations%set_vmr(active_gases(iGas), 0._kind_phys))
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

  end subroutine lw_gas_optics_init

  ! #########################################################################################
  ! SUBROUTINE sw_gas_optics_init
  ! #########################################################################################
  subroutine sw_gas_optics_init(Model,mpicomm, mpirank, mpiroot, sw_gas_props,       &
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
         sw_gas_props

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
    do iGas=1,nGases
       call check_error_msg('sw_gas_optics_init',gas_concentrations%set_vmr(active_gases(iGas), 0._kind_phys))
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

  end subroutine sw_gas_optics_init

  ! #########################################################################################
  ! SUBROUTINE sw_cloud_optics_init
  ! #########################################################################################
  subroutine sw_cloud_optics_init(Model,mpicomm, mpirank, mpiroot, sw_cloud_props,   &
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
    type(ty_cloud_optics),intent(inout) :: &
         sw_cloud_props

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer,          intent(out) :: &
         errflg     ! Error code

    ! Variables that will be passed to gas_optics%load()
    real(kind_phys) :: &
         radliq_lwr_sw,                      & ! used by RRTGMP cloud optics 
         radliq_upr_sw,                      & ! used by RRTGMP cloud optics 
         radliq_fac_sw,                      & ! used by RRTGMP cloud optics 
         radice_lwr_sw,                      & ! used by RRTGMP cloud optics 
         radice_upr_sw,                      & ! used by RRTGMP cloud optics 
         radice_fac_sw                         ! used by RRTGMP cloud optics 

    real(kind_phys), dimension(:), allocatable :: &
         pade_sizereg_extliq_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_ssaliq_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_asyliq_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_extice_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_ssaice_sw,             & ! used by RRTGMP cloud optics 
         pade_sizereg_asyice_sw                ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:,:), allocatable :: &
         lut_extliq_sw,                      & ! used by RRTGMP cloud optics 
         lut_ssaliq_sw,                      & ! used by RRTGMP cloud optics 
         lut_asyliq_sw,                      & ! used by RRTGMP cloud optics 
         band_lims_cldy_sw                     ! used by RRTGMP cloud optics                          

    real(kind_phys), dimension(:,:,:), allocatable :: &
         lut_extice_sw,                      & ! used by RRTGMP cloud optics 
         lut_ssaice_sw,                      & ! used by RRTGMP cloud optics 
         lut_asyice_sw,                      & ! used by RRTGMP cloud optics 
         pade_extliq_sw,                     & ! used by RRTGMP cloud optics 
         pade_ssaliq_sw,                     & ! used by RRTGMP cloud optics 
         pade_asyliq_sw                        ! used by RRTGMP cloud optics 
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         pade_extice_sw,                     & ! used by RRTGMP cloud optics 
         pade_ssaice_sw,                     & ! used by RRTGMP cloud optics 
         pade_asyice_sw                        ! used by RRTGMP cloud optics 
    ! Dimensions (to be broadcast across all processors)
    integer :: &
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
    integer,dimension(:),allocatable :: temp1,temp2,temp3,temp4,temp_log_array1,&
         temp_log_array2, temp_log_array3, temp_log_array4
    character(len=264) :: sw_cloud_props_file

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the gfs_physics_nml (scm/src/GFS_typedefs.F90)
    sw_cloud_props_file = trim(Model%rrtmgp_root)//trim(Model%sw_file_clouds)

    ! Read dimensions for k-distribution fields (only on master processor(0))
    if (mpirank .eq. mpiroot) then
       if(nf90_open(trim(sw_cloud_props_file), NF90_WRITE, ncid_sw_clds) == NF90_NOERR) then
          status = nf90_inq_dimid(ncid_sw_clds, 'nband', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=nbandSWcldy_sw)
          status = nf90_inq_dimid(ncid_sw_clds, 'nrghice', dimid)
          status = nf90_inquire_dimension(ncid_sw_clds, dimid, len=nrghice_sw)
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
    if (Model%rrtmgp_cld_optics .eq. 1 .or. Model%rrtmgp_cld_optics .eq. 2) then
       call MPI_BCAST(nbandSWcldy_sw,  1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nrghice_sw,         1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_liq_sw,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_ice_sw,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsizereg_sw,     1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ext_sw,   1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ssa_g_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nbound_sw,       1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(npairsSWcldy_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    endif
#endif

    if (Model%rrtmgp_cld_optics .eq. 1) then
       allocate(lut_extliq_sw(nsize_liq_sw, nBandSWcldy_sw))
       allocate(lut_ssaliq_sw(nsize_liq_sw, nBandSWcldy_sw))
       allocate(lut_asyliq_sw(nsize_liq_sw, nBandSWcldy_sw))
       allocate(lut_extice_sw(nsize_ice_sw, nBandSWcldy_sw, nrghice_sw))
       allocate(lut_ssaice_sw(nsize_ice_sw, nBandSWcldy_sw, nrghice_sw))
       allocate(lut_asyice_sw(nsize_ice_sw, nBandSWcldy_sw, nrghice_sw))
       allocate(band_lims_cldy_sw(2, nBandSWcldy_sw))
    endif
    if (Model%rrtmgp_cld_optics .eq. 2) then
       allocate(pade_extliq_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ext_sw ))
       allocate(pade_ssaliq_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ssa_g_sw))
       allocate(pade_asyliq_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ssa_g_sw))
       allocate(pade_extice_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ext_sw,   nrghice_sw))
       allocate(pade_ssaice_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ssa_g_sw, nrghice_sw))
       allocate(pade_asyice_sw(nbandSWcldy_sw, nsizereg_sw,  ncoeff_ssa_g_sw, nrghice_sw))
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
       if (Model%rrtmgp_cld_optics .eq. 1) then
          !
          if(nf90_open(trim(sw_cloud_props_file), NF90_WRITE, ncid_sw_clds) == NF90_NOERR) then
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
       if (Model%rrtmgp_cld_optics .eq. 2) then
          !
          if(nf90_open(trim(sw_cloud_props_file), NF90_WRITE, ncid_sw_clds) == NF90_NOERR) then
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
    if (Model%rrtmgp_cld_optics .eq. 1) then
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
    if (Model%rrtmgp_cld_optics .eq. 2) then
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
    if (Model%rrtmgp_cld_optics .eq. 1) then
       call check_error_msg('sw_cloud_optics_init',sw_cloud_props%set_ice_roughness(nrghice_sw))
       call check_error_msg('sw_cloud_optics_init',sw_cloud_props%load(band_lims_cldy_sw,   &
            radliq_lwr_sw, radliq_upr_sw, radliq_fac_sw, radice_lwr_sw, radice_upr_sw,      &
            radice_fac_sw, lut_extliq_sw, lut_ssaliq_sw, lut_asyliq_sw, lut_extice_sw,      &
            lut_ssaice_sw, lut_asyice_sw))
    endif
    if (Model%rrtmgp_cld_optics .eq. 2) then
       call check_error_msg('sw_cloud_optics_init',sw_cloud_props%set_ice_roughness(nrghice_sw))
       call check_error_msg('sw_cloud_optics_init', sw_cloud_props%load(band_lims_cldy_sw,  &
            pade_extliq_sw, pade_ssaliq_sw, pade_asyliq_sw, pade_extice_sw, pade_ssaice_sw, &
            pade_asyice_sw, pade_sizereg_extliq_sw, pade_sizereg_ssaliq_sw,                 &
            pade_sizereg_asyliq_sw, pade_sizereg_extice_sw, pade_sizereg_ssaice_sw,         &
            pade_sizereg_asyice_sw))
    endif

  end subroutine sw_cloud_optics_init

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
end module rrtmgp_aux
