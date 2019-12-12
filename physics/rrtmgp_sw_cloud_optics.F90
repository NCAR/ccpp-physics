module rrtmgp_sw_cloud_optics
  use machine,                  only: kind_phys
  use GFS_typedefs,             only: GFS_control_type
  use mo_rte_kind,              only: wl
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,          only: ty_cloud_optics
  use physparam,                only: isubcsw, iovrsw
  use mo_optical_props,         only: ty_optical_props_2str
  use mo_rrtmg_sw_cloud_optics, only: rrtmg_sw_cloud_optics   
  use rrtmgp_aux,               only: check_error_msg
  use netcdf

contains

  ! #########################################################################################
  ! SUBROUTINE sw_cloud_optics_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_optics_init
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_sw_cloud_optics_init(Model,mpicomm, mpirank, mpiroot, sw_cloud_props,   &
        errmsg, errflg)
    use netcdf
!#ifdef MPI
!    use mpi
!#endif

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model          ! DDT containing model control parameters
    integer,intent(in) :: &
         mpicomm,     & ! MPI communicator
         mpirank,     & ! Current MPI rank
         mpiroot        ! Master MPI rank

    ! Outputs
    type(ty_cloud_optics),intent(out) :: &
         sw_cloud_props ! DDT containing spectral information for RRTMGP SW radiation scheme
    character(len=*), intent(out) :: &
         errmsg         ! Error message
    integer,          intent(out) :: &
         errflg         ! Error code

    ! Variables that will be passed to cloud_optics%load()
    real(kind_phys) :: &
         radliq_lwr_sw,                      & !   
         radliq_upr_sw,                      & !   
         radliq_fac_sw,                      & !   
         radice_lwr_sw,                      & !   
         radice_upr_sw,                      & !   
         radice_fac_sw                         !   

    real(kind_phys), dimension(:), allocatable :: &
         pade_sizereg_extliq_sw,             & !   
         pade_sizereg_ssaliq_sw,             & !   
         pade_sizereg_asyliq_sw,             & !   
         pade_sizereg_extice_sw,             & !   
         pade_sizereg_ssaice_sw,             & !   
         pade_sizereg_asyice_sw                !   
    real(kind_phys), dimension(:,:), allocatable :: &
         lut_extliq_sw,                      & !   
         lut_ssaliq_sw,                      & !   
         lut_asyliq_sw,                      & !   
         band_lims_cldy_sw                     !                            

    real(kind_phys), dimension(:,:,:), allocatable :: &
         lut_extice_sw,                      & !   
         lut_ssaice_sw,                      & !   
         lut_asyice_sw,                      & !   
         pade_extliq_sw,                     & !   
         pade_ssaliq_sw,                     & !   
         pade_asyliq_sw                        !   
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         pade_extice_sw,                     & !   
         pade_ssaice_sw,                     & !   
         pade_asyice_sw                        !   
    ! Dimensions (to be broadcast across all processors)
    integer :: &
         nrghice_sw,                         & ! Number of ice-roughness categories in file
         nbandSWcldy_sw,                     & !   
         nsize_liq_sw,                       & !   
         nsize_ice_sw,                       & !   
         nsizereg_sw,                        & !   
         ncoeff_ext_sw,                      & !   
         ncoeff_ssa_g_sw,                    & !   
         nbound_sw,                          & !    
         npairsSWcldy_sw                       !   

    ! Local variables
    integer :: status,ncid_sw_clds,dimid,varID
    character(len=264) :: sw_cloud_props_file
!#ifdef MPI
!    integer :: ierr
!#endif
    ! Initialize
    errmsg = ''
    errflg = 0

    if (Model%rrtmgp_cld_optics .eq. 0) return

    ! Filenames are set in the gfs_physics_nml (scm/src/GFS_typedefs.F90)
    sw_cloud_props_file = trim(Model%rrtmgp_root)//trim(Model%sw_file_clouds)

    ! Read dimensions for k-distribution fields (only on master processor(0))
!    if (mpirank .eq. mpiroot) then
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
 
       ! Check to ensure that number of ice-roughness categories is feasible.
       if (Model%rrtmgp_nrghice .gt. nrghice_sw) then
          errmsg = 'Number of RRTMGP ice-roughness categories requested in namelist file is not allowed'
       endif
!    endif

!    ! Broadcast dimensions to all processors
!#ifdef MPI
!    if (Model%rrtmgp_cld_optics .eq. 1 .or. Model%rrtmgp_cld_optics .eq. 2) then
!       call MPI_BCAST(nbandSWcldy_sw,  1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!       call MPI_BCAST(nrghice_sw,      1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!       call MPI_BCAST(nsize_liq_sw,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!       call MPI_BCAST(nsize_ice_sw,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!       call MPI_BCAST(nsizereg_sw,     1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!       call MPI_BCAST(ncoeff_ext_sw,   1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!       call MPI_BCAST(ncoeff_ssa_g_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!       call MPI_BCAST(nbound_sw,       1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!       call MPI_BCAST(npairsSWcldy_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
!    endif
!#endif

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
!    if (mpirank .eq. mpiroot) then
       ! 
       if (Model%rrtmgp_cld_optics .eq. 1) then
          write (*,*) 'Reading RRTMGP shortwave cloud data (LUT) ... '
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
          write (*,*) 'Reading RRTMGP shortwave cloud data (PADE) ... '
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
    !endif

    ! Broadcast arrays to all processors
!#ifdef MPI
!    if (Model%rrtmgp_cld_optics .eq. 1) then
!       write (*,*) 'Broadcasting RRTMGP shortwave cloud data (LUT) ... '
!#ifndef SINGLE_PREC
!       call MPI_BCAST(radliq_lwr_sw,           1,                            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radliq_upr_sw,           1,                            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radliq_fac_sw,           1,                            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radice_lwr_sw,           1,                            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radice_upr_sw,           1,                            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radice_fac_sw,           1,                            MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_extliq_sw,           size(lut_extliq_sw),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_ssaliq_sw,           size(lut_ssaliq_sw),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_asyliq_sw,           size(lut_asyliq_sw),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_extice_sw,           size(lut_extice_sw),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_ssaice_sw,           size(lut_ssaice_sw),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_asyice_sw,           size(lut_asyice_sw),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(band_lims_cldy_sw,       size(band_lims_cldy_sw),      MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)    
!#else
!       call MPI_BCAST(radliq_lwr_sw,           1,                            MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radliq_upr_sw,           1,                            MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radliq_fac_sw,           1,                            MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radice_lwr_sw,           1,                            MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radice_upr_sw,           1,                            MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(radice_fac_sw,           1,                            MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_extliq_sw,           size(lut_extliq_sw),          MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_ssaliq_sw,           size(lut_ssaliq_sw),          MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_asyliq_sw,           size(lut_asyliq_sw),          MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_extice_sw,           size(lut_extice_sw),          MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_ssaice_sw,           size(lut_ssaice_sw),          MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(lut_asyice_sw,           size(lut_asyice_sw),          MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(band_lims_cldy_sw,       size(band_lims_cldy_sw),      MPI_REAL,   mpiroot, mpicomm, ierr)
!#endif 
!    endif
!    if (Model%rrtmgp_cld_optics .eq. 2) then
!       write (*,*) 'Broadcasting RRTMGP shortwave cloud data (PADE) ... '
!#ifndef SINGLE_PREC
!       call MPI_BCAST(pade_extliq_sw,          size(pade_extliq_sw),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_ssaliq_sw,          size(pade_ssaliq_sw),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_asyliq_sw,          size(pade_asyliq_sw),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_extice_sw,          size(pade_extice_sw),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_ssaice_sw,          size(pade_ssaice_sw),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_asyice_sw,          size(pade_asyice_sw),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_extliq_sw,  size(pade_sizereg_extliq_sw), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_ssaliq_sw,  size(pade_sizereg_ssaliq_sw), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_asyliq_sw,  size(pade_sizereg_asyliq_sw), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_extice_sw,  size(pade_sizereg_extice_sw), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_ssaice_sw,  size(pade_sizereg_ssaice_sw), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_asyice_sw,  size(pade_sizereg_asyice_sw), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(band_lims_cldy_sw,       size(band_lims_cldy_sw),      MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)    
!#else
!       call MPI_BCAST(pade_extliq_sw,          size(pade_extliq_sw),         MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_ssaliq_sw,          size(pade_ssaliq_sw),         MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_asyliq_sw,          size(pade_asyliq_sw),         MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_extice_sw,          size(pade_extice_sw),         MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_ssaice_sw,          size(pade_ssaice_sw),         MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_asyice_sw,          size(pade_asyice_sw),         MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_extliq_sw,  size(pade_sizereg_extliq_sw), MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_ssaliq_sw,  size(pade_sizereg_ssaliq_sw), MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_asyliq_sw,  size(pade_sizereg_asyliq_sw), MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_extice_sw,  size(pade_sizereg_extice_sw), MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_ssaice_sw,  size(pade_sizereg_ssaice_sw), MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(pade_sizereg_asyice_sw,  size(pade_sizereg_asyice_sw), MPI_REAL,   mpiroot, mpicomm, ierr)
!       call MPI_BCAST(band_lims_cldy_sw,       size(band_lims_cldy_sw),      MPI_REAL,   mpiroot, mpicomm, ierr)
!#endif
!    endif
!#endif

    ! Load tables data for RRTMGP cloud-optics  
    if (Model%rrtmgp_cld_optics .eq. 1) then
       call check_error_msg('sw_cloud_optics_init',sw_cloud_props%set_ice_roughness(Model%rrtmgp_nrghice))
       call check_error_msg('sw_cloud_optics_init',sw_cloud_props%load(band_lims_cldy_sw,   &
            radliq_lwr_sw, radliq_upr_sw, radliq_fac_sw, radice_lwr_sw, radice_upr_sw,      &
            radice_fac_sw, lut_extliq_sw, lut_ssaliq_sw, lut_asyliq_sw, lut_extice_sw,      &
            lut_ssaice_sw, lut_asyice_sw))
    endif
    if (Model%rrtmgp_cld_optics .eq. 2) then
       call check_error_msg('sw_cloud_optics_init',sw_cloud_props%set_ice_roughness(Model%rrtmgp_nrghice))
       call check_error_msg('sw_cloud_optics_init', sw_cloud_props%load(band_lims_cldy_sw,  &
            pade_extliq_sw, pade_ssaliq_sw, pade_asyliq_sw, pade_extice_sw, pade_ssaice_sw, &
            pade_asyice_sw, pade_sizereg_extliq_sw, pade_sizereg_ssaliq_sw,                 &
            pade_sizereg_asyliq_sw, pade_sizereg_extice_sw, pade_sizereg_ssaice_sw,         &
            pade_sizereg_asyice_sw))
    endif
  end subroutine rrtmgp_sw_cloud_optics_init

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_optics_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_optics_run
!! \htmlinclude rrtmgp_sw_cloud_optics.html
!!
  subroutine rrtmgp_sw_cloud_optics_run(Model, ncol, cld_frac, cld_lwp, cld_reliq, cld_iwp, &
       cld_reice, cld_swp, cld_resnow, cld_rwp, cld_rerain, sw_cloud_props, sw_gas_props,   &
       nday, idxday, sw_optical_props_cloudsByBand, cldtausw, errmsg, errflg)
    
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model
    integer, intent(in) :: &
         ncol,             & ! Number of horizontal gridpoints
         nday                ! Number of daylit points.
    integer,intent(in),dimension(ncol) :: &
         idxday              ! Indices for daylit points.
    real(kind_phys), dimension(ncol,model%levs),intent(in) :: &
         cld_frac,         & ! Total cloud fraction by layer
         cld_lwp,          & ! Cloud liquid water path
         cld_reliq,        & ! Cloud liquid effective radius
         cld_iwp,          & ! Cloud ice water path
         cld_reice,        & ! Cloud ice effective radius
         cld_swp,          & ! Cloud snow water path
         cld_resnow,       & ! Cloud snow effective radius
         cld_rwp,          & ! Cloud rain water path
         cld_rerain          ! Cloud rain effective radius
    type(ty_cloud_optics),intent(in) :: &
         sw_cloud_props      ! RRTMGP DDT: 
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props        ! RRTMGP DDT: K-distribution data

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                     ! Error message
    integer,          intent(out) :: &
         errflg                     ! Error code
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_cloudsByBand  ! RRTMGP DDT: Shortwave optical properties (cloudy atmosphere)
    real(kind_phys), dimension(ncol,Model%levs), intent(out) :: &
         cldtausw                   ! approx 10.mu band layer cloud optical depth  

    ! Local variables
    logical,dimension(nday,model%levs) :: liqmask, icemask
    real(kind_phys), dimension(nday,model%levs,sw_gas_props%get_nband()) :: &
         tau_cld, ssa_cld, asy_cld

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. Model%lsswr) return
    if (nDay .gt. 0) then
       
       ! #######################################################################################
       ! Compute ice/liquid cloud masks, needed by rrtmgp_cloud_optics
       ! #######################################################################################
       liqmask = (cld_frac(idxday(1:nday),:) .gt. 0 .and. cld_lwp(idxday(1:nday),:) .gt. 0)
       icemask = (cld_frac(idxday(1:nday),:) .gt. 0 .and. cld_iwp(idxday(1:nday),:) .gt. 0)
       
       ! #######################################################################################
       ! Allocate space for RRTMGP DDTs containing cloud radiative properties
       ! #######################################################################################
       ! Cloud optics [nday,model%levs,nBands]
       call check_error_msg('rrtmgp_sw_cloud_optics_run',sw_optical_props_cloudsByBand%alloc_2str(&
            nday, model%levs, sw_gas_props%get_band_lims_wavenumber()))
 
       ! #######################################################################################
       ! Compute cloud-optics for RTE.
       ! #######################################################################################
       if (Model%rrtmgp_cld_optics .gt. 0) then
          ! RRTMGP cloud-optics.
          call check_error_msg('rrtmgp_sw_cloud_optics_run',sw_cloud_props%cloud_optics(&
               nday,                         & ! IN  - Number of daylit gridpoints
               model%levs,                   & ! IN  - Number of vertical layers
               sw_cloud_props%get_nband(),   & ! IN  - Number of SW bands
               Model%rrtmgp_nrghice,         & ! IN  - Number of ice-roughness categories
               liqmask,                      & ! IN  - Liquid-cloud mask
               icemask,                      & ! IN  - Ice-cloud mask
               cld_lwp(idxday(1:nday),:),    & ! IN  - Cloud liquid water path
               cld_iwp(idxday(1:nday),:),    & ! IN  - Cloud ice water path
               cld_reliq(idxday(1:nday),:),  & ! IN  - Cloud liquid effective radius
               cld_reice(idxday(1:nday),:),  & ! IN  - Cloud ice effective radius
               sw_optical_props_cloudsByBand)) ! OUT - RRTMGP DDT: Shortwave optical properties, 
                                               !       in each band (tau,ssa,g)
       else
          ! RRTMG cloud-optics
          if (any(cld_frac .gt. 0)) then
             sw_optical_props_cloudsByBand%tau(:,:,:) = 0._kind_phys
             sw_optical_props_cloudsByBand%ssa(:,:,:) = 0._kind_phys
             sw_optical_props_cloudsByBand%g(:,:,:)   = 0._kind_phys
             call rrtmg_sw_cloud_optics(nday, model%levs, sw_gas_props%get_nband(), &
                  cld_lwp(idxday(1:nday),:), cld_reliq(idxday(1:nday),:),           &
                  cld_iwp(idxday(1:nday),:), cld_reice(idxday(1:nday),:),           &
                  cld_rwp(idxday(1:nday),:), cld_rerain(idxday(1:nday),:),          &
                  cld_swp(idxday(1:nday),:), cld_resnow(idxday(1:nday),:),          &
                  cld_frac(idxday(1:nday),:), tau_cld, ssa_cld, asy_cld)
             sw_optical_props_cloudsByBand%tau(:,:,:) = tau_cld
             sw_optical_props_cloudsByBand%ssa(:,:,:) = ssa_cld
             sw_optical_props_cloudsByBand%g(:,:,:)   = asy_cld
          endif
       endif

       ! GFS_RRTMGP_POST_RUN() requires the SW optical depth ~0.55microns
       cldtausw(idxday(1:nDay),:) = sw_optical_props_cloudsByBand%tau(:,:,11)    
    endif

  end subroutine rrtmgp_sw_cloud_optics_run

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_optics_finalize()
  ! #########################################################################################  
  subroutine rrtmgp_sw_cloud_optics_finalize()
  end subroutine rrtmgp_sw_cloud_optics_finalize 

end module rrtmgp_sw_cloud_optics
