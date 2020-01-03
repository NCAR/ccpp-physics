module rrtmgp_lw_cloud_optics
  use machine,                  only: kind_phys
  use mo_rte_kind,              only: wl
  use mo_cloud_optics,          only: ty_cloud_optics
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use mo_optical_props,         only: ty_optical_props_1scl
  use mo_rrtmg_lw_cloud_optics, only: rrtmg_lw_cloud_optics   
  use rrtmgp_aux,               only: check_error_msg
  use netcdf

  public rrtmgp_lw_cloud_optics_init, rrtmgp_lw_cloud_optics_run, rrtmgp_lw_cloud_optics_finalize
contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_init()
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_cloud_optics_init
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_lw_cloud_optics_init(cld_optics_scheme, nrghice, rrtmgp_root_dir,       &
       rrtmgp_lw_file_clouds, mpicomm, mpirank, mpiroot, lw_cloud_props, errmsg, errflg)
#ifdef MPI
    use mpi
#endif

    ! Inputs
    integer, intent(in) :: &
         nrghice,            & ! Number of ice-roughness categories
         cld_optics_scheme,  & ! Cloud-optics scheme
         mpicomm,            & ! MPI communicator
         mpirank,            & ! Current MPI rank
         mpiroot               ! Master MPI rank
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,    & ! RTE-RRTMGP root directory
         rrtmgp_lw_file_clouds ! RRTMGP file containing coefficients used to compute clouds optical properties
 
    ! Outputs
    type(ty_cloud_optics),intent(out) :: &
         lw_cloud_props        ! RRTMGP DDT: spectral information for RRTMGP LW radiation scheme
    character(len=*), intent(out) :: &
         errmsg                ! Error message
    integer,          intent(out) :: &
         errflg                ! Error code

    ! Variables that will be passed to cloud_optics%load()
    ! cld_optics_scheme = 1
    real(kind_phys) :: &
         radliq_lwr,          & ! Liquid particle size lower bound for LUT interpolation   
         radliq_upr,          & ! Liquid particle size upper bound for LUT interpolation
         radliq_fac,          & ! Factor for calculating LUT interpolation indices for liquid   
         radice_lwr,          & ! Ice particle size upper bound for LUT interpolation  
         radice_upr,          & ! Ice particle size lower bound for LUT interpolation
         radice_fac             ! Factor for calculating LUT interpolation indices for ice  
    real(kind_phys), dimension(:,:), allocatable :: &
         lut_extliq,          & ! LUT shortwave liquid extinction coefficient  
         lut_ssaliq,          & ! LUT shortwave liquid single scattering albedo   
         lut_asyliq,          & ! LUT shortwave liquid asymmetry parameter  
         band_lims_cldy         ! Beginning and ending wavenumber [cm -1] for each band                           
    real(kind_phys), dimension(:,:,:), allocatable :: &
         lut_extice,          & ! LUT shortwave ice extinction coefficient
         lut_ssaice,          & ! LUT shortwave ice single scattering albedo
         lut_asyice             ! LUT shortwave ice asymmetry parameter
    ! cld_optics_scheme = 2
    real(kind_phys), dimension(:), allocatable :: &
         pade_sizereg_extliq, & ! Particle size regime boundaries for shortwave liquid extinction 
                                   ! coefficient for Pade interpolation  
         pade_sizereg_ssaliq, & ! Particle size regime boundaries for shortwave liquid single 
                                   ! scattering albedo for Pade interpolation 
         pade_sizereg_asyliq, & ! Particle size regime boundaries for shortwave liquid asymmetry 
                                   ! parameter for Pade interpolation  
         pade_sizereg_extice, & ! Particle size regime boundaries for shortwave ice extinction 
                                   ! coefficient for Pade interpolation  
         pade_sizereg_ssaice, & ! Particle size regime boundaries for shortwave ice single 
                                   ! scattering albedo for Pade interpolation 
         pade_sizereg_asyice    ! Particle size regime boundaries for shortwave ice asymmetry 
                                   ! parameter for Pade interpolation  
    real(kind_phys), dimension(:,:,:), allocatable :: &
         pade_extliq,         & ! PADE coefficients for shortwave liquid extinction
         pade_ssaliq,         & ! PADE coefficients for shortwave liquid single scattering albedo
         pade_asyliq            ! PADE coefficients for shortwave liquid asymmetry parameter
    real(kind_phys), dimension(:,:,:,:), allocatable :: &
         pade_extice,         & ! PADE coefficients for shortwave ice extinction
         pade_ssaice,         & ! PADE coefficients for shortwave ice single scattering albedo
         pade_asyice            ! PADE coefficients for shortwave ice asymmetry parameter
    ! Dimensions
    integer :: &
         nrghice_lw, nbandLWcldy, nsize_liq, nsize_ice, nsizereg,&
         ncoeff_ext, ncoeff_ssa_g, nbound, npairsLWcldy

    ! Local variables
    integer :: dimID,varID,status,ncid_lw_clds
    character(len=264) :: lw_cloud_props_file
    integer,parameter :: max_strlen=256
#ifdef MPI
    integer :: ierr
#endif

    ! Initialize
    errmsg = ''
    errflg = 0

    if (cld_optics_scheme .eq. 0) return

    ! Filenames are set in the physics_nml
    lw_cloud_props_file = trim(rrtmgp_root_dir)//trim(rrtmgp_lw_file_clouds)

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

       ! Check to ensure that number of ice-roughness categories is feasible.
       if (nrghice .gt. nrghice_lw) then
          errmsg = 'Number of RRTMGP ice-roughness categories requested in namelist file is not allowed'
       endif
    endif

    ! Broadcast dimensions to all processors
#ifdef MPI
    call MPI_BCAST(nbandSWcldy_sw,     1, MPI_INTEGER, mpiroot, mpicomm, ierr)
    if (cld_optics_scheme .eq. 1) then
       call MPI_BARRIER(mpicomm, ierr)
       call MPI_BCAST(nrghice_sw,      1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_liq_sw,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsize_ice_sw,    1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BARRIER(mpicomm, ierr)
    endif
    if (cld_optics_scheme .eq. 2) then
       call MPI_BARRIER(mpicomm, ierr)
       call MPI_BCAST(nrghice_sw,      1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nsizereg_sw,     1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ext_sw,   1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(ncoeff_ssa_g_sw, 1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BCAST(nbound_sw,       1, MPI_INTEGER, mpiroot, mpicomm, ierr)
       call MPI_BARRIER(mpicomm, ierr)
    endif
#endif

    if (Cld_optics_scheme .eq. 1) then
       allocate(lut_extliq(nsize_liq, nBandLWcldy))
       allocate(lut_ssaliq(nsize_liq, nBandLWcldy))
       allocate(lut_asyliq(nsize_liq, nBandLWcldy))
       allocate(lut_extice(nsize_ice, nBandLWcldy, nrghice_lw))
       allocate(lut_ssaice(nsize_ice, nBandLWcldy, nrghice_lw))
       allocate(lut_asyice(nsize_ice, nBandLWcldy, nrghice_lw))
       allocate(band_lims_cldy(2, nBandLWcldy))
    endif
    if (Cld_optics_scheme .eq. 2) then
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
       if (Cld_optics_scheme .eq. 1) then
          write (*,*) 'Reading RRTMGP longwave cloud data (LUT) ... '
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
       if (Cld_optics_scheme .eq. 2) then
          write (*,*) 'Reading RRTMGP longwave cloud data (PADE) ... '
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
    if (cld_optics_scheme .eq. 1) then
       if (mpirank==mpiroot) write (*,*) 'Broadcasting RRTMGP shortwave cloud data (LUT) ... '
       call MPI_BARRIER(mpicomm, ierr)
#ifndef SINGLE_PREC
       call MPI_BCAST(radliq_lwr,           1,                         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radliq_upr,           1,                         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radliq_fac,           1,                         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_lwr,           1,                         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_upr,           1,                         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_fac,           1,                         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_extliq,           size(lut_extliq),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_ssaliq,           size(lut_ssaliq),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_asyliq,           size(lut_asyliq),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_extice,           size(lut_extice),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_ssaice,           size(lut_ssaice),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_asyice,           size(lut_asyice),          MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(band_lims_cldy,       size(band_lims_cldy),      MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)    
#else
       call MPI_BCAST(radliq_lwr,           1,                         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(radliq_upr,           1,                         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(radliq_fac,           1,                         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_lwr,           1,                         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_upr,           1,                         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(radice_fac,           1,                         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_extliq,           size(lut_extliq),          MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_ssaliq,           size(lut_ssaliq),          MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_asyliq,           size(lut_asyliq),          MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_extice,           size(lut_extice),          MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_ssaice,           size(lut_ssaice),          MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(lut_asyice,           size(lut_asyice),          MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(band_lims_cldy,       size(band_lims_cldy),      MPI_REAL,               mpiroot, mpicomm, ierr)
#endif 
       call MPI_BARRIER(mpicomm, ierr)
    endif
    if (cld_optics_scheme .eq. 2) then
       if (mpirank==mpiroot) write (*,*) 'Broadcasting RRTMGP shortwave cloud data (PADE) ... '
       call MPI_BARRIER(mpicomm, ierr)
#ifndef SINGLE_PREC
       call MPI_BCAST(pade_extliq,          size(pade_extliq),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_ssaliq,          size(pade_ssaliq),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_asyliq,          size(pade_asyliq),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_extice,          size(pade_extice),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_ssaice,          size(pade_ssaice),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_asyice,          size(pade_asyice),         MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_extliq,  size(pade_sizereg_extliq), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_ssaliq,  size(pade_sizereg_ssaliq), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_asyliq,  size(pade_sizereg_asyliq), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_extice,  size(pade_sizereg_extice), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_ssaice,  size(pade_sizereg_ssaice), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_asyice,  size(pade_sizereg_asyice), MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)
       call MPI_BCAST(band_lims_cldy,       size(band_lims_cldy),      MPI_DOUBLE_PRECISION,   mpiroot, mpicomm, ierr)    
#else
       call MPI_BCAST(pade_extliq,          size(pade_extliq),         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_ssaliq,          size(pade_ssaliq),         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_asyliq,          size(pade_asyliq),         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_extice,          size(pade_extice),         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_ssaice,          size(pade_ssaice),         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_asyice,          size(pade_asyice),         MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_extliq,  size(pade_sizereg_extliq), MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_ssaliq,  size(pade_sizereg_ssaliq), MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_asyliq,  size(pade_sizereg_asyliq), MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_extice,  size(pade_sizereg_extice), MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_ssaice,  size(pade_sizereg_ssaice), MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(pade_sizereg_asyice,  size(pade_sizereg_asyice), MPI_REAL,               mpiroot, mpicomm, ierr)
       call MPI_BCAST(band_lims_cldy,       size(band_lims_cldy),      MPI_REAL,               mpiroot, mpicomm, ierr)
#endif
       call MPI_BARRIER(mpicomm, ierr)
    endif
#endif

    ! Load tables data for RRTMGP cloud-optics  
    if (cld_optics_scheme .eq. 1) then
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%set_ice_roughness(nrghice))
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%load(band_lims_cldy, &
            radliq_lwr, radliq_upr, radliq_fac, radice_lwr, radice_upr, radice_fac,    &
            lut_extliq, lut_ssaliq, lut_asyliq, lut_extice, lut_ssaice, lut_asyice))
    endif
    if (cld_optics_scheme .eq. 2) then
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%set_ice_roughness(nrghice))
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%load(band_lims_cldy, &
            pade_extliq, pade_ssaliq, pade_asyliq, pade_extice, pade_ssaice,           &
            pade_asyice, pade_sizereg_extliq, pade_sizereg_ssaliq, pade_sizereg_asyliq,&
            pade_sizereg_extice, pade_sizereg_ssaice, pade_sizereg_asyice))
    endif
  end subroutine rrtmgp_lw_cloud_optics_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_cloud_optics_run
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_lw_cloud_optics_run(doLWrad, nCol, nLev, cld_optics_scheme, nrghice,    &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,      &
       cld_rerain, lw_cloud_props, lw_gas_props, cldtaulw, lw_optical_props_cloudsByBand,   &
       errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doLWrad            ! Logical flag for longwave radiation call
    integer, intent(in) :: &
         nCol,             & ! Number of horizontal gridpoints
         nLev,             & ! Number of vertical levels
         nrghice,          & ! Number of ice-roughness categories
         cld_optics_scheme   ! Cloud-optics scheme
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
         cld_frac,         & ! Total cloud fraction by layer
         cld_lwp,          & ! Cloud liquid water path
         cld_reliq,        & ! Cloud liquid effective radius
         cld_iwp,          & ! Cloud ice water path
         cld_reice,        & ! Cloud ice effective radius
         cld_swp,          & ! Cloud snow water path       (used only for RRTMG legacy scheme)
         cld_resnow,       & ! Cloud snow effective radius (used only for RRTMG legacy scheme)
         cld_rwp,          & ! Cloud rain water path       (used only for RRTMG legacy scheme)
         cld_rerain          ! Cloud rain effective radius (used only for RRTMG legacy scheme)
    type(ty_cloud_optics),intent(in) :: &
         lw_cloud_props      ! RRTMGP DDT: spectral information for RRTMGP LW radiation scheme
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props        ! RRTMGP DDT: spectral information for RRTMGP LW radiation scheme
 
    ! Outputs
    real(kind_phys), dimension(ncol,nLev), intent(out) :: &
         cldtaulw                      ! Approx. 10.mu band layer cloud optical depth  
    type(ty_optical_props_1scl),intent(out) :: &
         lw_optical_props_cloudsByBand ! RRTMGP DDT: longwave cloud optical properties in each band
    integer, intent(out) :: &
         errflg                        ! CCPP error flag
    character(len=*), intent(out) :: &
         errmsg                        ! CCPP error message

    ! Local variables
    logical,dimension(ncol,nLev) :: liqmask, icemask
    real(kind_phys), dimension(ncol,nLev,lw_gas_props%get_nband()) :: &
         tau_cld

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return
    
    ! Compute ice/liquid cloud masks, needed by rrtmgp_cloud_optics
    liqmask = (cld_frac .gt. 0 .and. cld_lwp .gt. 0)
    icemask = (cld_frac .gt. 0 .and. cld_iwp .gt. 0)

    ! Allocate space for RRTMGP DDTs containing cloud radiative properties
    ! Cloud optics [nCol,nLev,nBands]
    call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_optical_props_cloudsByBand%alloc_1scl(&
         ncol, nLev, lw_gas_props%get_band_lims_wavenumber()))

    ! Compute cloud-optics for RTE.
    if (rrtmgp_cld_optics .gt. 0) then
       ! i) RRTMGP cloud-optics.
       call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_cloud_props%cloud_optics(&
            ncol,                          & ! IN  - Number of horizontal gridpoints 
            nLev,                          & ! IN  - Number of vertical layers
            lw_cloud_props%get_nband(),    & ! IN  - Number of LW bands
            nrghice,                       & ! IN  - Number of ice-roughness categories
            liqmask,                       & ! IN  - Liquid-cloud mask
            icemask,                       & ! IN  - Ice-cloud mask
            cld_lwp,                       & ! IN  - Cloud liquid water path
            cld_iwp,                       & ! IN  - Cloud ice water path
            cld_reliq,                     & ! IN  - Cloud liquid effective radius
            cld_reice,                     & ! IN  - Cloud ice effective radius
            lw_optical_props_cloudsByBand))  ! OUT - RRTMGP DDT containing cloud radiative properties
                                             !       in each band
    else
       ! ii) RRTMG cloud-optics.
       if (any(cld_frac .gt. 0)) then
          call rrtmg_lw_cloud_optics(ncol, nLev, lw_gas_props%get_nband(), cld_lwp,     &
               cld_reliq, cld_iwp, cld_reice, cld_rwp, cld_rerain, cld_swp, cld_resnow, &
               cld_frac, tau_cld)
          lw_optical_props_cloudsByBand%tau = tau_cld
       endif
    endif

    ! All-sky LW optical depth ~10microns
    cldtaulw = lw_optical_props_cloudsByBand%tau(:,:,7)

  end subroutine rrtmgp_lw_cloud_optics_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_finalize()
  ! #########################################################################################
  subroutine rrtmgp_lw_cloud_optics_finalize()
  end subroutine rrtmgp_lw_cloud_optics_finalize
end module rrtmgp_lw_cloud_optics
