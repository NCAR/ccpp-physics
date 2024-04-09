!> \file rrtmgp_lw_cloud_optics.F90
!!
!> \defgroup rrtmgp_lw_cloud_optics rrtmgp_lw_cloud_optics.F90
!!
!! \brief This module contains two routines: The first initializes data and functions
!! needed to compute the longwave cloud radiative properteis in RRTMGP. The second routine
!! is a ccpp scheme within the "radiation loop", where the shortwave optical prperties
!! (optical-depth, single-scattering albedo, asymmetry parameter) are computed for ALL
!! cloud types visible to RRTMGP.
!!
module rrtmgp_lw_cloud_optics
  use machine,                  only: kind_phys
  use mo_rte_kind,              only: wl
  use mo_cloud_optics,          only: ty_cloud_optics
  use rrtmgp_lw_gas_optics,     only: lw_gas_props
  use radiation_tools,          only: check_error_msg
  use netcdf
#ifdef MPI
  use mpi
#endif

  implicit none

  type(ty_cloud_optics) :: lw_cloud_props   
  integer :: &
       nrghice_fromfileLW, nBandLW, nSize_liqLW, nSize_iceLW, nSizeRegLW, &
       nCoeff_extLW, nCoeff_ssa_gLW, nBoundLW, npairsLW
  real(kind_phys) :: &
       radliq_facLW,          & ! Factor for calculating LUT interpolation indices for liquid
       radice_facLW             ! Factor for calculating LUT interpolation indices for ice  
  real(kind_phys), dimension(:,:), allocatable :: &
       lut_extliqLW,          & ! LUT shortwave liquid extinction coefficient  
       lut_ssaliqLW,          & ! LUT shortwave liquid single scattering albedo   
       lut_asyliqLW,          & ! LUT shortwave liquid asymmetry parameter  
       band_limsCLDLW           ! Beginning and ending wavenumber [cm -1] for each band                           
  real(kind_phys), dimension(:,:,:), allocatable :: &
       lut_exticeLW,          & ! LUT shortwave ice extinction coefficient
       lut_ssaiceLW,          & ! LUT shortwave ice single scattering albedo
       lut_asyiceLW             ! LUT shortwave ice asymmetry parameter
  real(kind_phys), dimension(:), allocatable :: &
       pade_sizereg_extliqLW, & ! Particle size regime boundaries for shortwave liquid extinction 
                                ! coefficient for Pade interpolation  
       pade_sizereg_ssaliqLW, & ! Particle size regime boundaries for shortwave liquid single 
                                ! scattering albedo for Pade interpolation 
       pade_sizereg_asyliqLW, & ! Particle size regime boundaries for shortwave liquid asymmetry 
                                ! parameter for Pade interpolation  
       pade_sizereg_exticeLW, & ! Particle size regime boundaries for shortwave ice extinction 
                                ! coefficient for Pade interpolation  
       pade_sizereg_ssaiceLW, & ! Particle size regime boundaries for shortwave ice single 
                                ! scattering albedo for Pade interpolation 
       pade_sizereg_asyiceLW    ! Particle size regime boundaries for shortwave ice asymmetry 
                                ! parameter for Pade interpolation  
  real(kind_phys), dimension(:,:,:), allocatable :: &
       pade_extliqLW,         & ! PADE coefficients for shortwave liquid extinction
       pade_ssaliqLW,         & ! PADE coefficients for shortwave liquid single scattering albedo
       pade_asyliqLW            ! PADE coefficients for shortwave liquid asymmetry parameter
  real(kind_phys), dimension(:,:,:,:), allocatable :: &
       pade_exticeLW,         & ! PADE coefficients for shortwave ice extinction
       pade_ssaiceLW,         & ! PADE coefficients for shortwave ice single scattering albedo
       pade_asyiceLW            ! PADE coefficients for shortwave ice asymmetry parameter
  
  ! Parameters used for rain and snow(+groupel) RRTMGP cloud-optics
  real(kind_phys), parameter :: &
       absrain  = 0.33e-3, & ! Rain drop absorption coefficient \f$(m^{2}/g)\f$ .
       abssnow0 = 1.5,     & ! Snow flake absorption coefficient (micron), fu coeff
       abssnow1 = 2.34e-3    ! Snow flake absorption coefficient \f$(m^{2}/g)\f$, ncar coef
  real(kind_phys) :: &
       radliq_lwrLW,         & ! Liquid particle size lower bound for LUT interpolation   
       radliq_uprLW,         & ! Liquid particle size upper bound for LUT interpolation
       radice_lwrLW,         & ! Ice particle size upper bound for LUT interpolation  
       radice_uprLW            ! Ice particle size lower bound for LUT interpolation

contains

  ! ######################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_init()
  ! ######################################################################################
  subroutine rrtmgp_lw_cloud_optics_init(rrtmgp_root_dir, rrtmgp_lw_file_clouds,         &
       doGP_cldoptics_PADE, doGP_cldoptics_LUT, nrghice, mpicomm, mpirank, mpiroot,      &
       errmsg, errflg)

    ! Inputs
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,    & ! RTE-RRTMGP root directory
         rrtmgp_lw_file_clouds ! RRTMGP file containing clouds optics data

    logical, intent(in) :: &
         doGP_cldoptics_PADE,& ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT    ! Use RRTMGP cloud-optics: LUTs?
    integer, intent(inout) :: &
         nrghice               ! Number of ice-roughness categories
    integer, intent(in) :: & 
         mpicomm,            & ! MPI communicator
         mpirank,            & ! Current MPI rank
         mpiroot               ! Master MPI rank

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                ! Error message
    integer,          intent(out) :: &
         errflg                ! Error code

    ! Local variables
    integer :: dimID,varID,status,ncid,mpierr
    character(len=264) :: lw_cloud_props_file

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the physics_nml
    lw_cloud_props_file = trim(rrtmgp_root_dir)//trim(rrtmgp_lw_file_clouds)

    ! #######################################################################################
    !
    ! Read dimensions for longwave cloud-optics fields...
    ! (ONLY master processor(0), if MPI enabled)
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif
       write (*,*) 'Reading RRTMGP longwave cloud-optics metadata ... '

       ! Open file
       status = nf90_open(trim(lw_cloud_props_file), NF90_NOWRITE, ncid)
       
       ! Read dimensions
       status = nf90_inq_dimid(ncid, 'nband', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nBandLW)
       status = nf90_inq_dimid(ncid, 'nrghice', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nrghice_fromfileLW)
       status = nf90_inq_dimid(ncid, 'nsize_liq', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nSize_liqLW)
       status = nf90_inq_dimid(ncid, 'nsize_ice', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nSize_iceLW)
       status = nf90_inq_dimid(ncid, 'nsizereg', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nSizeRegLW)
       status = nf90_inq_dimid(ncid, 'ncoeff_ext', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nCoeff_extLW)
       status = nf90_inq_dimid(ncid, 'ncoeff_ssa_g', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nCoeff_ssa_gLW)
       status = nf90_inq_dimid(ncid, 'nbound', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nBoundLW)
       status = nf90_inq_dimid(ncid, 'pair', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=npairsLW)

#ifdef MPI
    endif ! On master processor

    ! Other processors waiting...
    call mpi_barrier(mpicomm, mpierr)

    ! #######################################################################################
    !
    ! Broadcast dimensions...
    ! (ALL processors)
    !
    ! #######################################################################################
    call mpi_bcast(nBandLW,            1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nSize_liqLW,        1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nSize_iceLW,        1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nSizeregLW,         1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nCoeff_extLW,       1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nCoeff_ssa_gLW,     1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nBoundLW,           1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nPairsLW,           1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
#endif

    ! Has the number of ice-roughnesses to use been provided from the namelist?
    ! If so, override nrghice from cloud-optics file
    if (nrghice .ne. 0) nrghice_fromfileLW = nrghice
#ifdef MPI
    call mpi_bcast(nrghice_fromfileLW, 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
#endif

    ! #######################################################################################
    !
    ! Allocate space for arrays...
    ! (ALL processors)
    !
    ! #######################################################################################
    if (doGP_cldoptics_LUT) then
       allocate(lut_extliqLW(nSize_liqLW, nBandLW))
       allocate(lut_ssaliqLW(nSize_liqLW, nBandLW))
       allocate(lut_asyliqLW(nSize_liqLW, nBandLW))
       allocate(lut_exticeLW(nSize_iceLW, nBandLW, nrghice_fromfileLW))
       allocate(lut_ssaiceLW(nSize_iceLW, nBandLW, nrghice_fromfileLW))
       allocate(lut_asyiceLW(nSize_iceLW, nBandLW, nrghice_fromfileLW))
    endif
    if (doGP_cldoptics_PADE) then
       allocate(pade_extliqLW(nBandLW, nSizeRegLW,  nCoeff_extLW ))
       allocate(pade_ssaliqLW(nBandLW, nSizeRegLW,  nCoeff_ssa_gLW))
       allocate(pade_asyliqLW(nBandLW, nSizeRegLW,  nCoeff_ssa_gLW))
       allocate(pade_exticeLW(nBandLW, nSizeRegLW,  nCoeff_extLW,   nrghice_fromfileLW))
       allocate(pade_ssaiceLW(nBandLW, nSizeRegLW,  nCoeff_ssa_gLW, nrghice_fromfileLW))
       allocate(pade_asyiceLW(nBandLW, nSizeRegLW,  nCoeff_ssa_gLW, nrghice_fromfileLW))
       allocate(pade_sizereg_extliqLW(nBoundLW))
       allocate(pade_sizereg_ssaliqLW(nBoundLW))
       allocate(pade_sizereg_asyliqLW(nBoundLW))
       allocate(pade_sizereg_exticeLW(nBoundLW))
       allocate(pade_sizereg_ssaiceLW(nBoundLW))
       allocate(pade_sizereg_asyiceLW(nBoundLW))
    endif
    allocate(band_limsCLDLW(2,nBandLW))
       
    ! #######################################################################################
    !
    ! Read in data ...
    ! (ONLY master processor(0), if MPI enabled) 
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif
       ! Read in fields from file
       if (doGP_cldoptics_LUT) then
          write (*,*) 'Reading RRTMGP longwave cloud data (LUT) ... '
          status = nf90_inq_varid(ncid,'radliq_lwr',varID)
          status = nf90_get_var(ncid,varID,radliq_lwrLW)
          status = nf90_inq_varid(ncid,'radliq_upr',varID)
          status = nf90_get_var(ncid,varID,radliq_uprLW)
          status = nf90_inq_varid(ncid,'radliq_fac',varID)
          status = nf90_get_var(ncid,varID,radliq_facLW)
          status = nf90_inq_varid(ncid,'radice_lwr',varID)
          status = nf90_get_var(ncid,varID,radice_lwrLW)
          status = nf90_inq_varid(ncid,'radice_upr',varID)
          status = nf90_get_var(ncid,varID,radice_uprLW)
          status = nf90_inq_varid(ncid,'radice_fac',varID)
          status = nf90_get_var(ncid,varID,radice_facLW)
          status = nf90_inq_varid(ncid,'lut_extliq',varID)
          status = nf90_get_var(ncid,varID,lut_extliqLW)
          status = nf90_inq_varid(ncid,'lut_ssaliq',varID)
          status = nf90_get_var(ncid,varID,lut_ssaliqLW)
          status = nf90_inq_varid(ncid,'lut_asyliq',varID)
          status = nf90_get_var(ncid,varID,lut_asyliqLW)
          status = nf90_inq_varid(ncid,'lut_extice',varID)
          status = nf90_get_var(ncid,varID,lut_exticeLW)
          status = nf90_inq_varid(ncid,'lut_ssaice',varID)
          status = nf90_get_var(ncid,varID,lut_ssaiceLW)
          status = nf90_inq_varid(ncid,'lut_asyice',varID)
          status = nf90_get_var(ncid,varID,lut_asyiceLW)
          status = nf90_inq_varid(ncid,'bnd_limits_wavenumber',varID)
          status = nf90_get_var(ncid,varID,band_limsCLDLW)
       endif
       if (doGP_cldoptics_PADE) then
          write (*,*) 'Reading RRTMGP longwave cloud data (PADE) ... '
          status = nf90_inq_varid(ncid,'radliq_lwr',varID)
          status = nf90_get_var(ncid,varID,radliq_lwrLW)
          status = nf90_inq_varid(ncid,'radliq_upr',varID)
          status = nf90_get_var(ncid,varID,radliq_uprLW)
          status = nf90_inq_varid(ncid,'radliq_fac',varID)
          status = nf90_get_var(ncid,varID,radliq_facLW)
          status = nf90_inq_varid(ncid,'radice_lwr',varID)
          status = nf90_get_var(ncid,varID,radice_lwrLW)
          status = nf90_inq_varid(ncid,'radice_upr',varID)
          status = nf90_get_var(ncid,varID,radice_uprLW)
          status = nf90_inq_varid(ncid,'radice_fac',varID)
          status = nf90_get_var(ncid,varID,radice_facLW)
          status = nf90_inq_varid(ncid,'pade_extliq',varID)
          status = nf90_get_var(ncid,varID,pade_extliqLW)
          status = nf90_inq_varid(ncid,'pade_ssaliq',varID)
          status = nf90_get_var(ncid,varID,pade_ssaliqLW)
          status = nf90_inq_varid(ncid,'pade_asyliq',varID)
          status = nf90_get_var(ncid,varID,pade_asyliqLW)
          status = nf90_inq_varid(ncid,'pade_extice',varID)
          status = nf90_get_var(ncid,varID,pade_exticeLW)
          status = nf90_inq_varid(ncid,'pade_ssaice',varID)
          status = nf90_get_var(ncid,varID,pade_ssaiceLW)
          status = nf90_inq_varid(ncid,'pade_asyice',varID)
          status = nf90_get_var(ncid,varID,pade_asyiceLW)
          status = nf90_inq_varid(ncid,'pade_sizreg_extliq',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_extliqLW)
          status = nf90_inq_varid(ncid,'pade_sizreg_ssaliq',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_ssaliqLW)
          status = nf90_inq_varid(ncid,'pade_sizreg_asyliq',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_asyliqLW)
          status = nf90_inq_varid(ncid,'pade_sizreg_extice',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_exticeLW)
          status = nf90_inq_varid(ncid,'pade_sizreg_ssaice',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_ssaiceLW)
          status = nf90_inq_varid(ncid,'pade_sizreg_asyice',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_asyiceLW)
          status = nf90_inq_varid(ncid,'bnd_limits_wavenumber',varID)
          status = nf90_get_var(ncid,varID,band_limsCLDLW)
       endif
          
       ! Close file
       status = nf90_close(ncid)       
#ifdef MPI
    endif ! Master process

    ! Other processors waiting...
    call mpi_barrier(mpicomm, mpierr)

    ! #######################################################################################
    !
    ! Broadcast data... 
    ! (ALL processors)
    !
    ! #######################################################################################

    ! Real scalars
    call mpi_bcast(radliq_facLW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(radice_facLW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(radliq_lwrLW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(radliq_uprLW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(radice_lwrLW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(radice_uprLW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)

    ! Real arrays
    call mpi_bcast(band_limsCLDLW, size(band_limsCLDLW),  &
         MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    if (doGP_cldoptics_LUT) then
       call mpi_bcast(lut_extliqLW, size(lut_extliqLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(lut_ssaliqLW, size(lut_ssaliqLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(lut_asyliqLW, size(lut_asyliqLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(lut_exticeLW, size(lut_exticeLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(lut_ssaiceLW, size(lut_ssaiceLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(lut_asyiceLW, size(lut_asyiceLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
    if (doGP_cldoptics_PADE) then
       call mpi_bcast(pade_extliqLW, size(pade_extliqLW),                   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_ssaliqLW, size(pade_ssaliqLW),                   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_asyliqLW, size(pade_asyliqLW),                   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_exticeLW, size(pade_exticeLW),                   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_ssaiceLW, size(pade_ssaiceLW),                   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_asyiceLW, size(pade_asyiceLW),                   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_sizereg_extliqLW, size(pade_sizereg_extliqLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_sizereg_ssaliqLW, size(pade_sizereg_ssaliqLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_sizereg_asyliqLW, size(pade_sizereg_asyliqLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_sizereg_exticeLW, size(pade_sizereg_exticeLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_sizereg_ssaiceLW, size(pade_sizereg_ssaiceLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
       call mpi_bcast(pade_sizereg_asyiceLW, size(pade_sizereg_asyiceLW),   &
            MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    endif
#endif

    ! #######################################################################################
    !   
    ! Initialize RRTMGP DDT's...
    !
    ! #######################################################################################
    if (doGP_cldoptics_LUT) then
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%load(band_limsCLDLW,      &
            radliq_lwrLW, radliq_uprLW, radliq_facLW, radice_lwrLW, radice_uprLW,           &
            radice_facLW, lut_extliqLW, lut_ssaliqLW, lut_asyliqLW, lut_exticeLW,           &
            lut_ssaiceLW, lut_asyiceLW))
    endif
 
    if (doGP_cldoptics_PADE) then
       call check_error_msg('lw_cloud_optics_init', lw_cloud_props%load(band_limsCLDLW,     &
            pade_extliqLW, pade_ssaliqLW, pade_asyliqLW, pade_exticeLW, pade_ssaiceLW,      &
            pade_asyiceLW, pade_sizereg_extliqLW, pade_sizereg_ssaliqLW,                    &
            pade_sizereg_asyliqLW, pade_sizereg_exticeLW, pade_sizereg_ssaiceLW,            &
            pade_sizereg_asyiceLW))
    endif

    call check_error_msg('lw_cloud_optics_init',lw_cloud_props%set_ice_roughness(nrghice))
 
  end subroutine rrtmgp_lw_cloud_optics_init
end module rrtmgp_lw_cloud_optics
