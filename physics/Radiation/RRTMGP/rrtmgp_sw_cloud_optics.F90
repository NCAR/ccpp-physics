!>\file rrtmgp_sw_cloud_optics.F90
!!

!> This module contains the cloud optics properties calculation for RRTMGP-SW
module rrtmgp_sw_cloud_optics
  use machine,                  only: kind_phys
  use mo_rte_kind,              only: wl
  use mo_cloud_optics_rrtmgp,   only: ty_cloud_optics => ty_cloud_optics_rrtmgp
  use rrtmgp_sw_gas_optics,     only: sw_gas_props
  use radiation_tools,          only: check_error_msg
  use netcdf
#ifdef MPI
  use mpi_f08
#endif

  implicit none
  
  type(ty_cloud_optics) :: sw_cloud_props 
  integer :: &
       nrghice_fromfileSW, nBandSW, nSize_liqSW, nSize_iceSW, nSizeregSW, &
       nCoeff_extSW, nCoeff_ssa_gSW, nBoundSW, nPairsSW
  real(kind_phys), dimension(:,:), allocatable :: &
       lut_extliqSW,          & !< LUT shortwave liquid extinction coefficient  
       lut_ssaliqSW,          & !< LUT shortwave liquid single scattering albedo   
       lut_asyliqSW,          & !< LUT shortwave liquid asymmetry parameter  
       band_limsCLDSW           !< Beginning and ending wavenumber [cm -1] for each band                           
  real(kind_phys), dimension(:,:,:), allocatable :: &
       lut_exticeSW,          & !< LUT shortwave ice extinction coefficient
       lut_ssaiceSW,          & !< LUT shortwave ice single scattering albedo
       lut_asyiceSW             !< LUT shortwave ice asymmetry parameter
  real(kind_phys) :: &
       radliq_lwrSW,          & !< Liquid particle size lower bound for LUT interpolation
       radliq_uprSW,          & !< Liquid particle size upper bound for LUT interpolation
       radice_lwrSW,          & !< Ice particle    size upper bound for LUT interpolation
       radice_uprSW             !< Ice particle    size lower bound for LUT interpolation

  ! Parameters used for rain and snow(+groupel) RRTMGP cloud-optics. *NOTE* Same as in RRTMG
  ! Need to document these magic numbers below.
  real(kind_phys),parameter :: &
       a0r = 3.07e-3,        & !
       a0s = 0.0,            & !
       a1s = 1.5               !  
  real(kind_phys),dimension(:),allocatable :: b0r,b0s,b1s,c0r,c0s

contains
  ! ######################################################################################
  ! SUBROUTINE sw_cloud_optics_init
  ! ######################################################################################
!>
  subroutine rrtmgp_sw_cloud_optics_init( rrtmgp_root_dir, rrtmgp_sw_file_clouds,        &
       nrghice, mpicomm, mpirank, mpiroot, errmsg, errflg)

    ! Inputs
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,    & !< RTE-RRTMGP root directory
         rrtmgp_sw_file_clouds !< RRTMGP file containing cloud-optic data
    integer, intent(inout) :: &
         nrghice               !< Number of ice-roughness categories
    type(MPI_Comm), intent(in) :: &
         mpicomm               !< MPI communicator
    integer, intent(in) :: &
         mpirank,            & !< Current MPI rank
         mpiroot               !< Master MPI rank

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                !< CCPP error message
    integer,          intent(out) :: &
         errflg                !< CCPP error code
    
    ! Local variables
    integer :: status,ncid,dimid,varID,mpierr
    character(len=264) :: sw_cloud_props_file

    ! Initialize
    errmsg = ''
    errflg = 0

    ! Filenames are set in the physics_nml
    sw_cloud_props_file = trim(rrtmgp_root_dir)//trim(rrtmgp_sw_file_clouds)

    ! #######################################################################################
    !
    ! Read dimensions for shortwave cloud-optics fields...
    ! (ONLY master processor(0), if MPI enabled)
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif
       write (*,*) 'Reading RRTMGP shortwave cloud-optics metadata ... '

       ! Open file
       status = nf90_open(trim(sw_cloud_props_file), NF90_NOWRITE, ncid)

       ! Read dimensions
       status = nf90_inq_dimid(ncid, 'nband', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nBandSW)
       status = nf90_inq_dimid(ncid, 'nrghice', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nrghice_fromfileSW)
       status = nf90_inq_dimid(ncid, 'nsize_liq', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nSize_liqSW)
       status = nf90_inq_dimid(ncid, 'nsize_ice', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nSize_iceSW)
       status = nf90_inq_dimid(ncid, 'nsizereg', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nSizeregSW)
       status = nf90_inq_dimid(ncid, 'ncoeff_ext', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nCoeff_extSW)
       status = nf90_inq_dimid(ncid, 'ncoeff_ssa_g', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nCoeff_ssa_gSW)
       status = nf90_inq_dimid(ncid, 'nbound', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nBoundSW)
       status = nf90_inq_dimid(ncid, 'pair', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nPairsSW) 
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
    call mpi_bcast(nBandSW,            1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nSize_liqSW,        1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nSize_iceSW,        1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nSizeregSW,         1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nCoeff_extSW,       1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nCoeff_ssa_gSW,     1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nBoundSW,           1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
    call mpi_bcast(nPairsSW,           1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
#endif
    
    ! Has the number of ice-roughnes categories been provided from the namelist?
    ! If so, override nrghice from cloud-optics file
    if (nrghice .ne. 0) nrghice_fromfileSW = nrghice
#ifdef MPI
    call mpi_bcast(nrghice_fromfileSW, 1, MPI_INTEGER, mpiroot, mpicomm, mpierr)
#endif

    ! #######################################################################################
    !
    ! Allocate space for arrays...
    ! (ALL processors)
    !
    ! #######################################################################################
    allocate(lut_extliqSW(nSize_liqSW, nBandSW))
    allocate(lut_ssaliqSW(nSize_liqSW, nBandSW))
    allocate(lut_asyliqSW(nSize_liqSW, nBandSW))
    allocate(lut_exticeSW(nSize_iceSW, nBandSW, nrghice_fromfileSW))
    allocate(lut_ssaiceSW(nSize_iceSW, nBandSW, nrghice_fromfileSW))
    allocate(lut_asyiceSW(nSize_iceSW, nBandSW, nrghice_fromfileSW))
    allocate(band_limsCLDSW(2,nBandSW))

    ! #######################################################################################
    !
    ! Read in data ...
    ! (ONLY master processor(0), if MPI enabled) 
    !
    ! #######################################################################################
#ifdef MPI
    if (mpirank .eq. mpiroot) then
#endif 
       write (*,*) 'Reading RRTMGP shortwave cloud data (LUT) ... '
       status = nf90_inq_varid(ncid,'radliq_lwr',varID)
       status = nf90_get_var(ncid,varID,radliq_lwrSW)
       status = nf90_inq_varid(ncid,'radliq_upr',varID)
       status = nf90_get_var(ncid,varID,radliq_uprSW)
       status = nf90_inq_varid(ncid,'radice_lwr',varID)
       status = nf90_get_var(ncid,varID,radice_lwrSW)
       status = nf90_inq_varid(ncid,'radice_upr',varID)
       status = nf90_get_var(ncid,varID,radice_uprSW)
       status = nf90_inq_varid(ncid,'lut_extliq',varID)
       status = nf90_get_var(ncid,varID,lut_extliqSW)
       status = nf90_inq_varid(ncid,'lut_ssaliq',varID)
       status = nf90_get_var(ncid,varID,lut_ssaliqSW)
       status = nf90_inq_varid(ncid,'lut_asyliq',varID)
       status = nf90_get_var(ncid,varID,lut_asyliqSW)
       status = nf90_inq_varid(ncid,'lut_extice',varID)
       status = nf90_get_var(ncid,varID,lut_exticeSW)
       status = nf90_inq_varid(ncid,'lut_ssaice',varID)
       status = nf90_get_var(ncid,varID,lut_ssaiceSW)
       status = nf90_inq_varid(ncid,'lut_asyice',varID)
       status = nf90_get_var(ncid,varID,lut_asyiceSW)
       status = nf90_inq_varid(ncid,'bnd_limits_wavenumber',varID)
       status = nf90_get_var(ncid,varID,band_limsCLDSW)

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
    call mpi_bcast(radliq_lwrSW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(radliq_uprSW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(radice_lwrSW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(radice_uprSW, 1, MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)

    ! Real arrays
    call mpi_bcast(band_limsCLDSW, size(band_limsCLDSW),  &
         MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(lut_extliqSW, size(lut_extliqSW),      &
         MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(lut_ssaliqSW, size(lut_ssaliqSW),      &
         MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(lut_asyliqSW, size(lut_asyliqSW),      &
         MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(lut_exticeSW, size(lut_exticeSW),      &
         MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(lut_ssaiceSW, size(lut_ssaiceSW),      &
         MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
    call mpi_bcast(lut_asyiceSW, size(lut_asyiceSW),      &
         MPI_DOUBLE_PRECISION, mpiroot, mpicomm, mpierr)
#endif

    ! #######################################################################################
    !   
    ! Initialize RRTMGP DDT's...
    !
    ! #######################################################################################
    call check_error_msg('sw_cloud_optics_init',sw_cloud_props%load(band_limsCLDSW,         &
            radliq_lwrSW, radliq_uprSW, radice_lwrSW, radice_uprSW,                         &
            lut_extliqSW, lut_ssaliqSW, lut_asyliqSW,                                       &
            lut_exticeSW, lut_ssaiceSW, lut_asyiceSW))

    call check_error_msg('sw_cloud_optics_init',sw_cloud_props%set_ice_roughness(nrghice_fromfileSW))

    ! Initialize coefficients for rain and snow(+groupel) cloud optics
    allocate(b0r(sw_cloud_props%get_nband()),b0s(sw_cloud_props%get_nband()), &
             b1s(sw_cloud_props%get_nband()),c0r(sw_cloud_props%get_nband()), &
             c0s(sw_cloud_props%get_nband()))
    b0r = (/0.496, 0.466,   0.437,   0.416, 0.391, 0.374, 0.352,    &
            0.183, 0.048,   0.012,   0.000, 0.000, 0.000, 0.000/)
    b0s = (/0.460, 0.460,   0.460,   0.460, 0.460, 0.460, 0.460,    &
            0.460, 0.000,   0.000,   0.000, 0.000, 0.000, 0.000/)
    b1s = (/0.000, 0.000,   0.000,   0.000, 0.000, 0.000, 0.000,    &
            0.000, 1.62e-5, 1.62e-5, 0.000, 0.000, 0.000, 0.000/)
    c0r = (/0.980, 0.975,   0.965,   0.960, 0.955, 0.952, 0.950,    &
            0.944, 0.894,   0.884,   0.883, 0.883, 0.883, 0.883/)
    c0s = (/0.970, 0.970,   0.970,   0.970, 0.970, 0.970, 0.970,    &
            0.970, 0.970,   0.970,   0.700, 0.700, 0.700, 0.700/)

  end subroutine rrtmgp_sw_cloud_optics_init
end module rrtmgp_sw_cloud_optics
