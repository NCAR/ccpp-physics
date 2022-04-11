module rrtmgp_lw_cloud_optics
  use machine,                  only: kind_phys
  use mo_rte_kind,              only: wl
  use mo_cloud_optics,          only: ty_cloud_optics
  use mo_optical_props,         only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_rrtmg_lw_cloud_optics, only: rrtmg_lw_cloud_optics   
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
!! \section arg_table_rrtmgp_lw_cloud_optics_init
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_lw_cloud_optics_init(nrghice, mpicomm, mpirank, mpiroot,             &
       doG_cldoptics, doGP_cldoptics_PADE, doGP_cldoptics_LUT, rrtmgp_root_dir,          &
       rrtmgp_lw_file_clouds, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doG_cldoptics,                 & ! Use legacy RRTMG cloud-optics?
         doGP_cldoptics_PADE,           & ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT               ! Use RRTMGP cloud-optics: LUTs?
    integer, intent(inout) :: &
         nrghice                          ! Number of ice-roughness categories
    integer, intent(in) :: & 
         mpicomm,                       & ! MPI communicator
         mpirank,                       & ! Current MPI rank
         mpiroot                          ! Master MPI rank
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,               & ! RTE-RRTMGP root directory
         rrtmgp_lw_file_clouds            ! RRTMGP file containing coefficients used to compute clouds optical properties

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                           ! Error message
    integer,          intent(out) :: &
         errflg                           ! Error code

    ! Local variables
    integer :: dimID,varID,status,ncid,mpierr
    character(len=264) :: lw_cloud_props_file
    integer,parameter :: max_strlen=256, nrghice_default=2

    ! Initialize
    errmsg = ''
    errflg = 0

    ! If not using RRTMGP cloud optics, return.
    if (doG_cldoptics) return
    
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

  ! ######################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_run()
  ! ######################################################################################
!! \section arg_table_rrtmgp_lw_cloud_optics_run
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_lw_cloud_optics_run(doLWrad, doG_cldoptics, icliq_lw, icice_lw,      &
       doGP_cldoptics_PADE, doGP_cldoptics_LUT, doGP_lwscat, do_mynnedmf, imfdeepcnv,    &
       imfdeepcnv_gf, imfdeepcnv_samf, nCol, nLev, nbndsGPlw , p_lay, cld_frac, cld_lwp, &
       cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp, cld_rerain,          &
       precip_frac, cld_cnv_lwp, cld_cnv_reliq, cld_cnv_iwp, cld_cnv_reice, cld_pbl_lwp, &
       cld_pbl_reliq, cld_pbl_iwp, cld_pbl_reice, lon, lat, cldtaulw,                    &
       lw_optical_props_cloudsByBand, lw_optical_props_cnvcloudsByBand,                  &
       lw_optical_props_MYNNcloudsByBand, lw_optical_props_precipByBand, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doLWrad,             & ! Logical flag for longwave radiation call
         doG_cldoptics,       & ! Use legacy RRTMG cloud-optics?
         doGP_cldoptics_PADE, & ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT,  & ! Use RRTMGP cloud-optics: LUTs?
         doGP_lwscat,         & ! Include scattering in LW cloud-optics?
         do_mynnedmf            !
    integer, intent(in) ::    &
         nbndsGPlw,           & !
         nCol,                & ! Number of horizontal gridpoints
         nLev,                & ! Number of vertical levels
         icliq_lw,            & ! Choice of treatment of liquid cloud optical properties (RRTMG legacy)
         icice_lw,            & ! Choice of treatment of ice cloud optical properties (RRTMG legacy) 
         imfdeepcnv,          & ! 
         imfdeepcnv_gf,       & ! 
         imfdeepcnv_samf        ! 
    real(kind_phys), dimension(:), intent(in) :: &
         lon,                 & ! Longitude
         lat                    ! Latitude
    real(kind_phys), dimension(:,:),intent(in) :: &
         p_lay,               & ! Layer pressure (Pa)
         cld_frac,            & ! Total cloud fraction by layer
         cld_lwp,             & ! Cloud liquid water path
         cld_reliq,           & ! Cloud liquid effective radius
         cld_iwp,             & ! Cloud ice water path
         cld_reice,           & ! Cloud ice effective radius
         cld_swp,             & ! Cloud snow water path       
         cld_resnow,          & ! Cloud snow effective radius 
         cld_rwp,             & ! Cloud rain water path      
         cld_rerain,          & ! Cloud rain effective radius 
         precip_frac,         & ! Precipitation fraction by layer.
         cld_cnv_lwp,         & ! Water path for       convective liquid cloud-particles (microns) 
         cld_cnv_reliq,       & ! Effective radius for convective liquid cloud-particles (microns)
         cld_cnv_iwp,         & ! Water path for       convective ice cloud-particles (microns)
         cld_cnv_reice,       & ! Effective radius for convective ice cloud-particles (microns) 
         cld_pbl_lwp,         & ! Water path for       SGS PBL liquid cloud-particles
         cld_pbl_reliq,       & ! Effective radius for SGS PBL liquid cloud-particles
         cld_pbl_iwp,         & ! Water path for       SGS PBL ice    cloud-particles
         cld_pbl_reice          ! Effective radius for SGS PBL ice    cloud-particles
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                             ! CCPP error message
    integer, intent(out) :: &
         errflg                             ! CCPP error flag
    type(ty_optical_props_2str),intent(inout) :: &
         lw_optical_props_cloudsByBand,     & ! RRTMGP DDT: Longwave optical properties in each band (clouds)
         lw_optical_props_cnvcloudsByBand,  & ! RRTMGP DDT: Longwave optical properties in each band (convective cloud)
         lw_optical_props_MYNNcloudsByBand, & ! RRTMGP DDT: Longwave optical properties in each band (MYNN-PBL cloud)
         lw_optical_props_precipByBand        ! RRTMGP DDT: Longwave optical properties in each band (precipitation)
    real(kind_phys), dimension(:,:), intent(inout) :: &
         cldtaulw                           ! Approx 10.mu band layer cloud optical depth  
         
    ! Local variables
    real(kind_phys) :: tau_rain, tau_snow
    real(kind_phys), dimension(ncol,nLev,nbndsGPlw) :: &
         tau_cld, tau_precip
    integer :: iCol, iLay, iBand

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    ! Initialize locals
    tau_cld    = 0._kind_phys
    tau_precip = 0._kind_phys

    if (.not. doLWrad) return

    ! Compute cloud-optics for RTE.
    if (doGP_cldoptics_PADE .or. doGP_cldoptics_LUT) then
       
       ! i) Cloud-optics.
       lw_optical_props_cloudsByBand%band_lims_wvn     = lw_gas_props%get_band_lims_wavenumber()
       do iBand=1,lw_gas_props%get_nband()
          lw_optical_props_cloudsByBand%band2gpt(1:2,iBand)     = iBand
          lw_optical_props_cloudsByBand%gpt2band(iBand)         = iBand
       end do
       call check_error_msg('rrtmgp_lw_cloud_optics_run - clouds',lw_cloud_props%cloud_optics(&
            cld_lwp,                           & ! IN  - Cloud liquid water path (g/m2)
            cld_iwp,                           & ! IN  - Cloud ice water path (g/m2)
            cld_reliq,                         & ! IN  - Cloud liquid effective radius (microns)
            cld_reice,                         & ! IN  - Cloud ice effective radius (microns)
            lw_optical_props_cloudsByBand))      ! OUT - RRTMGP DDT containing cloud radiative properties
                                                 !       in each band
       ! ii) Convective cloud-optics
       if (imfdeepcnv == imfdeepcnv_samf .or. imfdeepcnv == imfdeepcnv_gf) then
          lw_optical_props_cnvcloudsByBand%band_lims_wvn = lw_gas_props%get_band_lims_wavenumber()
          do iBand=1,lw_gas_props%get_nband()
             lw_optical_props_cnvcloudsByBand%band2gpt(1:2,iBand) = iBand
             lw_optical_props_cnvcloudsByBand%gpt2band(iBand)     = iBand
          end do
          call check_error_msg('rrtmgp_lw_cnvcloud_optics_run - convective cloud',lw_cloud_props%cloud_optics(&
               cld_cnv_lwp,                       & ! IN  - Convective cloud liquid water path (g/m2)
               cld_cnv_iwp,                       & ! IN  - Convective cloud ice water path (g/m2)
               cld_cnv_reliq,                     & ! IN  - Convective cloud liquid effective radius (microns)
               cld_cnv_reice,                     & ! IN  - Convective cloud ice effective radius (microns)
               lw_optical_props_cnvcloudsByBand))   ! OUT - RRTMGP DDT containing convective cloud radiative properties
                                                    !       in each band
       endif

       ! iii) MYNN cloud-optics
       if (do_mynnedmf) then
          lw_optical_props_MYNNcloudsByBand%band_lims_wvn = lw_gas_props%get_band_lims_wavenumber()
          do iBand=1,lw_gas_props%get_nband()
             lw_optical_props_MYNNcloudsByBand%band2gpt(1:2,iBand) = iBand
             lw_optical_props_MYNNcloudsByBand%gpt2band(iBand)     = iBand
          end do
          call check_error_msg('rrtmgp_lw_MYNNcloud_optics_run - MYNN-EDMF cloud',lw_cloud_props%cloud_optics(&
               cld_pbl_lwp,                       & ! IN  - MYNN-EDMF PBL cloud liquid water path (g/m2)
               cld_pbl_iwp,                       & ! IN  - MYNN-EDMF PBL cloud ice water path (g/m2)
               cld_pbl_reliq,                     & ! IN  - MYNN-EDMF PBL cloud liquid effective radius (microns)
               cld_pbl_reice,                     & ! IN  - MYNN-EDMF PBL cloud ice effective radius (microns)
               lw_optical_props_MYNNcloudsByBand))   ! OUT - RRTMGP DDT containing MYNN-EDMF PBL  cloud radiative properties
                                                     !       in each band
       endif

       ! iv) Cloud precipitation optics: rain and snow(+groupel)   
       lw_optical_props_precipByBand%band_lims_wvn     = lw_gas_props%get_band_lims_wavenumber()
       do iBand=1,lw_gas_props%get_nband()
          lw_optical_props_precipByBand%band2gpt(1:2,iBand)     = iBand
          lw_optical_props_precipByBand%gpt2band(iBand)         = iBand
       end do
       do iCol=1,nCol
          do iLay=1,nLev                                      
             if (cld_frac(iCol,iLay) .gt. 0.) then
                ! Rain optical-depth (No band dependence)
                tau_rain = absrain*cld_rwp(iCol,iLay)
                
                ! Snow (+groupel) optical-depth (No band dependence)
                if (cld_swp(iCol,iLay) .gt. 0. .and. cld_resnow(iCol,iLay) .gt. 10._kind_phys) then
                   tau_snow = abssnow0*1.05756*cld_swp(iCol,iLay)/cld_resnow(iCol,iLay)
                else
                   tau_snow = 0.0
                endif
                do iBand=1,nbndsGPlw
                   lw_optical_props_precipByBand%tau(iCol,iLay,iBand) = tau_rain + tau_snow
                enddo
             endif
          enddo
       enddo
    endif

    ! All-sky LW optical depth ~10microns (DJS asks: Same as SW, move to cloud-diagnostics?)
    cldtaulw = lw_optical_props_cloudsByBand%tau(:,:,7)
        
  end subroutine rrtmgp_lw_cloud_optics_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_finalize()
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_cloud_optics_finalize
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_lw_cloud_optics_finalize()
  end subroutine rrtmgp_lw_cloud_optics_finalize

end module rrtmgp_lw_cloud_optics
