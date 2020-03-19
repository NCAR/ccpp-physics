module rrtmgp_lw_cloud_optics
  use machine,                  only: kind_phys
  use mo_rte_kind,              only: wl
  use mo_cloud_optics,          only: ty_cloud_optics
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use mo_optical_props,         only: ty_optical_props_1scl
  use mo_rrtmg_lw_cloud_optics, only: rrtmg_lw_cloud_optics   
  use rrtmgp_aux,               only: check_error_msg
  use netcdf

  implicit none

  public rrtmgp_lw_cloud_optics_init, rrtmgp_lw_cloud_optics_run, rrtmgp_lw_cloud_optics_finalize

contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_init()
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_cloud_optics_init
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_lw_cloud_optics_init(cld_optics_scheme, nrghice, rrtmgp_root_dir, &
       rrtmgp_lw_file_clouds, mpicomm, mpirank, mpiroot, lw_cloud_props, errmsg, errflg)

    ! Inputs
    integer, intent(inout) :: &
         nrghice                ! Number of ice-roughness categories
   integer, intent(in) :: &
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
         band_lims         ! Beginning and ending wavenumber [cm -1] for each band                           
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
         nrghice_fromfile, nBand, nSize_liq, nSize_ice, nSizeReg,&
         nCoeff_ext, nCoeff_ssa_g, nBound, npairs

    ! Local variables
    integer :: dimID,varID,status,ncid
    character(len=264) :: lw_cloud_props_file
    integer,parameter :: max_strlen=256, nrghice_default=2

    ! Initialize
    errmsg = ''
    errflg = 0

    if (cld_optics_scheme .eq. 0) return

    ! Filenames are set in the physics_nml
    lw_cloud_props_file = trim(rrtmgp_root_dir)//trim(rrtmgp_lw_file_clouds)

    ! On master processor only...
!    if (mpirank .eq. mpiroot) then
       ! Open file
       status = nf90_open(trim(lw_cloud_props_file), NF90_WRITE, ncid)

       ! Read dimensions
       status = nf90_inq_dimid(ncid, 'nband', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nBand)
       status = nf90_inq_dimid(ncid, 'nrghice', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nrghice_fromfile)
       status = nf90_inq_dimid(ncid, 'nsize_liq', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nSize_liq)
       status = nf90_inq_dimid(ncid, 'nsize_ice', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nSize_ice)
       status = nf90_inq_dimid(ncid, 'nsizereg', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nSizeReg)
       status = nf90_inq_dimid(ncid, 'ncoeff_ext', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nCoeff_ext)
       status = nf90_inq_dimid(ncid, 'ncoeff_ssa_g', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nCoeff_ssa_g)
       status = nf90_inq_dimid(ncid, 'nbound', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nBound)
       status = nf90_inq_dimid(ncid, 'pair', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=npairs)
       status = nf90_close(ncid)

       ! Has the number of ice-roughnesses to use been provided from the namelist?
       ! If not provided, use default number of ice-roughness categories
       if (nrghice .eq. 0) then
          nrghice = nrghice_default
       else
          nrghice = nrghice_fromfile
          ! If provided in the namelist, check to ensure that number of ice-roughness categories is feasible.
          if (nrghice .gt. nrghice_fromfile) then
             errmsg  = 'Number of RRTMGP ice-roughness categories requested in namelist file is not allowed. Using default number of categories.'
             nrghice = nrghice_default
          endif
       endif

       ! Allocate space for arrays
       if (cld_optics_scheme .eq. 1) then
          allocate(lut_extliq(nSize_liq, nBand))
          allocate(lut_ssaliq(nSize_liq, nBand))
          allocate(lut_asyliq(nSize_liq, nBand))
          allocate(lut_extice(nSize_ice, nBand, nrghice_fromfile))
          allocate(lut_ssaice(nSize_ice, nBand, nrghice_fromfile))
          allocate(lut_asyice(nSize_ice, nBand, nrghice_fromfile))
       endif
       if (cld_optics_scheme .eq. 2) then
          allocate(pade_extliq(nBand, nSizeReg,  nCoeff_ext ))
          allocate(pade_ssaliq(nBand, nSizeReg,  nCoeff_ssa_g))
          allocate(pade_asyliq(nBand, nSizeReg,  nCoeff_ssa_g))
          allocate(pade_extice(nBand, nSizeReg,  nCoeff_ext,   nrghice_fromfile))
          allocate(pade_ssaice(nBand, nSizeReg,  nCoeff_ssa_g, nrghice_fromfile))
          allocate(pade_asyice(nBand, nSizeReg,  nCoeff_ssa_g, nrghice_fromfile))
          allocate(pade_sizereg_extliq(nBound))
          allocate(pade_sizereg_ssaliq(nBound))
          allocate(pade_sizereg_asyliq(nBound))
          allocate(pade_sizereg_extice(nBound))
          allocate(pade_sizereg_ssaice(nBound))
          allocate(pade_sizereg_asyice(nBound))
       endif
       allocate(band_lims(2,nBand))
       
       ! Read in fields from file
       if (cld_optics_scheme .eq. 1) then
          write (*,*) 'Reading RRTMGP longwave cloud data (LUT) ... '
          status = nf90_inq_varid(ncid,'radliq_lwr',varID)
          status = nf90_get_var(ncid,varID,radliq_lwr)
          status = nf90_inq_varid(ncid,'radliq_upr',varID)
          status = nf90_get_var(ncid,varID,radliq_upr)
          status = nf90_inq_varid(ncid,'radliq_fac',varID)
          status = nf90_get_var(ncid,varID,radliq_fac)
          status = nf90_inq_varid(ncid,'radice_lwr',varID)
          status = nf90_get_var(ncid,varID,radice_lwr)
          status = nf90_inq_varid(ncid,'radice_upr',varID)
          status = nf90_get_var(ncid,varID,radice_upr)
          status = nf90_inq_varid(ncid,'radice_fac',varID)
          status = nf90_get_var(ncid,varID,radice_fac)
          status = nf90_inq_varid(ncid,'lut_extliq',varID)
          status = nf90_get_var(ncid,varID,lut_extliq)
          status = nf90_inq_varid(ncid,'lut_ssaliq',varID)
          status = nf90_get_var(ncid,varID,lut_ssaliq)
          status = nf90_inq_varid(ncid,'lut_asyliq',varID)
          status = nf90_get_var(ncid,varID,lut_asyliq)
          status = nf90_inq_varid(ncid,'lut_extice',varID)
          status = nf90_get_var(ncid,varID,lut_extice)
          status = nf90_inq_varid(ncid,'lut_ssaice',varID)
          status = nf90_get_var(ncid,varID,lut_ssaice)
          status = nf90_inq_varid(ncid,'lut_asyice',varID)
          status = nf90_get_var(ncid,varID,lut_asyice)
          status = nf90_inq_varid(ncid,'bnd_limits_wavenumber',varID)
          status = nf90_get_var(ncid,varID,band_lims)
       endif
       if (cld_optics_scheme .eq. 2) then
          write (*,*) 'Reading RRTMGP longwave cloud data (PADE) ... '
          status = nf90_inq_varid(ncid,'radliq_lwr',varID)
          status = nf90_get_var(ncid,varID,radliq_lwr)
          status = nf90_inq_varid(ncid,'radliq_upr',varID)
          status = nf90_get_var(ncid,varID,radliq_upr)
          status = nf90_inq_varid(ncid,'radliq_fac',varID)
          status = nf90_get_var(ncid,varID,radliq_fac)
          status = nf90_inq_varid(ncid,'radice_lwr',varID)
          status = nf90_get_var(ncid,varID,radice_lwr)
          status = nf90_inq_varid(ncid,'radice_upr',varID)
          status = nf90_get_var(ncid,varID,radice_upr)
          status = nf90_inq_varid(ncid,'radice_fac',varID)
          status = nf90_get_var(ncid,varID,radice_fac)
          status = nf90_inq_varid(ncid,'pade_extliq',varID)
          status = nf90_get_var(ncid,varID,pade_extliq)
          status = nf90_inq_varid(ncid,'pade_ssaliq',varID)
          status = nf90_get_var(ncid,varID,pade_ssaliq)
          status = nf90_inq_varid(ncid,'pade_asyliq',varID)
          status = nf90_get_var(ncid,varID,pade_asyliq)
          status = nf90_inq_varid(ncid,'pade_extice',varID)
          status = nf90_get_var(ncid,varID,pade_extice)
          status = nf90_inq_varid(ncid,'pade_ssaice',varID)
          status = nf90_get_var(ncid,varID,pade_ssaice)
          status = nf90_inq_varid(ncid,'pade_asyice',varID)
          status = nf90_get_var(ncid,varID,pade_asyice)
          status = nf90_inq_varid(ncid,'pade_sizreg_extliq',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_extliq)
          status = nf90_inq_varid(ncid,'pade_sizreg_ssaliq',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_ssaliq)
          status = nf90_inq_varid(ncid,'pade_sizreg_asyliq',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_asyliq)
          status = nf90_inq_varid(ncid,'pade_sizreg_extice',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_extice)
          status = nf90_inq_varid(ncid,'pade_sizreg_ssaice',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_ssaice)
          status = nf90_inq_varid(ncid,'pade_sizreg_asyice',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_asyice)
          status = nf90_inq_varid(ncid,'bnd_limits_wavenumber',varID)
          status = nf90_get_var(ncid,varID,band_lims)
       endif
          
       ! Close file
       status = nf90_close(ncid)       
!    endif
 
    ! Load tables data for RRTMGP cloud-optics  
    if (cld_optics_scheme .eq. 1) then
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%load(band_lims, &
            radliq_lwr, radliq_upr, radliq_fac, radice_lwr, radice_upr, radice_fac,    &
            lut_extliq, lut_ssaliq, lut_asyliq, lut_extice, lut_ssaice, lut_asyice))
    endif
    if (cld_optics_scheme .eq. 2) then
       call check_error_msg('lw_cloud_optics_init', lw_cloud_props%load(band_lims,  &
            pade_extliq, pade_ssaliq, pade_asyliq, pade_extice, pade_ssaice, pade_asyice,&
            pade_sizereg_extliq, pade_sizereg_ssaliq, pade_sizereg_asyliq,               &
            pade_sizereg_extice, pade_sizereg_ssaice, pade_sizereg_asyice))
    endif
    call check_error_msg('lw_cloud_optics_init', lw_cloud_props%set_ice_roughness(nrghice))
 
  end subroutine rrtmgp_lw_cloud_optics_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_cloud_optics_run
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_lw_cloud_optics_run(doLWrad, nCol, nLev, cld_optics_scheme, nrghice,    &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,      &
       cld_rerain, p_lay, lw_cloud_props, lw_gas_props, lon, lat,                           &
       cldtaulw, lw_optical_props_cloudsByBand, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doLWrad             ! Logical flag for longwave radiation call
    integer, intent(in) :: &
         nCol,             & ! Number of horizontal gridpoints
         nLev,             & ! Number of vertical levels
         nrghice,          & ! Number of ice-roughness categories
         cld_optics_scheme   ! Cloud-optics scheme
    real(kind_phys), dimension(nCol), intent(in) :: &
         lon,              & ! Longitude
         lat                 ! Latitude
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
         p_lay,            & ! Layer pressure (Pa)
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
    integer :: iCol, iLay

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    tau_cld = 0.

    if (.not. doLWrad) return
    
    ! Compute ice/liquid cloud masks, needed by rrtmgp_cloud_optics
    liqmask = (cld_frac .gt. 0 .and. cld_lwp .gt. 0)
    icemask = (cld_frac .gt. 0 .and. cld_iwp .gt. 0)

    ! Allocate space for RRTMGP DDTs containing cloud radiative properties
    ! Cloud optics [nCol,nLev,nBands]
    call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_optical_props_cloudsByBand%alloc_1scl(&
         ncol, nLev, lw_gas_props%get_band_lims_wavenumber()))
    lw_optical_props_cloudsByBand%tau(:,:,:) = 0._kind_phys

    ! Compute cloud-optics for RTE.
    if (cld_optics_scheme .gt. 0) then
       ! i) RRTMGP cloud-optics.
       call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_cloud_props%cloud_optics(&
            !ncol,                          & ! IN  - Number of horizontal gridpoints 
            !nLev,                          & ! IN  - Number of vertical layers
            !lw_cloud_props%get_nband(),    & ! IN  - Number of LW bands
            !nrghice,                       & ! IN  - Number of ice-roughness categories
            !liqmask,                       & ! IN  - Liquid-cloud mask (1)
            !icemask,                       & ! IN  - Ice-cloud mask (1)
            cld_lwp,                       & ! IN  - Cloud liquid water path (g/m2)
            cld_iwp,                       & ! IN  - Cloud ice water path (g/m2)
            cld_reliq,                     & ! IN  - Cloud liquid effective radius (microns)
            cld_reice,                     & ! IN  - Cloud ice effective radius (microns)
            lw_optical_props_cloudsByBand))  ! OUT - RRTMGP DDT containing cloud radiative properties
                                             !       in each band
    else
       ! ii) RRTMG cloud-optics.
       if (any(cld_frac .gt. 0)) then
          call rrtmg_lw_cloud_optics(ncol, nLev, lw_gas_props%get_nband(), cld_lwp,     &
               cld_reliq, cld_iwp, cld_reice, cld_rwp, cld_rerain, cld_swp, cld_resnow, &
               cld_frac, tau_cld)
       endif
       lw_optical_props_cloudsByBand%tau = tau_cld
    endif    
 
    ! All-sky LW optical depth ~10microns
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
