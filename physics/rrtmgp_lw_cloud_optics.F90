module rrtmgp_lw_cloud_optics
  use machine,                  only: kind_phys
  use mo_rte_kind,              only: wl
  use mo_cloud_optics,          only: ty_cloud_optics
  use mo_optical_props,         only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_rrtmg_lw_cloud_optics, only: rrtmg_lw_cloud_optics   
  use rrtmgp_lw_gas_optics,     only: lw_gas_props
  use rrtmgp_aux,               only: check_error_msg
  use netcdf

  implicit none

  public rrtmgp_lw_cloud_optics_init, rrtmgp_lw_cloud_optics_run, rrtmgp_lw_cloud_optics_finalize

  type(ty_cloud_optics) :: lw_cloud_props   
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
  subroutine rrtmgp_lw_cloud_optics_init(nCol, nLev, nbndsGPlw, doG_cldoptics,        &
       doGP_cldoptics_PADE, doGP_cldoptics_LUT, nrghice, rrtmgp_root_dir,                &
       rrtmgp_lw_file_clouds, mpicomm, mpirank, mpiroot, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doG_cldoptics,       & ! Use legacy RRTMG cloud-optics?
         doGP_cldoptics_PADE, & ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT     ! Use RRTMGP cloud-optics: LUTs?
    integer, intent(inout) :: &
         nrghice                ! Number of ice-roughness categories
    integer, intent(in) :: &
         nbndsGPlw,           & ! Number of longwave bands
         nCol,                & ! Number of horizontal gridpoints
         nLev,                & ! Number of vertical levels 
         mpicomm,             & ! MPI communicator
         mpirank,             & ! Current MPI rank
         mpiroot                ! Master MPI rank
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,     & ! RTE-RRTMGP root directory
         rrtmgp_lw_file_clouds  ! RRTMGP file containing coefficients used to compute clouds optical properties

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                           ! Error message
    integer,          intent(out) :: &
         errflg                           ! Error code

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

    ! If not using RRTMGP cloud optics, return.
    if (doG_cldoptics) return

    !
    ! Otherwise, using RRTMGP cloud-optics, continue with initialization...
    !
    
    ! Filenames are set in the physics_nml
    lw_cloud_props_file = trim(rrtmgp_root_dir)//trim(rrtmgp_lw_file_clouds)

    ! On master processor only...
!    if (mpirank .eq. mpiroot) then
       ! Open file
       status = nf90_open(trim(lw_cloud_props_file), NF90_NOWRITE, ncid)

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

       ! Has the number of ice-roughnesses to use been provided from the namelist?
       ! If not, use nrghice from cloud-optics data file.
       if (nrghice .eq. 0) nrghice = nrghice_fromfile

       ! Allocate space for arrays
       if (doGP_cldoptics_LUT) then
          allocate(lut_extliqLW(nSize_liq, nBand))
          allocate(lut_ssaliqLW(nSize_liq, nBand))
          allocate(lut_asyliqLW(nSize_liq, nBand))
          allocate(lut_exticeLW(nSize_ice, nBand, nrghice))
          allocate(lut_ssaiceLW(nSize_ice, nBand, nrghice))
          allocate(lut_asyiceLW(nSize_ice, nBand, nrghice))
       endif
       if (doGP_cldoptics_PADE) then
          allocate(pade_extliqLW(nBand, nSizeReg,  nCoeff_ext ))
          allocate(pade_ssaliqLW(nBand, nSizeReg,  nCoeff_ssa_g))
          allocate(pade_asyliqLW(nBand, nSizeReg,  nCoeff_ssa_g))
          allocate(pade_exticeLW(nBand, nSizeReg,  nCoeff_ext,   nrghice))
          allocate(pade_ssaiceLW(nBand, nSizeReg,  nCoeff_ssa_g, nrghice))
          allocate(pade_asyiceLW(nBand, nSizeReg,  nCoeff_ssa_g, nrghice))
          allocate(pade_sizereg_extliqLW(nBound))
          allocate(pade_sizereg_ssaliqLW(nBound))
          allocate(pade_sizereg_asyliqLW(nBound))
          allocate(pade_sizereg_exticeLW(nBound))
          allocate(pade_sizereg_ssaiceLW(nBound))
          allocate(pade_sizereg_asyiceLW(nBound))
       endif
       allocate(band_limsCLDLW(2,nBand))
       
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
!    endif
 
    ! Load tables data for RRTMGP cloud-optics  
    if (doGP_cldoptics_LUT) then
!$omp critical (load_lw_cloud_props_LUTs)
       call check_error_msg('lw_cloud_optics_init',lw_cloud_props%load(band_limsCLDLW,        &
            radliq_lwrLW, radliq_uprLW, radliq_facLW, radice_lwrLW, radice_uprLW,  radice_facLW,     &
            lut_extliqLW, lut_ssaliqLW, lut_asyliqLW, lut_exticeLW, lut_ssaiceLW, lut_asyiceLW))
!$omp end critical (load_lw_cloud_props_LUTs)  
    endif
    if (doGP_cldoptics_PADE) then
!$omp critical (load_lw_cloud_props_PADE_approx)
       call check_error_msg('lw_cloud_optics_init', lw_cloud_props%load(band_limsCLDLW,       &
            pade_extliqLW, pade_ssaliqLW, pade_asyliqLW, pade_exticeLW, pade_ssaiceLW, pade_asyiceLW,&
            pade_sizereg_extliqLW, pade_sizereg_ssaliqLW, pade_sizereg_asyliqLW,               &
            pade_sizereg_exticeLW, pade_sizereg_ssaiceLW, pade_sizereg_asyiceLW))
!$omp endcritical (load_lw_cloud_props_PADE_approx)
    endif
!$omp critical (load_lw_cloud_props_nrghice)
    call check_error_msg('lw_cloud_optics_init',lw_cloud_props%set_ice_roughness(nrghice))
!$omp end critical (load_lw_cloud_props_nrghice) 
 
  end subroutine rrtmgp_lw_cloud_optics_init

  ! ######################################################################################
  ! SUBROUTINE rrtmgp_lw_cloud_optics_run()
  ! ######################################################################################
!! \section arg_table_rrtmgp_lw_cloud_optics_run
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_lw_cloud_optics_run(doLWrad, doG_cldoptics, icliq_lw, icice_lw,      &
       doGP_cldoptics_PADE, doGP_cldoptics_LUT, doGP_lwscat, nCol, nLev, nbndsGPlw, p_lay, &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,   &
       cld_rerain, precip_frac, lon, lat, cldtaulw,        &
       lw_optical_props_cloudsByBand, lw_optical_props_precipByBand, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doLWrad,             & ! Logical flag for longwave radiation call
         doG_cldoptics,       & ! Use legacy RRTMG cloud-optics?
         doGP_cldoptics_PADE, & ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT,  & ! Use RRTMGP cloud-optics: LUTs?
         doGP_lwscat            ! Include scattering in LW cloud-optics?
    integer, intent(in) ::    &
         nbndsGPlw,           & ! Number of longwave bands
         nCol,                & ! Number of horizontal gridpoints
         nLev,                & ! Number of vertical levels
         icliq_lw,            & ! Choice of treatment of liquid cloud optical properties (RRTMG legacy)
         icice_lw               ! Choice of treatment of ice cloud optical properties (RRTMG legacy) 
    real(kind_phys), dimension(nCol), intent(in) :: &
         lon,                 & ! Longitude
         lat                    ! Latitude
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
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
         precip_frac            ! Precipitation fraction by layer.
 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                             ! CCPP error message
    integer, intent(out) :: &
         errflg                             ! CCPP error flag
    type(ty_optical_props_2str),intent(out) :: &
         lw_optical_props_cloudsByBand,   & ! RRTMGP DDT: Longwave optical properties in each band (clouds)
         lw_optical_props_precipByBand      ! RRTMGP DDT: Longwave optical properties in each band (precipitation)
    real(kind_phys), dimension(ncol,nLev), intent(inout) :: &
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
       call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_optical_props_cloudsByBand%alloc_2str(&
            ncol, nLev, lw_cloud_props%get_band_lims_wavenumber()))
       lw_optical_props_cloudsByBand%tau(:,:,:) = 0._kind_phys
       call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_optical_props_precipByBand%alloc_2str(&
            ncol, nLev, lw_cloud_props%get_band_lims_wavenumber()))
       lw_optical_props_precipByBand%tau(:,:,:) = 0._kind_phys

       ! i) RRTMGP cloud-optics.
       call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_cloud_props%cloud_optics(&
            cld_lwp,                       & ! IN  - Cloud liquid water path (g/m2)
            cld_iwp,                       & ! IN  - Cloud ice water path (g/m2)
            cld_reliq,                     & ! IN  - Cloud liquid effective radius (microns)
            cld_reice,                     & ! IN  - Cloud ice effective radius (microns)
            lw_optical_props_cloudsByBand))  ! OUT - RRTMGP DDT containing cloud radiative properties
                                             !       in each band
       ! Add in rain and snow(+groupel) 
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
    if (doG_cldoptics) then
       call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_optical_props_cloudsByBand%alloc_2str(&
            ncol, nLev, lw_gas_props%get_band_lims_wavenumber()))
       lw_optical_props_cloudsByBand%tau(:,:,:) = 0._kind_phys
       call check_error_msg('rrtmgp_lw_cloud_optics_run',lw_optical_props_precipByBand%alloc_2str(&
            ncol, nLev, lw_gas_props%get_band_lims_wavenumber()))
       lw_optical_props_precipByBand%tau(:,:,:) = 0._kind_phys
       ! ii) RRTMG cloud-optics.
       if (any(cld_frac .gt. 0)) then
          call rrtmg_lw_cloud_optics(ncol, nLev, nbndsGPlw, cld_lwp, cld_reliq, cld_iwp,&
               cld_reice, cld_rwp, cld_rerain, cld_swp, cld_resnow, cld_frac, icliq_lw, &
               icice_lw, tau_cld, tau_precip)
          lw_optical_props_cloudsByBand%tau = tau_cld
          lw_optical_props_precipByBand%tau = tau_precip
        endif
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
