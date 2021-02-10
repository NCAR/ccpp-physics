module rrtmgp_sw_cloud_optics
  use machine,                  only: kind_phys
  use mo_rte_kind,              only: wl
  use mo_cloud_optics,          only: ty_cloud_optics
  use mo_optical_props,         only: ty_optical_props_2str
  use mo_rrtmg_sw_cloud_optics, only: rrtmg_sw_cloud_optics   
  use rrtmgp_sw_gas_optics,     only: sw_gas_props
  use rrtmgp_aux,               only: check_error_msg
  use netcdf

  implicit none

  public rrtmgp_sw_cloud_optics_init, rrtmgp_sw_cloud_optics_run, rrtmgp_sw_cloud_optics_finalize
  
  type(ty_cloud_optics) :: sw_cloud_props 
  real(kind_phys) :: &
       radliq_facSW,          & ! Factor for calculating LUT interpolation indices for liquid   
       radice_facSW             ! Factor for calculating LUT interpolation indices for ice  
  real(kind_phys), dimension(:,:), allocatable :: &
       lut_extliqSW,          & ! LUT shortwave liquid extinction coefficient  
       lut_ssaliqSW,          & ! LUT shortwave liquid single scattering albedo   
       lut_asyliqSW,          & ! LUT shortwave liquid asymmetry parameter  
       band_limsCLDSW           ! Beginning and ending wavenumber [cm -1] for each band                           
  real(kind_phys), dimension(:,:,:), allocatable :: &
       lut_exticeSW,          & ! LUT shortwave ice extinction coefficient
       lut_ssaiceSW,          & ! LUT shortwave ice single scattering albedo
       lut_asyiceSW             ! LUT shortwave ice asymmetry parameter
  real(kind_phys), dimension(:), allocatable :: &
       pade_sizereg_extliqSW, & ! Particle size regime boundaries for shortwave liquid extinction 
                                ! coefficient for Pade interpolation  
       pade_sizereg_ssaliqSW, & ! Particle size regime boundaries for shortwave liquid single 
                                ! scattering albedo for Pade interpolation 
       pade_sizereg_asyliqSW, & ! Particle size regime boundaries for shortwave liquid asymmetry 
                                ! parameter for Pade interpolation  
       pade_sizereg_exticeSW, & ! Particle size regime boundaries for shortwave ice extinction 
                                ! coefficient for Pade interpolation  
       pade_sizereg_ssaiceSW, & ! Particle size regime boundaries for shortwave ice single 
                                ! scattering albedo for Pade interpolation 
       pade_sizereg_asyiceSW    ! Particle size regime boundaries for shortwave ice asymmetry 
                                ! parameter for Pade interpolation  
  real(kind_phys), dimension(:,:,:), allocatable :: &
       pade_extliqSW,         & ! PADE coefficients for shortwave liquid extinction
       pade_ssaliqSW,         & ! PADE coefficients for shortwave liquid single scattering albedo
       pade_asyliqSW            ! PADE coefficients for shortwave liquid asymmetry parameter
  real(kind_phys), dimension(:,:,:,:), allocatable :: &
       pade_exticeSW,         & ! PADE coefficients for shortwave ice extinction
       pade_ssaiceSW,         & ! PADE coefficients for shortwave ice single scattering albedo
       pade_asyiceSW            ! PADE coefficients for shortwave ice asymmetry parameter

  ! Parameters used for rain and snow(+groupel) RRTMGP cloud-optics
  real(kind_phys),parameter :: &
       a0r = 3.07e-3, & !
       a0s = 0.0,     & !
       a1s = 1.5        !  
  real(kind_phys),dimension(:),allocatable :: b0r,b0s,b1s,c0r,c0s
  real(kind_phys) :: &
       radliq_lwrSW,         & ! Liquid particle size lower bound for LUT interpolation   
       radliq_uprSW,         & ! Liquid particle size upper bound for LUT interpolation
       radice_lwrSW,         & ! Ice particle size upper bound for LUT interpolation  
       radice_uprSW            ! Ice particle size lower bound for LUT interpolation

contains
  ! ######################################################################################
  ! SUBROUTINE sw_cloud_optics_init
  ! ######################################################################################
!! \section arg_table_rrtmgp_sw_cloud_optics_init
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_sw_cloud_optics_init(nCol, nLev, nbndsGPsw, doG_cldoptics, doGP_cldoptics_PADE,             &
       doGP_cldoptics_LUT, nrghice, rrtmgp_root_dir, rrtmgp_sw_file_clouds, mpicomm,     &
       mpirank, mpiroot, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doG_cldoptics,       & ! Use legacy RRTMG cloud-optics?
         doGP_cldoptics_PADE, & ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT     ! Use RRTMGP cloud-optics: LUTs?    
    integer, intent(inout) :: &
         nrghice               ! Number of ice-roughness categories
    integer, intent(in) :: &
         nbndsGPsw,          & ! Number of bands used in shortwave.
         mpicomm,            & ! MPI communicator
         mpirank,            & ! Current MPI rank
         mpiroot,            & ! Master MPI rank
         nCol,               & ! Number of horizontal gridpoints
         nLev                  ! Number of vertical levels 
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,    & ! RTE-RRTMGP root directory
         rrtmgp_sw_file_clouds ! RRTMGP file containing coefficients used to compute clouds optical properties

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                ! CCPP error message
    integer,          intent(out) :: &
         errflg                ! CCPP error code
    
    ! Dimensions
    integer :: &
         nrghice_fromfile, nBand, nSize_liq, nSize_ice, nSizereg,&
         nCoeff_ext, nCoeff_ssa_g, nBound, nPairs

    ! Local variables
    integer :: status,ncid,dimid,varID
    character(len=264) :: sw_cloud_props_file
    integer,parameter ::  nrghice_default=2

    ! Initialize
    errmsg = ''
    errflg = 0

    if (doG_cldoptics) return

    ! Filenames are set in the physics_nml
    sw_cloud_props_file = trim(rrtmgp_root_dir)//trim(rrtmgp_sw_file_clouds)

    ! On master processor only...
!    if (mpirank .eq. mpiroot) then
       ! Open file
       status = nf90_open(trim(sw_cloud_props_file), NF90_NOWRITE, ncid)

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
       status = nf90_inquire_dimension(ncid, dimid, len=nSizereg)
       status = nf90_inq_dimid(ncid, 'ncoeff_ext', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nCoeff_ext)
       status = nf90_inq_dimid(ncid, 'ncoeff_ssa_g', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nCoeff_ssa_g)
       status = nf90_inq_dimid(ncid, 'nbound', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nBound)
       status = nf90_inq_dimid(ncid, 'pair', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nPairs)
 
       ! Has the number of ice-roughnesses provided from the namelist?
       ! If not, use nrghice from cloud-optics file
       if (nrghice .eq. 0) nrghice = nrghice_fromfile

       ! Allocate space for arrays
       if (doGP_cldoptics_LUT) then
          allocate(lut_extliqSW(nSize_liq, nBand))
          allocate(lut_ssaliqSW(nSize_liq, nBand))
          allocate(lut_asyliqSW(nSize_liq, nBand))
          allocate(lut_exticeSW(nSize_ice, nBand, nrghice))
          allocate(lut_ssaiceSW(nSize_ice, nBand, nrghice))
          allocate(lut_asyiceSW(nSize_ice, nBand, nrghice))
       endif
       if (doGP_cldoptics_PADE) then
          allocate(pade_extliqSW(nBand, nSizeReg,  nCoeff_ext ))
          allocate(pade_ssaliqSW(nBand, nSizeReg,  nCoeff_ssa_g))
          allocate(pade_asyliqSW(nBand, nSizeReg,  nCoeff_ssa_g))
          allocate(pade_exticeSW(nBand, nSizeReg,  nCoeff_ext,   nrghice))
          allocate(pade_ssaiceSW(nBand, nSizeReg,  nCoeff_ssa_g, nrghice))
          allocate(pade_asyiceSW(nBand, nSizeReg,  nCoeff_ssa_g, nrghice))
          allocate(pade_sizereg_extliqSW(nBound))
          allocate(pade_sizereg_ssaliqSW(nBound))
          allocate(pade_sizereg_asyliqSW(nBound))
          allocate(pade_sizereg_exticeSW(nBound))
          allocate(pade_sizereg_ssaiceSW(nBound))
          allocate(pade_sizereg_asyiceSW(nBound))
       endif
       allocate(band_limsCLDSW(2,nBand))

       ! Read in fields from file
       if (doGP_cldoptics_LUT) then
          write (*,*) 'Reading RRTMGP shortwave cloud data (LUT) ... '
          status = nf90_inq_varid(ncid,'radliq_lwr',varID)
          status = nf90_get_var(ncid,varID,radliq_lwrSW)
          status = nf90_inq_varid(ncid,'radliq_upr',varID)
          status = nf90_get_var(ncid,varID,radliq_uprSW)
          status = nf90_inq_varid(ncid,'radliq_fac',varID)
          status = nf90_get_var(ncid,varID,radliq_facSW)
          status = nf90_inq_varid(ncid,'radice_lwr',varID)
          status = nf90_get_var(ncid,varID,radice_lwrSW)
          status = nf90_inq_varid(ncid,'radice_upr',varID)
          status = nf90_get_var(ncid,varID,radice_uprSW)
          status = nf90_inq_varid(ncid,'radice_fac',varID)
          status = nf90_get_var(ncid,varID,radice_facSW)
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
       endif
       if (doGP_cldoptics_PADE) then
          write (*,*) 'Reading RRTMGP shortwave cloud data (PADE) ... '
          status = nf90_inq_varid(ncid,'radliq_lwr',varID)
          status = nf90_get_var(ncid,varID,radliq_lwrSW)
          status = nf90_inq_varid(ncid,'radliq_upr',varID)
          status = nf90_get_var(ncid,varID,radliq_uprSW)
          status = nf90_inq_varid(ncid,'radliq_fac',varID)
          status = nf90_get_var(ncid,varID,radliq_facSW)
          status = nf90_inq_varid(ncid,'radice_lwr',varID)
          status = nf90_get_var(ncid,varID,radice_lwrSW)
          status = nf90_inq_varid(ncid,'radice_upr',varID)
          status = nf90_get_var(ncid,varID,radice_uprSW)
          status = nf90_inq_varid(ncid,'radice_fac',varID)
          status = nf90_get_var(ncid,varID,radice_facSW)
          status = nf90_inq_varid(ncid,'pade_extliq',varID)
          status = nf90_get_var(ncid,varID,pade_extliqSW)
          status = nf90_inq_varid(ncid,'pade_ssaliq',varID)
          status = nf90_get_var(ncid,varID,pade_ssaliqSW)
          status = nf90_inq_varid(ncid,'pade_asyliq',varID)
          status = nf90_get_var(ncid,varID,pade_asyliqSW)
          status = nf90_inq_varid(ncid,'pade_extice',varID)
          status = nf90_get_var(ncid,varID,pade_exticeSW)
          status = nf90_inq_varid(ncid,'pade_ssaice',varID)
          status = nf90_get_var(ncid,varID,pade_ssaiceSW)
          status = nf90_inq_varid(ncid,'pade_asyice',varID)
          status = nf90_get_var(ncid,varID,pade_asyiceSW)
          status = nf90_inq_varid(ncid,'pade_sizreg_extliq',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_extliqSW)
          status = nf90_inq_varid(ncid,'pade_sizreg_ssaliq',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_ssaliqSW)
          status = nf90_inq_varid(ncid,'pade_sizreg_asyliq',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_asyliqSW)
          status = nf90_inq_varid(ncid,'pade_sizreg_extice',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_exticeSW)
          status = nf90_inq_varid(ncid,'pade_sizreg_ssaice',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_ssaiceSW)
          status = nf90_inq_varid(ncid,'pade_sizreg_asyice',varID)
          status = nf90_get_var(ncid,varID,pade_sizereg_asyiceSW)
          status = nf90_inq_varid(ncid,'bnd_limits_wavenumber',varID)
          status = nf90_get_var(ncid,varID,band_limsCLDSW)
       endif

       ! Close file
       status = nf90_close(ncid)       
!    endif

    ! Load tables data for RRTMGP cloud-optics  
    if (doGP_cldoptics_LUT) then
!$omp critical (load_sw_cloud_props_LUTs)
       call check_error_msg('sw_cloud_optics_init',sw_cloud_props%load(band_limsCLDSW,        &
            radliq_lwrSW, radliq_uprSW, radliq_facSW, radice_lwrSW, radice_uprSW,  radice_facSW,     &
            lut_extliqSW, lut_ssaliqSW, lut_asyliqSW, lut_exticeSW, lut_ssaiceSW, lut_asyiceSW))
!$omp end critical (load_sw_cloud_props_LUTs)
    endif
    if (doGP_cldoptics_PADE) then
!$omp critical (load_sw_cloud_props_PADE_approx)
       call check_error_msg('sw_cloud_optics_init', sw_cloud_props%load(band_limsCLDSW,       &
            pade_extliqSW, pade_ssaliqSW, pade_asyliqSW, pade_exticeSW, pade_ssaiceSW, pade_asyiceSW,&
            pade_sizereg_extliqSW, pade_sizereg_ssaliqSW, pade_sizereg_asyliqSW,               &
            pade_sizereg_exticeSW, pade_sizereg_ssaiceSW, pade_sizereg_asyiceSW))
!$omp end critical (load_sw_cloud_props_PADE_approx) 
    endif
!$omp critical (load_sw_cloud_props_nrghice)
    call check_error_msg('sw_cloud_optics_init',sw_cloud_props%set_ice_roughness(nrghice))
!$omp end critical (load_sw_cloud_props_nrghice)

    ! Initialize coefficients for rain and snow(+groupel) cloud optics
!$omp critical (load_sw_precip_props) 
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
!$omp end critical (load_sw_precip_props)

  end subroutine rrtmgp_sw_cloud_optics_init

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_optics_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_optics_run
!! \htmlinclude rrtmgp_sw_cloud_optics.html
!!
  subroutine rrtmgp_sw_cloud_optics_run(doSWrad, doG_cldoptics, icliq_sw, icice_sw,         &
       doGP_cldoptics_PADE, doGP_cldoptics_LUT, nCol, nLev, nDay, nbndsGPsw, idxday, nrghice, cld_frac,&
       cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp, cld_rerain,    &
       precip_frac, sw_optical_props_cloudsByBand,            &
       sw_optical_props_precipByBand, cldtausw, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doSWrad,             & ! Logical flag for shortwave radiation call
         doG_cldoptics,       & ! Use legacy RRTMG cloud-optics?
         doGP_cldoptics_PADE, & ! Use RRTMGP cloud-optics: PADE approximation?
         doGP_cldoptics_LUT     ! Use RRTMGP cloud-optics: LUTs?
    integer, intent(in) :: &
         nbndsGPsw,           & ! Number of shortwave bands
         nCol,                & ! Number of horizontal gridpoints
         nLev,                & ! Number of vertical levels
         nday,                & ! Number of daylit points.
         nrghice,             & ! Number of ice-roughness categories
         icliq_sw,            & ! Choice of treatment of liquid cloud optical properties (RRTMG legacy)
         icice_sw               ! Choice of treatment of ice cloud optical properties (RRTMG legacy) 
    integer,intent(in),dimension(ncol) :: &
         idxday                 ! Indices for daylit points.
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
         cld_frac,            & ! Total cloud fraction by layer
         cld_lwp,             & ! Cloud liquid water path
         cld_reliq,           & ! Cloud liquid effective radius
         cld_iwp,             & ! Cloud ice water path
         cld_reice,           & ! Cloud ice effective radius
         cld_swp,             & ! Cloud snow water path
         cld_resnow,          & ! Cloud snow effective radius
         cld_rwp,             & ! Cloud rain water path
         cld_rerain,          & ! Cloud rain effective radius
         precip_frac            ! Precipitation fraction by layer

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                             ! CCPP error message
    integer,          intent(out) :: &
         errflg                             ! CCPP error flag
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_cloudsByBand,   & ! RRTMGP DDT: Shortwave optical properties in each band (clouds)
         sw_optical_props_precipByBand      ! RRTMGP DDT: Shortwave optical properties in each band (cloud precipitation)
    real(kind_phys), dimension(ncol,NLev), intent(out) :: &
         cldtausw                           ! Approx 10.mu band layer cloud optical depth  
    
    ! Local variables
    integer :: iDay, iLay, iBand
    real(kind_phys) :: tau_rain, tau_snow, ssa_rain, ssa_snow, asy_rain, asy_snow, &
         tau_prec, asy_prec, ssa_prec, asyw, ssaw, za1, za2
    real(kind_phys), dimension(nday,nLev,nbndsGPsw) :: &
         tau_cld, ssa_cld, asy_cld, tau_precip, ssa_precip, asy_precip
    type(ty_optical_props_2str) :: sw_optical_props_cloudsByBand_daylit

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. doSWrad) return
    
    ! Only process sunlit points...
    if (nDay .gt. 0) then
       
       ! Compute cloud/precipitation optics.
       if (doGP_cldoptics_PADE .or. doGP_cldoptics_LUT) then
          call check_error_msg('rrtmgp_sw_cloud_optics_run',sw_optical_props_cloudsByBand%alloc_2str(&
               nday, nLev, sw_cloud_props%get_band_lims_wavenumber()))
          sw_optical_props_cloudsByBand%tau(:,:,:) = 0._kind_phys
          sw_optical_props_cloudsByBand%ssa(:,:,:) = 1._kind_phys
          sw_optical_props_cloudsByBand%g(:,:,:)   = 0._kind_phys
          call check_error_msg('rrtmgp_sw_cloud_optics_run',sw_optical_props_precipByBand%alloc_2str(&
               nday, nLev, sw_cloud_props%get_band_lims_wavenumber()))
          sw_optical_props_precipByBand%tau(:,:,:) = 0._kind_phys
          sw_optical_props_precipByBand%ssa(:,:,:) = 1._kind_phys
          sw_optical_props_precipByBand%g(:,:,:)   = 0._kind_phys

          ! RRTMGP cloud-optics.
          call check_error_msg('rrtmgp_sw_cloud_optics_run',sw_cloud_props%cloud_optics(&
               cld_lwp(idxday(1:nday),:),    & ! IN  - Cloud liquid water path
               cld_iwp(idxday(1:nday),:),    & ! IN  - Cloud ice water path
               cld_reliq(idxday(1:nday),:),  & ! IN  - Cloud liquid effective radius
               cld_reice(idxday(1:nday),:),  & ! IN  - Cloud ice effective radius
               sw_optical_props_cloudsByBand)) ! OUT - RRTMGP DDT: Shortwave optical properties, 
                                               !       in each band (tau,ssa,g)
          ! Cloud precipitation optics: rain and snow(+groupel) 
          do iDay=1,nDay
             do iLay=1,nLev                                      
                if (cld_frac(idxday(iDay),iLay) .gt. 1.e-12_kind_phys) then
                   ! Rain/Snow optical depth (No band dependence)
                   tau_rain = cld_rwp(idxday(iDay),iLay)*a0r                   
                   if (cld_swp(idxday(iDay),iLay) .gt. 0. .and. cld_resnow(idxday(iDay),iLay) .gt. 10._kind_phys) then
                      tau_snow = cld_swp(idxday(iDay),iLay)*1.09087*(a0s + a1s/(1.0315*cld_resnow(idxday(iDay),iLay)))     ! fu's formula 
                   else
                      tau_snow = 0._kind_phys
                   endif
                   
                   ! Rain/Snow single-scattering albedo and asymmetry (Band dependent)
                   do iBand=1,nbndsGPsw
                      ! By species
                      ssa_rain = tau_rain*(1.-b0r(iBand))
                      asy_rain = ssa_rain*c0r(iBand)
                      ssa_snow = tau_snow*(1.-(b0s(iBand)+b1s(iBand)*1.0315*cld_resnow(idxday(iDay),iLay)))
                      asy_snow = ssa_snow*c0s(iBand)
                      ! Combine
                      tau_prec = max(1.e-12_kind_phys, tau_rain + tau_snow)
                      ssa_prec = max(1.e-12_kind_phys, ssa_rain + ssa_snow)
                      asy_prec = max(1.e-12_kind_phys, asy_rain + asy_snow)
                      asyw     = asy_prec/max(1.e-12_kind_phys, ssa_prec)
                      ssaw     = min(1._kind_phys-0.000001, ssa_prec/tau_prec)
                      za1      = asyw * asyw
                      za2      = ssaw * za1                      
                      sw_optical_props_precipByBand%tau(idxday(iDay),iLay,iBand) = (1._kind_phys - za2) * tau_prec
                      sw_optical_props_precipByBand%ssa(idxday(iDay),iLay,iBand) = (ssaw - za2) / (1._kind_phys - za2)
                      sw_optical_props_precipByBand%g(idxday(iDay),iLay,iBand)   = asyw/(1+asyw)
                   enddo
                endif
             enddo
          enddo
       endif
       if (doG_cldoptics) then
          call check_error_msg('rrtmgp_sw_cloud_optics_run',sw_optical_props_cloudsByBand%alloc_2str(&
               nday, nLev, sw_gas_props%get_band_lims_wavenumber()))
          sw_optical_props_cloudsByBand%tau(:,:,:) = 0._kind_phys
          sw_optical_props_cloudsByBand%ssa(:,:,:) = 1._kind_phys
          sw_optical_props_cloudsByBand%g(:,:,:)   = 0._kind_phys
          call check_error_msg('rrtmgp_sw_cloud_optics_run',sw_optical_props_precipByBand%alloc_2str(&
               nday, nLev, sw_gas_props%get_band_lims_wavenumber()))
          sw_optical_props_precipByBand%tau(:,:,:) = 0._kind_phys
          sw_optical_props_precipByBand%ssa(:,:,:) = 1._kind_phys
          sw_optical_props_precipByBand%g(:,:,:)   = 0._kind_phys

          ! RRTMG cloud(+precipitation) optics
          if (any(cld_frac .gt. 0)) then
             call rrtmg_sw_cloud_optics(nday, nLev, sw_gas_props%get_nband(), &
                  cld_lwp(idxday(1:nday),:), cld_reliq(idxday(1:nday),:),  &
                  cld_iwp(idxday(1:nday),:), cld_reice(idxday(1:nday),:),  &
                  cld_rwp(idxday(1:nday),:), cld_rerain(idxday(1:nday),:), &
                  cld_swp(idxday(1:nday),:), cld_resnow(idxday(1:nday),:), &
                  cld_frac(idxday(1:nday),:), icliq_sw, icice_sw,          &
                  tau_cld,    ssa_cld,    asy_cld,                         &
                  tau_precip, ssa_precip, asy_precip)

            ! Cloud-optics (Need to reorder from G->GP band conventions)
             sw_optical_props_cloudsByBand%tau(:,:,1) = tau_cld(:,:,sw_gas_props%get_nband())
             sw_optical_props_cloudsByBand%ssa(:,:,1) = ssa_cld(:,:,sw_gas_props%get_nband())
             sw_optical_props_cloudsByBand%g(:,:,1)   = asy_cld(:,:,sw_gas_props%get_nband())
             sw_optical_props_cloudsByBand%tau(:,:,2:sw_gas_props%get_nband()) = tau_cld(:,:,1:sw_gas_props%get_nband()-1)
             sw_optical_props_cloudsByBand%ssa(:,:,2:sw_gas_props%get_nband()) = ssa_cld(:,:,1:sw_gas_props%get_nband()-1)
             sw_optical_props_cloudsByBand%g(:,:,2:sw_gas_props%get_nband())   = asy_cld(:,:,1:sw_gas_props%get_nband()-1)
             ! Precipitation-optics (Need to reorder from G->GP band conventions)
             sw_optical_props_precipByBand%tau(:,:,1) = tau_precip(:,:,sw_gas_props%get_nband())
             sw_optical_props_precipByBand%ssa(:,:,1) = ssa_precip(:,:,sw_gas_props%get_nband())
             sw_optical_props_precipByBand%g(:,:,1)   = asy_precip(:,:,sw_gas_props%get_nband())
             sw_optical_props_precipByBand%tau(:,:,2:sw_gas_props%get_nband()) = tau_precip(:,:,1:sw_gas_props%get_nband()-1)
             sw_optical_props_precipByBand%ssa(:,:,2:sw_gas_props%get_nband()) = ssa_precip(:,:,1:sw_gas_props%get_nband()-1)
             sw_optical_props_precipByBand%g(:,:,2:sw_gas_props%get_nband())   = asy_precip(:,:,1:sw_gas_props%get_nband()-1)
          
          endif
       endif
       
       ! All-sky SW optical depth ~0.55microns (DJS asks: Move to cloud diagnostics?)
       cldtausw(idxday(1:nDay),:) = sw_optical_props_cloudsByBand%tau(:,:,11)    
    endif
 
  end subroutine rrtmgp_sw_cloud_optics_run

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_optics_finalize()
  ! #########################################################################################  
  subroutine rrtmgp_sw_cloud_optics_finalize()
  end subroutine rrtmgp_sw_cloud_optics_finalize 

end module rrtmgp_sw_cloud_optics
