module rrtmgp_sw_cloud_optics
  use machine,                  only: kind_phys
  use mo_rte_kind,              only: wl
  use mo_gas_optics_rrtmgp,     only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,          only: ty_cloud_optics
  use physparam,                only: isubcsw, iovrsw
  use mo_optical_props,         only: ty_optical_props_2str
  use mo_rrtmg_sw_cloud_optics, only: rrtmg_sw_cloud_optics   
  use rrtmgp_aux,               only: check_error_msg
  use netcdf

  implicit none

  public rrtmgp_sw_cloud_optics_init, rrtmgp_sw_cloud_optics_run, rrtmgp_sw_cloud_optics_finalize

contains
  ! #########################################################################################
  ! SUBROUTINE sw_cloud_optics_init
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_optics_init
!! \htmlinclude rrtmgp_lw_cloud_optics.html
!!
  subroutine rrtmgp_sw_cloud_optics_init(cld_optics_scheme, nrghice, rrtmgp_root_dir,       &
       rrtmgp_sw_file_clouds, mpicomm, mpirank, mpiroot, sw_cloud_props, errmsg, errflg)

    ! Inputs
    integer, intent(inout) :: &
         nrghice               ! Number of ice-roughness categories
    integer, intent(in) :: &
         cld_optics_scheme,  & ! Cloud-optics scheme
         mpicomm,            & ! MPI communicator
         mpirank,            & ! Current MPI rank
         mpiroot               ! Master MPI rank
    character(len=128),intent(in) :: &
         rrtmgp_root_dir,    & ! RTE-RRTMGP root directory
         rrtmgp_sw_file_clouds ! RRTMGP file containing coefficients used to compute clouds optical properties

    ! Outputs
    type(ty_cloud_optics),intent(out) :: &
         sw_cloud_props        ! RRTMGP DDT: shortwave spectral information
    character(len=*), intent(out) :: &
         errmsg                ! CCPP error message
    integer,          intent(out) :: &
         errflg                ! CCPP error code

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
         nrghice_fromfile, nBand, nSize_liq, nSize_ice, nSizereg,&
         nCoeff_ext, nCoeff_ssa_g, nBound, nPairs

    ! Local variables
    integer :: status,ncid,dimid,varID
    character(len=264) :: sw_cloud_props_file
    integer,parameter ::  nrghice_default=2

    ! Initialize
    errmsg = ''
    errflg = 0

    if (cld_optics_scheme .eq. 0) return

    ! Filenames are set in the physics_nml
    sw_cloud_props_file = trim(rrtmgp_root_dir)//trim(rrtmgp_sw_file_clouds)

    ! On master processor only...
!    if (mpirank .eq. mpiroot) then
       ! Open file
       status = nf90_open(trim(sw_cloud_props_file), NF90_WRITE, ncid)

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
          write (*,*) 'Reading RRTMGP shortwave cloud data (LUT) ... '
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
          write (*,*) 'Reading RRTMGP shortwave cloud data (PADE) ... '
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
       call check_error_msg('sw_cloud_optics_init',sw_cloud_props%load(band_lims,        &
            radliq_lwr, radliq_upr, radliq_fac, radice_lwr, radice_upr,  radice_fac,     &
            lut_extliq, lut_ssaliq, lut_asyliq, lut_extice, lut_ssaice, lut_asyice))
    endif
    if (cld_optics_scheme .eq. 2) then
       call check_error_msg('sw_cloud_optics_init', sw_cloud_props%load(band_lims,       &
            pade_extliq, pade_ssaliq, pade_asyliq, pade_extice, pade_ssaice, pade_asyice,&
            pade_sizereg_extliq, pade_sizereg_ssaliq, pade_sizereg_asyliq,               &
            pade_sizereg_extice, pade_sizereg_ssaice, pade_sizereg_asyice))
    endif
    call check_error_msg('sw_cloud_optics_init',sw_cloud_props%set_ice_roughness(nrghice))
  end subroutine rrtmgp_sw_cloud_optics_init

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_optics_run()
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_cloud_optics_run
!! \htmlinclude rrtmgp_sw_cloud_optics.html
!!
  subroutine rrtmgp_sw_cloud_optics_run(doSWrad, nCol, nLev, nDay, idxday, nrghice,         &
       cld_optics_scheme, cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp,        &
       cld_resnow, cld_rwp, cld_rerain, sw_cloud_props, sw_gas_props,                       &
       sw_optical_props_cloudsByBand, cldtausw, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doSWrad                       ! Logical flag for shortwave radiation call
    integer, intent(in) :: &
         nCol,                        & ! Number of horizontal gridpoints
         nLev,                        & ! Number of vertical levels
         nday,                        & ! Number of daylit points.
         nrghice,                     & ! Number of ice-roughness categories
         cld_optics_scheme              ! Cloud-optics scheme
    integer,intent(in),dimension(ncol) :: &
         idxday                         ! Indices for daylit points.
    real(kind_phys), dimension(ncol,nLev),intent(in) :: &
         cld_frac,                    & ! Total cloud fraction by layer
         cld_lwp,                     & ! Cloud liquid water path
         cld_reliq,                   & ! Cloud liquid effective radius
         cld_iwp,                     & ! Cloud ice water path
         cld_reice,                   & ! Cloud ice effective radius
         cld_swp,                     & ! Cloud snow water path
         cld_resnow,                  & ! Cloud snow effective radius
         cld_rwp,                     & ! Cloud rain water path
         cld_rerain                     ! Cloud rain effective radius
    type(ty_cloud_optics),intent(in) :: &
         sw_cloud_props                 ! RRTMGP DDT: shortwave cloud properties
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props                   ! RRTMGP DDT: shortwave K-distribution data

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                         ! CCPP error message
    integer,          intent(out) :: &
         errflg                         ! CCPP error code
    type(ty_optical_props_2str),intent(out) :: &
         sw_optical_props_cloudsByBand  ! RRTMGP DDT: Shortwave optical properties (cloudy atmosphere)
    real(kind_phys), dimension(ncol,NLev), intent(out) :: &
         cldtausw                       ! approx 10.mu band layer cloud optical depth  

    ! Local variables
    logical,dimension(nday,nLev) :: liqmask, icemask
    real(kind_phys), dimension(nday,nLev,sw_gas_props%get_nband()) :: &
         tau_cld, ssa_cld, asy_cld

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doSWrad) return
    if (nDay .gt. 0) then
       
       ! Compute ice/liquid cloud masks, needed by rrtmgp_cloud_optics
       liqmask = (cld_frac(idxday(1:nday),:) .gt. 0 .and. cld_lwp(idxday(1:nday),:) .gt. 0)
       icemask = (cld_frac(idxday(1:nday),:) .gt. 0 .and. cld_iwp(idxday(1:nday),:) .gt. 0)
       
       ! Allocate space for RRTMGP DDTs containing cloud radiative properties
       ! Cloud optics [nday,nLev,nBands]
       call check_error_msg('rrtmgp_sw_cloud_optics_run',sw_optical_props_cloudsByBand%alloc_2str(&
            nday, nLev, sw_gas_props%get_band_lims_wavenumber()))
       sw_optical_props_cloudsByBand%tau(:,:,:) = 0._kind_phys
       sw_optical_props_cloudsByBand%ssa(:,:,:) = 0._kind_phys
       sw_optical_props_cloudsByBand%g(:,:,:)   = 0._kind_phys 

       ! Compute cloud-optics for RTE.
       if (cld_optics_scheme .gt. 0) then
          ! RRTMGP cloud-optics.
          call check_error_msg('rrtmgp_sw_cloud_optics_run',sw_cloud_props%cloud_optics(&
               cld_lwp(idxday(1:nday),:),    & ! IN  - Cloud liquid water path
               cld_iwp(idxday(1:nday),:),    & ! IN  - Cloud ice water path
               cld_reliq(idxday(1:nday),:),  & ! IN  - Cloud liquid effective radius
               cld_reice(idxday(1:nday),:),  & ! IN  - Cloud ice effective radius
               sw_optical_props_cloudsByBand)) ! OUT - RRTMGP DDT: Shortwave optical properties, 
                                               !       in each band (tau,ssa,g)
       else
          ! RRTMG cloud-optics
          tau_cld(:,:,:) = 0._kind_phys
          ssa_cld(:,:,:) = 0._kind_phys
          asy_cld(:,:,:) = 0._kind_phys
          if (any(cld_frac .gt. 0)) then
             call rrtmg_sw_cloud_optics(nday, nLev, sw_gas_props%get_nband(),       &
                  cld_lwp(idxday(1:nday),:), cld_reliq(idxday(1:nday),:),           &
                  cld_iwp(idxday(1:nday),:), cld_reice(idxday(1:nday),:),           &
                  cld_rwp(idxday(1:nday),:), cld_rerain(idxday(1:nday),:),          &
                  cld_swp(idxday(1:nday),:), cld_resnow(idxday(1:nday),:),          &
                  cld_frac(idxday(1:nday),:), tau_cld, ssa_cld, asy_cld)
          endif
          sw_optical_props_cloudsByBand%tau(:,:,:) = tau_cld
          sw_optical_props_cloudsByBand%ssa(:,:,:) = ssa_cld
          sw_optical_props_cloudsByBand%g(:,:,:)   = asy_cld
       endif

       ! All-sky SW optical depth ~0.55microns
       cldtausw(idxday(1:nDay),:) = sw_optical_props_cloudsByBand%tau(:,:,11)    
    endif

  end subroutine rrtmgp_sw_cloud_optics_run

  ! #########################################################################################
  ! SUBROTUINE rrtmgp_sw_cloud_optics_finalize()
  ! #########################################################################################  
  subroutine rrtmgp_sw_cloud_optics_finalize()
  end subroutine rrtmgp_sw_cloud_optics_finalize 

end module rrtmgp_sw_cloud_optics
