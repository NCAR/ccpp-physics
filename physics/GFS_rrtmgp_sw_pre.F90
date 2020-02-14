module GFS_rrtmgp_sw_pre
  use physparam
  use machine, only: &
       kind_phys                   ! Working type
  use module_radiation_astronomy,only: &
       coszmn                      ! Function to compute cos(SZA)
  use module_radiation_surface,  only: &
       NF_ALBD,                  & ! Number of surface albedo categories (4; nir-direct, nir-diffuse, uvvis-direct, uvvis-diffuse)
       setalb                      ! Routine to compute surface albedo
  use surface_perturbation, only: & 
       cdfnor                      ! Routine to compute CDF (used to compute percentiles)
  use mo_gas_optics_rrtmgp,  only: &
       ty_gas_optics_rrtmgp
  public GFS_rrtmgp_sw_pre_run,GFS_rrtmgp_sw_pre_init,GFS_rrtmgp_sw_pre_finalize
  
contains

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_pre_init
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_pre_init ()
  end subroutine GFS_rrtmgp_sw_pre_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_pre_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_sw_pre_run
!! \htmlinclude GFS_rrtmgp_sw_pre.html
!!
  subroutine GFS_rrtmgp_sw_pre_run(doLWrad, do_sfcperts, ncol, nlev, ntrac, nsfcpert, nmtvr,mpi_rank, solhr, &
       pertalb, sfc_wts, xlon, coslat, sinlat, slmsk, snowd, sncovr, snoalb, zorl, coszen, coszdg, tsfc,&
       hprime, alvsf, alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc, p_lay, p_lev, tv_lay,  &
       relhum, tracer, sw_gas_props, nday, idxday, alb1d, sfalb, sfc_alb_nir_dir, sfc_alb_nir_dif,  &
       sfc_alb_uvvis_dir, sfc_alb_uvvis_dif,             errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         doLWrad ,          &  ! Flag for longwave radiation call
         do_sfcperts           ! Flag for stochastic surface perturbations option
 
    integer, intent(in)    :: &
         ncol,             & ! Number of horizontal grid points
         nlev,             & ! Number of vertical levels
         ntrac,            & ! Number of tracers
         nsfcpert,         & ! number of surface perturbations
         nmtvr,            & ! number of topographic variables in GWD
         mpi_rank            ! Current MPI-rank
    real(kind_phys), intent(in) :: &
         solhr               ! Time after 00z at the current timestep (hours)
    real(kind_phys), dimension(nsfcpert), intent(in) :: &
         pertalb             ! Magnitude of surface albedo perturbation
    real(kind_phys), dimension(ncol,nsfcpert), intent(in) :: &
         sfc_wts             ! Magnitude of surface albedo perturbation
    real(kind_phys), dimension(ncol), intent(in) :: &
         xlon,             & ! Longitude
         coslat,           & ! Cosine of latitude
         sinlat,           & ! Sine of latitude
         slmsk,            & ! Lank/sea mask
         snowd,            & ! Water equivalent snow depth (mm)
         sncovr,           & ! Surface snow area fraction
         snoalb,           & ! Maximum snow albedo
         zorl,             & ! Surface roughness length
         tsfc,             & ! Surface skin temperature (K)
         alvsf,            & ! Mean vis albedo with strong cosz dependency
         alnsf,            & ! Mean nIR albedo with strong cosz dependency
         alvwf,            & ! Mean vis albedo with weak   cosz dependency
         alnwf,            & ! Mean nIR albedo with weak   cosz dependency
         facsf,            & ! Fractional coverage with strong cosz dependency
         facwf,            & ! Fractional coverage with weak cosz dependency
         fice,             & ! Ice fraction over open water
         tisfc               ! Sea ice surface skin temperature
    real(kind_phys), dimension(ncol,nmtvr), intent(in) :: &
         hprime               ! orographic metrics 
    real(kind_phys), dimension(ncol,nlev),intent(in) :: &
         p_lay,             & ! Layer pressure
         tv_lay,            & ! Layer virtual-temperature
         relhum               ! Layer relative-humidity
    real(kind_phys), dimension(ncol, nlev, 2:ntrac),intent(in) :: &
         tracer               ! Chemical tracers (g/g)
    real(kind_phys), dimension(ncol,nlev+1),intent(in) :: &
         p_lev                ! Pressure @ layer interfaces (Pa)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props         ! RRTMGP DDT: spectral information for SW calculation

    ! Outputs
    integer, intent(out)   :: &
         nday                 ! Number of daylit points
    integer, dimension(ncol), intent(out) :: &
         idxday               ! Indices for daylit points
    real(kind_phys), dimension(ncol), intent(out) :: &
         coszen,           & ! mean cos of zenith angle over rad call period
         coszdg,           & ! daytime mean cosz over rad call period
         sfalb,            & ! mean surface diffused SW albedo
         alb1d               ! Surface albedo pertubation
    real(kind_phys), dimension(sw_gas_props%get_nband(),ncol), intent(out) :: &
         sfc_alb_nir_dir,   & ! Surface albedo (direct) 
         sfc_alb_nir_dif,   & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir, & ! Surface albedo (direct)
         sfc_alb_uvvis_dif    ! Surface albedo (diffuse)
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag

    ! Local variables
    integer :: i, j, iCol, iBand, iLay
    real(kind_phys), dimension(ncol, NF_ALBD) :: sfcalb

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. doLWrad) return
    
    ! #######################################################################################
    ! Compute cosine of zenith angle (only when SW is called)
    ! #######################################################################################
    call coszmn (xlon, sinlat, coslat, solhr, NCOL, mpi_rank, &
         coszen, coszdg)

    ! #######################################################################################
    ! For SW gather daylit points
    ! #######################################################################################
    nday   = 0
    idxday = 0
    do i = 1, NCOL
       if (coszen(i) >= 0.0001) then
          nday = nday + 1
          idxday(nday) = i
       endif
    enddo

    ! #######################################################################################
    ! mg, sfc-perts
    !  ---  scale random patterns for surface perturbations with perturbation size
    !  ---  turn vegetation fraction pattern into percentile pattern
    ! #######################################################################################
    alb1d(:) = 0.
    if (do_sfcperts) then
       if (pertalb(1) > 0.) then
          do i=1,ncol
             call cdfnor(sfc_wts(i,5),alb1d(i))
          enddo
       endif
    endif  
    
    ! #######################################################################################
    ! Call module_radiation_surface::setalb() to setup surface albedo.
    ! #######################################################################################
    call setalb (slmsk, snowd, sncovr, snoalb, zorl, coszen, tsfc, tsfc, hprime(:,1), alvsf,     &
         alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc, NCOL, alb1d, pertalb, sfcalb)
       
    ! Approximate mean surface albedo from vis- and nir-  diffuse values.
    sfalb(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
  
    ! Spread across all SW bands
    do iBand=1,sw_gas_props%get_nband()
       sfc_alb_nir_dir(iBand,1:NCOL)   = sfcalb(1:NCOL,1)
       sfc_alb_nir_dif(iBand,1:NCOL)   = sfcalb(1:NCOL,2)
       sfc_alb_uvvis_dir(iBand,1:NCOL) = sfcalb(1:NCOL,3)
       sfc_alb_uvvis_dif(iBand,1:NCOL) = sfcalb(1:NCOL,4)
    enddo 

  end subroutine GFS_rrtmgp_sw_pre_run
  
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_pre_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_pre_finalize ()
  end subroutine GFS_rrtmgp_sw_pre_finalize

end module GFS_rrtmgp_sw_pre
