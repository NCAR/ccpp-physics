!> \file GFS_rrtmgp_sw_pre.f90
!! This file contains
module GFS_rrtmgp_sw_pre
  use physparam
  use machine, only: &
       kind_phys                   ! Working type
  use GFS_typedefs, only:        &
       GFS_sfcprop_type,         & ! Surface fields
       GFS_control_type,         & ! Model control parameters
       GFS_grid_type,            & ! Grid and interpolation related data
       GFS_coupling_type,        & !
       GFS_statein_type,         & !
       GFS_radtend_type,         & ! Radiation tendencies needed in physics
       GFS_interstitial_type
  use module_radiation_astronomy,only: &
       coszmn                      ! Function to compute cos(SZA)
  use module_radiation_surface,  only: &
       NF_ALBD,                  & ! Number of surface albedo categories (4; nir-direct, nir-diffuse, uvvis-direct, uvvis-diffuse)
       setalb                      ! Routine to compute surface albedo
  ! DJS2019: This radiation_aerosols_module is a whole-lotta mess that needs some love. As it stands now, it's
  !          entirely dependent on RRTMG legacy code.
  use module_radiation_aerosols, only: &
       NF_AESW,                  & ! Number of optical-fields in SW output (3=tau+g+omega)
       NF_AELW,                  & ! Number of optical-fields in LW output (3=tau+g+omega)
       setaer,                   & ! Routine to compute aerosol radiative properties (tau,g,omega)
       NSPC1                        ! Number of species for vertically integrated aerosol optical-depth
  use surface_perturbation, only: & 
       cdfnor                      ! Routine to compute CDF (used to compute percentiles)
  use mo_gas_optics_rrtmgp,  only: &
       ty_gas_optics_rrtmgp
  public GFS_rrtmgp_sw_pre_run,GFS_rrtmgp_sw_pre_init,GFS_rrtmgp_sw_pre_finalize
  
contains

  subroutine GFS_rrtmgp_sw_pre_init ()
  end subroutine GFS_rrtmgp_sw_pre_init

!> \section arg_table_GFS_rrtmgp_sw_pre_run
!! \htmlinclude GFS_rrtmgp_sw_pre.html
!!
  subroutine GFS_rrtmgp_sw_pre_run (Model, Interstitial, Grid,   Sfcprop, Statein, ncol, p_lay, &
       p_lev, tv_lay, relhum, tracer, sw_gas_props, nday, idxday, alb1d, RadTend, &
       Coupling, aerosolssw, aerodp, errmsg, errflg)
    
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                ! Fortran DDT containing FV3-GFS model control parameters
    type(GFS_Interstitial_type),intent(inout) :: &
         Interstitial
    type(GFS_grid_type), intent(in) :: &
         Grid                 ! Fortran DDT containing FV3-GFS grid and interpolation related data 
    type(GFS_sfcprop_type), intent(in) :: &
         Sfcprop              ! Fortran DDT containing FV3-GFS surface fields
    type(GFS_statein_type), intent(in) :: &
         Statein              ! Fortran DDT containing FV3-GFS prognostic state data in from dycore    
    integer, intent(in)    :: &
         ncol                 ! Number of horizontal grid points
    real(kind_phys), dimension(ncol,Model%levs),intent(in) :: &
         p_lay,             & ! Layer pressure
         tv_lay,            & ! Layer virtual-temperature
         relhum               ! Layer relative-humidity
    real(kind_phys), dimension(ncol, Model%levs, 2:Model%ntrac),intent(in) :: &
         tracer
    real(kind_phys), dimension(ncol,Model%levs+1),intent(in) :: &
         p_lev                ! Interface (level) pressure
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props         ! RRTMGP DDT containing spectral information for SW calculation

    ! Outputs
    integer, intent(out)   :: &
         nday                 ! Number of daylit points
    integer, dimension(ncol), intent(out) :: &
         idxday               ! Indices for daylit points
    real(kind_phys), dimension(ncol), intent(out) :: &
         alb1d                ! Surface albedo pertubation
    type(GFS_radtend_type), intent(inout) :: &
         Radtend              ! Fortran DDT containing FV3-GFS radiation tendencies 
    type(GFS_coupling_type), intent(inout) :: &
         Coupling
    real(kind_phys), dimension(ncol,Model%levs,sw_gas_props%get_nband(),NF_AESW), intent(out) ::&
         aerosolssw               ! Aerosol radiative properties in each SW band.
    real(kind_phys), dimension(ncol,NSPC1), intent(inout) :: &
         aerodp               ! Vertical integrated optical depth for various aerosol species  
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag

    ! Local variables
    integer :: i, j, iCol, iBand, iSFC, iTOA, iLay
    real(kind_phys), dimension(ncol, NF_ALBD) :: sfcalb
    real(kind_phys), dimension(ncol, Model%levs, sw_gas_props%get_nband(), NF_AESW) :: &
         aerosolssw2
    real(kind_phys), dimension(ncol, Model%levs, Model%rrtmgp_nBandsLW, NF_AELW) :: &
         aerosolslw
    logical :: top_at_1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. Model%lsswr) return

    ! #######################################################################################
    ! What is vertical ordering?
    ! #######################################################################################
    top_at_1 = (Statein%prsi(1,1) .lt.  Statein%prsi(1, Model%levs))
    if (top_at_1) then 
       iSFC = Model%levs
       iTOA = 1
    else
       iSFC = 1
       iTOA = Model%levs
    endif
    
    ! #######################################################################################
    ! Compute cosine of zenith angle (only when SW is called)
    ! #######################################################################################
    call coszmn (Grid%xlon, Grid%sinlat, Grid%coslat, Model%solhr, NCOL, Model%me, &
         Radtend%coszen, Radtend%coszdg)

    ! #######################################################################################
    ! For SW, gather daylit points, compute surface albedo in each band,
    ! #######################################################################################
    ! Check for daytime points for SW radiation.
    nday = 0
    idxday = 0
    do i = 1, NCOL
       if (Radtend%coszen(i) >= 0.0001) then
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
    if (Model%do_sfcperts) then
       if (Model%pertalb(1) > 0.) then
          do i=1,ncol
             call cdfnor(Coupling%sfc_wts(i,5),alb1d(i))
          enddo
       endif
    endif  
    
    ! Call module_radiation_surface::setalb() to setup surface albedo.
    call setalb (Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr, Sfcprop%snoalb, Sfcprop%zorl, &
         Radtend%coszen, Sfcprop%tsfc, Sfcprop%tsfc, Sfcprop%hprime(:,1), Sfcprop%alvsf,           &
         Sfcprop%alnsf, Sfcprop%alvwf, Sfcprop%alnwf, Sfcprop%facsf, Sfcprop%facwf,          &
         Sfcprop%fice, Sfcprop%tisfc, NCOL, alb1d, Model%pertalb, sfcalb)
       
    ! Approximate mean surface albedo from vis- and nir-  diffuse values.
    Radtend%sfalb(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
  
    ! Spread across all SW bands
    do iBand=1,sw_gas_props%get_nband()
       Interstitial%sfc_alb_nir_dir(iBand,1:NCOL)   = sfcalb(1:NCOL,1)
       Interstitial%sfc_alb_nir_dif(iBand,1:NCOL)   = sfcalb(1:NCOL,2)
       Interstitial%sfc_alb_uvvis_dir(iBand,1:NCOL) = sfcalb(1:NCOL,3)
       Interstitial%sfc_alb_uvvis_dif(iBand,1:NCOL) = sfcalb(1:NCOL,4)
    enddo 

    ! #######################################################################################
    ! Call module_radiation_aerosols::setaer(),to setup aerosols property profile
    ! #######################################################################################
    call setaer(p_lev, p_lay, Statein%prslk(1:NCOL,iSFC:iTOA), tv_lay, relhum,              &
         Sfcprop%slmsk,  tracer, Grid%xlon, Grid%xlat, NCOL, Model%levs, Model%levs+1,      &
         Model%lsswr, .true., aerosolssw2, aerosolslw, aerodp)
    
    ! Store aerosol optical properties
    ! SW. 
    ! For RRTMGP SW the bands are now ordered from [IR(band) -> nIR -> UV], in RRTMG the 
    ! band ordering was [nIR -> UV -> IR(band)]
    aerosolssw(1:NCOL,1:Model%levs,1,1)                          = aerosolssw2(1:NCOL,1:Model%levs,sw_gas_props%get_nband(),1)
    aerosolssw(1:NCOL,1:Model%levs,1,2)                          = aerosolssw2(1:NCOL,1:Model%levs,sw_gas_props%get_nband(),2)
    aerosolssw(1:NCOL,1:Model%levs,1,3)                          = aerosolssw2(1:NCOL,1:Model%levs,sw_gas_props%get_nband(),3)
    aerosolssw(1:NCOL,1:Model%levs,2:sw_gas_props%get_nband(),1) = aerosolssw2(1:NCOL,1:Model%levs,1:sw_gas_props%get_nband()-1,1)
    aerosolssw(1:NCOL,1:Model%levs,2:sw_gas_props%get_nband(),2) = aerosolssw2(1:NCOL,1:Model%levs,1:sw_gas_props%get_nband()-1,2)
    aerosolssw(1:NCOL,1:Model%levs,2:sw_gas_props%get_nband(),3) = aerosolssw2(1:NCOL,1:Model%levs,1:sw_gas_props%get_nband()-1,3)

  end subroutine GFS_rrtmgp_sw_pre_run
  
!> \section arg_table_GFS_rrtmgp_sw_pre_finalize Argument Table
!!
  subroutine GFS_rrtmgp_sw_pre_finalize ()
  end subroutine GFS_rrtmgp_sw_pre_finalize

end module GFS_rrtmgp_sw_pre
