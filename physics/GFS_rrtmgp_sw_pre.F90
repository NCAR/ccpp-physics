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
  subroutine GFS_rrtmgp_sw_pre_run(Model, Grid, Sfcprop, Statein, ncol, p_lay,  p_lev,      &
       tv_lay, relhum, tracer, sw_gas_props, nday, idxday, alb1d, sfc_alb_nir_dir,          &
       sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, RadTend, Coupling,            &
       errmsg, errflg)
    
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                ! DDT: FV3-GFS model control parameters
    type(GFS_grid_type), intent(in) :: &
         Grid                 ! DDT: FV3-GFS grid and interpolation related data 
    type(GFS_sfcprop_type), intent(in) :: &
         Sfcprop              ! DDT: FV3-GFS surface fields
    type(GFS_statein_type), intent(in) :: &
         Statein              ! DDT: FV3-GFS prognostic state data in from dycore    
    integer, intent(in)    :: &
         ncol                 ! Number of horizontal grid points
    real(kind_phys), dimension(ncol,Model%levs),intent(in) :: &
         p_lay,             & ! Layer pressure
         tv_lay,            & ! Layer virtual-temperature
         relhum               ! Layer relative-humidity
    real(kind_phys), dimension(ncol, Model%levs, 2:Model%ntrac),intent(in) :: &
         tracer
    real(kind_phys), dimension(ncol,Model%levs+1),intent(in) :: &
         p_lev                ! Pressure @ layer interfaces (Pa)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props         ! RRTMGP DDT: spectral information for SW calculation

    ! Outputs
    integer, intent(out)   :: &
         nday                 ! Number of daylit points
    integer, dimension(ncol), intent(out) :: &
         idxday               ! Indices for daylit points
    real(kind_phys), dimension(ncol), intent(out) :: &
         alb1d                ! Surface albedo pertubation
    real(kind_phys), dimension(sw_gas_props%get_nband(),ncol), intent(out) :: &
         sfc_alb_nir_dir,   & ! Surface albedo (direct) 
         sfc_alb_nir_dif,   & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir, & ! Surface albedo (direct)
         sfc_alb_uvvis_dif    ! Surface albedo (diffuse)
    type(GFS_radtend_type), intent(inout) :: &
         Radtend              ! DDT: FV3-GFS radiation tendencies 
    type(GFS_coupling_type), intent(inout) :: &
         Coupling             ! DDT: FV3-GFS coupling arrays
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
    
    if (.not. Model%lsswr) return
    
    ! #######################################################################################
    ! Compute cosine of zenith angle (only when SW is called)
    ! #######################################################################################
    call coszmn (Grid%xlon, Grid%sinlat, Grid%coslat, Model%solhr, NCOL, Model%me, &
         Radtend%coszen, Radtend%coszdg)

    ! #######################################################################################
    ! For SW gather daylit points
    ! #######################################################################################
    nday   = 0
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
    
    ! #######################################################################################
    ! Call module_radiation_surface::setalb() to setup surface albedo.
    ! #######################################################################################
    call setalb (Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr, Sfcprop%snoalb, Sfcprop%zorl, &
         Radtend%coszen, Sfcprop%tsfc, Sfcprop%tsfc, Sfcprop%hprime(:,1), Sfcprop%alvsf,     &
         Sfcprop%alnsf, Sfcprop%alvwf, Sfcprop%alnwf, Sfcprop%facsf, Sfcprop%facwf,          &
         Sfcprop%fice, Sfcprop%tisfc, NCOL, alb1d, Model%pertalb, sfcalb)
       
    ! Approximate mean surface albedo from vis- and nir-  diffuse values.
    Radtend%sfalb(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
  
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
