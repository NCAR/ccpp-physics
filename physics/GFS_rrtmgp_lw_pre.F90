module GFS_rrtmgp_lw_pre
  use physparam
  use machine, only: &
       kind_phys                   ! Working type
  use GFS_typedefs, only:        &
       GFS_control_type,         & !
       GFS_sfcprop_type,         & ! Surface fields
       GFS_grid_type,            & ! Grid and interpolation related data
       GFS_statein_type,         & !
       GFS_radtend_type            ! Radiation tendencies needed in physics
  use module_radiation_surface,  only: &
       setemis                     ! Routine to compute surface-emissivity
  use mo_gas_optics_rrtmgp,  only: &
       ty_gas_optics_rrtmgp

  public GFS_rrtmgp_lw_pre_run,GFS_rrtmgp_lw_pre_init,GFS_rrtmgp_lw_pre_finalize
  
contains

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_lw_pre_init
  ! #########################################################################################
  subroutine GFS_rrtmgp_lw_pre_init ()
  end subroutine GFS_rrtmgp_lw_pre_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_lw_pre_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_lw_pre_run
!! \htmlinclude GFS_rrtmgp_lw_pre.html
!!
  subroutine GFS_rrtmgp_lw_pre_run (Model, Grid,  Sfcprop, Statein, ncol,  p_lay, p_lev,    &
       tv_lay, relhum, tracer, lw_gas_props, Radtend, sfc_emiss_byband, errmsg, errflg)
    
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
         tracer               ! trace gas concentrations
    real(kind_phys), dimension(ncol,Model%levs+1),intent(in) :: &
         p_lev                ! Interface (level) pressure
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props         ! RRTMGP DDT: spectral information for LW calculation

    ! Outputs
    type(GFS_radtend_type), intent(inout) :: &
         Radtend              ! DDT: FV3-GFS radiation tendencies 
    real(kind_phys), dimension(lw_gas_props%get_nband(),ncol), intent(out) :: &
         sfc_emiss_byband     ! Surface emissivity in each band
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. Model%lslwr) return

    ! #######################################################################################
    ! Call module_radiation_surface::setemis(),to setup surface emissivity for LW radiation.
    ! #######################################################################################
    call setemis (Grid%xlon, Grid%xlat, Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr,       &
         Sfcprop%zorl, Sfcprop%tsfc,Sfcprop%tsfc, Sfcprop%hprime(:,1), NCOL, Radtend%semis)

    ! Assign same emissivity to all bands
    do iBand=1,lw_gas_props%get_nband()
       sfc_emiss_byband(iBand,1:NCOL) = Radtend%semis(1:NCOL)
    enddo

  end subroutine GFS_rrtmgp_lw_pre_run
  
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_lw_pre_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_lw_pre_finalize ()
  end subroutine GFS_rrtmgp_lw_pre_finalize

end module GFS_rrtmgp_lw_pre
