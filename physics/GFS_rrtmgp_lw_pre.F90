!> \file GFS_rrtmgp_lw_pre.f90
!! This file contains
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
  use module_radiation_aerosols, only: &
       NF_AESW,                  & ! Number of optical-fields in SW output (3=tau+g+omega)
       NF_AELW,                  & ! Number of optical-fields in LW output (3=tau+g+omega)
       setaer,                   & ! Routine to compute aerosol radiative properties (tau,g,omega)
       NSPC1                       ! Number of species for vertically integrated aerosol optical-depth
  use mo_gas_optics_rrtmgp,  only: &
       ty_gas_optics_rrtmgp

  public GFS_rrtmgp_lw_pre_run,GFS_rrtmgp_lw_pre_init,GFS_rrtmgp_lw_pre_finalize
  
contains

  subroutine GFS_rrtmgp_lw_pre_init ()
  end subroutine GFS_rrtmgp_lw_pre_init

!> \section arg_table_GFS_rrtmgp_lw_pre_run
!! \htmlinclude GFS_rrtmgp_lw_pre.html
!!
  subroutine GFS_rrtmgp_lw_pre_run (Model, Grid,  Sfcprop, Statein, ncol,  p_lay, p_lev, &
       tv_lay, relhum, tracer, lw_gas_props, Radtend, aerosolslw, aerodp, errmsg, errflg)
    
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                ! Fortran DDT containing FV3-GFS model control parameters
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
         lw_gas_props         ! RRTMGP DDT containing spectral information for LW calculation

    ! Outputs
    type(GFS_radtend_type), intent(inout) :: &
         Radtend              ! Fortran DDT containing FV3-GFS radiation tendencies 
    real(kind_phys), dimension(ncol,Model%levs,lw_gas_props%get_nband(),NF_AELW), intent(out) ::&
         aerosolslw               ! Aerosol radiative properties in each SW band.
    real(kind_phys), dimension(ncol,NSPC1), intent(inout) :: &
         aerodp               ! Vertical integrated optical depth for various aerosol species  
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag

    ! Local
    integer :: iSFC, iTOA
    logical :: top_at_1
    real(kind_phys), dimension(ncol, Model%levs, Model%rrtmgp_nBandsSW, NF_AESW) :: &
         aerosolssw2
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. Model%lslwr) return

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
    ! Call module_radiation_surface::setemis(),to setup surface emissivity for LW radiation.
    ! #######################################################################################
    call setemis (Grid%xlon, Grid%xlat, Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr,    &
         Sfcprop%zorl, Sfcprop%tsfc,Sfcprop%tsfc, Sfcprop%hprim, NCOL,   &
         Radtend%semis)
    do iBand=1,lw_gas_props%get_nband()
       Radtend%sfc_emiss_byband(iBand,1:NCOL) = Radtend%semis(1:NCOL)
    enddo

    ! #######################################################################################
    ! Call module_radiation_aerosols::setaer(),to setup aerosols property profile
    ! #######################################################################################
    call setaer(p_lev, p_lay, Statein%prslk(1:NCOL,iSFC:iTOA), tv_lay, relhum,              &
         Sfcprop%slmsk,  tracer, Grid%xlon, Grid%xlat, ncol, Model%levs, Model%levs+1,      &
         .false., Model%lslwr, aerosolssw2, aerosolslw, aerodp)
    

  end subroutine GFS_rrtmgp_lw_pre_run
  
!> \section arg_table_GFS_rrtmgp_lw_pre_finalize Argument Table
!!
  subroutine GFS_rrtmgp_lw_pre_finalize ()
  end subroutine GFS_rrtmgp_lw_pre_finalize

end module GFS_rrtmgp_lw_pre
