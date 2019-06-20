!>\file GFS_rrtmgp_lw_post
!!This file contains
module GFS_rrtmgp_lw_post 
  use machine,                    only: kind_phys
  use GFS_typedefs,               only: GFS_coupling_type,  &
                                        GFS_control_type,   &
                                        GFS_grid_type,      &
                                        GFS_radtend_type
  use module_radiation_aerosols, only: NSPC1
  use module_radlw_parameters,   only: topflw_type, sfcflw_type, proflw_type
  ! RRTMGP DDT's
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_fluxes_byband,          only: ty_fluxes_byband
  use mo_heating_rates,          only: compute_heating_rate
  use rrtmgp_aux,                only: check_error_msg
  implicit none
  
  public GFS_rrtmgp_lw_post_init,GFS_rrtmgp_lw_post_run,GFS_rrtmgp_lw_post_finalize

contains

  subroutine GFS_rrtmgp_lw_post_init()
  end subroutine GFS_rrtmgp_lw_post_init

  ! PGI compiler does not accept lines longer than 264 characters, remove during pre-processing
#ifndef __PGI
!> \section arg_table_GFS_rrtmgp_lw_post_run Argument Table
!! | local_name        | standard_name                                                                                  | long_name                                                                    | units    | rank |  type                |   kind    | intent | optional |
!! |-------------------|------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------|----------|------|----------------------|-----------|--------|----------|
!! | Model             | GFS_control_type_instance                                                                      | Fortran DDT containing FV3-GFS model control parameters                      | DDT      |    0 | GFS_control_type     |           | in     | F        |
!! | Grid              | GFS_grid_type_instance                                                                         | Fortran DDT containing FV3-GFS grid and interpolation related data           | DDT      |    0 | GFS_grid_type        |           | in     | F        |
!! | Radtend           | GFS_radtend_type_instance                                                                      | Fortran DDT containing FV3-GFS radiation tendencies                          | DDT      |    0 | GFS_radtend_type     |           | inout  | F        |
!! | Coupling          | GFS_coupling_type_instance                                                                     | Fortran DDT containing FV3-GFS fields to/from coupling with other components | DDT      |    0 | GFS_coupling_type    |           | inout  | F        |
!! | im                | horizontal_loop_extent                                                                         | horizontal loop extent                                                       | count    |    0 | integer              |           | in     | F        |
!! | tsfa              | surface_air_temperature_for_radiation                                                          | lowest model layer air temperature for radiation                             | K        |    1 | real                 | kind_phys | in     | F        |
!! | p_lev             | air_pressure_at_interface_for_RRTMGP_in_hPa                                                    | air pressure level                                                           | hPa      |    2 | real                 | kind_phys | in     | F        |
!! | fluxlwUP_allsky   | lw_flux_profile_upward_allsky                                                                  | RRTMGP upward longwave all-sky flux profile                                  | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxlwDOWN_allsky | lw_flux_profile_downward_allsky                                                                | RRTMGP downward longwave all-sky flux profile                                | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxlwUP_clrsky   | lw_flux_profile_upward_clrsky                                                                  | RRTMGP upward longwave clr-sky flux profile                                  | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxlwDOWN_clrsky | lw_flux_profile_downward_clrsky                                                                | RRTMGP downward longwave clr-sky flux profile                                | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | hlwc              | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step                     | longwave total sky heating rate                                              | K s-1    |    2 | real                 | kind_phys | out    | F        |
!! | topflx_lw         | lw_fluxes_top_atmosphere                                                                       | longwave total sky fluxes at the top of the atm                              | W m-2    |    1 | topflw_type          |           | inout  | F        |
!! | sfcflx_lw         | lw_fluxes_sfc                                                                                  | longwave total sky fluxes at the Earth surface                               | W m-2    |    1 | sfcflw_type          |           | inout  | F        |
!! | flxprf_lw         | lw_fluxes                                                                                      | lw fluxes total sky / csk and up / down at levels                            | W m-2    |    2 | proflw_type          |           | inout  | T        |
!! | hlw0              | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step  | longwave clear sky heating rate                                              | K s-1    |    2 | real                 | kind_phys | inout  | T        |
!! | errmsg            | ccpp_error_message                                                                             | error message for error handling in CCPP                                     | none     |    0 | character            | len=*     | out    | F        |
!! | errflg            | ccpp_error_flag                                                                                | error flag for error handling in CCPP                                        | flag     |    0 | integer              |           | out    | F        |
!!
#endif
  subroutine GFS_rrtmgp_lw_post_run (Model, Grid, Radtend,  &
              Coupling, im, p_lev,          &
              tsfa, fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky, fluxlwDOWN_clrsky, &
              hlwc, topflx_lw, sfcflx_lw, flxprf_lw, hlw0, errmsg, errflg)

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model             ! Fortran DDT containing FV3-GFS model control parameters
    type(GFS_grid_type), intent(in) :: &
         Grid              ! Fortran DDT containing FV3-GFS grid and interpolation related data 
   type(GFS_coupling_type), intent(inout) :: &
         Coupling          ! Fortran DDT containing FV3-GFS fields to/from coupling with other components 
    type(GFS_radtend_type), intent(inout) :: &
         Radtend           ! Fortran DDT containing FV3-GFS radiation tendencies 
    integer, intent(in) :: &
         im                ! Horizontal loop extent 
    real(kind_phys), dimension(size(Grid%xlon,1)), intent(in) ::  &
         tsfa              ! Lowest model layer air temperature for radiation 
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         p_lev             ! Pressure @ model layer-interfaces    (hPa)
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         fluxlwUP_allsky,   & ! LW All-sky flux                    (W/m2)
         fluxlwDOWN_allsky, & ! LW All-sky flux                    (W/m2)
         fluxlwUP_clrsky,   & ! LW Clear-sky flux                  (W/m2)
         fluxlwDOWN_clrsky    ! LW All-sky flux                    (W/m2)

    ! Outputs (mandatory)
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg
    real(kind_phys),dimension(size(Grid%xlon,1), Model%levs),intent(out) :: &
         hlwc             ! Longwave all-sky heating-rate          (K/sec)
    type(topflw_type), dimension(size(Grid%xlon,1)), intent(inout) :: &
         topflx_lw        ! radiation fluxes at top, components:
                          ! upfxc - total sky upward flux at top   (w/m2)
                          ! upfx0 - clear sky upward flux at top   (w/m2)
    type(sfcflw_type), dimension(size(Grid%xlon,1)), intent(inout) :: &
         sfcflx_lw        ! radiation fluxes at sfc, components:
                          ! upfxc - total sky upward flux at sfc   (w/m2)  
                          ! upfx0 - clear sky upward flux at sfc   (w/m2)
                          ! dnfxc - total sky downward flux at sfc (w/m2)
                          ! dnfx0 - clear sky downward flux at sfc (w/m2)
 
    ! Outputs (optional)
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs), optional, intent(inout) :: &
         hlw0             ! Longwave clear-sky heating rate          (K/sec)
    type(proflw_type), dimension(size(Grid%xlon,1), Model%levs+1), optional, intent(inout) :: &
         flxprf_lw        ! 2D radiative fluxes, components:
                          ! upfxc - total sky upward flux            (W/m2)
                          ! dnfxc - total sky dnward flux            (W/m2)
                          ! upfx0 - clear sky upward flux            (W/m2)
                          ! dnfx0 - clear sky dnward flux            (W/m2)

    ! Local variables
    integer :: k, iSFC, iTOA
    logical :: l_clrskylw_hr, l_fluxeslw2d, top_at_1

   ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. Model%lslwr) return

    ! Are any optional outputs requested?
    l_clrskylw_hr   = present(hlw0)
    l_fluxeslw2d    = present(flxprf_lw)

    ! #######################################################################################
    ! What is vertical ordering?
    ! #######################################################################################
    top_at_1 = (p_lev(1,1) .lt. p_lev(1, Model%levs))
    if (top_at_1) then 
       iSFC = Model%levs+1
       iTOA = 1
    else
       iSFC = 1
       iTOA = Model%levs+1
    endif

    ! #######################################################################################
    ! Compute LW heating-rates. 
    ! #######################################################################################
    if (Model%lslwr) then
       ! Clear-sky heating-rate (optional)
       if (l_clrskylw_hr) then
          call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(  &
               fluxlwUP_clrsky,                 &
               fluxlwDOWN_clrsky,               &
               p_lev,                           &
               hlw0))
       endif
       ! All-sky heating-rate (mandatory)
       call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(     &
            fluxlwUP_allsky,                    &
            fluxlwDOWN_allsky,                  &
            p_lev,                              &
            hlwc))
       
       ! Copy fluxes from RRTGMP types into model radiation types.
       ! Mandatory outputs
       topflx_lw%upfxc = fluxlwUP_allsky(:,iTOA)
       topflx_lw%upfx0 = fluxlwUP_clrsky(:,iTOA)
       sfcflx_lw%upfxc = fluxlwUP_allsky(:,iSFC)
       sfcflx_lw%upfx0 = fluxlwUP_clrsky(:,iSFC)
       sfcflx_lw%dnfxc = fluxlwDOWN_allsky(:,iSFC)
       sfcflx_lw%dnfx0 = fluxlwDOWN_clrsky(:,iSFC)
       
       ! Optional outputs
       if(l_fluxeslw2d) then
          flxprf_lw%upfxc = fluxlwUP_allsky
          flxprf_lw%dnfxc = fluxlwDOWN_allsky
          flxprf_lw%upfx0 = fluxlwUP_clrsky
          flxprf_lw%dnfx0 = fluxlwDOWN_clrsky
       endif
    endif
    
    ! #######################################################################################
    !  Save LW outputs.
    ! #######################################################################################
    if (Model%lslwr) then
       ! Save surface air temp for diurnal adjustment at model t-steps
       Radtend%tsflw (:) = tsfa(:)
       
       ! All-sky heating rate profile
       do k = 1, model%levs
          Radtend%htrlw(1:im,k) = hlwc(1:im,k)
       enddo
       if (Model%lwhtr) then
          do k = 1, model%levs
             Radtend%lwhc(1:im,k) = hlw0(1:im,k)
          enddo
       endif
       
       ! Radiation fluxes for other physics processes
       Coupling%sfcdlw(:) = Radtend%sfcflw(:)%dnfxc
    endif 

  end subroutine GFS_rrtmgp_lw_post_run

  subroutine GFS_rrtmgp_lw_post_finalize ()
  end subroutine GFS_rrtmgp_lw_post_finalize

end module GFS_rrtmgp_lw_post
