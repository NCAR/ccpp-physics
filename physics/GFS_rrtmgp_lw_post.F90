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
!> \section arg_table_GFS_rrtmgp_lw_post
!! \htmlinclude GFS_rrtmgp_lw_post.html
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
