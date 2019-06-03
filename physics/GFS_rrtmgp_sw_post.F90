!>\file GFS_rrtmgp_sw_post
!!This file contains
module GFS_rrtmgp_sw_post 
  use machine,                    only: kind_phys
  use GFS_typedefs,               only: GFS_statein_type,   &
                                        GFS_coupling_type,  &
                                        GFS_control_type,   &
                                        GFS_grid_type,      &
                                        GFS_radtend_type,   &
                                        GFS_diag_type
  use module_radiation_aerosols, only: NSPC1
  use module_radsw_parameters,   only: topfsw_type, sfcfsw_type, profsw_type, cmpfsw_type
  ! RRTMGP DDT's
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_fluxes_byband,          only: ty_fluxes_byband
  use mo_heating_rates,          only: compute_heating_rate
  use rrtmgp_aux,                only: check_error_msg
  implicit none
  
  public GFS_rrtmgp_sw_post_init,GFS_rrtmgp_sw_post_run,GFS_rrtmgp_sw_post_finalize

contains

  subroutine GFS_rrtmgp_sw_post_init()
  end subroutine GFS_rrtmgp_sw_post_init

  ! PGI compiler does not accept lines longer than 264 characters, remove during pre-processing
#ifndef __PGI
!> \section arg_table_GFS_rrtmgp_sw_post_run Argument Table
!! | local_name        | standard_name                                                                                  | long_name                                                                    | units    | rank |  type                |   kind    | intent | optional |
!! |-------------------|------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------|----------|------|----------------------|-----------|--------|----------|
!! | Model             | GFS_control_type_instance                                                                      | Fortran DDT containing FV3-GFS model control parameters                      | DDT      |    0 | GFS_control_type     |           | in     | F        |
!! | Grid              | GFS_grid_type_instance                                                                         | Fortran DDT containing FV3-GFS grid and interpolation related data           | DDT      |    0 | GFS_grid_type        |           | in     | F        |
!! | Diag              | GFS_diag_type_instance                                                                         | Fortran DDT containing FV3-GFS diagnotics data                               | DDT      |    0 | GFS_diag_type        |           | inout  | F        |
!! | Radtend           | GFS_radtend_type_instance                                                                      | Fortran DDT containing FV3-GFS radiation tendencies                          | DDT      |    0 | GFS_radtend_type     |           | inout  | F        |
!! | Statein           | GFS_statein_type_instance                                                                      | Fortran DDT containing FV3-GFS prognostic state data in from dycore          | DDT      |    0 | GFS_statein_type     |           | in     | F        |
!! | Coupling          | GFS_coupling_type_instance                                                                     | Fortran DDT containing FV3-GFS fields to/from coupling with other components | DDT      |    0 | GFS_coupling_type    |           | inout  | F        |
!! | scmpsw            | components_of_surface_downward_shortwave_fluxes                                                | derived type for special components of surface downward shortwave fluxes     | W m-2    |    1 | cmpfsw_type          |           | inout  | T        |
!! | im                | horizontal_loop_extent                                                                         | horizontal loop extent                                                       | count    |    0 | integer              |           | in     | F        |
!! | tsfa              | surface_air_temperature_for_radiation                                                          | lowest model layer air temperature for radiation                             | K        |    1 | real                 | kind_phys | in     | F        |
!! | p_lev             | air_pressure_at_interface_for_RRTMGP_in_hPa                                                    | air pressure level                                                           | hPa      |    2 | real                 | kind_phys | in     | F        |
!! | nday              | daytime_points_dimension                                                                       | daytime points dimension                                                     | count    |    0 | integer              |           | in     | F        |
!! | idxday            | daytime_points                                                                                 | daytime points                                                               | index    |    1 | integer              |           | in     | F        |
!! | fluxswUP_allsky   | sw_flux_profile_upward_allsky                                                                  | RRTMGP upward shortwave all-sky flux profile                                 | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxswDOWN_allsky | sw_flux_profile_downward_allsky                                                                | RRTMGP downward shortwave all-sky flux profile                               | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxswUP_clrsky   | sw_flux_profile_upward_clrsky                                                                  | RRTMGP upward shortwave clr-sky flux profile                                 | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxswDOWN_clrsky | sw_flux_profile_downward_clrsky                                                                | RRTMGP downward shortwave clr-sky flux profile                               | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | sw_gas_props      | coefficients_for_sw_gas_optics                                                                 | DDT containing spectral information for RRTMGP SW radiation scheme           | DDT      |    0 | ty_gas_optics_rrtmgp |           | in     | F        |
!! | sfc_alb_nir_dir   | surface_shortwave_albedo_near_infrared_direct_in_each_band                                     | surface sw near-infrared direct albedo in each SW band                       | frac     |    2 | real                 | kind_phys | in     | F        |
!! | sfc_alb_nir_dif   | surface_shortwave_albedo_near_infrared_diffuse_in_each_band                                    | surface sw near-infrared diffuse albedo in each SW band                      | frac     |    2 | real                 | kind_phys | in     | F        |
!! | sfc_alb_uvvis_dir | surface_shortwave_albedo_uv_visible_direct_in_each_band                                        | surface sw uv-visible direct albedo in each SW band                          | frac     |    2 | real                 | kind_phys | in     | F        |
!! | sfc_alb_uvvis_dif | surface_shortwave_albedo_uv_visible_diffuse_in_each_band                                       | surface sw uv-visible diffuse albedo in each SW band                         | frac     |    2 | real                 | kind_phys | in     | F        |
!! | hswc              | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step                    | shortwave total sky heating rate                                             | K s-1    |    2 | real                 | kind_phys | out    | F        |
!! | topflx_sw         | sw_fluxes_top_atmosphere                                                                       | shortwave total sky fluxes at the top of the atm                             | W m-2    |    1 | topfsw_type          |           | inout  | F        |
!! | sfcflx_sw         | sw_fluxes_sfc                                                                                  | shortwave total sky fluxes at the Earth surface                              | W m-2    |    1 | sfcfsw_type          |           | inout  | F        |
!! | flxprf_sw         | sw_fluxes                                                                                      | sw fluxes total sky / csk and up / down at levels                            | W m-2    |    2 | profsw_type          |           | inout  | T        |
!! | hsw0              | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | shortwave clear sky heating rate                                             | K s-1    |    2 | real                 | kind_phys | inout  | T        |
!! | errmsg            | ccpp_error_message                                                                             | error message for error handling in CCPP                                     | none     |    0 | character            | len=*     | out    | F        |
!! | errflg            | ccpp_error_flag                                                                                | error flag for error handling in CCPP                                        | flag     |    0 | integer              |           | out    | F        |
!!
#endif
  subroutine GFS_rrtmgp_sw_post_run (Model, Grid, Diag, Radtend, Statein, Coupling,       & 
       scmpsw, im, p_lev, sw_gas_props, sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir,  &
       sfc_alb_uvvis_dif, tsfa, nday, idxday, fluxswUP_allsky, fluxswDOWN_allsky,         &
       fluxswUP_clrsky, fluxswDOWN_clrsky, hswc, topflx_sw, sfcflx_sw, flxprf_sw, hsw0,   &
       errmsg, errflg)

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model             ! Fortran DDT containing FV3-GFS model control parameters
    type(GFS_grid_type), intent(in) :: &
         Grid              ! Fortran DDT containing FV3-GFS grid and interpolation related data 
    type(GFS_statein_type), intent(in) :: &
         Statein           ! Fortran DDT containing FV3-GFS prognostic state data in from dycore    
    type(GFS_coupling_type), intent(inout) :: &
         Coupling          ! Fortran DDT containing FV3-GFS fields to/from coupling with other components 
    type(GFS_radtend_type), intent(inout) :: &
         Radtend           ! Fortran DDT containing FV3-GFS radiation tendencies 
    type(GFS_diag_type), intent(inout) :: &
         Diag              ! Fortran DDT containing FV3-GFS diagnotics data  
    integer, intent(in) :: &
         im,             & ! Horizontal loop extent 
         nDay              ! Number of daylit columns
    integer, intent(in), dimension(nday) :: &
         idxday            ! Index array for daytime points
    real(kind_phys), dimension(size(Grid%xlon,1)), intent(in) ::  &
         tsfa              ! Lowest model layer air temperature for radiation 
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props          ! DDT containing SW spectral information
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         p_lev             ! Pressure @ model layer-interfaces    (hPa)
    real(kind_phys),dimension(sw_gas_props%get_nband(),size(Grid%xlon,1)),intent(in) :: &
         sfc_alb_nir_dir,   & ! Shortwave surface albedo (nIR-direct) 
         sfc_alb_nir_dif,   & ! Shortwave surface albedo (nIR-diffuse)
         sfc_alb_uvvis_dir, & ! Shortwave surface albedo (uvvis-direct)
         sfc_alb_uvvis_dif    ! Shortwave surface albedo (uvvis-diffuse)    
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         fluxswUP_allsky,   & ! SW All-sky flux                    (W/m2)
         fluxswDOWN_allsky, & ! SW All-sky flux                    (W/m2)
         fluxswUP_clrsky,   & ! SW Clear-sky flux                  (W/m2)
         fluxswDOWN_clrsky    ! SW All-sky flux                    (W/m2)

    ! Outputs (mandatory)
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg
    real(kind_phys),dimension(size(Grid%xlon,1), Model%levs),intent(out) :: &
         hswc             ! Shortwave all-sky heating-rate         (K/sec)
    type(topfsw_type), dimension(size(Grid%xlon,1)), intent(inout) :: &
         topflx_sw        ! radiation fluxes at top, components:
                          ! upfxc - total sky upward flux at top   (w/m2)
                          ! upfx0 - clear sky upward flux at top   (w/m2)
    type(sfcfsw_type), dimension(size(Grid%xlon,1)), intent(inout) :: &
         sfcflx_sw        ! radiation fluxes at sfc, components:
                          ! upfxc - total sky upward flux at sfc   (w/m2)  
                          ! upfx0 - clear sky upward flux at sfc   (w/m2)
                          ! dnfxc - total sky downward flux at sfc (w/m2)
                          ! dnfx0 - clear sky downward flux at sfc (w/m2)
    
    ! Outputs (optional)
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs), optional, intent(inout) :: &
         hsw0             ! Shortwave clear-sky heating-rate         (K/sec)
    type(profsw_type), dimension(size(Grid%xlon,1), Model%levs+1), intent(inout), optional :: &
         flxprf_sw        ! 2D radiative fluxes, components:
                          ! upfxc - total sky upward flux            (W/m2)
                          ! dnfxc - total sky dnward flux            (W/m2)
                          ! upfx0 - clear sky upward flux            (W/m2)
                          ! dnfx0 - clear sky dnward flux            (W/m2)
    type(cmpfsw_type), dimension(size(Grid%xlon,1)), intent(inout), optional :: &
         scmpsw           ! 2D surface fluxes, components:
                          ! uvbfc - total sky downward uv-b flux at  (W/m2)
                          ! uvbf0 - clear sky downward uv-b flux at  (W/m2)
                          ! nirbm - downward nir direct beam flux    (W/m2)
                          ! nirdf - downward nir diffused flux       (W/m2)
                          ! visbm - downward uv+vis direct beam flux (W/m2)
                          ! visdf - downward uv+vis diffused flux    (W/m2)
    ! Local variables
    integer :: i, j, k, k1, itop, ibtc, iBand, iSFC, iTOA
    real(kind_phys) :: tem0d, tem1, tem2
    real(kind_phys), dimension(nDay, Model%levs) :: thetaTendClrSky, thetaTendAllSky
    logical :: l_clrskysw_hr, l_fluxessw2d, top_at_1, l_sfcFluxessw1D

    ! Are any optional outputs requested?
    l_clrskysw_hr   = present(hsw0)
    l_fluxessw2d    = present(flxprf_sw)
    l_sfcfluxessw1D = present(scmpsw)

    ! #######################################################################################
    ! What is vertical ordering?
    ! #######################################################################################
    top_at_1 = (p_lev(1,1) .lt. p_lev(1, Model%levs))
    if (top_at_1) then 
       iSFC = Model%levs
       iTOA = 1
    else
       iSFC = 1
       iTOA = Model%levs
    endif
    ! #######################################################################################
    ! Compute SW heating-rates
    ! #######################################################################################
    ! Initialize outputs
    hswc(:,:) = 0.
    topflx_sw = topfsw_type ( 0., 0., 0. )
    sfcflx_sw = sfcfsw_type ( 0., 0., 0., 0. )
    if (l_clrskysw_hr) then
       hsw0(:,:) = 0.
    endif
    if (l_fluxessw2D) then
       flxprf_sw = profsw_type ( 0., 0., 0., 0. )
    endif
    if (l_sfcfluxessw1D) then
       scmpsw = cmpfsw_type (0.,0.,0.,0.,0.,0.)
    endif

    if (Model%lsswr .and. nDay .gt. 0) then
       ! Clear-sky heating-rate (optional)
       if (l_clrskysw_HR) then
          call check_error_msg('GFS_rrtmgp_post',compute_heating_rate( &
               fluxswUP_clrsky,                &
               fluxswDOWN_clrsky,                &
               p_lev(idxday,1:Model%levs+1),     &
               thetaTendClrSky))
          hsw0(idxday,:)=thetaTendClrSky
       endif
       ! All-sky heating-rate (mandatory)
       call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(    &
            fluxswUP_allsky,                   &
            fluxswDOWN_allsky,                   &
            p_lev(idxday,1:Model%levs+1),        &
            thetaTendAllSky))
       hswc(idxday,:) = thetaTendAllSky
       
       ! Copy fluxes from RRTGMP types into model radiation types.
       ! Mandatory outputs
       topflx_sw%upfxc = fluxswUP_allsky(:,iTOA)
       topflx_sw%upfx0 = fluxswUP_clrsky(:,iTOA)
       sfcflx_sw%upfxc = fluxswUP_allsky(:,iSFC)
       sfcflx_sw%upfx0 = fluxswUP_clrsky(:,iSFC)
       sfcflx_sw%dnfxc = fluxswDOWN_allsky(:,iSFC)
       sfcflx_sw%dnfx0 = fluxswDOWN_clrsky(:,iSFC)

       ! Optional output
       if(l_fluxessw2D) then
          flxprf_sw%upfxc = fluxswUP_allsky
          flxprf_sw%dnfxc = fluxswDOWN_allsky
          flxprf_sw%upfx0 = fluxswUP_clrsky
          flxprf_sw%dnfx0 = fluxswDOWN_clrsky
       endif
    endif

    ! #######################################################################################
    ! Save SW outputs
    ! #######################################################################################
    if (Model%lsswr) then
       if (nday > 0) then
          ! All-sky heating rate
          do k = 1, Model%levs
             Radtend%htrsw(1:im,k) = hswc(1:im,k)
          enddo
          ! Clear-sk heating rate
          if (Model%swhtr) then
             do k = 1, Model%levs
                Radtend%swhc(1:im,k) = hsw0(1:im,k)
             enddo
          endif
          
          ! Surface down and up spectral component fluxes
          ! - Save two spectral bands' surface downward and upward fluxes for output.
          do i=1,im
             Coupling%nirbmdi(i) = scmpsw(i)%nirbm
             Coupling%nirdfdi(i) = scmpsw(i)%nirdf
             Coupling%visbmdi(i) = scmpsw(i)%visbm
             Coupling%visdfdi(i) = scmpsw(i)%visdf
             
             Coupling%nirbmui(i) = scmpsw(i)%nirbm * sfc_alb_nir_dir(1,i)
             Coupling%nirdfui(i) = scmpsw(i)%nirdf * sfc_alb_nir_dif(1,i)
             Coupling%visbmui(i) = scmpsw(i)%visbm * sfc_alb_uvvis_dir(1,i)
             Coupling%visdfui(i) = scmpsw(i)%visdf * sfc_alb_uvvis_dif(1,i)
          enddo
       else                   ! if_nday_block
          Radtend%htrsw(:,:) = 0.0
          Radtend%sfcfsw     = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          Diag%topfsw        = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw             = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
          
          do i=1,im
             Coupling%nirbmdi(i) = 0.0
             Coupling%nirdfdi(i) = 0.0
             Coupling%visbmdi(i) = 0.0
             Coupling%visdfdi(i) = 0.0
             
             Coupling%nirbmui(i) = 0.0
             Coupling%nirdfui(i) = 0.0
             Coupling%visbmui(i) = 0.0
             Coupling%visdfui(i) = 0.0
          enddo
          
          if (Model%swhtr) then
             Radtend%swhc(:,:) = 0
          endif
       endif                  ! end_if_nday
       
       ! Radiation fluxes for other physics processes
       do i=1,im
          Coupling%sfcnsw(i) = Radtend%sfcfsw(i)%dnfxc - Radtend%sfcfsw(i)%upfxc
          Coupling%sfcdsw(i) = Radtend%sfcfsw(i)%dnfxc
       enddo  
    endif                                ! end_if_lsswr

 
  end subroutine GFS_rrtmgp_sw_post_run

  subroutine GFS_rrtmgp_sw_post_finalize ()
  end subroutine GFS_rrtmgp_sw_post_finalize

end module GFS_rrtmgp_sw_post
