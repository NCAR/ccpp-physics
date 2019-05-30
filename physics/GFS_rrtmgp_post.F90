!>\file GFS_rrtmgp_post.f90
!! This file contains
module GFS_rrtmgp_post
  use machine,                    only: kind_phys
  use GFS_typedefs,               only: GFS_statein_type,   &
                                        GFS_coupling_type,  &
                                        GFS_control_type,   &
                                        GFS_grid_type,      &
                                        GFS_radtend_type,   &
                                        GFS_diag_type
  use module_radiation_aerosols, only: NSPC1
  use module_radlw_parameters,   only: topflw_type, sfcflw_type, proflw_type
  use module_radsw_parameters,   only: topfsw_type, sfcfsw_type, profsw_type, cmpfsw_type
  ! RRTMGP DDT's
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_fluxes_byband,          only: ty_fluxes_byband
  use mo_heating_rates,          only: compute_heating_rate

  implicit none
contains

!>\defgroup GFS_rrtmgp_post GFS RRTMGP Scheme Post
!! @{
!> \section arg_table_GFS_rrtmgp_post_init Argument Table
!!
  subroutine GFS_rrtmgp_post_init ()
  end subroutine GFS_rrtmgp_post_init

!> \section arg_table_GFS_rrtmgp_post_run Argument Table
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
!! | raddt             | time_step_for_radiation                                                                        | radiation time step                                                          | s        |    0 | real                 | kind_phys | in     | F        |
!! | aerodp            | atmosphere_optical_thickness_due_to_ambient_aerosol_particles                                  | vertical integrated optical depth for various aerosol species                | none     |    2 | real                 | kind_phys | in     | F        |
!! | cldsa             | cloud_area_fraction_for_radiation                                                              | fraction of clouds for low, middle, high, total and BL                       | frac     |    2 | real                 | kind_phys | in     | F        |
!! | mtopa             | model_layer_number_at_cloud_top                                                                | vertical indices for low, middle and high cloud tops                         | index    |    2 | integer              |           | in     | F        |
!! | mbota             | model_layer_number_at_cloud_base                                                               | vertical indices for low, middle and high cloud bases                        | index    |    2 | integer              |           | in     | F        |
!! | cloud_fraction    | total_cloud_fraction                                                                           | layer total cloud fraction                                                   | frac     |    2 | real                 | kind_phys | in     | F        |
!! | cldtaulw          | cloud_optical_depth_layers_at_10mu_band                                                        | approx 10mu band layer cloud optical depth                                   | none     |    2 | real                 | kind_phys | in     | F        |
!! | cldtausw          | cloud_optical_depth_layers_at_0.55mu_band                                                      | approx .55mu band layer cloud optical depth                                  | none     |    2 | real                 | kind_phys | in     | F        |
!! | tsfa              | surface_air_temperature_for_radiation                                                          | lowest model layer air temperature for radiation                             | K        |    1 | real                 | kind_phys | in     | F        |
!! | p_lev             | air_pressure_at_interface_for_RRTMGP_in_hPa                                                    | air pressure level                                                           | hPa      |    2 | real                 | kind_phys | in     | F        |
!! | nday              | daytime_points_dimension                                                                       | daytime points dimension                                                     | count    |    0 | integer              |           | in     | F        |
!! | idxday            | daytime_points                                                                                 | daytime points                                                               | index    |    1 | integer              |           | in     | F        |
!! | fluxswUP_allsky   | sw_flux_profile_upward_allsky                                                                  | RRTMGP upward shortwave all-sky flux profile                                 | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxswDOWN_allsky | sw_flux_profile_downward_allsky                                                                | RRTMGP downward shortwave all-sky flux profile                               | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxswUP_clrsky   | sw_flux_profile_upward_clrsky                                                                  | RRTMGP upward shortwave clr-sky flux profile                                 | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxswDOWN_clrsky | sw_flux_profile_downward_clrsky                                                                | RRTMGP downward shortwave clr-sky flux profile                               | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxlwUP_allsky   | lw_flux_profile_upward_allsky                                                                  | RRTMGP upward longwave all-sky flux profile                                  | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxlwDOWN_allsky | lw_flux_profile_downward_allsky                                                                | RRTMGP downward longwave all-sky flux profile                                | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxlwUP_clrsky   | lw_flux_profile_upward_clrsky                                                                  | RRTMGP upward longwave clr-sky flux profile                                  | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | fluxlwDOWN_clrsky | lw_flux_profile_downward_clrsky                                                                | RRTMGP downward longwave clr-sky flux profile                                | W m-2    |    2 | real                 | kind_phys | in     | F        |
!! | kdist_lw          | K_distribution_file_for_RRTMGP_LW_scheme                                                       | DDT containing spectral information for RRTMGP LW radiation scheme           | DDT      |    0 | ty_gas_optics_rrtmgp |           | in     | F        |
!! | kdist_sw          | K_distribution_file_for_RRTMGP_SW_scheme                                                       | DDT containing spectral information for RRTMGP SW radiation scheme           | DDT      |    0 | ty_gas_optics_rrtmgp |           | in     | F        |
!! | sfc_alb_nir_dir   | surface_shortwave_albedo_near_infrared_direct_in_each_band                                     | surface sw near-infrared direct albedo in each SW band                       | frac     |    2 | real                 | kind_phys | in     | F        |
!! | sfc_alb_nir_dif   | surface_shortwave_albedo_near_infrared_diffuse_in_each_band                                    | surface sw near-infrared diffuse albedo in each SW band                      | frac     |    2 | real                 | kind_phys | in     | F        |
!! | sfc_alb_uvvis_dir | surface_shortwave_albedo_uv_visible_direct_in_each_band                                        | surface sw uv-visible direct albedo in each SW band                          | frac     |    2 | real                 | kind_phys | in     | F        |
!! | sfc_alb_uvvis_dif | surface_shortwave_albedo_uv_visible_diffuse_in_each_band                                       | surface sw uv-visible diffuse albedo in each SW band                         | frac     |    2 | real                 | kind_phys | in     | F        |
!! | hlwc              | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step                     | longwave total sky heating rate                                              | K s-1    |    2 | real                 | kind_phys | out    | F        |
!! | hswc              | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step                    | shortwave total sky heating rate                                             | K s-1    |    2 | real                 | kind_phys | out    | F        |
!! | topflx_lw         | lw_fluxes_top_atmosphere                                                                       | longwave total sky fluxes at the top of the atm                              | W m-2    |    1 | topflw_type          |           | inout  | F        |
!! | sfcflx_lw         | lw_fluxes_sfc                                                                                  | longwave total sky fluxes at the Earth surface                               | W m-2    |    1 | sfcflw_type          |           | inout  | F        |
!! | topflx_sw         | sw_fluxes_top_atmosphere                                                                       | shortwave total sky fluxes at the top of the atm                             | W m-2    |    1 | topfsw_type          |           | inout  | F        |
!! | sfcflx_sw         | sw_fluxes_sfc                                                                                  | shortwave total sky fluxes at the Earth surface                              | W m-2    |    1 | sfcfsw_type          |           | inout  | F        |
!! | flxprf_lw         | lw_fluxes                                                                                      | lw fluxes total sky / csk and up / down at levels                            | W m-2    |    2 | proflw_type          |           | inout  | T        |
!! | flxprf_sw         | sw_fluxes                                                                                      | sw fluxes total sky / csk and up / down at levels                            | W m-2    |    2 | profsw_type          |           | inout  | T        |
!! | hlw0              | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step  | longwave clear sky heating rate                                              | K s-1    |    2 | real                 | kind_phys | inout  | T        |
!! | hsw0              | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | shortwave clear sky heating rate                                             | K s-1    |    2 | real                 | kind_phys | inout  | T        |
!! | errmsg            | ccpp_error_message                                                                             | error message for error handling in CCPP                                     | none     |    0 | character            | len=*     | out    | F        |
!! | errflg            | ccpp_error_flag                                                                                | error flag for error handling in CCPP                                        | flag     |    0 | integer              |           | out    | F        |
!!
  subroutine GFS_rrtmgp_post_run (Model, Grid, Diag, Radtend, Statein, &
              Coupling, scmpsw, im, raddt, aerodp,   &
              cldsa, mtopa, mbota, cloud_fraction, cldtaulw, cldtausw, p_lev, kdist_lw, kdist_sw,         &
              sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir,               &
              sfc_alb_uvvis_dif, &
              tsfa, nday, idxday, fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky, fluxlwDOWN_clrsky, &
              fluxswUP_allsky, fluxswDOWN_allsky, fluxswUP_clrsky, fluxswDOWN_clrsky, &
              hlwc, hswc, topflx_sw, sfcflx_sw, flxprf_sw, topflx_lw, sfcflx_lw, flxprf_lw, hlw0, hsw0, errmsg, errflg)

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
    real(kind_phys), intent(in) :: &
         raddt             ! Radiation time step
    real(kind_phys), dimension(size(Grid%xlon,1)), intent(in) ::  &
         tsfa              ! Lowest model layer air temperature for radiation 
    real(kind_phys), dimension(size(Grid%xlon,1),NSPC1), intent(in) :: &
         aerodp            ! Vertical integrated optical depth for various aerosol species  
    real(kind_phys), dimension(size(Grid%xlon,1),5), intent(in) :: &
         cldsa             ! Fraction of clouds for low, middle, high, total and BL 
    integer,         dimension(size(Grid%xlon,1),3), intent(in) ::&
         mbota,          & ! vertical indices for low, middle and high cloud tops 
         mtopa             ! vertical indices for low, middle and high cloud bases
    real(kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(in) :: &
         cloud_fraction, & ! Total cloud fraction in each layer
         cldtausw,       & ! approx .55mu band layer cloud optical depth  
         cldtaulw          ! approx 10mu band layer cloud optical depth  
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         kdist_lw,       & ! DDT containing LW spectral information
         kdist_sw          ! DDT containing SW spectral information
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         p_lev             ! Pressure @ model layer-interfaces    (hPa)
    real(kind_phys),dimension(kdist_sw%get_nband(),size(Grid%xlon,1)),intent(in) :: &
         sfc_alb_nir_dir,   & ! Shortwave surface albedo (nIR-direct) 
         sfc_alb_nir_dif,   & ! Shortwave surface albedo (nIR-diffuse)
         sfc_alb_uvvis_dir, & ! Shortwave surface albedo (uvvis-direct)
         sfc_alb_uvvis_dif    ! Shortwave surface albedo (uvvis-diffuse)    
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         fluxswUP_allsky,   & ! SW All-sky flux                    (W/m2)
         fluxswDOWN_allsky, & ! SW All-sky flux                    (W/m2)
         fluxswUP_clrsky,   & ! SW Clear-sky flux                  (W/m2)
         fluxswDOWN_clrsky, & ! SW All-sky flux                    (W/m2)
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
         hlwc,          & ! Longwave all-sky heating-rate          (K/sec)
         hswc             ! Shortwave all-sky heating-rate         (K/sec)
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
         hlw0,          & ! Longwave clear-sky heating rate          (K/sec)
         hsw0             ! Shortwave clear-sky heating-rate         (K/sec)
    type(proflw_type), dimension(size(Grid%xlon,1), Model%levs+1), optional, intent(inout) :: &
         flxprf_lw        ! 2D radiative fluxes, components:
                          ! upfxc - total sky upward flux            (W/m2)
                          ! dnfxc - total sky dnward flux            (W/m2)
                          ! upfx0 - clear sky upward flux            (W/m2)
                          ! dnfx0 - clear sky dnward flux            (W/m2)
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
    logical :: l_clrskylw_hr,l_clrskysw_hr, l_fluxeslw2d, l_fluxessw2d, top_at_1, l_sfcFluxessw1D

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. (Model%lsswr .or. Model%lslwr)) return

    ! Are any optional outputs requested?
    l_clrskylw_hr   = present(hlw0)
    l_fluxeslw2d    = present(flxprf_lw)
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

   
    ! #######################################################################################
    ! Compute LW heating-rates. (Note. This piece was originally in rrtmg_lw.F90:_run())
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


    ! #######################################################################################
    ! #######################################################################################
    !>  - For time averaged output quantities (including total-sky and
    !!    clear-sky SW and LW fluxes at TOA and surface; conventional
    !!    3-domain cloud amount, cloud top and base pressure, and cloud top
    !!    temperature; aerosols AOD, etc.), store computed results in
    !!    corresponding slots of array fluxr with appropriate time weights.
    !  --- ...  collect the fluxr data for wrtsfc
    ! #######################################################################################
    ! #######################################################################################

    if (Model%lssav) then
       if (Model%lsswr) then
          do i=1,im
             Diag%fluxr(i,34) = Diag%fluxr(i,34) + Model%fhswr*aerodp(i,1)  ! total aod at 550nm
             Diag%fluxr(i,35) = Diag%fluxr(i,35) + Model%fhswr*aerodp(i,2)  ! DU aod at 550nm
             Diag%fluxr(i,36) = Diag%fluxr(i,36) + Model%fhswr*aerodp(i,3)  ! BC aod at 550nm
             Diag%fluxr(i,37) = Diag%fluxr(i,37) + Model%fhswr*aerodp(i,4)  ! OC aod at 550nm
             Diag%fluxr(i,38) = Diag%fluxr(i,38) + Model%fhswr*aerodp(i,5)  ! SU aod at 550nm
             Diag%fluxr(i,39) = Diag%fluxr(i,39) + Model%fhswr*aerodp(i,6)  ! SS aod at 550nm
          enddo
       endif
       
       !  Save LW TOA and SFC fluxes
       if (Model%lslwr) then
          do i=1,im
             ! LW all-sky fluxes
             Diag%fluxr(i,1 ) = Diag%fluxr(i,1 ) + Model%fhlwr *    Diag%topflw(i)%upfxc   ! total sky top lw up
             Diag%fluxr(i,19) = Diag%fluxr(i,19) + Model%fhlwr * Radtend%sfcflw(i)%dnfxc   ! total sky sfc lw dn
             Diag%fluxr(i,20) = Diag%fluxr(i,20) + Model%fhlwr * Radtend%sfcflw(i)%upfxc   ! total sky sfc lw up
             ! LW clear-sky fluxes
             Diag%fluxr(i,28) = Diag%fluxr(i,28) + Model%fhlwr *    Diag%topflw(i)%upfx0   ! clear sky top lw up
             Diag%fluxr(i,30) = Diag%fluxr(i,30) + Model%fhlwr * Radtend%sfcflw(i)%dnfx0   ! clear sky sfc lw dn
             Diag%fluxr(i,33) = Diag%fluxr(i,33) + Model%fhlwr * Radtend%sfcflw(i)%upfx0   ! clear sky sfc lw up
          enddo
       endif
       
       ! Save sw toa and sfc fluxes with proper diurnal sw wgt. coszen=mean cosz over daylight
       ! part of sw calling interval, while coszdg= mean cosz over entire interval
       if (Model%lsswr) then
          do i = 1, IM
             if (Radtend%coszen(i) > 0.) then
                ! SW all-sky fluxes
                tem0d = Model%fhswr * Radtend%coszdg(i) / Radtend%coszen(i)
                Diag%fluxr(i,2 ) = Diag%fluxr(i,2)  +    Diag%topfsw(i)%upfxc * tem0d  ! total sky top sw up
                Diag%fluxr(i,3 ) = Diag%fluxr(i,3)  + Radtend%sfcfsw(i)%upfxc * tem0d  ! total sky sfc sw up
                Diag%fluxr(i,4 ) = Diag%fluxr(i,4)  + Radtend%sfcfsw(i)%dnfxc * tem0d  ! total sky sfc sw dn
                ! SW uv-b fluxes
                Diag%fluxr(i,21) = Diag%fluxr(i,21) + scmpsw(i)%uvbfc * tem0d          ! total sky uv-b sw dn
                Diag%fluxr(i,22) = Diag%fluxr(i,22) + scmpsw(i)%uvbf0 * tem0d          ! clear sky uv-b sw dn
                ! SW TOA incoming fluxes
                Diag%fluxr(i,23) = Diag%fluxr(i,23) + Diag%topfsw(i)%dnfxc * tem0d     ! top sw dn
                ! SW SFC flux components
                Diag%fluxr(i,24) = Diag%fluxr(i,24) + scmpsw(i)%visbm * tem0d          ! uv/vis beam sw dn
                Diag%fluxr(i,25) = Diag%fluxr(i,25) + scmpsw(i)%visdf * tem0d          ! uv/vis diff sw dn
                Diag%fluxr(i,26) = Diag%fluxr(i,26) + scmpsw(i)%nirbm * tem0d          ! nir beam sw dn
                Diag%fluxr(i,27) = Diag%fluxr(i,27) + scmpsw(i)%nirdf * tem0d          ! nir diff sw dn
                ! SW clear-sky fluxes
                Diag%fluxr(i,29) = Diag%fluxr(i,29) + Diag%topfsw(i)%upfx0 * tem0d     ! clear sky top sw up
                Diag%fluxr(i,31) = Diag%fluxr(i,31) + Radtend%sfcfsw(i)%upfx0 * tem0d  ! clear sky sfc sw up
                Diag%fluxr(i,32) = Diag%fluxr(i,32) + Radtend%sfcfsw(i)%dnfx0 * tem0d  ! clear sky sfc sw dn
             endif
          enddo
       endif
       
       !  ---  save total and boundary layer clouds
       
       if (Model%lsswr .or. Model%lslwr) then
          do i=1,im
             Diag%fluxr(i,17) = Diag%fluxr(i,17) + raddt * cldsa(i,4)
             Diag%fluxr(i,18) = Diag%fluxr(i,18) + raddt * cldsa(i,5)
          enddo
          
          !  ---  save cld frac,toplyr,botlyr and top temp, note that the order
          !       of h,m,l cloud is reversed for the fluxr output.
          !  ---  save interface pressure (pa) of top/bot
          
          do j = 1, 3
             do i = 1, IM
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j)
                ibtc  = mbota(i,j)
                Diag%fluxr(i, 8-j) = Diag%fluxr(i, 8-j) + tem0d
                Diag%fluxr(i,11-j) = Diag%fluxr(i,11-j) + tem0d * Statein%prsi(i,itop)
                Diag%fluxr(i,14-j) = Diag%fluxr(i,14-j) + tem0d * Statein%prsi(i,ibtc)
                Diag%fluxr(i,17-j) = Diag%fluxr(i,17-j) + tem0d * Statein%tgrs(i,itop)
                
                !       Anning adds optical depth and emissivity output
                tem1 = 0.
                tem2 = 0.
                do k=ibtc,itop
                   tem1 = tem1 + cldtausw(i,k)      ! approx .55 mu channel
                   tem2 = tem2 + cldtaulw(i,k)      ! approx 10. mu channel
                enddo
                Diag%fluxr(i,43-j) = Diag%fluxr(i,43-j) + tem0d * tem1
                Diag%fluxr(i,46-j) = Diag%fluxr(i,46-j) + tem0d * (1.0-exp(-tem2))
             enddo
          enddo
       endif
       
       !       if (.not. Model%uni_cld) then
       if (Model%lgocart .or. Model%ldiag3d) then
          do k = 1, Model%levs
             Coupling%cldcovi(1:im,k) = cloud_fraction(1:im,k)
          enddo
       endif
    endif                                ! end_if_lssav
    !
  end subroutine GFS_rrtmgp_post_run

!> \section arg_table_GFS_rrtmgp_post_finalize Argument Table
!!
  subroutine GFS_rrtmgp_post_finalize ()

  end subroutine GFS_rrtmgp_post_finalize
  subroutine check_error_msg(routine_name, error_msg)
    character(len=*), intent(in) :: &
         error_msg, routine_name
    
    if(error_msg /= "") then
       print*,"ERROR("//trim(routine_name)//"): "
       print*,trim(error_msg)
       return
    end if
  end subroutine check_error_msg   
end module GFS_rrtmgp_post
