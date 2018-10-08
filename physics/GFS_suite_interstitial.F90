!> \file GFS_suite_interstitial.f90
!!  Contains code related to more than one scheme in the GFS physics suite.

    module GFS_suite_interstitial_rad_reset

    contains

    subroutine GFS_suite_interstitial_rad_reset_init ()
    end subroutine GFS_suite_interstitial_rad_reset_init

    subroutine GFS_suite_interstitial_rad_reset_finalize()
    end subroutine GFS_suite_interstitial_rad_reset_finalize

!> \section arg_table_GFS_suite_interstitial_rad_reset_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Interstitial   | FV3-GFS_Interstitial_type                              | derived type GFS_interstitial_type in FV3               | DDT           |    0 | GFS_interstitial_type |           | inout  | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
    subroutine GFS_suite_interstitial_rad_reset_run (Interstitial, errmsg, errflg)

      use GFS_typedefs, only: GFS_interstitial_type

      implicit none

      ! interface variables
      type(GFS_interstitial_type), intent(inout) :: Interstitial
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call Interstitial%rad_reset()

    end subroutine GFS_suite_interstitial_rad_reset_run

    end module GFS_suite_interstitial_rad_reset


    module GFS_suite_interstitial_phys_reset

    contains

    subroutine GFS_suite_interstitial_phys_reset_init ()
    end subroutine GFS_suite_interstitial_phys_reset_init

    subroutine GFS_suite_interstitial_phys_reset_finalize()
    end subroutine GFS_suite_interstitial_phys_reset_finalize

!> \section arg_table_GFS_suite_interstitial_phys_reset_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Interstitial   | FV3-GFS_Interstitial_type                              | derived type GFS_interstitial_type in FV3               | DDT           |    0 | GFS_interstitial_type |           | inout  | F        |
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters | DDT           |    0 | GFS_control_type      |           | in     | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
    subroutine GFS_suite_interstitial_phys_reset_run (Interstitial, Model, errmsg, errflg)

      use GFS_typedefs, only: GFS_control_type, GFS_interstitial_type

      implicit none

      ! interface variables
      type(GFS_interstitial_type), intent(inout) :: Interstitial
      type(GFS_control_type),      intent(in)    :: Model
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call Interstitial%phys_reset(Model)

    end subroutine GFS_suite_interstitial_phys_reset_run

    end module GFS_suite_interstitial_phys_reset


    module GFS_suite_interstitial_1

    contains

    subroutine GFS_suite_interstitial_1_init ()
    end subroutine GFS_suite_interstitial_1_init

    subroutine GFS_suite_interstitial_1_finalize()
    end subroutine GFS_suite_interstitial_1_finalize

!> \section arg_table_GFS_suite_interstitial_1_run Argument Table
!! | local_name     | standard_name                                                             | long_name                                                               | units         | rank | type       |    kind   | intent | optional |
!! |----------------|---------------------------------------------------------------------------|-------------------------------------------------------------------------|---------------|------|------------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                    | horizontal loop extent                                                  | count         |    0 | integer    |           | in     | F        |
!! | levs           | vertical_dimension                                                        | vertical layer dimension                                                | count         |    0 | integer    |           | in     | F        |
!! | ntrac          | number_of_tracers                                                         | number of tracers                                                       | count         |    0 | integer    |           | in     | F        |
!! | crtrh          | critical_relative_humidity_at_sfc_pbltop_toa                              | critical relative humidity at SFC, PBL top and TOA                      | frac          |    1 | real       | kind_phys | in     | F        |
!! | dtf            | time_step_for_dynamics                                                    | dynamics timestep                                                       | s             |    0 | real       | kind_phys | in     | F        |
!! | dtp            | time_step_for_physics                                                     | physics timestep                                                        | s             |    0 | real       | kind_phys | in     | F        |
!! | slmsk          | sea_land_ice_mask_real                                                    | landmask: sea/land/ice=0/1/2                                            | flag          |    1 | real       | kind_phys | in     | F        |
!! | area           | cell_area                                                                 | area of the grid cell                                                   | m2            |    1 | real       | kind_phys | in     | F        |
!! | dxmin          | minimum_scaling_factor_for_critical_relative_humidity                     | minimum scaling factor for critical relative humidity                   | m2 rad-2      |    0 | real       | kind_phys | in     | F        |
!! | dxinv          | inverse_scaling_factor_for_critical_relative_humidity                     | inverse scaling factor for critical relative humidity                   | rad2 m-2      |    0 | real       | kind_phys | in     | F        |
!! | pgr            | surface_air_pressure                                                      | surface pressure                                                        | Pa            |    1 | real       | kind_phys | in     | F        |
!! | rhbbot         | critical_relative_humidity_at_surface                                     | critical relative humidity at the surface                               | frac          |    0 | real       | kind_phys | out    | F        |
!! | rhpbl          | critical_relative_humidity_at_PBL_top                                     | critical relative humidity at the PBL top                               | frac          |    0 | real       | kind_phys | out    | F        |
!! | rhbtop         | critical_relative_humidity_at_top_of_atmosphere                           | critical relative humidity at the top of atmosphere                     | frac          |    0 | real       | kind_phys | out    | F        |
!! | frain          | dynamics_to_physics_timestep_ratio                                        | ratio of dynamics timestep to physics timestep                          | none          |    0 | real       | kind_phys | out    | F        |
!! | islmsk         | sea_land_ice_mask                                                         | landmask: sea/land/ice=0/1/2                                            | flag          |    1 | integer    |           | out    | F        |
!! | frland         | land_area_fraction                                                        | land area fraction                                                      | frac          |    1 | real       | kind_phys | out    | F        |
!! | work1          | grid_size_related_coefficient_used_in_scale-sensitive_schemes             | grid size related coefficient used in scale-sensitive schemes           | none          |    1 | real       | kind_phys | out    | F        |
!! | work2          | grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement  | complement to work1                                                     | none          |    1 | real       | kind_phys | out    | F        |
!! | psurf          | surface_air_pressure_diag                                                 | surface air pressure diagnostic                                         | Pa            |    1 | real       | kind_phys | out    | F        |
!! | dudt           | tendency_of_x_wind_due_to_model_physics                                   | updated tendency of the x wind                                          | m s-2         |    2 | real       | kind_phys | out    | F        |
!! | dvdt           | tendency_of_y_wind_due_to_model_physics                                   | updated tendency of the y wind                                          | m s-2         |    2 | real       | kind_phys | out    | F        |
!! | dtdt           | tendency_of_air_temperature_due_to_model_physics                          | updated tendency of the temperature                                     | K s-1         |    2 | real       | kind_phys | out    | F        |
!! | dtdtc          | tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky   | clear sky radiative (shortwave + longwave) heating rate at current time | K s-1         |    2 | real       | kind_phys | out    | F        |
!! | dqdt           | tendency_of_tracers_due_to_model_physics                                  | updated tendency of the tracers                                         | kg kg-1 s-1   |    3 | real       | kind_phys | out    | F        |
!! | errmsg         | ccpp_error_message                                                        | error message for error handling in CCPP                                | none          |    0 | character  | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                           | error flag for error handling in CCPP                                   | flag          |    0 | integer    |           | out    | F        |
!!
    subroutine GFS_suite_interstitial_1_run (im, levs, ntrac, crtrh, dtf, dtp, slmsk, area, dxmin, dxinv, pgr, &
      rhbbot, rhpbl, rhbtop, frain, islmsk, frland, work1, work2, psurf, dudt, dvdt, dtdt, dtdtc, dqdt, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! interface variables
      integer,              intent(in) :: im, levs, ntrac
      real(kind=kind_phys), intent(in) :: dtf, dtp, dxmin, dxinv
      real(kind=kind_phys), intent(in), dimension(3) :: crtrh
      real(kind=kind_phys), intent(in), dimension(im) :: slmsk, area, pgr

      real(kind=kind_phys), intent(out) :: rhbbot, rhpbl, rhbtop, frain
      integer,              intent(out), dimension(im) :: islmsk
      real(kind=kind_phys), intent(out), dimension(im) :: frland, work1, work2, psurf
      real(kind=kind_phys), intent(out), dimension(im,levs) :: dudt, dvdt, dtdt, dtdtc
      real(kind=kind_phys), intent(out), dimension(im,levs,ntrac) ::  dqdt
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! local variables
      integer :: i, k, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      rhbbot = crtrh(1)
      rhpbl  = crtrh(2)
      rhbtop = crtrh(3)

      frain = dtf / dtp

      do i = 1, im
        islmsk(i)   = nint(slmsk(i))
        if (islmsk(i) == 1) then
          frland(i) = 1.0
        else
          frland(i) = 0.0
        endif
        work1(i) = (log(area(i)) - dxmin) * dxinv
        work1(i) = max(0.0, min(1.0,work1(i)))
        work2(i) = 1.0 - work1(i)
        psurf(i) = pgr(i)
      end do

      do k=1,levs
        do i=1,im
          dudt(i,k)  = 0.
          dvdt(i,k)  = 0.
          dtdt(i,k)  = 0.
          dtdtc(i,k) = 0.
        enddo
      enddo
      do n=1,ntrac
        do k=1,levs
          do i=1,im
            dqdt(i,k,n) = 0.
          enddo
        enddo
      enddo

    end subroutine GFS_suite_interstitial_1_run

  end module GFS_suite_interstitial_1


  module GFS_suite_interstitial_2

  contains

    subroutine GFS_suite_interstitial_2_init ()
    end subroutine GFS_suite_interstitial_2_init

    subroutine GFS_suite_interstitial_2_finalize()
    end subroutine GFS_suite_interstitial_2_finalize
#if 0
!> \section arg_table_GFS_suite_interstitial_2_run Argument Table
!! | local_name       | standard_name                                                                                 | long_name                                                                   | units         | rank | type       |    kind   | intent | optional |
!! |------------------|-----------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------|---------------|------|------------|-----------|--------|----------|
!! | im               | horizontal_loop_extent                                                                        | horizontal loop extent                                                      | count         |    0 | integer    |           | in     | F        |
!! | levs             | vertical_dimension                                                                            | vertical layer dimension                                                    | count         |    0 | integer    |           | in     | F        |
!! | lssav            | flag_diagnostics                                                                              | logical flag for storing diagnostics                                        | flag          |    0 | logical    |           | in     | F        |
!! | ldiag3d          | flag_diagnostics_3D                                                                           | flag for 3d diagnostic fields                                               | flag          |    0 | logical    |           | in     | F        |
!! | lsidea           | flag_idealized_physics                                                                        | flag for idealized physics                                                  | flag          |    0 | logical    |           | in     | F        |
!! | cplflx           | flag_for_flux_coupling                                                                        | flag controlling cplflx collection (default off)                            | flag          |    0 | logical    |           | in     | F        |
!! | flag_cice        | flag_for_cice                                                                                 | flag for cice                                                               | flag          |    1 | logical    |           | in     | F        |
!! | shal_cnv         | flag_for_shallow_convection                                                                   | flag for calling shallow convection                                         | flag          |    0 | logical    |           | in     | F        |
!! | old_monin        | flag_for_old_PBL_scheme                                                                       | flag for using old PBL schemes                                              | flag          |    0 | logical    |           | in     | F        |
!! | mstrat           | flag_for_moorthi_stratus                                                                      | flag for moorthi approach for stratus                                       | flag          |    0 | logical    |           | in     | F        |
!! | do_shoc          | flag_for_shoc                                                                                 | flag for SHOC                                                               | flag          |    0 | logical    |           | in     | F        |
!! | imfshalcnv       | flag_for_mass_flux_shallow_convection_scheme                                                  | flag for mass-flux shallow convection scheme                                | flag          |    0 | integer    |           | in     | F        |
!! | dtf              | time_step_for_dynamics                                                                        | dynamics timestep                                                           | s             |    0 | real       | kind_phys | in     | F        |
!! | xcosz            | instantaneous_cosine_of_zenith_angle                                                          | cosine of zenith angle at current time                                      | none          |    1 | real       | kind_phys | in     | F        |
!! | adjsfcdsw        | surface_downwelling_shortwave_flux                                                            | surface downwelling shortwave flux at current time                          | W m-2         |    1 | real       | kind_phys | in     | F        |
!! | adjsfcdlw        | surface_downwelling_longwave_flux                                                             | surface downwelling longwave flux at current time                           | W m-2         |    1 | real       | kind_phys | in     | F        |
!! | pgr              | surface_air_pressure                                                                          | surface pressure                                                            | Pa            |    1 | real       | kind_phys | in     | F        |
!! | ulwsfc_cice      | surface_upwelling_longwave_flux_for_cice                                                      | surface upwelling longwave flux for cice                                    | W m-2         |    1 | real       | kind_phys | in     | F        |
!! | lwhd             | tendency_of_air_temperature_due_to_longwave_heating_for_idea                                  | idea sky lw heating rates                                                   | K s-1         |    3 | real       | kind_phys | in     | F        |
!! | htrsw            | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep                    | total sky sw heating rate                                                   | K s-1         |    2 | real       | kind_phys | in     | F        |
!! | htrlw            | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep                     | total sky lw heating rate                                                   | K s-1         |    2 | real       | kind_phys | in     | F        |
!! | xmu              | zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes                                  | zenith angle temporal adjustment factor for shortwave fluxes                | none          |    1 | real       | kind_phys | in     | F        |
!! | ctei_rm          | critical_cloud_top_entrainment_instability_criteria                                           | critical cloud top entrainment instability criteria                         | none          |    1 | real       | kind_phys | in     | F        |
!! | work1            | grid_size_related_coefficient_used_in_scale-sensitive_schemes                                 | grid size related coefficient used in scale-sensitive schemes               | none          |    1 | real       | kind_phys | in     | F        |
!! | work2            | grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement                      | complement to work1                                                         | none          |    1 | real       | kind_phys | in     | F        |
!! | prsi             | air_pressure_at_interface                                                                     | air pressure at model layer interfaces                                      | Pa            |    2 | real       | kind_phys | in     | F        |
!! | tgrs             | air_temperature                                                                               | model layer mean temperature                                                | K             |    2 | real       | kind_phys | in     | F        |
!! | prsl             | air_pressure                                                                                  | mean layer pressure                                                         | Pa            |    2 | real       | kind_phys | in     | F        |
!! | qgrs_water_vapor | water_vapor_specific_humidity                                                                 | water vapor specific humidity                                               | kg kg-1       |    2 | real       | kind_phys | in     | F        |
!! | qgrs_cloud_water | cloud_condensed_water_mixing_ratio                                                            | moist (dry+vapor, no condensates) mixing ratio of cloud water (condensate)  | kg kg-1       |    2 | real       | kind_phys | in     | F        |
!! | cp               | specific_heat_of_dry_air_at_constant_pressure                                                 | specific heat of dry air at constant pressure                               | J kg-1 K-1    |    0 | real       | kind_phys | in     | F        |
!! | hvap             | latent_heat_of_vaporization_of_water_at_0C                                                    | latent heat of evaporation/sublimation                                      | J kg-1        |    0 | real       | kind_phys | in     | F        |
!! | prslk            | dimensionless_exner_function_at_model_layers                                                  | dimensionless Exner function at model layer centers                         | none          |    2 | real       | kind_phys | in     | F        |
!! | suntim           | duration_of_sunshine                                                                          | sunshine duration time                                                      | s             |    1 | real       | kind_phys | inout  | F        |
!! | adjsfculw        | surface_upwelling_longwave_flux                                                               | surface upwelling longwave flux at current time                             | W m-2         |    1 | real       | kind_phys | inout  | F        |
!! | dlwsfc           | cumulative_surface_downwelling_longwave_flux_multiplied_by_timestep                           | cumulative surface downwelling LW flux multiplied by timestep               | W m-2 s       |    1 | real       | kind_phys | inout  | F        |
!! | ulwsfc           | cumulative_surface_upwelling_longwave_flux_multiplied_by_timestep                             | cumulative surface upwelling LW flux multiplied by timestep                 | W m-2 s       |    1 | real       | kind_phys | inout  | F        |
!! | psmean           | cumulative_surface_pressure_multiplied_by_timestep                                            | cumulative surface pressure multiplied by timestep                          | Pa s          |    1 | real       | kind_phys | inout  | F        |
!! | dt3dt_lw         | cumulative_change_in_temperature_due_to_longwave_radiation                                    | cumulative change in temperature due to longwave radiation                  | K             |    2 | real       | kind_phys | inout  | F        |
!! | dt3dt_sw         | cumulative_change_in_temperature_due_to_shortwave_radiation_and_orographic_gravity_wave_drag  | cumulative change in temperature due to SW rad and oro. GWD                 | K             |    2 | real       | kind_phys | inout  | F        |
!! | dt3dt_pbl        | cumulative_change_in_temperature_due_to_PBL                                                   | cumulative change in temperature due to PBL                                 | K             |    2 | real       | kind_phys | inout  | F        |
!! | dt3dt_dcnv       | cumulative_change_in_temperature_due_to_deep_convection                                       | cumulative change in temperature due to deep conv.                          | K             |    2 | real       | kind_phys | inout  | F        |
!! | dt3dt_scnv       | cumulative_change_in_temperature_due_to_shal_convection                                       | cumulative change in temperature due to shal conv.                          | K             |    2 | real       | kind_phys | inout  | F        |
!! | dt3dt_mp         | cumulative_change_in_temperature_due_to_microphysics                                          | cumulative change in temperature due to microphysics                        | K             |    2 | real       | kind_phys | inout  | F        |
!! | ctei_rml         | grid_sensitive_critical_cloud_top_entrainment_instability_criteria                            | grid sensitive critical cloud top entrainment instability criteria          | none          |    1 | real       | kind_phys | inout  | F        |
!! | ctei_r           | cloud_top_entrainment_instability_value                                                       | cloud top entrainment instability value                                     | none          |    1 | real       | kind_phys | inout  | F        |
!! | kinver           | index_of_highest_temperature_inversion                                                        | index of highest temperature inversion                                      | index         |    1 | integer    |           | inout  | F        |
!! | errmsg           | ccpp_error_message                                                                            | error message for error handling in CCPP                                    | none          |    0 | character  | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                                               | error flag for error handling in CCPP                                       | flag          |    0 | integer    |           | out    | F        |
!!
#endif
    subroutine GFS_suite_interstitial_2_run (im, levs, lssav, ldiag3d, lsidea, cplflx, flag_cice, shal_cnv, old_monin, mstrat,  &
      do_shoc, imfshalcnv, dtf, xcosz, adjsfcdsw, adjsfcdlw, pgr, ulwsfc_cice, lwhd, htrsw, htrlw, xmu, ctei_rm, work1, work2,  &
      prsi, tgrs, prsl, qgrs_water_vapor, qgrs_cloud_water, cp, hvap, prslk,                                                    &
      suntim, adjsfculw, dlwsfc, ulwsfc, psmean, dt3dt_lw, dt3dt_sw, dt3dt_pbl, dt3dt_dcnv, dt3dt_scnv, dt3dt_mp, ctei_rml,     &
      ctei_r, kinver, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! interface variables
      integer,              intent(in) :: im, levs, imfshalcnv
      logical,              intent(in) :: lssav, ldiag3d, lsidea, cplflx, shal_cnv, old_monin, mstrat, do_shoc
      real(kind=kind_phys), intent(in) :: dtf, cp, hvap

      logical,              intent(in), dimension(im) :: flag_cice
      real(kind=kind_phys), intent(in), dimension(2) :: ctei_rm
      real(kind=kind_phys), intent(in), dimension(im) :: xcosz, adjsfcdsw, adjsfcdlw, pgr, xmu, ulwsfc_cice, work1, work2
      real(kind=kind_phys), intent(in), dimension(im, levs) :: htrsw, htrlw, tgrs, prsl, qgrs_water_vapor, qgrs_cloud_water, prslk
      real(kind=kind_phys), intent(in), dimension(im, levs+1) :: prsi
      real(kind=kind_phys), intent(in), dimension(im, levs, 6) :: lwhd

      integer,              intent(inout), dimension(im) :: kinver
      real(kind=kind_phys), intent(inout), dimension(im) :: suntim, dlwsfc, ulwsfc, psmean, adjsfculw, ctei_rml, ctei_r
      real(kind=kind_phys), intent(inout), dimension(im, levs) :: dt3dt_lw, dt3dt_sw, dt3dt_pbl, dt3dt_dcnv, dt3dt_scnv, dt3dt_mp

      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! local variables
      real(kind=kind_phys), parameter :: czmin   = 0.0001      ! cos(89.994)
      integer :: i, k
      real(kind=kind_phys) :: tem1, tem2, tem, hocp
      logical, dimension(im) :: invrsn
      real(kind=kind_phys), dimension(im) :: tx1, tx2

      real(kind=kind_phys), parameter :: qmin    = 1.0e-10

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      hocp = hvap/cp

      if (lssav) then      !  --- ...  accumulate/save output variables

!  --- ...  sunshine duration time is defined as the length of time (in mdl output
!           interval) that solar radiation falling on a plane perpendicular to the
!           direction of the sun >= 120 w/m2

        do i = 1, im
          if ( xcosz(i) >= czmin ) then   ! zenth angle > 89.994 deg
            tem1 = adjsfcdsw(i) / xcosz(i)
            if ( tem1 >= 120.0 ) then
              suntim(i) = suntim(i) + dtf
            endif
          endif
        enddo

!  --- ...  sfc lw fluxes used by atmospheric model are saved for output
        if (cplflx) then
          do i=1,im
            if (flag_cice(i)) adjsfculw(i) = ulwsfc_cice(i)
          enddo
        endif
        do i=1,im
          dlwsfc(i) = dlwsfc(i) +   adjsfcdlw(i)*dtf
          ulwsfc(i) = ulwsfc(i) +   adjsfculw(i)*dtf
          psmean(i) = psmean(i) +   pgr(i)*dtf        ! mean surface pressure
        end do

        if (ldiag3d) then
          if (lsidea) then
            do k=1,levs
              do i=1,im
                dt3dt_lw(i,k) = dt3dt_lw(i,k) + lwhd(i,k,1)*dtf
                dt3dt_sw(i,k) = dt3dt_sw(i,k) + lwhd(i,k,2)*dtf
                dt3dt_pbl(i,k) = dt3dt_pbl(i,k) + lwhd(i,k,3)*dtf
                dt3dt_dcnv(i,k) = dt3dt_dcnv(i,k) + lwhd(i,k,4)*dtf
                dt3dt_scnv(i,k) = dt3dt_scnv(i,k) + lwhd(i,k,5)*dtf
                dt3dt_mp(i,k) = dt3dt_mp(i,k) + lwhd(i,k,6)*dtf
              end do
            end do
          else
            do k=1,levs
              do i=1,im
                dt3dt_lw(i,k) = dt3dt_lw(i,k) + htrlw(i,k)*dtf
                dt3dt_sw(i,k) = dt3dt_sw(i,k) + htrsw(i,k)*dtf*xmu(i)
              enddo
            enddo
          endif
        endif
      endif    ! end if_lssav_block

      do i=1, im
        invrsn(i) = .false.
        tx1(i) = 0.0
        tx2(i) = 10.0
        ctei_r(i) = 10.0
      end do

      if ((((imfshalcnv == 0 .and. shal_cnv) .or. old_monin) .and. mstrat) &
         .or. do_shoc) then
        ctei_rml(:) = ctei_rm(1)*work1(:) + ctei_rm(2)*work2(:)
        do k=1,levs/2
          do i=1,im
            if (prsi(i,1)-prsi(i,k+1) < 0.35*prsi(i,1)       &
                .and. (.not. invrsn(i))) then
              tem = (tgrs(i,k+1) - tgrs(i,k))  &
                  / (prsl(i,k)   - prsl(i,k+1))

              if (((tem > 0.00010) .and. (tx1(i) < 0.0)) .or.  &
                  ((tem-abs(tx1(i)) > 0.0) .and. (tx2(i) < 0.0))) then
                invrsn(i) = .true.

                if (qgrs_water_vapor(i,k) > qgrs_water_vapor(i,k+1)) then
                  tem1 = tgrs(i,k+1) + hocp*max(qgrs_water_vapor(i,k+1),qmin)
                  tem2 = tgrs(i,k)   + hocp*max(qgrs_water_vapor(i,k),qmin)

                  tem1 = tem1 / prslk(i,k+1) - tem2 / prslk(i,k)

!  --- ...  (cp/l)(deltathetae)/(deltatwater) > ctei_rm -> conditon for CTEI
                  ctei_r(i) = (1.0/hocp)*tem1/(qgrs_water_vapor(i,k+1)-qgrs_water_vapor(i,k)  &
                            + qgrs_cloud_water(i,k+1)-qgrs_cloud_water(i,k))
                else
                  ctei_r(i) = 10
                endif

                if ( ctei_rml(i) > ctei_r(i) ) then
                  kinver(i) = k
                else
                  kinver(i) = levs
                endif
              endif

              tx2(i) = tx1(i)
              tx1(i) = tem
            endif
          enddo
        enddo
      endif

    end subroutine GFS_suite_interstitial_2_run

  end module GFS_suite_interstitial_2


  module GFS_suite_stateout_reset

  contains

    subroutine GFS_suite_stateout_reset_init ()
    end subroutine GFS_suite_stateout_reset_init

    subroutine GFS_suite_stateout_reset_finalize()
    end subroutine GFS_suite_stateout_reset_finalize

!> \section arg_table_GFS_suite_stateout_reset_run Argument Table
!! | local_name     | standard_name                                                | long_name                                                             | units         | rank | type       |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|------------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                       | horizontal loop extent                                                | count         |    0 | integer    |           | in     | F        |
!! | levs           | vertical_dimension                                           | vertical layer dimension                                              | count         |    0 | integer    |           | in     | F        |
!! | ntrac          | number_of_tracers                                            | number of tracers                                                     | count         |    0 | integer    |           | in     | F        |
!! | tgrs           | air_temperature                                              | model layer mean temperature                                          | K             |    2 | real       | kind_phys | in     | F        |
!! | ugrs           | x_wind                                                       | zonal wind                                                            | m s-1         |    2 | real       | kind_phys | in     | F        |
!! | vgrs           | y_wind                                                       | meridional wind                                                       | m s-1         |    2 | real       | kind_phys | in     | F        |
!! | qgrs           | tracer_concentration                                         | model layer mean tracer concentration                                 | kg kg-1       |    3 | real       | kind_phys | in     | F        |
!! | gt0            | air_temperature_updated_by_physics                           | temperature updated by physics                                        | K             |    2 | real       | kind_phys | out    | F        |
!! | gu0            | x_wind_updated_by_physics                                    | zonal wind updated by physics                                         | m s-1         |    2 | real       | kind_phys | out    | F        |
!! | gv0            | y_wind_updated_by_physics                                    | meridional wind updated by physics                                    | m s-1         |    2 | real       | kind_phys | out    | F        |
!! | gq0            | tracer_concentration_updated_by_physics                      | tracer concentration updated by physics                               | kg kg-1       |    3 | real       | kind_phys | out    | F        |
!! | errmsg         | ccpp_error_message                                           | error message for error handling in CCPP                              | none          |    0 | character  | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                              | error flag for error handling in CCPP                                 | flag          |    0 | integer    |           | out    | F        |
!!
    subroutine GFS_suite_stateout_reset_run (im, levs, ntrac,        &
                                             tgrs, ugrs, vgrs, qgrs, &
                                             gt0 , gu0 , gv0 , gq0 , &
                                             errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! interface variables
      integer, intent(in) :: im
      integer, intent(in) :: levs
      integer, intent(in) :: ntrac
      real(kind=kind_phys), dimension(im,levs),       intent(in)  :: tgrs, ugrs, vgrs
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(in)  :: qgrs
      real(kind=kind_phys), dimension(im,levs),       intent(out) :: gt0, gu0, gv0
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(out) :: gq0

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      gt0(:,:)   = tgrs(:,:)
      gu0(:,:)   = ugrs(:,:)
      gv0(:,:)   = vgrs(:,:)
      gq0(:,:,:) = qgrs(:,:,:)

    end subroutine GFS_suite_stateout_reset_run

  end module GFS_suite_stateout_reset


  module GFS_suite_stateout_update

  contains

    subroutine GFS_suite_stateout_update_init ()
    end subroutine GFS_suite_stateout_update_init

    subroutine GFS_suite_stateout_update_finalize()
    end subroutine GFS_suite_stateout_update_finalize

!> \section arg_table_GFS_suite_stateout_update_run Argument Table
!! | local_name     | standard_name                                                | long_name                                                             | units         | rank | type       |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|------------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                       | horizontal loop extent                                                | count         |    0 | integer    |           | in     | F        |
!! | levs           | vertical_dimension                                           | vertical layer dimension                                              | count         |    0 | integer    |           | in     | F        |
!! | ntrac          | number_of_tracers                                            | number of tracers                                                     | count         |    0 | integer    |           | in     | F        |
!! | dtp            | time_step_for_physics                                        | physics timestep                                                      | s             |    0 | real       | kind_phys | in     | F        |
!! | tgrs           | air_temperature                                              | model layer mean temperature                                          | K             |    2 | real       | kind_phys | in     | F        |
!! | ugrs           | x_wind                                                       | zonal wind                                                            | m s-1         |    2 | real       | kind_phys | in     | F        |
!! | vgrs           | y_wind                                                       | meridional wind                                                       | m s-1         |    2 | real       | kind_phys | in     | F        |
!! | qgrs           | tracer_concentration                                         | model layer mean tracer concentration                                 | kg kg-1       |    3 | real       | kind_phys | in     | F        |
!! | dudt           | tendency_of_x_wind_due_to_model_physics                      | updated tendency of the x wind                                        | m s-2         |    2 | real       | kind_phys | in     | F        |
!! | dvdt           | tendency_of_y_wind_due_to_model_physics                      | updated tendency of the y wind                                        | m s-2         |    2 | real       | kind_phys | in     | F        |
!! | dtdt           | tendency_of_air_temperature_due_to_model_physics             | updated tendency of the temperature                                   | K s-1         |    2 | real       | kind_phys | in     | F        |
!! | dqdt           | tendency_of_tracers_due_to_model_physics                     | updated tendency of the tracers                                       | kg kg-1 s-1   |    3 | real       | kind_phys | in     | F        |
!! | gt0            | air_temperature_updated_by_physics                           | temperature updated by physics                                        | K             |    2 | real       | kind_phys | out    | F        |
!! | gu0            | x_wind_updated_by_physics                                    | zonal wind updated by physics                                         | m s-1         |    2 | real       | kind_phys | out    | F        |
!! | gv0            | y_wind_updated_by_physics                                    | meridional wind updated by physics                                    | m s-1         |    2 | real       | kind_phys | out    | F        |
!! | gq0            | tracer_concentration_updated_by_physics                      | tracer concentration updated by physics                               | kg kg-1       |    3 | real       | kind_phys | out    | F        |
!! | errmsg         | ccpp_error_message                                           | error message for error handling in CCPP                              | none          |    0 | character  | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                              | error flag for error handling in CCPP                                 | flag          |    0 | integer    |           | out    | F        |
!!
    subroutine GFS_suite_stateout_update_run (im, levs, ntrac, dtp,  &
                     tgrs, ugrs, vgrs, qgrs, dudt, dvdt, dtdt, dqdt, &
                     gt0, gu0, gv0, gq0, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! Interface variables
      integer,              intent(in) :: im
      integer,              intent(in) :: levs
      integer,              intent(in) :: ntrac
      real(kind=kind_phys), intent(in) :: dtp

      real(kind=kind_phys), dimension(im,levs),       intent(in)  :: tgrs, ugrs, vgrs
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(in)  :: qgrs
      real(kind=kind_phys), dimension(im,levs),       intent(in)  :: dudt, dvdt, dtdt
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(in)  :: dqdt
      real(kind=kind_phys), dimension(im,levs),       intent(out) :: gt0, gu0, gv0
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(out) :: gq0

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      gt0(:,:)   = tgrs(:,:)   + dtdt(:,:)   * dtp
      gu0(:,:)   = ugrs(:,:)   + dudt(:,:)   * dtp
      gv0(:,:)   = vgrs(:,:)   + dvdt(:,:)   * dtp
      gq0(:,:,:) = qgrs(:,:,:) + dqdt(:,:,:) * dtp

    end subroutine GFS_suite_stateout_update_run

  end module GFS_suite_stateout_update


  module GFS_suite_interstitial_3

  contains

    subroutine GFS_suite_interstitial_3_init ()
    end subroutine GFS_suite_interstitial_3_init

    subroutine GFS_suite_interstitial_3_finalize()
    end subroutine GFS_suite_interstitial_3_finalize

#if 0
!> \section arg_table_GFS_suite_interstitial_3_run Argument Table
!! | local_name                 | standard_name                                                                                 | long_name                                                         | units         | rank | type       |    kind   | intent | optional |
!! |----------------------------|-----------------------------------------------------------------------------------------------|-------------------------------------------------------------------|---------------|------|------------|-----------|--------|----------|
!! | im                         | horizontal_loop_extent                                                                        | horizontal loop extent                                            | count         |    0 | integer    |           | in     | F        |
!! | levs                       | vertical_dimension                                                                            | vertical layer dimension                                          | count         |    0 | integer    |           | in     | F        |
!! | nn                         | number_of_tracers_for_convective_transport                                                    | number of tracers for convective transport                        | count         |    0 | integer    |           | in     | F        |
!! | cscnv                      | flag_for_Chikira_Sugiyama_deep_convection                                                     | flag for Chikira-Sugiyama convection                              | flag          |    0 | logical    |           | in     | F        |
!! | satmedmf                   | flag_for_scale_aware_TKE_moist_EDMF_PBL                                                       | flag for scale-aware TKE moist EDMF PBL scheme                    | flag          |    0 | logical    |           | in     | F        |
!! | trans_trac                 | flag_for_convective_transport_of_tracers                                                      | flag for convective transport of tracers                          | flag          |    0 | logical    |           | in     | F        |
!! | do_shoc                    | flag_for_shoc                                                                                 | flag for SHOC                                                     | flag          |    0 | logical    |           | in     | F        |
!! | ltaerosol                  | flag_for_aerosol_physics                                                                      | flag for aerosol physics                                          | flag          |    0 | logical    |           | in     | F        |
!! | ntrac                      | number_of_tracers                                                                             | number of tracers                                                 | count         |    0 | integer    |           | in     | F        |
!! | ntcw                       | index_for_liquid_cloud_condensate                                                             | tracer index for cloud condensate (or liquid water)               | index         |    0 | integer    |           | in     | F        |
!! | ntiw                       | index_for_ice_cloud_condensate                                                                | tracer index for  ice water                                       | index         |    0 | integer    |           | in     | F        |
!! | ntclamt                    | index_for_cloud_amount                                                                        | tracer index for cloud amount integer                             | index         |    0 | integer    |           | in     | F        |
!! | ntrw                       | index_for_rain_water                                                                          | tracer index for rain water                                       | index         |    0 | integer    |           | in     | F        |
!! | ntsw                       | index_for_snow_water                                                                          | tracer index for snow water                                       | index         |    0 | integer    |           | in     | F        |
!! | ntrnc                      | index_for_rain_number_concentration                                                           | tracer index for rain   number concentration                      | index         |    0 | integer    |           | in     | F        |
!! | ntsnc                      | index_for_snow_number_concentration                                                           | tracer index for snow   number concentration                      | index         |    0 | integer    |           | in     | F        |
!! | ntgl                       | index_for_graupel                                                                             | tracer index for graupel                                          | index         |    0 | integer    |           | in     | F        |
!! | ntgnc                      | index_for_graupel_number_concentration                                                        | tracer index for graupel number concentration                     | index         |    0 | integer    |           | in     | F        |
!! | gq0                        | tracer_concentration_updated_by_physics                                                       | tracer concentration updated by physics                           | kg kg-1       |    3 | real       | kind_phys | in     | F        |
!! | imp_physics                | flag_for_microphysics_scheme                                                                  | choice of microphysics scheme                                     | flag          |    0 | integer    |           | in     | F        |
!! | imp_physics_mg             | flag_for_morrison_gettelman_microphysics_scheme                                               | choice of Morrison-Gettelman rmicrophysics scheme                 | flag          |    0 | integer    |           | in     | F        |
!! | imp_physics_zhao_carr      | flag_for_zhao_carr_microphysics_scheme                                                        | choice of Zhao-Carr microphysics scheme                           | flag          |    0 | integer    |           | in     | F        |
!! | imp_physics_zhao_carr_pdf  | flag_for_zhao_carr_pdf_microphysics_scheme                                                    | choice of Zhao-Carr microphysics scheme with PDF clouds           | flag          |    0 | integer    |           | in     | F        |
!! | imp_physics_gfdl           | flag_for_gfdl_microphysics_scheme                                                             | choice of GFDL microphysics scheme                                | flag          |    0 | integer    |           | in     | F        |
!! | imp_physics_thompson       | flag_for_thompson_microphysics_scheme                                                         | choice of Thompson microphysics scheme                            | flag          |    0 | integer    |           | in     | F        |
!! | imp_physics_wsm6           | flag_for_wsm6_microphysics_scheme                                                             | choice of WSM6 microphysics scheme                                | flag          |    0 | integer    |           | in     | F        |
!! | prsi                       | air_pressure_at_interface                                                                     | air pressure at model layer interfaces                            | Pa            |    2 | real       | kind_phys | in     | F        |
!! | prsl                       | air_pressure                                                                                  | mean layer pressure                                               | Pa            |    2 | real       | kind_phys | in     | F        |
!! | prslk                      | dimensionless_exner_function_at_model_layers                                                  | dimensionless Exner function at model layer centers               | none          |    2 | real       | kind_phys | in     | F        |
!! | rhcbot                     | critical_relative_humidity_at_surface                                                         | critical relative humidity at the surface                         | frac          |    0 | real       | kind_phys | in     | F        |
!! | rhcpbl                     | critical_relative_humidity_at_PBL_top                                                         | critical relative humidity at the PBL top                         | frac          |    0 | real       | kind_phys | in     | F        |
!! | rhctop                     | critical_relative_humidity_at_top_of_atmosphere                                               | critical relative humidity at the top of atmosphere               | frac          |    0 | real       | kind_phys | in     | F        |
!! | rhcmax                     | maximum_critical_relative_humidity                                                            | maximum critical relative humidity                                | frac          |    0 | real       | kind_phys | in     | F        |
!! | islmsk                     | sea_land_ice_mask                                                                             | sea/land/ice mask (=0/1/2)                                        | flag          |    1 | integer    |           | in     | F        |
!! | work1                      | grid_size_related_coefficient_used_in_scale-sensitive_schemes                                 | grid size related coefficient used in scale-sensitive schemes     | none          |    1 | real       | kind_phys | in     | F        |
!! | work2                      | grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement                      | complement to work1                                               | none          |    1 | real       | kind_phys | in     | F        |
!! | kpbl                       | vertical_index_at_top_of_atmosphere_boundary_layer                                            | vertical index at top atmospheric boundary layer                  | index         |    1 | integer    |           | in     | F        |
!! | clw                        | convective_transportable_tracers                                                              | array to contain cloud water and other convective trans. tracers  | kg kg-1       |    3 | real       | kind_phys | inout  | F        |
!! | rhc                        | critical_relative_humidity                                                                    | critical relative humidity                                        | frac          |    2 | real       | kind_phys | inout  | F        |
!! | save_qc                    | cloud_liquid_water_mixing_ratio_save                                                          | cloud liquid water mixing ratio before entering a physics scheme  | kg kg-1       |    2 | real       | kind_phys | inout  | F        |
!! | save_qi                    | cloud_ice_water_mixing_ratio_save                                                             | cloud ice water mixing ratio before entering a physics scheme     | kg kg-1       |    2 | real       | kind_phys | inout  | F        |
!! | errmsg                     | ccpp_error_message                                                                            | error message for error handling in CCPP                          | none          |    0 | character  | len=*     | out    | F        |
!! | errflg                     | ccpp_error_flag                                                                               | error flag for error handling in CCPP                             | flag          |    0 | integer    |           | out    | F        |
!!
#endif
    subroutine GFS_suite_interstitial_3_run (im, levs, nn, cscnv, satmedmf, trans_trac, do_shoc, ltaerosol, ntrac, ntcw,  &
      ntiw, ntclamt, ntrw, ntsw, ntrnc, ntsnc, ntgl, ntgnc, gq0, imp_physics, imp_physics_mg, imp_physics_zhao_carr,      &
      imp_physics_zhao_carr_pdf, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6, prsi, prsl, prslk, rhcbot,     &
      rhcpbl, rhctop, rhcmax, islmsk, work1, work2, kpbl,                                                                 &
      clw, rhc, save_qc, save_qi, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! interface variables
      integer,                                          intent(in) :: im, levs, nn, ntrac, ntcw, ntiw, ntclamt, ntrw,     &
        ntsw, ntrnc, ntsnc, ntgl, ntgnc, imp_physics, imp_physics_mg, imp_physics_zhao_carr, imp_physics_zhao_carr_pdf,   &
        imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      integer, dimension(im),                           intent(in) :: islmsk, kpbl
      logical,                                          intent(in) :: cscnv, satmedmf, trans_trac, do_shoc, ltaerosol

      real(kind=kind_phys),                             intent(in) :: rhcbot, rhcmax, rhcpbl, rhctop
      real(kind=kind_phys), dimension(im),              intent(in) :: work1, work2
      real(kind=kind_phys), dimension(im, levs),        intent(in) :: prsl, prslk
      real(kind=kind_phys), dimension(im, levs+1),      intent(in) :: prsi
      real(kind=kind_phys), dimension(im, levs, ntrac), intent(in) :: gq0

      real(kind=kind_phys), dimension(im, levs),      intent(inout) :: rhc, save_qc, save_qi
      real(kind=kind_phys), dimension(im, levs, nn),  intent(inout) :: clw

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! local variables
      integer :: i,k,n,tracers,kk
      real(kind=kind_phys) :: tem, tem1, tem2
      real(kind=kind_phys), dimension(im) :: tx1, tx2

      real(kind=kind_phys),parameter :: slope_mg = 0.02, slope_upmg = 0.02,  &
                         turnrhcrit = 0.900, turnrhcrit_upper = 0.150

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !GF* The following section (initializing convective variables) is already executed in GFS_typedefs%interstitial_phys_reset
      ! do k=1,levs
      !   do i=1,im
      !     clw(i,k,1) = 0.0
      !     clw(i,k,2) = -999.9
      !   enddo
      ! enddo
      ! if ((imfdeepcnv >= 0) .or. imfshalcnv > 0)) then
      !   do k=1,levs
      !     do i=1,im
      !       cnvc(i,k)  = 0.0
      !       cnvw(i,k)  = 0.0
      !     enddo
      !   enddo
      ! endif
      ! if(imp_physics == 8) then
      !   if(Model%ltaerosol) then
      !     ice00 (:,:) = 0.0
      !     liq0  (:,:) = 0.0
      !   else
      !     ice00 (:,:) = 0.0
      !   endif
      ! endif
      !*GF

      if (cscnv .or. satmedmf .or. trans_trac ) then
        tracers = 2
        do n=2,ntrac
          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc) then
            tracers = tracers + 1
            do k=1,levs
              do i=1,im
                clw(i,k,tracers) = gq0(i,k,n)
              enddo
            enddo
          endif
        enddo
      endif ! end if_ras or cfscnv or samf

      if (ntcw > 0) then
        if (imp_physics == imp_physics_mg .and. .not. do_shoc) then ! compute rhc for GMAO macro physics cloud pdf
          do i=1,im
            tx1(i) = 1.0 / prsi(i,1)
            tx2(i) = 1.0 - rhcbot
          enddo
          do k = 1, levs
            do i = 1, im
              tem  = prsl(i,k) * tx1(i)
              tem1 = min(max((tem-turnrhcrit)/slope_mg, -20.0), 20.0)
!             tem2 = min(max((0.3-0.2*abs(cos(Grid%xlat(i)))-tem)/slope_upmg, -20.0), 20.0) ! Anning
              tem2 = min(max((turnrhcrit_upper-tem)/slope_upmg, -20.0), 20.0)
              if (islmsk(i) > 0) then
                tem1 = 1.0 / (1.0+exp(tem1+tem1))
              else
                tem1 = 2.0 / (1.0+exp(tem1+tem1))
              endif
              tem2 = 1.0 / (1.0+exp(tem2))

!             rhc(i,k) = min(rhcmax, max(0.7, 1.0-tx2(i)*tem1*tem2))
              rhc(i,k) = min(rhcmax, rhcmax*work1(i) + (1.0-tx2(i)*tem1*tem2)*work2(i))
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              kk = max(10,kpbl(i))
              if (k < kk) then
                tem    = rhcbot - (rhcbot-rhcpbl) * (1.0-prslk(i,k))                &
                                                 / (1.0-prslk(i,kk))
              else
                tem    = rhcpbl - (rhcpbl-rhctop) * (prslk(i,kk)-prslk(i,k)) &
                                                / prslk(i,kk)
              endif
              tem      = rhcmax * work1(i) + tem * work2(i)
              rhc(i,k) = max(0.0, min(1.0,tem))
            enddo
          enddo
        endif
      endif

      if (imp_physics == imp_physics_zhao_carr .or. imp_physics == imp_physics_zhao_carr_pdf) then   ! zhao-carr microphysics
        !GF* move to GFS_MP_generic_pre (from gscond/precpd)
        ! do i=1,im
        !   psautco_l(i) = Model%psautco(1)*work1(i) + Model%psautco(2)*work2(i)
        !   prautco_l(i) = Model%prautco(1)*work1(i) + Model%prautco(2)*work2(i)
        ! enddo
        !*GF
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntcw)
          enddo
        enddo
      elseif (imp_physics == imp_physics_gfdl) then
        clw(1:im,:,1) = gq0(1:im,:,ntcw)
      elseif (imp_physics == imp_physics_thompson) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
            clw(i,k,2) = gq0(i,k,ntcw)                    ! water
          enddo
        enddo
        if(ltaerosol) then
          save_qi(:,:) = clw(:,:,1)
          save_qc(:,:)  = clw(:,:,2)
        else
          save_qi(:,:) = clw(:,:,1)
        endif
      elseif (imp_physics == imp_physics_wsm6 .or. imp_physics == imp_physics_mg) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
            clw(i,k,2) = gq0(i,k,ntcw)                    ! water
          enddo
        enddo
      else       ! if_ntcw
        !GF* never executed unless imp_physics = imp_physics_zhao_carr or imp_physics_zhao_carr_pdf
        ! do i=1,im
        !   psautco_l(i) = Model%psautco(1)*work1(i) + Model%psautco(2)*work2(i)
        !   prautco_l(i) = Model%prautco(1)*work1(i) + Model%prautco(2)*work2(i)
        ! enddo
        !*GF
        rhc(:,:) = 1.0
      endif   ! end if_ntcw

    end subroutine GFS_suite_interstitial_3_run

  end module GFS_suite_interstitial_3

  module GFS_suite_interstitial_4

  contains

    subroutine GFS_suite_interstitial_4_init ()
    end subroutine GFS_suite_interstitial_4_init

    subroutine GFS_suite_interstitial_4_finalize()
    end subroutine GFS_suite_interstitial_4_finalize

!> \section arg_table_GFS_suite_interstitial_4_run Argument Table
!! | local_name                 | standard_name                                                                 | long_name                                                         | units         | rank | type             |    kind   | intent | optional |
!! |----------------------------|-------------------------------------------------------------------------------|-------------------------------------------------------------------|---------------|------|------------------|-----------|--------|----------|
!! | im                         | horizontal_loop_extent                                                        | horizontal loop extent                                            | count         |    0 | integer          |           | in     | F        |
!! | levs                       | vertical_dimension                                                            | vertical layer dimension                                          | count         |    0 | integer          |           | in     | F        |
!! | ltaerosol                  | flag_for_aerosol_physics                                                      | flag for aerosol physics                                          | flag          |    0 | logical          |           | in     | F        |
!! | lgocart                    | flag_gocart                                                                   | flag for 3d diagnostic fields for gocart 1                        | flag          |    0 | logical          |           | in     | F        |
!! | tracers_total              | number_of_total_tracers                                                       | total number of tracers                                           | count         |    0 | integer          |           | in     | F        |
!! | ntrac                      | number_of_tracers                                                             | number of tracers                                                 | count         |    0 | integer          |           | in     | F        |
!! | ntcw                       | index_for_liquid_cloud_condensate                                             | tracer index for cloud condensate (or liquid water)               | index         |    0 | integer          |           | in     | F        |
!! | ntiw                       | index_for_ice_cloud_condensate                                                | tracer index for  ice water                                       | index         |    0 | integer          |           | in     | F        |
!! | ntclamt                    | index_for_cloud_amount                                                        | tracer index for cloud amount integer                             | index         |    0 | integer          |           | in     | F        |
!! | ntrw                       | index_for_rain_water                                                          | tracer index for rain water                                       | index         |    0 | integer          |           | in     | F        |
!! | ntsw                       | index_for_snow_water                                                          | tracer index for snow water                                       | index         |    0 | integer          |           | in     | F        |
!! | ntrnc                      | index_for_rain_number_concentration                                           | tracer index for rain   number concentration                      | index         |    0 | integer          |           | in     | F        |
!! | ntsnc                      | index_for_snow_number_concentration                                           | tracer index for snow   number concentration                      | index         |    0 | integer          |           | in     | F        |
!! | ntgl                       | index_for_graupel                                                             | tracer index for graupel                                          | index         |    0 | integer          |           | in     | F        |
!! | ntgnc                      | index_for_graupel_number_concentration                                        | tracer index for graupel number concentration                     | index         |    0 | integer          |           | in     | F        |
!! | ntlnc                      | index_for_liquid_cloud_number_concentration                                   | tracer index for liquid number concentration                      | index         |    0 | integer          |           | in     | F        |
!! | ntinc                      | index_for_ice_cloud_number_concentration                                      | tracer index for ice    number concentration                      | index         |    0 | integer          |           | in     | F        |
!! | nn                         | number_of_tracers_for_convective_transport                                    | number of tracers for convective transport                        | count         |    0 | integer          |           | in     | F        |
!! | imp_physics                | flag_for_microphysics_scheme                                                  | choice of microphysics scheme                                     | flag          |    0 | integer          |           | in     | F        |
!! | imp_physics_gfdl           | flag_for_gfdl_microphysics_scheme                                             | choice of GFDL microphysics scheme                                | flag          |    0 | integer          |           | in     | F        |
!! | imp_physics_thompson       | flag_for_thompson_microphysics_scheme                                         | choice of Thompson microphysics scheme                            | flag          |    0 | integer          |           | in     | F        |
!! | imp_physics_zhao_carr      | flag_for_zhao_carr_microphysics_scheme                                        | choice of Zhao-Carr microphysics scheme                           | flag          |    0 | integer          |           | in     | F        |
!! | imp_physics_zhao_carr_pdf  | flag_for_zhao_carr_pdf_microphysics_scheme                                    | choice of Zhao-Carr microphysics scheme with PDF clouds           | flag          |    0 | integer          |           | in     | F        |
!! | dtf                        | time_step_for_dynamics                                                        | dynamics timestep                                                 | s             |    0 | real             | kind_phys | in     | F        |
!! | save_qc                    | cloud_liquid_water_mixing_ratio_save                                          | cloud liquid water mixing ratio before entering a physics scheme  | kg kg-1       |    2 | real             | kind_phys | in     | F        |
!! | save_qi                    | cloud_ice_water_mixing_ratio_save                                             | cloud ice water mixing ratio before entering a physics scheme     | kg kg-1       |    2 | real             | kind_phys | in     | F        |
!! | con_pi                     | pi                                                                            | ratio of a circle's circumference to its diameter                 | radians       |    0 | real             | kind_phys | in     | F        |
!! | gq0                        | tracer_concentration_updated_by_physics                                       | tracer concentration updated by physics                           | kg kg-1       |    3 | real             | kind_phys | inout  | F        |
!! | clw                        | convective_transportable_tracers                                              | array to contain cloud water and other convective trans. tracers  | kg kg-1       |    3 | real             | kind_phys | inout  | F        |
!! | dqdti                      | instantaneous_water_vapor_specific_humidity_tendency_due_to_convection        | instantaneous moisture tendency due to convection                 | kg kg-1 s-1   |    2 | real             | kind_phys | inout  | F        |
!! | errmsg                     | ccpp_error_message                                                            | error message for error handling in CCPP                          | none          |    0 | character        | len=*     | out    | F        |
!! | errflg                     | ccpp_error_flag                                                               | error flag for error handling in CCPP                             | flag          |    0 | integer          |           | out    | F        |
!!
    subroutine GFS_suite_interstitial_4_run (im, levs, ltaerosol, lgocart, tracers_total, ntrac, ntcw, ntiw, ntclamt,      &
      ntrw, ntsw, ntrnc, ntsnc, ntgl, ntgnc, ntlnc, ntinc, nn, imp_physics, imp_physics_gfdl, imp_physics_thompson,       &
      imp_physics_zhao_carr, imp_physics_zhao_carr_pdf, dtf, save_qc, save_qi, con_pi,                                    &
      gq0, clw, dqdti, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! interface variables

      integer,                                  intent(in) :: im, levs, tracers_total, ntrac, ntcw, ntiw, ntclamt, ntrw,  &
        ntsw, ntrnc, ntsnc, ntgl, ntgnc, ntlnc, ntinc, nn, imp_physics, imp_physics_gfdl, imp_physics_thompson,           &
        imp_physics_zhao_carr, imp_physics_zhao_carr_pdf

      logical,                                  intent(in) :: ltaerosol, lgocart

      real(kind=kind_phys),                     intent(in) :: con_pi, dtf
      real(kind=kind_phys), dimension(im,levs), intent(in) :: save_qc, save_qi

      real(kind=kind_phys), dimension(im,levs,ntrac), intent(inout) :: gq0
      real(kind=kind_phys), dimension(im,levs,nn),    intent(inout) :: clw
      real(kind=kind_phys), dimension(im,levs),       intent(inout) :: dqdti

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! local variables
      integer :: i,k,n,tracers

      real(kind=kind_phys) :: liqm, icem

      liqm = 4./3.*con_pi*1.e-12
      icem = 4./3.*con_pi*3.2768*1.e-14*890.

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  --- update the tracers due to deep & shallow cumulus convective transport
!           (except for suspended water and ice)

      if (tracers_total > 0) then
        tracers = 2
        do n=2,ntrac
!         if ( n /= ntcw .and. n /= ntiw .and. n /= ntclamt) then
          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc ) then
              tracers = tracers + 1
            do k=1,levs
              do i=1,im
                gq0(i,k,n) = clw(i,k,tracers)
              enddo
            enddo
          endif
        enddo
      endif

      if (ntcw > 0) then

!  for microphysics
        if (imp_physics == imp_physics_zhao_carr_pdf .or. imp_physics == imp_physics_zhao_carr    &
                               .or. imp_physics == imp_physics_gfdl) then
           gq0(1:im,:,ntcw) = clw(1:im,:,1) + clw(1:im,:,2)
        elseif (ntiw > 0) then
          do k=1,levs
            do i=1,im
              gq0(i,k,ntiw) = clw(i,k,1)                     ! ice
              gq0(i,k,ntcw) = clw(i,k,2)                     ! water
            enddo
          enddo
          if (imp_physics == imp_physics_thompson) then
            if (ltaerosol) then
              do k=1,levs
                do i=1,im
                  gq0(i,k,ntlnc) = gq0(i,k,ntlnc)  &
                           +  max(0.0, (clw(i,k,2)-save_qc(i,k))) / liqm
                  gq0(i,k,ntinc) = gq0(i,k,ntinc)  &
                           +  max(0.0, (clw(i,k,1)-save_qi(i,k))) / icem
                enddo
              enddo
            else
              do k=1,levs
                do i=1,im
                  gq0(i,k,ntinc) = gq0(i,k,ntinc)  &
                           +  max(0.0, (clw(i,k,1)-save_qi(i,k))) / icem
                enddo
              enddo
            endif
          endif
        else
          do k=1,levs
            do i=1,im
              gq0(i,k,ntcw) = clw(i,k,1) + clw(i,k,2)
            enddo
          enddo
        endif   ! end if_ntiw
      else
        do k=1,levs
          do i=1,im
            clw(i,k,1) = clw(i,k,1) + clw(i,k,2)
          enddo
        enddo
      endif   ! end if_ntcw

! dqdt_v : instaneous moisture tendency (kg/kg/sec)
      if (lgocart) then
        do k=1,levs
          do i=1,im
            dqdti(i,k) = dqdti(i,k) * (1.0 / dtf)
          enddo
        enddo
      endif

    end subroutine GFS_suite_interstitial_4_run

  end module GFS_suite_interstitial_4
