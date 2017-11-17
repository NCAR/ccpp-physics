!> \file GFS_suite_interstitial.f90
!!  Contains code related to more than one scheme in the GFS physics suite.

      module GFS_suite_interstitial_1

      contains

      subroutine GFS_suite_interstitial_1_init ()
      end subroutine GFS_suite_interstitial_1_init

      subroutine GFS_suite_interstitial_1_finalize()
      end subroutine GFS_suite_interstitial_1_finalize

!> \section arg_table_GFS_suite_interstitial_1_run Argument Table
!! | local var name | longname                                               | description                                                           | units         | rank | type                          |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_typedefs%GFS_control_type |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_typedefs%GFS_grid_type    |           | in     | F        |
!! | tottracer      | number_of_total_tracers                                | total number of tracers                                               | count         |    0 | integer                       |           |   out  | F        |
!! | trc_shft       | start_index_of_other_tracers                           | beginning index of the non-water tracer species                       | index         |    0 | integer                       |           |   out  | F        |
!! | tracers        | number_of_water_tracers                                | number of water-related tracers                                       | index         |    0 | integer                       |           |   out  | F        |
!! | ntk            | index_of_TKE                                           | index of TKE in the tracer array                                      | index         |    0 | integer                       |           |   out  | F        |
!! | skip_macro     | flag_skip_macro                                        | flag to skip cloud macrophysics in Morrison scheme                    | flag          |    1 | logical                       |           |   out  | F        |
!! | clw            | convective_transportable_tracers                       | array to contain cloud water and other convective trans. tracers      | kg kg-1       |    3 | real                          | kind_phys |   out  | F        |
!! | cnvw           | convective_cloud_water_specific_humidity               | convective cloud water specific humidity                              | kg kg-1       |    2 | real                          | kind_phys |   out  | F        |
!! | cnvc           | convective_cloud_cover                                 | convective cloud cover                                                | frac          |    2 | real                          | kind_phys |   out  | F        |
!!
      subroutine GFS_suite_interstitial_1_run (Model, Grid, tottracer, trc_shft, tracers, ntk, skip_macro, clw, cnvc, cnvw)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_grid_type

        type(GFS_control_type),           intent(in) :: Model
        type(GFS_grid_type),              intent(in) :: Grid
        integer,                          intent(out) :: tottracer, trc_shft, tracers, ntk
        logical, dimension(size(Grid%xlon,1)), intent(out) :: skip_macro
        real(kind=kind_phys), allocatable, intent(out) :: clw(:,:,:), cnvc(:,:), cnvw(:,:)

        tottracer = 0            ! no convective transport of tracers
        if (Model%trans_trac .or. Model%cscnv) then
          if (Model%ntcw > 0) then
            if (Model%ntoz < Model%ntcw) then
              trc_shft = Model%ntcw + Model%ncld - 1
            else
              trc_shft = Model%ntoz
            endif
          elseif (Model%ntoz > 0) then
            trc_shft = Model%ntoz
          else
            trc_shft = 1
          endif

          tracers   = Model%ntrac - trc_shft
          tottracer = tracers
          if (Model%ntoz > 0) tottracer = tottracer + 1  ! ozone is added separately
        endif
        if (Model%ntke > 0) ntk = Model%ntke - trc_shft + 3

        skip_macro = .false.

        allocate ( clw(size(Grid%xlon,1),Model%levs,tottracer+2) )
        if (Model%imfdeepcnv >= 0 .or. Model%imfshalcnv > 0) then
          allocate (cnvc(size(Grid%xlon,1),Model%levs), cnvw(size(Grid%xlon,1),Model%levs))
        endif

      end subroutine GFS_suite_interstitial_1_run

    end module

    module GFS_suite_interstitial_2

    contains

    subroutine GFS_suite_interstitial_2_init ()
    end subroutine GFS_suite_interstitial_2_init

    subroutine GFS_suite_interstitial_2_finalize()
    end subroutine GFS_suite_interstitial_2_finalize

!> \section arg_table_GFS_suite_interstitial_2_run Argument Table
!! | local var name | longname                                                                | description                                                             | units         | rank | type                          |    kind   | intent | optional |
!! |----------------|-------------------------------------------------------------------------|-------------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                                    | Fortran DDT containing FV3-GFS model control parameters                 | DDT           |    0 | GFS_typedefs%GFS_control_type |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                                       | Fortran DDT containing FV3-GFS grid and interpolation related data      | DDT           |    0 | GFS_typedefs%GFS_grid_type    |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                                    | Fortran DDT containing FV3-GFS surface fields                           | DDT           |    0 | GFS_typedefs%GFS_sfcprop_type |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                                    | Fortran DDT containing FV3-GFS prognostic state data in from dycore     | DDT           |    0 | GFS_typedefs%GFS_statein_type |           | in     | F        |
!! | Diag           | FV3-GFS_diag_type                                                       | Fortran DDT containing FV3-GFS fields targeted for diagnostic output    | DDT           |    0 | GFS_typedefs%GFS_diag_type    |           | inout  | F        |
!! | rhbbot         | critical_relative_humidity_at_surface                                   | critical relative humidity at the surface                               | frac          |    0 | real                          | kind_phys |   out  | F        |
!! | rhpbl          | critical_relative_humidity_at_PBL_top                                   | critical relative humidity at the PBL top                               | frac          |    0 | real                          | kind_phys |   out  | F        |
!! | rhbtop         | critical_relative_humidity_at_top_of_atmosphere                         | critical relative humidity at the top of atmosphere                     | frac          |    0 | real                          | kind_phys |   out  | F        |
!! | frain          | dynamics_to_physics_timestep_ratio                                      | ratio of dynamics timestep to physics timestep                          | none          |    0 | real                          | kind_phys |   out  | F        |
!! | islmsk         | sea_land_ice_mask                                                       | landmask: sea/land/ice=0/1/2                                            | flag          |    1 | integer                       |           |   out  | F        |
!! | work1          | grid_related_coefficient                                                | grid size related coefficient used in scale-sensitive schemes           | none          |    1 | real                          | kind_phys |   out  | F        |
!! | work2          | grid_related_coefficient_complement                                     | complement to work1                                                     | none          |    1 | real                          | kind_phys |   out  | F        |
!! | dudt           | tendency_of_x_wind_due_to_model_physics                                 | updated tendency of the x wind                                          | m s-2         |    2 | real                          | kind_phys |   out  | F        |
!! | dvdt           | tendency_of_y_wind_due_to_model_physics                                 | updated tendency of the y wind                                          | m s-2         |    2 | real                          | kind_phys |   out  | F        |
!! | dtdt           | tendency_of_air_temperature_due_to_model_physics                        | updated tendency of the temperature                                     | K s-1         |    2 | real                          | kind_phys |   out  | F        |
!! | dtdtc          | tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky | clear sky radiative (shortwave + longwave) heating rate at current time | K s-1         |    2 | real                          | kind_phys |   out  | F        |
!! | dqdt           | tendency_of_tracers_due_to_model_physics                                | updated tendency of the tracers                                         | kg kg-1 s-1   |    3 | real                          | kind_phys |   out  | F        |
!!
    subroutine GFS_suite_interstitial_2_run (Model, Grid, Sfcprop, Statein, Diag, rhbbot, rhpbl, rhbtop, frain, islmsk, work1, work2, dudt, dvdt, dtdt, dtdtc, dqdt)

      use machine,               only: kind_phys
      use physcons,              only: dxmin, dxinv
      use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_sfcprop_type, GFS_statein_type, GFS_diag_type

      type(GFS_control_type),           intent(in) :: Model
      type(GFS_grid_type),              intent(in) :: Grid
      type(GFS_sfcprop_type),           intent(in) :: Sfcprop
      type(GFS_statein_type),           intent(in) :: Statein
      type(GFS_diag_type),              intent(inout) :: Diag

      real(kind=kind_phys), intent(out) :: rhbbot, rhpbl, rhbtop, frain
      integer, dimension(size(Grid%xlon,1)), intent(out) :: islmsk
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(out)  :: work1, work2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(out) :: dudt, dvdt, dtdt, dtdtc
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,Model%ntrac), intent(out) ::  dqdt

      integer :: i

      rhbbot = Model%crtrh(1)
      rhpbl  = Model%crtrh(2)
      rhbtop = Model%crtrh(3)

      frain = Model%dtf / Model%dtp

      do i = 1, size(Grid%xlon,1)
        islmsk(i)   = nint(Sfcprop%slmsk(i))
        work1(i)   = (log(Grid%area(i)) - dxmin) * dxinv
        work1(i)   = max(0.0, min(1.0,work1(i)))
        work2(i)   = 1.0 - work1(i)
        Diag%psurf(i)   = Statein%pgr(i)
      end do

      dudt(:,:)  = 0.
      dvdt(:,:)  = 0.
      dtdt(:,:)  = 0.
      dtdtc(:,:) = 0.
      dqdt(:,:,:) = 0.

    end subroutine GFS_suite_interstitial_2_run

  end module

  module GFS_suite_interstitial_3

  contains

  subroutine GFS_suite_interstitial_3_init ()
  end subroutine GFS_suite_interstitial_3_init

  subroutine GFS_suite_interstitial_3_finalize()
  end subroutine GFS_suite_interstitial_3_finalize

!> \section arg_table_GFS_suite_interstitial_3_run Argument Table
!! | local var name | longname                                                     | description                                                           | units         | rank | type                          |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                         | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_typedefs%GFS_control_type |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                            | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_typedefs%GFS_grid_type    |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                         | Fortran DDT containing FV3-GFS prognostic state data in from dycore   | DDT           |    0 | GFS_typedefs%GFS_statein_type |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                         | Fortran DDT containing FV3-GFS radiation tendencies needed in physics | DDT           |    0 | GFS_typedefs%GFS_radtend_type |           | in     | F        |
!! | xcosz          | instantaneous_cosine_of_zenith_angle                         | cosine of zenith angle at current time                                | none          | 1    | real                          | kind_phys | in     | F        |
!! | adjsfcdsw      | surface_downwelling_shortwave_flux                           | surface downwelling shortwave flux at current time                    | W m-2         | 1    | real                          | kind_phys | in     | F        |
!! | adjsfcdlw      | surface_downwelling_longwave_flux                            | surface downwelling longwave flux at current time                     | W m-2         | 1    | real                          | kind_phys | in     | F        |
!! | adjsfculw      | surface_upwelling_longwave_flux                              | surface upwelling longwave flux at current time                       | W m-2         | 1    | real                          | kind_phys | in     | F        |
!! | xmu            | zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes | zenith angle temporal adjustment factor for shortwave fluxes          | none          | 1    | real                          | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_diag_type                                            | Fortran DDT containing FV3-GFS fields targeted for diagnostic output  | DDT           |    0 | GFS_typedefs%GFS_diag_type    |           | inout  | F        |
!! | kcnv           | flag_deep_convection                                         | flag indicating whether convection occurs in column (0 or 1)          | index         | 1    | integer                       |           |   out  | F        |
!! | heat           | kinematic_surface_upward_sensible_heat_flux                  | kinematic surface upward sensible heat flux                           | K m s-1       |    1 | real                          | kind_phys |   out  | F        |
!! | evap           | kinematic_surface_upward_latent_heat_flux                    | kinematic surface upward latent heat flux                             | kg kg-1 m s-1 |    1 | real                          | kind_phys |   out  | F        |
!!
  subroutine GFS_suite_interstitial_3_run (Model, Grid, Statein, Radtend, xcosz, adjsfcdsw, adjsfcdlw, adjsfculw, xmu, Diag, kcnv, hflx, evap)

    use machine,               only: kind_phys
    use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_statein_type, GFS_radtend_type, GFS_diag_type

    type(GFS_control_type),           intent(in)    :: Model
    type(GFS_grid_type),              intent(in)    :: Grid
    type(GFS_statein_type),           intent(in)    :: Statein
    type(GFS_radtend_type),           intent(in)    :: Radtend
    type(GFS_diag_type),              intent(inout) :: Diag

    integer, dimension(size(Grid%xlon,1)), intent(out) :: kcnv
    real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: xcosz, adjsfcdsw, adjsfcdlw, adjsfculw, xmu
    real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(out) :: hflx, evap

    real(kind=kind_phys), parameter :: czmin   = 0.0001      ! cos(89.994)

    integer :: i

    real(kind=kind_phys) :: tem1

    if (Model%lssav) then      !  --- ...  accumulate/save output variables

  !  --- ...  sunshine duration time is defined as the length of time (in mdl output
  !           interval) that solar radiation falling on a plane perpendicular to the
  !           direction of the sun >= 120 w/m2

      do i = 1, size(Grid%xlon,1)
        if ( xcosz(i) >= czmin ) then   ! zenth angle > 89.994 deg
          tem1 = adjsfcdsw(i) / xcosz(i)
          if ( tem1 >= 120.0 ) then
            Diag%suntim(i) = Diag%suntim(i) + Model%dtf
          endif
        endif
      enddo

  !  --- ...  sfc lw fluxes used by atmospheric model are saved for output

      Diag%dlwsfc(:) = Diag%dlwsfc(:) +   adjsfcdlw(:)*Model%dtf
      Diag%ulwsfc(:) = Diag%ulwsfc(:) +   adjsfculw(:)*Model%dtf
      Diag%psmean(:) = Diag%psmean(:) + Statein%pgr(:)*Model%dtf        ! mean surface pressure

      if (Model%ldiag3d) then
        if (Model%lsidea) then
          Diag%dt3dt(:,:,1) = Diag%dt3dt(:,:,1) + Radtend%lwhd(:,:,1)*Model%dtf
          Diag%dt3dt(:,:,2) = Diag%dt3dt(:,:,2) + Radtend%lwhd(:,:,2)*Model%dtf
          Diag%dt3dt(:,:,3) = Diag%dt3dt(:,:,3) + Radtend%lwhd(:,:,3)*Model%dtf
          Diag%dt3dt(:,:,4) = Diag%dt3dt(:,:,4) + Radtend%lwhd(:,:,4)*Model%dtf
          Diag%dt3dt(:,:,5) = Diag%dt3dt(:,:,5) + Radtend%lwhd(:,:,5)*Model%dtf
          Diag%dt3dt(:,:,6) = Diag%dt3dt(:,:,6) + Radtend%lwhd(:,:,6)*Model%dtf
        else
          do k = 1, Model%levs
            Diag%dt3dt(:,k,1) = Diag%dt3dt(:,k,1) + Radtend%htrlw(:,k)*Model%dtf
            Diag%dt3dt(:,k,2) = Diag%dt3dt(:,k,2) + Radtend%htrsw(:,k)*Model%dtf*xmu(:)
          enddo
        endif
      endif
    endif    ! end if_lssav_block

    kcnv(:)   = 0

    hflx(:)       = 0.0
    evap(:)       = 0.0

  end subroutine GFS_suite_interstitial_3_run

end module
