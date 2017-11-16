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
!! | tottracer      | number_of_total_tracers                                | total number of tracers                                               | count         |    0 | integer                       |           |   out  | F        |
!! | trc_shft       | start_index_of_other_tracers                           | beginning index of the non-water tracer species                       | index         |    0 | integer                       |           |   out  | F        |
!! | tracers        | number_of_water_tracers                                | number of water-related tracers                                       | index         |    0 | integer                       |           |   out  | F        |
!! | ntk            | index_of_TKE                                           | index of TKE in the tracer array                                      | index         |    0 | integer                       |           |   out  | F        |
!!
      subroutine GFS_suite_interstitial_1_run (Model, Grid, tottracer, trc_shft, tracers, ntk, skip_macro, clw, cnvc, cnvw)

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
