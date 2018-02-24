!>\file GFS_radlw_post
!!This file contains
      module GFS_radlw_post 
      contains

!>\defgroup GFS_radlw_post GFS RRTMG/RADLW Scheme Post
!! @{
!> \section arg_table_GFS_radlw_post_init Argument Table
!!
      subroutine GFS_radlw_post_init()
      end subroutine GFS_radlw_post_init

! PGI compiler does not accept lines longer than 264 characters, remove during pre-processing
#ifndef __PGI
!> \section arg_table_GFS_radlw_post_run Argument Table
!! | local var name    | longname                                                                                      | description                                                                   | units    | rank |  type                         |   kind    | intent    | optional |
!! |-------------------|-----------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------|----------|------|-------------------------------|-----------|-----------|----------|
!! |   Model           | FV3-GFS_Control_type                                                                          | Fortran DDT containing FV3-GFS model control parameters                       | DDT      |  0   | GFS_control_type              |           | in        | F        |
!! |   Grid            | FV3-GFS_Grid_type                                                                             | Fortran DDT containing FV3-GFS grid and interpolation related data            | DDT      |  0   | GFS_grid_type                 |           | in        | F        |
!! |   Radtend         | FV3-GFS_Radtend_type                                                                          | Fortran DDT containing FV3-GFS fields targetted for diagnostic output         | DDT      |  0   | GFS_radtend_type              |           | inout     | F        |
!! |   Coupling        | FV3-GFS_Coupling_type                                                                         | Fortran DDT containing FV3-GFS fields to/from coupling with other components  | DDT      |  0   | GFS_coupling_type             |           | inout     | F        |
!! |   ltp             | extra_top_layer                                                                               | extra top layers                                                              | none     |  0   | integer                       |           | in        | F        |
!! |   lm              | vertical_layer_dimension_for_radiation                                                        | number of vertical layers for radiation calculation                           | count    |  0   | integer                       |           | in        | F        |
!! |   kd              | vertical_index_difference_between_inout_and_local                                             | vertical index difference between in/out and local                            | index    |  0   | integer                       |           | in        | F        |
!! |   tsfa            | surface_air_temperature_for_radiation                                                         | lowest model layer air temperature for radiation                              | K        |  1   | real                          | kind_phys | in        | F        |
!! |   htlwc           | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step                    | total sky heating rate due to longwave radiation                              | K s-1    |  2   | real                          | kind_phys | in        | F        |
!! |   htlw0           | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step | clear sky heating rate due to longwave radiation                              | K s-1    |  2   | real                          | kind_phys | in        | F        |
!!
#endif
      subroutine GFS_radlw_post_run (Model, Grid, Radtend, Coupling,   &
                 ltp, lm, kd, tsfa, htlwc, htlw0)
    
      use machine,                   only: kind_phys
      use GFS_typedefs,              only: GFS_coupling_type,          &
                                           GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type
      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_grid_type),            intent(in)    :: Grid
      type(GFS_radtend_type),         intent(inout) :: Radtend
      integer,                        intent(in)    :: ltp, LM, kd
      real(kind=kind_phys), dimension(size(Grid%xlon,1), Model%levr+LTP), intent(in) ::  htlwc
      real(kind=kind_phys), dimension(size(Grid%xlon,1), Model%levr+LTP), intent(in) ::  htlw0
      real(kind=kind_phys), dimension(size(Grid%xlon,1)),                 intent(in) ::  tsfa
      ! local variables
      integer :: k1, k

      if (Model%lslwr) then
!> -# Save calculation results
!>  - Save surface air temp for diurnal adjustment at model t-steps

        Radtend%tsflw (:) = tsfa(:)

        do k = 1, LM
          k1 = k + kd
            Radtend%htrlw(:,k) = htlwc(:,k1)
        enddo
        ! --- repopulate the points above levr
        if (Model%levr < Model%levs) then
          do k = LM,Model%levs
            Radtend%htrlw (:,k) = Radtend%htrlw (:,LM)
          enddo
        endif

        if (Model%lwhtr) then
          do k = 1, lm
            k1 = k + kd
            Radtend%lwhc(:,k) = htlw0(:,k1)
          enddo
          ! --- repopulate the points above levr
          if (Model%levr < Model%levs) then
            do k = LM,Model%levs
              Radtend%lwhc(:,k) = Radtend%lwhc(:,LM)
            enddo
           endif
         endif

! --- radiation fluxes for other physics processes
        Coupling%sfcdlw(:) = Radtend%sfcflw(:)%dnfxc

      endif                                ! end_if_lslwr
       
      end subroutine GFS_radlw_post_run

!> \section arg_table_GFS_radlw_post_finalize Argument Table
!!
      subroutine GFS_radlw_post_finalize ()
      end subroutine GFS_radlw_post_finalize

!! @}
      end module GFS_radlw_post
