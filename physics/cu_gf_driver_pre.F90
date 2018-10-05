!> \file cu_gf_driver_pre.F90
!!  Contains code related to GF convective schemes to be used within the GFS physics suite.

module cu_gf_driver_pre

   implicit none

   private

   public :: cu_gf_driver_pre_init, cu_gf_driver_pre_run, cu_gf_driver_pre_finalize

   contains

   subroutine cu_gf_driver_pre_init ()
   end subroutine cu_gf_driver_pre_init

   subroutine cu_gf_driver_pre_finalize()
   end subroutine cu_gf_driver_pre_finalize

!> \section arg_table_cu_gf_driver_pre_run Argument Table
!! | local_name     | standard_name                                          | long_name                                        | units         | rank | type      |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|--------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | kdt            | index_of_time_step                                     | current forecast iteration                       | index         |    0 | integer   |           | in     | F        |
!! | fhour          | forecast_time                                          | curent forecast time                             | h             |    0 | real      | kind_phys | in     | F        |
!! | dtp            | time_step_for_physics                                  | physics timestep                                 | s             |    0 | real      | kind_phys | in     | F        |
!! | t              | air_temperature                                        | model layer mean temperature                     | K             |    2 | real      | kind_phys | in     | F        |
!! | q              | water_vapor_specific_humidity                          | water vapor specific humidity                    | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | prevst         | temperature_from_previous_timestep                     | temperature from previous time step              | K             |    2 | real      | kind_phys | in     | F        |
!! | prevsq         | moisture_from_previous_timestep                        | moisture from previous time step                 | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | forcet         | temperature_tendency_due_to_dynamics                   | temperature tendency due to dynamics only        | K s-1         |    2 | real      | kind_phys | out    | F        |
!! | forceq         | moisture_tendency_due_to_dynamics                      | moisture tendency due to dynamics only           | kg kg-1 s-1   |    2 | real      | kind_phys | out    | F        |
!! | cactiv         | conv_activity_counter                                  | convective activity memory                       | none          |    1 | integer   |           | out    | F        |
!! | conv_act       | gf_memory_counter                                      | Memory counter for GF                            | none          |    1 | real      | kind_phys | in     | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP         | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP            | flag          |    0 | integer   |           | out    | F        |
!!
   subroutine cu_gf_driver_pre_run (kdt, fhour, dtp, t, q, prevst, prevsq, forcet, forceq, cactiv, conv_act, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      integer,          intent(in)  :: kdt
      real(kind_phys),  intent(in)  :: fhour
      real(kind_phys),  intent(in)  :: dtp
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys),  intent(in)  :: prevst(:,:)
      real(kind_phys),  intent(in)  :: prevsq(:,:)
      real(kind_phys),  intent(out) :: forcet(:,:)
      real(kind_phys),  intent(out) :: forceq(:,:)
      integer,          intent(out) :: cactiv(:)
      real(kind_phys),  intent(in)  :: conv_act(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! local variables
      real(kind=kind_phys) :: dtdyn

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if(kdt.gt.1) then
        dtdyn=3600.0*(fhour)/kdt
        if(dtp > dtdyn) then
          forcet(:,:)=(t(:,:) - prevst(:,:))/dtp
          forceq(:,:)=(q(:,:) - prevsq(:,:))/dtp
        else
          forcet(:,:)=(t(:,:) - prevst(:,:))/dtdyn
          forceq(:,:)=(q(:,:) - prevsq(:,:))/dtdyn
        endif
      else
        forcet(:,:)=0.0
        forceq(:,:)=0.0
      endif

      cactiv(:)=nint(conv_act(:))

   end subroutine cu_gf_driver_pre_run

end module cu_gf_driver_pre
