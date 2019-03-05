!> \file cu_ntiedtke_pre.F90
!!  Contains code related to New Tiedtke convective scheme

module cu_ntiedtke_pre

   implicit none

   private

   public :: cu_ntiedtke_pre_init, cu_ntiedtke_pre_run, cu_ntiedtke_pre_finalize

   contains

   subroutine cu_ntiedtke_pre_init ()
   end subroutine cu_ntiedtke_pre_init

   subroutine cu_ntiedtke_pre_finalize()
   end subroutine cu_ntiedtke_pre_finalize

!> \section arg_table_cu_ntiedtke_pre_run Argument Table
!! | local_name     | standard_name                                          | long_name                                        | units         | rank | type      |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|--------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | flag_init      | flag_for_first_time_step                       | flag signaling first time step for time integration loop | flag          |    0 | logical   |           | in     | F        |
!! | flag_restart   | flag_for_restart                               | flag for restart (warmstart) or coldstart                | flag          |    0 | logical   |           | in     | F        |
!! | kdt            | index_of_time_step                                     | current forecast iteration                       | index         |    0 | integer   |           | in     | F        |
!! | fhour          | forecast_time                                          | curent forecast time                             | h             |    0 | real      | kind_phys | in     | F        |
!! | dtp            | time_step_for_physics                                  | physics timestep                                 | s             |    0 | real      | kind_phys | in     | F        |
!! | t              | air_temperature                                        | model layer mean temperature                     | K             |    2 | real      | kind_phys | in     | F        |
!! | q              | water_vapor_specific_humidity                          | water vapor specific humidity                    | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | prevst         | temperature_from_previous_timestep                     | temperature from previous time step              | K             |    2 | real      | kind_phys | in     | F        |
!! | prevsq         | moisture_from_previous_timestep                        | moisture from previous time step                 | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | forcet         | temperature_tendency_due_to_dynamics                   | temperature tendency due to dynamics only        | K s-1         |    2 | real      | kind_phys | out    | F        |
!! | forceq         | moisture_tendency_due_to_dynamics                      | moisture tendency due to dynamics only           | kg kg-1 s-1   |    2 | real      | kind_phys | out    | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP         | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP            | flag          |    0 | integer   |           | out    | F        |
!!
   subroutine cu_ntiedtke_pre_run (flag_init, flag_restart, kdt, fhour, dtp, t, q, prevst, prevsq, &
                                   forcet, forceq, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      logical,          intent(in)  :: flag_init
      logical,          intent(in)  :: flag_restart
      integer,          intent(in)  :: kdt
      real(kind_phys),  intent(in)  :: fhour
      real(kind_phys),  intent(in)  :: dtp
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys),  intent(in)  :: prevst(:,:)
      real(kind_phys),  intent(in)  :: prevsq(:,:)
      real(kind_phys),  intent(out) :: forcet(:,:)
      real(kind_phys),  intent(out) :: forceq(:,:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! local variables
      real(kind=kind_phys) :: dtdyn

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! For restart runs, can assume that prevst and prevsq
      ! are read from the restart files beforehand, same
      ! for conv_act.
      if(flag_init .and. .not.flag_restart) then
        forcet(:,:)=0.0
        forceq(:,:)=0.0
      else
        dtdyn=3600.0*(fhour)/kdt
        if(dtp > dtdyn) then
          forcet(:,:)=(t(:,:) - prevst(:,:))/dtp
          forceq(:,:)=(q(:,:) - prevsq(:,:))/dtp
        else
          forcet(:,:)=(t(:,:) - prevst(:,:))/dtdyn
          forceq(:,:)=(q(:,:) - prevsq(:,:))/dtdyn
        endif
      endif

   end subroutine cu_ntiedtke_pre_run

end module cu_ntiedtke_pre
