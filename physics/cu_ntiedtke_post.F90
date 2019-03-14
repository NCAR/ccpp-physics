!> \file cu_ntiedtke_post.F90
!!  Contains code related to New Tiedtke convective scheme

module cu_ntiedtke_post

   implicit none

   private

   public :: cu_ntiedtke_post_init, cu_ntiedtke_post_run, cu_ntiedtke_post_finalize

   contains

   subroutine cu_ntiedtke_post_init ()
   end subroutine cu_ntiedtke_post_init

   subroutine cu_ntiedtke_post_finalize()
   end subroutine cu_ntiedtke_post_finalize

!> \section arg_table_cu_ntiedtke_post_run Argument Table
!! | local_name     | standard_name                                          | long_name                                        | units   | rank | type      |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|--------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | t              | air_temperature_updated_by_physics                     | temperature updated by physics                   | K       |    2 | real      | kind_phys | in     | F        |
!! | q              | water_vapor_specific_humidity_updated_by_physics       | water vapor specific humidity updated by physics | kg kg-1 |    2 | real      | kind_phys | in     | F        |
!! | prevst         | temperature_from_previous_timestep                     | temperature from previous time step              | K       |    2 | real      | kind_phys | out    | F        |
!! | prevsq         | moisture_from_previous_timestep                        | moisture from previous time step                 | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP         | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP            | flag    |    0 | integer   |           | out    | F        |
!!
   subroutine cu_ntiedtke_post_run (t, q, prevst, prevsq, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! Interface variables
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys),  intent(out) :: prevst(:,:)
      real(kind_phys),  intent(out) :: prevsq(:,:)
      character(len=*), intent(out) :: errmsg
      integer, intent(out)          :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      prevst(:,:) = t(:,:)
      prevsq(:,:) = q(:,:)

   end subroutine cu_ntiedtke_post_run

end module cu_ntiedtke_post
