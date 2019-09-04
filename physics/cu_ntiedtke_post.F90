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
!! \htmlinclude cu_ntiedtke_post_run.html
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
