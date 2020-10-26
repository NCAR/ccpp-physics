!> \file cu_gf_driver_post.F90
!!  Contains code related to GF convective schemes to be used within the GFS physics suite.

module cu_gf_driver_post

   implicit none

   private

   public :: cu_gf_driver_post_init, cu_gf_driver_post_run, cu_gf_driver_post_finalize

   contains

   subroutine cu_gf_driver_post_init ()
   end subroutine cu_gf_driver_post_init

   subroutine cu_gf_driver_post_finalize()
   end subroutine cu_gf_driver_post_finalize

!> \section arg_table_cu_gf_driver_post_run Argument Table
!! \htmlinclude cu_gf_driver_post_run.html
!!
   subroutine cu_gf_driver_post_run (im, t, q, prevst, prevsq, cactiv, conv_act, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! Interface variables
      integer,          intent(in)  :: im
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys),  intent(out) :: prevst(:,:)
      real(kind_phys),  intent(out) :: prevsq(:,:)
      integer,          intent(in)  :: cactiv(:)
      real(kind_phys),  intent(out) :: conv_act(:)
      character(len=*), intent(out) :: errmsg
      integer, intent(out)          :: errflg

      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      prevst(:,:) = t(:,:)
      prevsq(:,:) = q(:,:)

      do i = 1, im
        if (cactiv(i).gt.0) then
          conv_act(i) = conv_act(i)+1.0
        else
          conv_act(i)=0.0
        endif
      enddo

   end subroutine cu_gf_driver_post_run

end module cu_gf_driver_post
