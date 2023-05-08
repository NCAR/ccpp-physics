!> \file cu_unified_driver_post.F90
!!  Contains code related to unified convective schemes to be used within the GFS physics suite.

module cu_unified_driver_post

   implicit none

   private

   public :: cu_unified_driver_post_run

   contains

!>\ingroup cu_unified_group
!> \section arg_table_cu_unified_driver_post_run Argument Table
!! \htmlinclude cu_unified_driver_post_run.html
!!
   subroutine cu_unified_driver_post_run (im, t, q, prevst, prevsq, cactiv, cactiv_m, conv_act, conv_act_m, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! Interface variables
      integer,          intent(in)  :: im
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys),  intent(out) :: prevst(:,:)
      real(kind_phys),  intent(out) :: prevsq(:,:)
      integer,          intent(in)  :: cactiv(:)
      integer,          intent(in)  :: cactiv_m(:)
      real(kind_phys),  intent(out) :: conv_act(:)
      real(kind_phys),  intent(out) :: conv_act_m(:)
      character(len=*), intent(out) :: errmsg
!$acc declare copyin(t,q,cactiv,cactiv_m) copyout(prevst,prevsq,conv_act,conv_act_m)
      integer, intent(out)          :: errflg

      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!$acc kernels
      prevst(:,:) = t(:,:)
      prevsq(:,:) = q(:,:)

      do i = 1, im
        if (cactiv(i).gt.0) then
          conv_act(i) = conv_act(i)+1.0
        else
          conv_act(i)=0.0
        endif
        if (cactiv_m(i).gt.0) then
          conv_act_m(i) = conv_act_m(i)+1.0
        else
          conv_act_m(i)=0.0
        endif
      enddo
!$acc end kernels

   end subroutine cu_unified_driver_post_run

end module cu_unified_driver_post
