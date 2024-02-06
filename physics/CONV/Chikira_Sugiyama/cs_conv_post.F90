!>  \file cs_conv_post.F90
!!  This file contains code to execute after the Chikira-Sugiyama Convection scheme.

module cs_conv_post
  contains

!> \section arg_table_cs_conv_post_run Argument Table
!! \htmlinclude cs_conv_post_run.html
!!
  subroutine cs_conv_post_run(im, kmax, do_aw, sigmatot, sigmafrac, errmsg, errflg)

  use machine ,   only : kind_phys

  implicit none

! --- inputs
  integer, intent(in)  :: im, kmax
  logical,  intent(in) :: do_aw
  real(kind_phys), dimension(:,:), intent(in) :: sigmatot

! --- input/output
  real(kind_phys), dimension(:,:), intent(out)  :: sigmafrac

  character(len=*), intent(out) :: errmsg
  integer,          intent(out) :: errflg

! --- locals
  integer :: i, k, kk

  ! Initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  if (do_aw) then
    do k=1,kmax
      kk = min(k+1,kmax)  ! assuming no cloud top reaches the model top
      do i=1,im                                               !DD
        sigmafrac(i,k) = 0.5 * (sigmatot(i,k)+sigmatot(i,kk))
      enddo
    enddo
  endif

  return
  end subroutine cs_conv_post_run

end module cs_conv_post