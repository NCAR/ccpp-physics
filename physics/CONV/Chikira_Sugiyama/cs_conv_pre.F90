!>  \file cs_conv_pre.F90
!!  This file contains preparation for the Chikira-Sugiyama Convection scheme.

module cs_conv_pre
  contains

!! \section arg_table_cs_conv_pre_run Argument Table
!! \htmlinclude cs_conv_pre_run.html
!!
  subroutine cs_conv_pre_run(im, levs, work1, work2, cs_parm1, cs_parm2, wcbmax,  &
     &                       fswtr, fscav, errmsg, errflg)


  use machine ,   only : kind_phys

  implicit none

! --- inputs
  integer, intent(in) :: im, levs
  real(kind_phys), dimension(:),   intent(in) :: work1, work2
  real(kind_phys), intent(in) :: cs_parm1, cs_parm2

! --- input/output
  real(kind_phys), dimension(:), intent(out) :: fswtr, fscav
  real(kind_phys), dimension(:), intent(out) :: wcbmax

  character(len=*), intent(out) :: errmsg
  integer,          intent(out) :: errflg

! --- locals
  integer :: i, k

  ! Initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  do i =1,im
   wcbmax(i) = cs_parm1 * work1(i) + cs_parm2 * work2(i)
  enddo

  fswtr(:) = 0.0
  fscav(:) = 0.0

  return
  end subroutine cs_conv_pre_run

end module cs_conv_pre