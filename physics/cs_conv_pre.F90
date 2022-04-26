!>  \file cs_conv_pre.F90
!!  This file contains preparation for the Chikira-Sugiyama Convection scheme.

module cs_conv_pre
  contains

!! \section arg_table_cs_conv_pre_run Argument Table
!! \htmlinclude cs_conv_pre_run.html
!!
  subroutine cs_conv_pre_run(im, levs, ntrac, q, clw1, clw2,            &
     &                       work1, work2, cs_parm1, cs_parm2, wcbmax,  &
     &                       fswtr, fscav, save_q1, save_q2, save_q3,   &
     &                       errmsg, errflg)


  use machine ,   only : kind_phys

  implicit none

! --- inputs
  integer, intent(in) :: im, levs, ntrac
  real(kind_phys), dimension(:,:), intent(in) :: q
  real(kind_phys), dimension(:,:), intent(in) :: clw1,clw2
  real(kind_phys), dimension(:),   intent(in) :: work1, work2
  real(kind_phys), intent(in) :: cs_parm1, cs_parm2

! --- input/output
  real(kind_phys), dimension(:), intent(out) :: fswtr, fscav
  real(kind_phys), dimension(:), intent(out) :: wcbmax
  real(kind_phys), dimension(:,:), intent(out) :: save_q1,save_q2
  ! save_q3 is not allocated for Zhao-Carr MP
  real(kind_phys), dimension(:,:), intent(out)     :: save_q3

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
  do k=1,levs
    do i=1,im
      ! DH* note - save_q1 assignment may be redundant,
      ! because already done in GFS_DCNV_generic_pre?
      ! Keep for using cs_conv w/o GFS_DCNV_generic_pre?
      save_q1(i,k) = q(i,k)
      save_q2(i,k) = max(0.0,clw2(i,k))
      save_q3(i,k) = max(0.0,clw1(i,k))
    enddo
  enddo

  return
  end subroutine cs_conv_pre_run

end module cs_conv_pre