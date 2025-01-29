! ########################################################################################
!>\file mp_pumas_post.F90
!!

!> This module contains the post-processing step after calling the PUMAS microphysics.
!
! ########################################################################################
module mp_pumas_post
  use machine, only: kind_phys, kind_dbl_postc

  implicit none
  public mp_pumas_post_init, mp_pumas_post_run, mp_pumas_post_finalize
contains
  ! ######################################################################################
  !> \section arg_table_mp_pumas_post_init Argument Table
  !! \htmlinclude mp_pumas_post_init.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_post_init(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_post_init

  ! ######################################################################################
  !> \section arg_table_mp_pumas_post_run Argument Table
  !! \htmlinclude mp_pumas_post_run.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_post_run(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_post_run

  ! ######################################################################################
  !> \section arg_table_mp_pumas_post_finalize Argument Table
  !! \htmlinclude mp_pumas_post_finalize.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_post_finalize(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_post_finalize

end module mp_pumas_post
