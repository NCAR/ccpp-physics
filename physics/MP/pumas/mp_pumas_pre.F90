! ########################################################################################
!>\file mp_pumas_pre.F90
!!

!> This module contains the pre-processing step prior to calling the PUMAS microphysics.
!
! ########################################################################################
module mp_pumas_pre
  use machine, only: kind_phys, kind_dbl_prec

  implicit none
  public mp_pumas_pre_init, mp_pumas_pre_run, mp_pumas_pre_finalize
contains
  ! ######################################################################################
  !> \section arg_table_mp_pumas_pre_init Argument Table
  !! \htmlinclude mp_pumas_pre_init.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_pre_init(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_pre_init

  ! ######################################################################################
  !> \section arg_table_mp_pumas_pre_run Argument Table
  !! \htmlinclude mp_pumas_pre_run.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_pre_run(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_pre_run

  ! ######################################################################################
  !> \section arg_table_mp_pumas_pre_finalize Argument Table
  !! \htmlinclude mp_pumas_pre_finalize.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_pre_finalize(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_pre_finalize

end module mp_pumas_pre
