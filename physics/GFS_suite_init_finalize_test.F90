    module GFS_suite_ini_fini_test

    contains

!> \section arg_table_GFS_suite_ini_fini_test_init Argument Table
!! \htmlinclude GFS_suite_ini_fini_test_init.html
!!
    subroutine GFS_suite_ini_fini_test_init (errmsg, errflg)

        implicit none

        ! interface variables
        character(len=*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        write(0,*) "DH DEBUG: IN GFS_suite_ini_fini_test_init"

    end subroutine GFS_suite_ini_fini_test_init

!> \section arg_table_GFS_suite_ini_fini_test_finalize Argument Table
!! \htmlinclude GFS_suite_ini_fini_test_finalize.html
!!
    subroutine GFS_suite_ini_fini_test_finalize(errmsg, errflg)

        implicit none

        ! interface variables
        character(len=*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        write(0,*) "DH DEBUG: IN GFS_suite_ini_fini_test_finalize"

    end subroutine GFS_suite_ini_fini_test_finalize

!> \section arg_table_GFS_suite_ini_fini_test_run Argument Table
!! \htmlinclude GFS_suite_ini_fini_test_run.html
!!
    subroutine GFS_suite_ini_fini_test_run (errmsg, errflg)

      use GFS_typedefs, only: GFS_interstitial_type

      implicit none

      ! interface variables
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      write(errmsg,'(a)') "DH ERROR: GFS_suite_ini_fini_test_run should not be called"
      errflg = 1

    end subroutine GFS_suite_ini_fini_test_run

    end module GFS_suite_ini_fini_test
