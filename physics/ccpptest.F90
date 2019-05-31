    module ccpptest

        use machine, only : kind_dyn

    contains

!> \section arg_table_ccpptest_init Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type      |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                   | flag          |    0 | integer   |           | out    | F        |
!!
    subroutine ccpptest_init (errmsg, errflg)

        implicit none

        ! interface variables
        character(len=*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        write(0,*) "DH DEBUG: IN ccpptest_init"

    end subroutine ccpptest_init

!> \section arg_table_ccpptest_finalize Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type      |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                   | flag          |    0 | integer   |           | out    | F        |
!!
    subroutine ccpptest_finalize(errmsg, errflg)

        implicit none

        ! interface variables
        character(len=*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        write(0,*) "DH DEBUG: IN ccpptest_finalize"

    end subroutine ccpptest_finalize

!> \section arg_table_ccpptest_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type      |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | blkno          | ccpp_block_number                                      | number of block for explicit data blocking in CCPP      | index         |    0 | integer   |           | in     | F        |
!! | convtest1      | dummy_var_unit_conversion_test_1                       | dummy variable 1 to test unit conversions (length)      | m             |    1 | real      | kind_dyn  | inout  | F        |
!! | convtest2      | dummy_var_unit_conversion_test_2                       | dummy variable 2 to test unit conversions (composed)    | Pa            |    2 | real      | kind_dyn  | in     | F        |
!! | convtest3      | dummy_var_unit_conversion_test_3                       | dummy variable 3 to test unit conversions (time)        | s             |    0 | integer   |           | out    | F        |
!! | convtest4      | dummy_var_unit_conversion_test_4                       | dummy variable 4 to test unit conversions (speed)       | km h-1        |    2 | real      | kind_dyn  | inout  | F        |
!! | convtest5      | dummy_var_unit_conversion_test_5                       | dummy variable 4 to test unit conversions (rad. flux)   | erg cm-2 s-1  |    3 | real      | kind_dyn  | inout  | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                   | flag          |    0 | integer   |           | out    | F        |
!!
    subroutine ccpptest_run (blkno, convtest1, convtest2, convtest3, convtest4, convtest5, errmsg, errflg)

      implicit none

      ! interface variables
      integer,          intent(in   ) :: blkno
      real(kind_dyn),   intent(inout) :: convtest1(:)
      real(kind_dyn),   intent(in   ) :: convtest2(:,:)
      integer,          intent(  out) :: convtest3
      real(kind_dyn),   intent(inout) :: convtest4(:,:)
      real(kind_dyn),   intent(inout) :: convtest5(:,:,:)
      character(len=*), intent(  out) :: errmsg
      integer,          intent(  out) :: errflg

      errmsg = ''
      errflg = 0

      if (.not. blkno==1) return

      write(0,*) "DH DEBUG: ccpptest_run  IN, convtest1 in            m is ", convtest1
      write(0,*) "DH DEBUG: ccpptest_run  IN, convtest2 in           Pa is ", convtest2
      !write(0,*) "DH DEBUG: ccpptest_run  IN, convtest3 in            s is ", convtest3
      write(0,*) "DH DEBUG: ccpptest_run  IN, convtest4 in       km h-1 is ", convtest4
      write(0,*) "DH DEBUG: ccpptest_run  IN, convtest5 in erg cm-2 s-1 is ", convtest5

      convtest1 = 10*convtest1
      convtest3 = 42
      convtest4 = 0.1*convtest4
      convtest5 = convtest5 + 1

      write(0,*) "DH DEBUG: ccpptest_run OUT, convtest1 in            m is ", convtest1
      write(0,*) "DH DEBUG: ccpptest_run OUT, convtest2 in           Pa is ", convtest2
      write(0,*) "DH DEBUG: ccpptest_run OUT, convtest3 in            s is ", convtest3
      write(0,*) "DH DEBUG: ccpptest_run OUT, convtest4 in       km h-1 is ", convtest4
      write(0,*) "DH DEBUG: ccpptest_run OUT, convtest5 in erg cm-2 s-1 is ", convtest5

    end subroutine ccpptest_run

    end module ccpptest
