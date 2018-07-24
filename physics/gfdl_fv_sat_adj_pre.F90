!> \file gfdl_fv_sat_adj_pre.F90
!!  Contains code related to preparing the gfdl_fv_sat_adj runs.

    module fv_sat_adj_pre

    contains

    subroutine fv_sat_adj_pre_init ()
    end subroutine fv_sat_adj_pre_init

    subroutine fv_sat_adj_pre_finalize()
    end subroutine fv_sat_adj_pre_finalize

!> \section arg_table_fv_sat_adj_pre_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                   |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|------------------------|-----------|--------|----------|
!! | Interstitial   | CCPP_Interstitial_type                                 | derived type CCPP_interstitial_type                     | DDT           |    0 | CCPP_interstitial_type |           | inout  | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                | none          |    0 | character              | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                   | flag          |    0 | integer                |           | out    | F        |
!!
    subroutine fv_sat_adj_pre_run (Interstitial, errmsg, errflg)

      use CCPP_typedefs, only: CCPP_interstitial_type

      implicit none

      ! interface variables
      type(CCPP_interstitial_type), intent(inout) :: Interstitial
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call Interstitial%reset()

    end subroutine fv_sat_adj_pre_run

    end module fv_sat_adj_pre
