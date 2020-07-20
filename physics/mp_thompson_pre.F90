!>\file mp_thompson_pre.F90
!!

! CCPP license goes here, as well as further documentation
!>\ingroup aathompson
module mp_thompson_pre

      use machine, only : kind_phys

      implicit none

      public :: mp_thompson_pre_init, mp_thompson_pre_run, mp_thompson_pre_finalize

      private

   contains

      subroutine mp_thompson_pre_init()
      end subroutine mp_thompson_pre_init

!! \section arg_table_mp_thompson_pre_run Argument Table
!! \htmlinclude mp_thompson_pre_run.html
!!
      subroutine mp_thompson_pre_run(ncol, nlev, tgrs, tgrs_save, errmsg, errflg)

         implicit none

         ! Interface variables
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         real(kind_phys),           intent(in   ) :: tgrs(1:ncol,1:nlev)
         real(kind_phys),           intent(  out) :: tgrs_save(1:ncol,1:nlev)

         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Save current air temperature for tendency limiters in mp_thompson_post
         tgrs_save = tgrs

      end subroutine mp_thompson_pre_run

      subroutine mp_thompson_pre_finalize()
      end subroutine mp_thompson_pre_finalize

end module mp_thompson_pre
