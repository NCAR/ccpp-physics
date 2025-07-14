!>\file mp_tempo_pre.F90
!!

! CCPP license goes here, as well as further documentation
!>\ingroup aatempo
module mp_tempo_pre

      use machine, only : kind_phys

      implicit none

      public :: mp_tempo_pre_run

      private

   contains

!> \section arg_table_mp_tempo_pre_run Argument Table
!! \htmlinclude mp_tempo_pre_run.html
!!
      subroutine mp_tempo_pre_run(ncol, nlev, tgrs, tgrs_save, errmsg, errflg)

         implicit none

         ! Interface variables
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         real(kind_phys),           intent(in   ) :: tgrs(:,:)
         real(kind_phys),           intent(  out) :: tgrs_save(:,:)

         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Save current air temperature for tendency limiters in mp_tempo_post
         tgrs_save = tgrs

      end subroutine mp_tempo_pre_run

end module mp_tempo_pre
