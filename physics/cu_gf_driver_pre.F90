!> \file cu_gf_driver_pre.F90
!!  Contains code related to GF convective schemes to be used within the GFS physics suite.

module cu_gf_driver_pre

   implicit none

   private

   public :: cu_gf_driver_pre_init, cu_gf_driver_pre_run, cu_gf_driver_pre_finalize

   contains

   subroutine cu_gf_driver_pre_init ()
   end subroutine cu_gf_driver_pre_init

   subroutine cu_gf_driver_pre_finalize()
   end subroutine cu_gf_driver_pre_finalize

!> \section arg_table_cu_gf_driver_pre_run Argument Table
!! \htmlinclude cu_gf_driver_pre_run.html
!!
   subroutine cu_gf_driver_pre_run (flag_init, flag_restart, kdt, fhour, dtp, t, q, prevst, prevsq, &
                                    forcet, forceq, cactiv, conv_act, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      logical,          intent(in)  :: flag_init
      logical,          intent(in)  :: flag_restart
      integer,          intent(in)  :: kdt
      real(kind_phys),  intent(in)  :: fhour
      real(kind_phys),  intent(in)  :: dtp
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys),  intent(in)  :: prevst(:,:)
      real(kind_phys),  intent(in)  :: prevsq(:,:)
      real(kind_phys),  intent(out) :: forcet(:,:)
      real(kind_phys),  intent(out) :: forceq(:,:)
      integer,          intent(out) :: cactiv(:)
      real(kind_phys),  intent(in)  :: conv_act(:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! local variables
      real(kind=kind_phys) :: dtdyn

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! For restart runs, can assume that prevst and prevsq
      ! are read from the restart files beforehand, same
      ! for conv_act.
      if(flag_init .and. .not.flag_restart) then
        forcet(:,:)=0.0
        forceq(:,:)=0.0
      else
        dtdyn=3600.0*(fhour)/kdt
        if(dtp > dtdyn) then
          forcet(:,:)=(t(:,:) - prevst(:,:))/dtp
          forceq(:,:)=(q(:,:) - prevsq(:,:))/dtp
        else
          forcet(:,:)=(t(:,:) - prevst(:,:))/dtdyn
          forceq(:,:)=(q(:,:) - prevsq(:,:))/dtdyn
        endif
      endif

      cactiv(:)=nint(conv_act(:))

   end subroutine cu_gf_driver_pre_run

end module cu_gf_driver_pre
