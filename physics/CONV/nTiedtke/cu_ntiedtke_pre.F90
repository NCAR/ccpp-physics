!> \file cu_ntiedtke_pre.F90
!!  Contains code related to New Tiedtke convective scheme

module cu_ntiedtke_pre

   implicit none

   private

   public ::  cu_ntiedtke_pre_run

   contains

!> \section arg_table_cu_ntiedtke_pre_run Argument Table
!! \htmlinclude cu_ntiedtke_pre_run.html
!!
   subroutine cu_ntiedtke_pre_run (flag_init, flag_restart, kdt, fhour, dtp, t, q, prevst, prevsq, &
                                   forcet, forceq, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      logical,          intent(in)  :: flag_init
      logical,          intent(in)  :: flag_restart
      integer,          intent(in)  :: kdt
      real(kind_phys),  intent(in)  :: fhour
      real(kind_phys),  intent(in)  :: dtp
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys),  intent(in),  optional :: prevst(:,:)
      real(kind_phys),  intent(in),  optional :: prevsq(:,:)
      real(kind_phys),  intent(out), optional :: forcet(:,:)
      real(kind_phys),  intent(out), optional :: forceq(:,:)
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

   end subroutine cu_ntiedtke_pre_run

end module cu_ntiedtke_pre
