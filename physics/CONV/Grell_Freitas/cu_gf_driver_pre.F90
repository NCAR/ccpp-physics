!> \file cu_gf_driver_pre.F90
!!  Contains code related to GF convective schemes to be used within the GFS physics suite.

module cu_gf_driver_pre

   implicit none

   private

   public :: cu_gf_driver_pre_run

   contains

!>\ingroup cu_gf_group
!> \section arg_table_cu_gf_driver_pre_run Argument Table
!! \htmlinclude cu_gf_driver_pre_run.html
!!
   subroutine cu_gf_driver_pre_run (flag_init, flag_restart, kdt, fhour, dtp, t, q, prevst, prevsq, &
                                    forcet, forceq, cactiv, cactiv_m, conv_act, conv_act_m,         &
                                    rrfs_sd, ntsmoke, ntdust, ntcoarsepm, chem3d, gq0,              &
                                    errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      logical,          intent(in)  :: flag_init
      logical,          intent(in)  :: flag_restart
      logical,          intent(in)  :: rrfs_sd
      integer,          intent(in)  :: kdt
      real(kind_phys),  intent(in)  :: fhour
      real(kind_phys),  intent(in)  :: dtp
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys),  intent(in),    optional :: prevst(:,:)
      real(kind_phys),  intent(in),    optional :: prevsq(:,:)
!$acc declare copyin(t,q,prevst,prevsq)
      real(kind_phys),  intent(out),   optional :: forcet(:,:)
      real(kind_phys),  intent(out),   optional :: forceq(:,:)
      integer,          intent(out),   optional :: cactiv(:)
      integer,          intent(out),   optional :: cactiv_m(:)
      integer,          intent(in)  :: ntsmoke, ntdust, ntcoarsepm
!$acc declare copyout(forcet,forceq,cactiv,cactiv_m)
      real(kind_phys),  intent(in),    optional :: conv_act(:)
      real(kind_phys),  intent(in),    optional :: conv_act_m(:)
      real(kind_phys),  intent(inout), optional :: chem3d(:,:,:)
      real(kind_phys),  intent(inout) :: gq0(:,:,:)
!$acc declare copyin(conv_act,conv_act_m) copy(chem3d,gq0)
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
!$acc kernels
        forcet(:,:)=0.0
        forceq(:,:)=0.0
!$acc end kernels
      else
        dtdyn=3600.0*(fhour)/kdt
        if(dtp > dtdyn) then
!$acc kernels
          forcet(:,:)=(t(:,:) - prevst(:,:))/dtp
          forceq(:,:)=(q(:,:) - prevsq(:,:))/dtp
!$acc end kernels
        else
!$acc kernels
          forcet(:,:)=(t(:,:) - prevst(:,:))/dtdyn
          forceq(:,:)=(q(:,:) - prevsq(:,:))/dtdyn
!$acc end kernels
        endif
      endif

!$acc kernels
      cactiv(:)=nint(conv_act(:))
      cactiv_m(:)=nint(conv_act_m(:))

      if (rrfs_sd) then
       chem3d(:,:,1) = gq0(:,:,ntsmoke)
       chem3d(:,:,2) = gq0(:,:,ntdust)
       chem3d(:,:,3) = gq0(:,:,ntcoarsepm)
      endif
!$acc end kernels

   end subroutine cu_gf_driver_pre_run

end module cu_gf_driver_pre
