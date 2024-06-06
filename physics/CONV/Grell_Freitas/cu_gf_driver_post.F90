!> \file cu_gf_driver_post.F90
!!  Contains code related to GF convective schemes to be used within the GFS physics suite.

module cu_gf_driver_post

   implicit none

   private

   public :: cu_gf_driver_post_run

   contains

!>\ingroup cu_gf_group
!> \section arg_table_cu_gf_driver_post_run Argument Table
!! \htmlinclude cu_gf_driver_post_run.html
!!
   subroutine cu_gf_driver_post_run (im, km, t, q, prevst, prevsq, cactiv, cactiv_m, conv_act, conv_act_m, rrfs_sd, ntsmoke, ntdust, ntcoarsepm, chem3d, gq0, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! Interface variables
      integer,          intent(in)  :: im, km
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys),  intent(out), optional :: prevst(:,:)
      real(kind_phys),  intent(out), optional :: prevsq(:,:)
      integer,          intent(in),  optional :: cactiv(:)
      integer,          intent(in),  optional :: cactiv_m(:)
      real(kind_phys),  intent(out), optional :: conv_act(:)
      real(kind_phys),  intent(out), optional :: conv_act_m(:)
      logical,          intent(in)  :: rrfs_sd
      integer,          intent(in)  :: ntsmoke, ntdust, ntcoarsepm
      real(kind_phys),  intent(inout), optional :: chem3d(:,:,:)
      real(kind_phys),  intent(inout) :: gq0(:,:,:)
      character(len=*), intent(out) :: errmsg
!$acc declare copyin(t,q,cactiv,cactiv_m) copyout(prevst,prevsq,conv_act,conv_act_m,chem3d,gq0)
      integer, intent(out)          :: errflg

      ! Local variables
      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!$acc kernels
      prevst(:,:) = t(:,:)
      prevsq(:,:) = q(:,:)

      do i = 1, im
        if (cactiv(i).gt.0) then
          conv_act(i) = conv_act(i)+1.0
        else
          conv_act(i)=0.0
        endif
        if (cactiv_m(i).gt.0) then
          conv_act_m(i) = conv_act_m(i)+1.0
        else
          conv_act_m(i)=0.0
        endif
      enddo

      if (rrfs_sd) then
       gq0(:,:,ntsmoke   ) = chem3d(:,:,1)
       gq0(:,:,ntdust    ) = chem3d(:,:,2)
       gq0(:,:,ntcoarsepm) = chem3d(:,:,3)
      endif
!$acc end kernels

   end subroutine cu_gf_driver_post_run

end module cu_gf_driver_post
