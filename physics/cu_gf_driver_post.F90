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
   subroutine cu_gf_driver_post_run (im, km, t, q, prevst, prevsq, cactiv, cactiv_m, conv_act, conv_act_m,dt, garea, raincv, maxupmf, refl_10cm, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! Interface variables
      integer,          intent(in)  :: im, km
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys), dimension(:),intent(in) :: garea
      real(kind_phys),  intent(out) :: prevst(:,:)
      real(kind_phys),  intent(out) :: prevsq(:,:)
      integer,          intent(in)  :: cactiv(:)
      integer,          intent(in)  :: cactiv_m(:)
      real(kind_phys),  intent(out) :: conv_act(:)
      real(kind_phys),  intent(out) :: conv_act_m(:)
      ! for Radar reflectivity
      real(kind_phys),  intent(in)  :: dt
      real(kind_phys),  intent(in)  :: raincv(:), maxupmf(:)
      real(kind_phys),  intent(inout) :: refl_10cm(:,:)
      character(len=*), intent(out) :: errmsg
!$acc declare copyin(t,q,cactiv,cactiv_m) copyout(prevst,prevsq,conv_act,conv_act_m)
      integer, intent(out)          :: errflg

      ! Local variables
      real(kind_phys), parameter :: dbzmin=-10.0
      real(kind_phys) :: cuprate
      real(kind_phys) :: ze, ze_conv, dbz_sum
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
        ! reflectivity parameterization for parameterized convection (reference:Unipost MDLFLD.f)
        if(sqrt(garea(i)).lt.6500.)then
        ze      = 0.0
        ze_conv = 0.0
        dbz_sum = 0.0
        cuprate = raincv(i) * 3600.0 / dt          ! cu precip rate (mm/h)
        ze_conv = 300.0 * cuprate**1.4
        if (maxupmf(i).gt.0.05) then
         do k = 1, km
          ze = 10._kind_phys ** (0.1 * refl_10cm(i,k))
          dbz_sum = max(dbzmin, 10.0 * log10(ze + ze_conv))
          refl_10cm(i,k) = dbz_sum
         enddo
        endif
        endif
      enddo
!$acc end kernels

   end subroutine cu_gf_driver_post_run

end module cu_gf_driver_post
