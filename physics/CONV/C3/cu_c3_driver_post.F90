!> \file cu_c3_driver_post.F90
!!  Contains code related to C3 convective schemes to be used within the GFS physics suite.

module cu_c3_driver_post

   implicit none

   private

   public :: cu_c3_driver_post_run

   contains

!>\ingroup cu_c3_group
!> \section arg_table_cu_c3_driver_post_run Argument Table
!! \htmlinclude cu_c3_driver_post_run.html
!!
   subroutine cu_c3_driver_post_run (im, km, t, q, prevst, prevsq, cactiv, cactiv_m, conv_act, conv_act_m, dt, garea, raincv, maxupmf, refl_10cm, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! Interface variables
      integer,          intent(in)  :: im, km
      real(kind_phys),  intent(in)  :: t(:,:)
      real(kind_phys),  intent(in)  :: q(:,:)
      real(kind_phys), dimension(:),intent(in) :: garea
      real(kind_phys),  intent(out), optional :: prevst(:,:)
      real(kind_phys),  intent(out), optional :: prevsq(:,:)
      integer,          intent(in),  optional :: cactiv(:)
      integer,          intent(in),  optional :: cactiv_m(:)
      real(kind_phys),  intent(out), optional :: conv_act(:)
      real(kind_phys),  intent(out), optional :: conv_act_m(:)
      ! for Radar reflectivity
      real(kind_phys),  intent(in)  :: dt
      real(kind_phys),  intent(in)  :: raincv(:)
      real(kind_phys),  intent(in), optional :: maxupmf(:)
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
        ze      = 0.0
        ze_conv = 0.0
        dbz_sum = 0.0
        cuprate = 1.e3*raincv(i) * 3600.0 / dt          ! cu precip rate (mm/h)
        if(cuprate .lt. 0.05) cuprate=0.
        ze_conv = 300.0 * cuprate**1.5
        if (maxupmf(i).gt.0.1 .and. cuprate.gt.0.) then
         do k = 1, km
          ze = 10._kind_phys ** (0.1 * refl_10cm(i,k))
          dbz_sum = max(dbzmin, 10.0 * log10(ze + ze_conv))
          refl_10cm(i,k) = dbz_sum
         enddo
        endif
      enddo
!$acc end kernels

   end subroutine cu_c3_driver_post_run

end module cu_c3_driver_post
