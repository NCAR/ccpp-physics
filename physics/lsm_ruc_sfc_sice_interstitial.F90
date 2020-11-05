module lsm_ruc_sfc_sice_pre

   use machine, only: kind_phys

   implicit none

   private

   public :: lsm_ruc_sfc_sice_pre_init, lsm_ruc_sfc_sice_pre_run, lsm_ruc_sfc_sice_pre_finalize

contains

   subroutine lsm_ruc_sfc_sice_pre_init ()
   end subroutine lsm_ruc_sfc_sice_pre_init

   subroutine lsm_ruc_sfc_sice_pre_finalize ()
   end subroutine lsm_ruc_sfc_sice_pre_finalize

#if 0
!> \section arg_table_lsm_ruc_sfc_sice_pre_run Argument Table
!! \htmlinclude lsm_ruc_sfc_sice_pre_run.html
!!
#endif
   subroutine lsm_ruc_sfc_sice_pre_run(im, lsoil_ruc, lsoil, kice, land, icy, stc, tslb, tiice, errmsg, errflg)

      implicit none

      ! Interface variables
      integer, intent(in) :: im, lsoil_ruc, lsoil, kice
      logical, dimension(im), intent(in) :: land, icy
!  --- on Noah levels
      real (kind=kind_phys), dimension(im,lsoil), intent(inout) :: stc
!  --- on RUC levels
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(in) :: tslb
      real (kind=kind_phys), dimension(im,kice), intent(inout) :: tiice

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        if (icy(i)) then
          do k=1,kice
            tiice(i,k) = tslb(i,k)
          end do
        else if (.not.land(i)) then
          do k=1,min(lsoil,lsoil_ruc)
            stc(i,k) = tslb(i,k)
          end do
        end if
      end do

      end subroutine lsm_ruc_sfc_sice_pre_run

end module lsm_ruc_sfc_sice_pre

module lsm_ruc_sfc_sice_post

   use machine, only: kind_phys

   implicit none

   private

   public :: lsm_ruc_sfc_sice_post_init, lsm_ruc_sfc_sice_post_run, lsm_ruc_sfc_sice_post_finalize

contains

   subroutine lsm_ruc_sfc_sice_post_init ()
   end subroutine lsm_ruc_sfc_sice_post_init

   subroutine lsm_ruc_sfc_sice_post_finalize ()
   end subroutine lsm_ruc_sfc_sice_post_finalize

#if 0
!> \section arg_table_lsm_ruc_sfc_sice_post_run Argument Table
!! \htmlinclude lsm_ruc_sfc_sice_post_run.html
!!
#endif
   subroutine lsm_ruc_sfc_sice_post_run(im, lsoil_ruc, lsoil, kice, land, icy, stc, tslb, tiice, errmsg, errflg)

      implicit none

      ! Interface variables
      integer, intent(in) :: im, lsoil_ruc, lsoil, kice
      logical, dimension(im), intent(in) :: land, icy
!  --- on Noah levels
      real (kind=kind_phys), dimension(im,lsoil), intent(in) :: stc
      real (kind=kind_phys), dimension(im,kice),  intent(in) :: tiice
!  --- on RUC levels
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: tslb

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        if (icy(i)) then
          do k=1,kice
            tslb(i,k) = tiice(i,k)
          end do
        else if (.not.land(i)) then
          do k=1,min(lsoil,lsoil_ruc)
            tslb(i,k) = stc(i,k)
          end do
        end if
      end do

      end subroutine lsm_ruc_sfc_sice_post_run

end module lsm_ruc_sfc_sice_post
