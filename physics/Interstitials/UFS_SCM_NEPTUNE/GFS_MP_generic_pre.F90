!> \file GFS_MP_generic_pre.F90
!! This file contains the subroutines that calculate diagnotics variables
!! before calling any microphysics scheme:

!> This module contains the CCPP-compliant MP generic pre interstitial codes.
      module GFS_MP_generic_pre
      contains

!> \section arg_table_GFS_MP_generic_pre_run Argument Table
!! \htmlinclude GFS_MP_generic_pre_run.html
!!
      subroutine GFS_MP_generic_pre_run(im, levs, ldiag3d, qdiag3d, do_aw, progsigma, ntcw, nncl, &
                                        ntrac, gt0, gq0, save_t, save_q, num_dfi_radar, errmsg, errflg)
!
      use machine,               only: kind_phys

      implicit none
      integer,                                intent(in) :: im, levs, ntcw, nncl, ntrac, num_dfi_radar
      logical,                                intent(in) :: ldiag3d, qdiag3d, do_aw, progsigma
      real(kind=kind_phys), dimension(:,:),   intent(in) :: gt0
      real(kind=kind_phys), dimension(:,:,:), intent(in) :: gq0

      real(kind=kind_phys), dimension(:,:),   intent(inout) :: save_t
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: save_q

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (ldiag3d .or. do_aw .or. num_dfi_radar>0) then
        do k=1,levs
          do i=1,im
            save_t(i,k) = gt0(i,k)
          enddo
        enddo
      endif
      if (ldiag3d .or. do_aw .or. progsigma) then
        if(qdiag3d) then
           do n=1,ntrac
              do k=1,levs
                 do i=1,im
                    save_q(i,k,n) = gq0(i,k,n)
                 enddo
              enddo
           enddo
        else if(do_aw .or. progsigma) then
           ! if qdiag3d, all q are saved already
           save_q(1:im,:,1) = gq0(1:im,:,1)
           do n=ntcw,ntcw+nncl-1
              save_q(1:im,:,n) = gq0(1:im,:,n)
           enddo
        endif
      endif

      end subroutine GFS_MP_generic_pre_run

      end module GFS_MP_generic_pre
