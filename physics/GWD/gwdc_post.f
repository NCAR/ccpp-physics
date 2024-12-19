!> \file gwdc_post.f This file contains code to execute after the original code for parameterization of
!! stationary convection forced gravity wave drag based on
!! Chun and Baik (1998) \cite chun_and_baik_1998.

      module gwdc_post

      contains

!> \section arg_table_gwdc_post_run Argument Table
!! \htmlinclude gwdc_post_run.html
!!
      subroutine gwdc_post_run(                                         &
     &  im, levs, lssav, ldiag3d, dtf, dtp, con_cp,                     &
     &  tauctx, taucty, gwdcu, gwdcv,                                   &
     &  dugwd, dvgwd, dtend, dtidx, index_of_x_wind, index_of_y_wind,   &
     &  index_of_process_nonorographic_gwd, gu0, gv0, gt0,              &
     &  errmsg, errflg)

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs
      logical, intent(in) :: lssav, ldiag3d
      real(kind=kind_phys), intent(in) :: dtf, dtp, con_cp
      real(kind=kind_phys), intent(in) ::                               &
     &  tauctx(:), taucty(:), gwdcu(:,:), gwdcv(:,:)

      real(kind=kind_phys), intent(inout) :: dugwd(:), dvgwd(:),        &
     &  gu0(:,:), gv0(:,:), gt0(:,:)
      real(kind=kind_phys), intent(inout), optional :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:)
      integer, intent(in) :: index_of_process_nonorographic_gwd
      integer, intent(in) :: index_of_x_wind, index_of_y_wind

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k, idtend
      real(kind=kind_phys) :: eng0, eng1

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  --- ...  write out cloud top stress and wind tendencies

      if (lssav) then
        dugwd(:) = dugwd(:) + tauctx(:)*dtf
        dvgwd(:) = dvgwd(:) + taucty(:)*dtf
      endif   ! end if_lssav

      if (ldiag3d) then
         idtend = dtidx(index_of_x_wind,index_of_process_nonorographic_g&
     &                  wd)
         if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + gwdcu(:,:)  * dtf
         endif
         idtend = dtidx(index_of_y_wind,index_of_process_nonorographic_g&
     &                  wd)
         if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + gwdcv(:,:)  * dtf
         endif
      endif

!  --- ...  update the wind components with  gwdc tendencies

      do k = 1, levs
        do i = 1, im
          eng0     = 0.5*(gu0(i,k)*gu0(i,k) + gv0(i,k)*gv0(i,k))
          gu0(i,k) = gu0(i,k) + gwdcu(i,k) * dtp
          gv0(i,k) = gv0(i,k) + gwdcv(i,k) * dtp
          eng1     = 0.5*(gu0(i,k)*gu0(i,k) + gv0(i,k)*gv0(i,k))
          gt0(i,k) = gt0(i,k) + (eng0-eng1)/(dtp*con_cp)
        enddo
!         if (lprnt) write(7000,*)' gu0=',gu0(ipr,k),' gwdcu=',
!    &gwdcu(ipr,k), ' gv0=', gv0(ipr,k),' gwdcv=',gwdcv(ipr,k)
!    &,' k=',k
      enddo

      end subroutine gwdc_post_run

      end module gwdc_post
