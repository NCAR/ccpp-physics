!> \file GFS_gwd_generic_post.F90 
!! This file contains the CCPP-compliant orographic gravity wave drag post
!! interstitial codes.
module GFS_GWD_generic_post

contains

!> \section arg_table_GFS_GWD_generic_post_run Argument Table
!! \htmlinclude GFS_GWD_generic_post_run.html
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!>  @{
      subroutine GFS_GWD_generic_post_run(im, levs, ntrac, tend_opt_gwd, lssav, ldiag3d, dtf, dtp, dusfcg, dvsfcg, ten_t, ten_u, ten_v, ten_q, dudt, dvdt, dtdt, dqdt, &
      &  dugwd, dvgwd, gt0, gq0, gu0, gv0, flag_for_gwd_generic_tend, dtend, dtidx, index_of_temperature, index_of_x_wind,  &
      &  index_of_y_wind, index_of_process_orographic_gwd, errmsg, errflg)

      use machine, only : kind_phys
      implicit none
      
      integer, intent(in) :: im, levs, ntrac, tend_opt_gwd
      logical, intent(in) :: lssav, ldiag3d, flag_for_gwd_generic_tend
      
      real(kind=kind_phys), intent(in) :: dusfcg(:), dvsfcg(:)
      real(kind=kind_phys), intent(in) :: ten_t(:,:), ten_u(:,:), ten_v(:,:)
      real(kind=kind_phys), intent(in) :: ten_q(:,:,:)
      real(kind=kind_phys), intent(inout) :: dudt(:,:), dvdt(:,:), dtdt(:,:)
      real(kind=kind_phys), intent(inout) :: dqdt(:,:,:)
      real(kind=kind_phys), intent(inout) :: gt0(:,:), gv0(:,:), gu0(:,:)
      real(kind=kind_phys), intent(inout) :: gq0(:,:,:)
      real(kind=kind_phys), intent(in) :: dtf, dtp
      
      real(kind=kind_phys), intent(inout) :: dugwd(:), dvgwd(:)

      ! dtend only allocated only if ldiag3d is .true.
      real(kind=kind_phys), intent(inout), optional :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:), index_of_temperature,          &
     &  index_of_x_wind, index_of_y_wind, index_of_process_orographic_gwd
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i,k,n,idtend

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      case_GWD_ten: select case (tend_opt_gwd)
        case (1) !immediately apply tendencies
                  !Current state = current state + dt*current tendency
                  !Accumulated tendency unchanged
          do k=1,levs
            do i=1,im
              gt0(i,k) = gt0(i,k) + dtp*ten_t(i,k)
              gu0(i,k) = gu0(i,k) + dtp*ten_u(i,k)
              gv0(i,k) = gv0(i,k) + dtp*ten_v(i,k)
              do n = 1, ntrac
                gq0(i,k,n) = gq0(i,k,n) + dtp*ten_q(i,k,n)
              end do
            end do
          end do
        case (2) !add tendencies to sum
                  !Accumulated tendency = accumulated tendency + current tendency
                  !Current state unchanged
          do k=1,levs
            do i=1,im
              dtdt(i,k) = dtdt(i,k) + ten_t(i,k)
              dudt(i,k) = dudt(i,k) + ten_u(i,k)
              dvdt(i,k) = dvdt(i,k) + ten_v(i,k)
              do n = 1, ntrac
                dqdt(i,k,n) = dqdt(i,k,n) + ten_q(i,k,n)
              end do
            end do
          end do
        case (3) !add tendencies to sum and apply
                  !Current state = current state + dt*(accumulated tendency + current tendency)
                  !Accumulated tendency = 0
          do k=1,levs
            do i=1,im
              gt0(i,k) = gt0(i,k) + dtp*(dtdt(i,k) + ten_t(i,k))
              dtdt(i,k) = 0.0
              gu0(i,k) = gu0(i,k) + dtp*(dudt(i,k) + ten_u(i,k))
              dudt(i,k) = 0.0
              gv0(i,k) = gv0(i,k) + dtp*(dvdt(i,k) + ten_v(i,k))
              dvdt(i,k) = 0.0
              do n = 1, ntrac
                gq0(i,k,n) = gq0(i,k,n) + dtp*(dqdt(i,k,n) + ten_q(i,k,n))
                dqdt(i,k,n) = 0.0
              end do
            end do
          end do
        case (4) !Current state unchanged
                  !Accumulated tendency unchanged
                  !Current tendency unchanged (but will be overwritten during next primary scheme)
          exit case_GWD_ten
        case default
          errflg = 1
          errmsg = 'A tendency application control was outside of the acceptable range (1-4)'
          return
      end select case_GWD_ten
      
      if (lssav) then
        dugwd(:) = dugwd(:) + dusfcg(:)*dtf
        dvgwd(:) = dvgwd(:) + dvsfcg(:)*dtf

        if (ldiag3d .and. flag_for_gwd_generic_tend) then
          idtend = dtidx(index_of_temperature, index_of_process_orographic_gwd)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + ten_t*dtf
          endif

          idtend = dtidx(index_of_x_wind, index_of_process_orographic_gwd)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + ten_u*dtf
          endif

          idtend = dtidx(index_of_y_wind, index_of_process_orographic_gwd)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + ten_v*dtf
          endif
        endif
      endif

    end subroutine GFS_GWD_generic_post_run
!> @}

end module GFS_GWD_generic_post
