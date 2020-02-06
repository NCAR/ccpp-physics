!>\file model_tend_post.F90
!! Calculates tendencies from all processes outside of CPPP

module model_tend_post

contains

  subroutine model_tend_post_init()
  end subroutine model_tend_post_init

  subroutine model_tend_post_finalize()
  end subroutine model_tend_post_finalize

  !> \section arg_table_model_tend_post_run Argument Table
  !! \htmlinclude model_tend_post_run.html
  !!
  subroutine model_tend_post_run(kdt,                          &
       gt0,gu0,gv0, gq0_water_vapor,                           &
       t_start,u_start,v_start,q_start,  &
       t_end, u_end, v_end, q_end,                             &
       dt3dt_ccpp, du3dt_ccpp, dv3dt_ccpp, dq3dt_ccpp,         &
!       dt3dt_total, du3dt_total, dv3dt_total, dq3dt_total,     &
       im, levs, ntrac, index_for_water_vapor,                 &
       lssav, ldiag3d, qdiag3d, errmsg,errflg)
    use machine,               only: kind_phys
    implicit none

    real(kind=kind_phys), dimension(:,:),       intent(in) :: gt0, gu0, gv0, gq0_water_vapor
    real(kind=kind_phys), dimension(:,:),       intent(in) :: t_start, u_start, v_start
    real(kind=kind_phys), dimension(:,:),       intent(in) :: q_start
    real(kind=kind_phys), dimension(:,:),       intent(inout) :: t_end, u_end, v_end
    real(kind=kind_phys), dimension(:,:),       intent(inout) :: q_end
    real(kind=kind_phys), dimension(:,:),       intent(inout) :: du3dt_ccpp, dv3dt_ccpp
    real(kind=kind_phys), dimension(:,:),       intent(inout) :: dt3dt_ccpp, dq3dt_ccpp
    ! real(kind=kind_phys), dimension(:,:),       intent(inout) :: du3dt_total, dv3dt_total
    ! real(kind=kind_phys), dimension(:,:),       intent(inout) :: dt3dt_total, dq3dt_total

    integer, intent(in) :: im, levs, ntrac, kdt
    integer, intent(in) :: index_for_water_vapor

    logical, intent(in) :: lssav, qdiag3d, ldiag3d

    character(len=*),     intent(out) :: errmsg
    integer,              intent(out) :: errflg

    real(kind=kind_phys) :: dt
    integer :: i,k

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    diag_enabled: if(lssav .and. ldiag3d) then
      if(any(gt0(1:im,1:levs)<1e-3)) then
        print *,'error: temperatures less than 1e-3'
      endif
      if(all(abs(gu0(1:im,1:levs))<1e-3)) then
        print *,'error: all u wind is near zero'
      endif
      if(all(abs(gv0(1:im,1:levs))<1e-3)) then
        print *,'error: all v wind is near zero'
      endif

      if(any(t_start(1:im,1:levs)<1e-3)) then
        print *,'error: start temperatures less than 1e-3'
      endif
      if(all(abs(u_start(1:im,1:levs))<1e-3)) then
        print *,'error: all start u wind is near zero'
      endif
      if(all(abs(v_start(1:im,1:levs))<1e-3)) then
        print *,'error: all start v wind is near zero'
      endif

      do k=1,levs
        do i=1,im
          ! if(t_end(i,k)>1e-3 .and. gt0(i,k)>1e-3) then
          !     dt3dt_total(i,k) = dt3dt_total(i,k) + gt0(i,k)-t_end(i,k)
          !     du3dt_total(i,k) = du3dt_total(i,k) + gu0(i,k)-u_end(i,k)
          !     dv3dt_total(i,k) = dv3dt_total(i,k) + gv0(i,k)-v_end(i,k)
          !     if(qdiag3d) then
          !       dq3dt_total(i,k) = dq3dt_total(i,k) + gq0_water_vapor(i,k)-q_end(i,k)
          !     endif
          ! endif
          t_end(i,k) = gt0(i,k)
          u_end(i,k) = gu0(i,k)
          v_end(i,k) = gv0(i,k)
          if(qdiag3d) then
            q_end(i,k) = gq0_water_vapor(i,k)
          endif
          if(t_end(i,k)>1e-3 .and. t_start(i,k)>1e-3) then
            dt3dt_ccpp(i,k) = dt3dt_ccpp(i,k) + t_end(i,k)-t_start(i,k)
            du3dt_ccpp(i,k) = du3dt_ccpp(i,k) + u_end(i,k)-u_start(i,k)
            dv3dt_ccpp(i,k) = dv3dt_ccpp(i,k) + v_end(i,k)-v_start(i,k)
            if(qdiag3d) then
              dq3dt_ccpp(i,k) = dq3dt_ccpp(i,k) + q_end(i,k)-q_start(i,k)
            endif
          endif
        enddo
      enddo

    endif diag_enabled

  end subroutine model_tend_post_run

end module model_tend_post
