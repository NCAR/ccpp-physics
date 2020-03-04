!>\file total_tend.F90
!!  Calculates tendencies from all processes outside of CPPP

module total_tend

contains

!> \section arg_table_total_tend_init Argument Table
!!
subroutine total_tend_init()
end subroutine total_tend_init

!> \section arg_table_total_tend_finalize Argument Table
!!
subroutine total_tend_finalize()
end subroutine total_tend_finalize

!> \section arg_table_total_tend_run Argument Table
!! \htmlinclude total_tend_run.html
!!
subroutine total_tend_run(dtp, kdt,                     &
     tgrs,ugrs,vgrs,qvgrs, t_start,u_start,v_start,q_start, &
     dt3dt_total,du3dt_total,dv3dt_total,dq3dt_total,       &
     gt0,gu0,gv0, gq0_water_vapor,                          &
     im, levs, ntrac,                                       &
     lssav, ldiag3d, qdiag3d, errmsg,errflg)
  use machine,               only: kind_phys
  implicit none

  real(kind=kind_phys), dimension(:,:),       intent(in)  :: gt0, gu0, gv0, gq0_water_vapor
  real(kind=kind_phys), dimension(:,:),       intent(in)  :: tgrs, ugrs, vgrs, qvgrs
  real(kind=kind_phys), dimension(:,:),       intent(out) :: t_start, u_start, v_start
  real(kind=kind_phys), dimension(:,:),       intent(out) :: q_start
  real(kind=kind_phys), dimension(:,:),       intent(inout) :: &
       dt3dt_total,du3dt_total,dv3dt_total,dq3dt_total
  
  integer, intent(in) :: im, levs, ntrac, kdt

  logical, intent(in) :: lssav, qdiag3d, ldiag3d

  real(kind=kind_phys) :: dtp
  
  character(len=*),     intent(out) :: errmsg
  integer,              intent(out) :: errflg

  integer :: i, k, good

  ! Initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  good=0

  if(Lssav .and. ldiag3d) then
    do k=1,levs
      do i=1,im
        if(t_start(i,k)>1e-3 .and. gt0(i,k)>1e-3) then
          good=good+1
          dt3dt_total(i,k) = dt3dt_total(i,k) + (gt0(i,k)-t_start(i,k))
          du3dt_total(i,k) = du3dt_total(i,k) + (gu0(i,k)-u_start(i,k))
          dv3dt_total(i,k) = dv3dt_total(i,k) + (gv0(i,k)-v_start(i,k))
          if(qdiag3d) then
            dq3dt_total(i,k) = dq3dt_total(i,k) + (gq0_water_vapor(i,k)-q_start(i,k))
          endif
        endif
        t_start(i,k)=gt0(i,k)
        u_start(i,k)=gu0(i,k)
        v_start(i,k)=gv0(i,k)
        if(qdiag3d) then
          q_start(i,k)=gq0_water_vapor(i,k)
        endif
        ! Alternative is to use the state in:
        ! if(t_start(i,k)>1e-3 .and. tgrs(i,k)>1e-3) then
        !   good=good+1
        !   dt3dt_total(i,k) = dt3dt_total(i,k) + (tgrs(i,k)-t_start(i,k))
        !   du3dt_total(i,k) = du3dt_total(i,k) + (ugrs(i,k)-u_start(i,k))
        !   dv3dt_total(i,k) = dv3dt_total(i,k) + (vgrs(i,k)-v_start(i,k))
        !   if(qdiag3d) then
        !     dq3dt_total(i,k) = dq3dt_total(i,k) + (qvgrs(i,k)-q_start(i,k))
        !   endif
        ! endif
        ! t_start(i,k)=tgrs(i,k)
        ! u_start(i,k)=ugrs(i,k)
        ! v_start(i,k)=vgrs(i,k)
        ! if(qdiag3d) then
        !   q_start(i,k)=qvgrs(i,k)
        ! endif
      enddo
    enddo
  endif
end subroutine total_tend_run

end module total_tend
