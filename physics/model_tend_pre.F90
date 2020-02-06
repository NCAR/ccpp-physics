!>\file model_tend_pre.F90
!!  Calculates tendencies from all processes outside of CPPP

module model_tend_pre

contains

!> \section arg_table_model_tend_pre_init Argument Table
!!
subroutine model_tend_pre_init()
end subroutine model_tend_pre_init

!> \section arg_table_model_tend_pre_finalize Argument Table
!!
subroutine model_tend_pre_finalize()
end subroutine model_tend_pre_finalize

!> \section arg_table_model_tend_pre_run Argument Table
!! \htmlinclude model_tend_pre_run.html
!!
subroutine model_tend_pre_run(dtp, kdt,                     &
     tgrs,ugrs,vgrs,qvgrs, t_start,u_start,v_start,q_start, &
     dt3dt_model,du3dt_model,dv3dt_model,dq3dt_model,       &
     t_end,u_end,v_end,q_end,                               &
     im, levs, ntrac,                                       &
     lssav, ldiag3d, qdiag3d, errmsg,errflg)
  use machine,               only: kind_phys
  implicit none

  real(kind=kind_phys), dimension(:,:),       intent(in)  :: tgrs, ugrs, vgrs, qvgrs
  real(kind=kind_phys), dimension(:,:),       intent(out) :: t_start, u_start, v_start
  real(kind=kind_phys), dimension(:,:),       intent(out) :: q_start
  real(kind=kind_phys), dimension(:,:),       intent(out) :: t_end, u_end, v_end
  real(kind=kind_phys), dimension(:,:),       intent(out) :: q_end
  real(kind=kind_phys), dimension(:,:),       intent(inout) :: &
       dt3dt_model,du3dt_model,dv3dt_model,dq3dt_model
  
  integer, intent(in) :: im, levs, ntrac, kdt

  logical, intent(in) :: lssav, qdiag3d, ldiag3d

  real(kind=kind_phys) :: dtp
  
  character(len=*),     intent(out) :: errmsg
  integer,              intent(out) :: errflg

  integer :: i, k

  ! Initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  if(Lssav .and. ldiag3d) then
    do k=1,levs
      do i=1,im
        t_start(i,k) = tgrs(i,k)
        u_start(i,k) = ugrs(i,k)
        v_start(i,k) = vgrs(i,k)
        if(qdiag3d) then
          q_start(i,k) = qvgrs(i,k)
        endif
        if(t_start(i,k)>1e-3 .and. t_end(i,k)>1e-3) then
          dt3dt_model(i,k) = dt3dt_model(i,k) + (t_start(i,k)-t_end(i,k))
          du3dt_model(i,k) = du3dt_model(i,k) + (u_start(i,k)-u_end(i,k))
          dv3dt_model(i,k) = dv3dt_model(i,k) + (v_start(i,k)-v_end(i,k))
          if(qdiag3d) then
            dq3dt_model(i,k) = dq3dt_model(i,k) + (q_start(i,k)-q_end(i,k))
          endif
        endif
      enddo
    enddo
  endif
end subroutine model_tend_pre_run

end module model_tend_pre
