! #########################################################################################
!> \file mynnedmf_wrapper_post.f90
!!
! #########################################################################################
module mynnedmf_wrapper_post
  use machine,        only: kind_phys
  implicit none
contains

! #########################################################################################
!> \section arg_table_mynnedmf_wrapper_post_run Argument Table
!! \htmlinclude mynnedmf_wrapper_post_run.html
!!
! #########################################################################################
  subroutine mynnedmf_wrapper_post_run (tend_opt_pbl, im, levs, ntrac, &
       dtp, ten_t, ten_u, ten_v, ten_q, ten_t_pbl, ten_q_pbl, gt0, gu0, gv0, gq0, dtdt, dudt, dvdt, dqdt, &
       errmsg, errflg)
    
    ! Inputs
    integer, intent(in) :: tend_opt_pbl, im, levs, ntrac
    real(kind=kind_phys), intent(in) :: dtp
    real(kind=kind_phys), intent(in), dimension(:,:) :: ten_u, ten_v, ten_t
    real(kind=kind_phys), intent(out),dimension(:,:) :: ten_t_pbl,ten_q_pbl
    real(kind=kind_phys), intent(in), dimension(:,:,:) :: ten_q
    real(kind=kind_phys), intent(inout), dimension(:,:) :: gt0, gu0, gv0
    real(kind=kind_phys), intent(inout), dimension(:,:,:) :: gq0
    real(kind=kind_phys), intent(inout), dimension(:,:) :: dtdt, dudt, dvdt
    real(kind=kind_phys), intent(inout), dimension(:,:,:) :: dqdt

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg          ! CCPP Error message.
    integer,  intent(out) :: &
         errflg          ! CCPP Error flag.
    
    ! Locals
    integer :: i,k,n
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    case_pbl_ten: select case (tend_opt_pbl)
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
        exit case_pbl_ten
      case default
        errflg = 1
        errmsg = 'A tendency application control was outside of the acceptable range (1-4)'
        return
    end select case_pbl_ten

    ten_t_pbl(:,:)=0.
    ten_q_pbl(:,:)=0.
    
    !Output t and q tenedncies for PBL only to be used
    !as input in other schemes
    do k=1,levs
       do i=1,im
          ten_t_pbl(i,k)=ten_t(i,k)
          ten_q_pbl(i,k)=ten_q(i,k,1)
       end do
    end do
    
  end subroutine mynnedmf_wrapper_post_run

end module mynnedmf_wrapper_post
