! #########################################################################################
!> \file GFS_photochemistry_post.f90
!!
! #########################################################################################
module GFS_photochemistry_post
  use machine,        only: kind_phys
  implicit none
contains

! #########################################################################################
!> \section arg_table_GFS_photochemistry_post_run Argument Table
!! \htmlinclude GFS_photochemistry_post_run.html
!!
! #########################################################################################
  subroutine GFS_photochemistry_post_run (tend_opt_photochem, im, levs, ntrac, &
       dtp, ten_t, ten_u, ten_v, ten_q, gt0, gu0, gv0, gq0, dtdt, dudt, dvdt, dqdt, &
       errmsg, errflg)
    
    ! Inputs
    integer, intent(in) :: tend_opt_photochem, im, levs, ntrac
    real(kind=kind_phys), intent(in) :: dtp
    real(kind=kind_phys), intent(in), dimension(:,:) :: ten_u, ten_v, ten_t
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
    
    case_photochemistry_ten: select case (tend_opt_photochem)
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
        exit case_photochemistry_ten
      case default
        errflg = 1
        errmsg = 'A tendency application control was outside of the acceptable range (1-4)'
        return
    end select case_photochemistry_ten

  end subroutine GFS_photochemistry_post_run

end module GFS_photochemistry_post
