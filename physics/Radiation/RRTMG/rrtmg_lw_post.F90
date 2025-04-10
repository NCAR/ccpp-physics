!>\file rrtmg_lw_post.F90

!> This module contains code executed after RRTMG-LW scheme
      module rrtmg_lw_post 
      contains

!>\defgroup rrtmg_lw_post GFS RRTMG scheme post
!! This module saves RRTMG-LW fluxes results.
!> @{

!> \section arg_table_rrtmg_lw_post_run Argument Table
!! \htmlinclude rrtmg_lw_post_run.html
!!
      subroutine rrtmg_lw_post_run (im, levs, ntrac, ltp, lm, kd, tend_opt_lwrad, lslwr, lwhtr, delt,      &
                 tsfa, htlwc, htlw0, sfcflw, tsflw, sfcdlw, htrlw, lwhc, gu0, gv0, gt0, gq0, dudt, dvdt, dtdt, dqdt, ten_t, ten_u, ten_v, ten_q,      &
                 errmsg, errflg)
    
      use machine,                   only: kind_phys
      use module_radlw_parameters,   only: sfcflw_type
      
      implicit none
      
      integer,                                     intent(in) :: im, levs, ntrac, ltp, lm, kd, tend_opt_lwrad
      logical,                                     intent(in) :: lslwr, lwhtr
      real(kind=kind_phys),                        intent(in) :: delt
      real(kind=kind_phys), dimension(im),         intent(in) ::  tsfa
      real(kind=kind_phys), dimension(im, LM+LTP), intent(in) ::  htlwc
      real(kind=kind_phys), dimension(im, LM+LTP), intent(in) ::  htlw0
      
      type(sfcflw_type), dimension(im),            intent(in) :: sfcflw
      
      real(kind=kind_phys), dimension(im),         intent(inout) ::  tsflw, sfcdlw
      real(kind=kind_phys), dimension(im, levs),   intent(inout) ::  htrlw, lwhc
      real(kind=kind_phys), dimension(:,:),        intent(inout) :: gu0, gv0, gt0
      real(kind=kind_phys), dimension(:,:,:),      intent(inout) :: gq0
      real(kind=kind_phys), dimension(:,:),        intent(inout) :: dudt, dvdt, dtdt
      real(kind=kind_phys), dimension(:,:,:),      intent(inout) :: dqdt
      real(kind=kind_phys), dimension(:,:),        intent(out) :: ten_t, ten_u, ten_v
      real(kind=kind_phys), dimension(:,:,:),      intent(out) :: ten_q
      character(len=*),                            intent(out) :: errmsg
      integer,                                     intent(out) :: errflg
      
      ! local variables
      integer :: i, k1, k, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      ten_t = 0.0
      ten_u = 0.0
      ten_v = 0.0
      ten_q = 0.0
      
      if (lslwr) then
! Save calculation results
! Save surface air temp for diurnal adjustment at model t-steps

        tsflw (:) = tsfa(:)

        do k = 1, LM
          k1 = k + kd
            htrlw(1:im,k) = htlwc(1:im,k1)
        enddo
        ! --- repopulate the points above levr
        if (lm < levs) then
          do k = lm+1, levs
            htrlw (1:im,k) = htrlw (1:im,LM)
          enddo
        endif

        if (lwhtr) then
          do k = 1, lm
            k1 = k + kd
            lwhc(1:im,k) = htlw0(1:im,k1)
          enddo
          ! --- repopulate the points above levr
          if (lm < levs) then
            do k = lm+1, levs
              lwhc(1:im,k) = lwhc(1:im,LM)
            enddo
          endif
        endif

! --- radiation fluxes for other physics processes
        sfcdlw(:) = sfcflw(:)%dnfxc

      endif                                ! end_if_lslwr
      
      ten_t = htrlw
      
      !This may belong in a separate GFS_radsw_post routine rather than here, although it would need to be created
      case_LWRAD_ten: select case (tend_opt_lwrad)
        case (1) !immediately apply tendencies
                  !Current state = current state + dt*current tendency
                  !Accumulated tendency unchanged
          do k=1,levs
            do i=1,im
              gt0(i,k) = gt0(i,k) + delt*ten_t(i,k)
              gu0(i,k) = gu0(i,k) + delt*ten_u(i,k)
              gv0(i,k) = gv0(i,k) + delt*ten_v(i,k)
              do n = 1, ntrac
                gq0(i,k,n) = gq0(i,k,n) + delt*ten_q(i,k,n)
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
              gt0(i,k) = gt0(i,k) + delt*(dtdt(i,k) + ten_t(i,k))
              dtdt(i,k) = 0.0
              gu0(i,k) = gu0(i,k) + delt*(dudt(i,k) + ten_u(i,k))
              dudt(i,k) = 0.0
              gv0(i,k) = gv0(i,k) + delt*(dvdt(i,k) + ten_v(i,k))
              dvdt(i,k) = 0.0
              do n = 1, ntrac
                gq0(i,k,n) = gq0(i,k,n) + delt*(dqdt(i,k,n) + ten_q(i,k,n))
                dqdt(i,k,n) = 0.0
              end do
            end do
          end do
        case (4) !Current state unchanged
                  !Accumulated tendency unchanged
                  !Current tendency unchanged (but will be overwritten during next primary scheme)
          exit case_LWRAD_ten
        case default
          errflg = 1
          errmsg = 'A tendency application control was outside of the acceptable range (1-4)'
          return
      end select case_LWRAD_ten
      
      end subroutine rrtmg_lw_post_run

!> @}
      end module rrtmg_lw_post
