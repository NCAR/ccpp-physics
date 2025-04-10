!>\file rrtmg_sw_post.F90

!> This module contains RRTMG-SW scheme post
      module rrtmg_sw_post
      contains

!>\defgroup rrtmg_sw_post GFS RRTMG-SW scheme post
!! This module  saves two spectral bands' surface downward and upward fluxes for
!!  output.
!> @{
!> \section arg_table_rrtmg_sw_post_run Argument Table
!! \htmlinclude rrtmg_sw_post_run.html
!!
      subroutine rrtmg_sw_post_run (im, levr, levs, ntrac, ltp, nday, lm, kd, tend_opt_swrad, lsswr,  &
                 swhtr, delt, sfcalb1, sfcalb2, sfcalb3, sfcalb4, htswc, htsw0,      &
                 nirbmdi, nirdfdi, visbmdi, visdfdi, nirbmui, nirdfui, visbmui,&
                 visdfui, sfcdsw, sfcnsw, htrsw, swhc, scmpsw, sfcfsw, topfsw, &
                 sfcdswc, gu0, gv0, gt0, gq0, dudt, dvdt, dtdt, dqdt, ten_t, ten_u, ten_v, ten_q, errmsg, errflg)

      use machine,                   only: kind_phys
      use module_radsw_parameters,   only: topfsw_type, sfcfsw_type,   &
                                           cmpfsw_type

      implicit none

      integer,                              intent(in)    :: im, levr, levs,   &
                                                             ltp, nday, lm, kd,&
                                                             ntrac, tend_opt_swrad
      logical,                              intent(in)    :: lsswr, swhtr 
      real(kind=kind_phys),                 intent(in)    :: delt
      real(kind=kind_phys), dimension(:),   intent(in)    :: sfcalb1, sfcalb2, &
                                                             sfcalb3, sfcalb4
      real(kind=kind_phys), dimension(:,:), intent(in)    :: htswc, htsw0
      
      real(kind=kind_phys), dimension(:),   intent(inout) :: nirbmdi, nirdfdi, &
                                                             visbmdi, visdfdi, &
                                                             nirbmui, nirdfui, &
                                                             visbmui, visdfui, &
                                                             sfcdsw,  sfcnsw,  &
                                                             sfcdswc
      real(kind=kind_phys), dimension(:,:), intent(inout) :: htrsw, swhc

      type(cmpfsw_type), dimension(:),     intent(inout) :: scmpsw
      type(sfcfsw_type), dimension(:),     intent(inout) :: sfcfsw
      type(topfsw_type), dimension(:),     intent(inout) :: topfsw
      
      real(kind=kind_phys), dimension(:,:),   intent(inout) :: gu0, gv0, gt0
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: gq0
      real(kind=kind_phys), dimension(:,:),   intent(inout) :: dudt, dvdt, dtdt
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: dqdt
      real(kind=kind_phys), dimension(:,:),   intent(out)   :: ten_t, ten_u, ten_v
      real(kind=kind_phys), dimension(:,:,:), intent(out)   :: ten_q

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer :: i, k1, k, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      ten_t = 0.0
      ten_u = 0.0
      ten_v = 0.0
      ten_q = 0.0
      
      if (lsswr) then
        if (nday > 0) then
          do k = 1, LM
            k1 = k + kd
            htrsw(1:im,k) = htswc(1:im,k1)
          enddo
          ! We are assuming that radiative tendencies are from bottom to top 
          ! --- repopulate the points above levr i.e. LM
          if (lm < levs) then
            do k = lm+1, levs
              htrsw (1:im,k) = htrsw (1:im,LM)
            enddo
          endif

          if (swhtr) then
            do k = 1, lm
               k1 = k + kd
               swhc(1:im,k) = htsw0(1:im,k1)
             enddo
             ! --- repopulate the points above levr i.e. LM
             if (lm < levs) then
               do k = lm+1, levs
                 swhc(1:im,k) = swhc(1:im,LM)
               enddo
             endif
          endif

!  --- surface down and up spectral component fluxes
!   Save two spectral bands' surface downward and upward fluxes for
!    output.

          do i=1,im
            nirbmdi(i) = scmpsw(i)%nirbm
            nirdfdi(i) = scmpsw(i)%nirdf
            visbmdi(i) = scmpsw(i)%visbm
            visdfdi(i) = scmpsw(i)%visdf

            nirbmui(i) = scmpsw(i)%nirbm * sfcalb1(i)
            nirdfui(i) = scmpsw(i)%nirdf * sfcalb2(i)
            visbmui(i) = scmpsw(i)%visbm * sfcalb3(i)
            visdfui(i) = scmpsw(i)%visdf * sfcalb4(i)
          enddo

        else                   ! if_nday_block

          htrsw(:,:) = 0.0

          sfcfsw = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          topfsw = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )

          do i=1,im
            nirbmdi(i) = 0.0
            nirdfdi(i) = 0.0
            visbmdi(i) = 0.0
            visdfdi(i) = 0.0

            nirbmui(i) = 0.0
            nirdfui(i) = 0.0
            visbmui(i) = 0.0
            visdfui(i) = 0.0
          enddo

          if (swhtr) then
            swhc(:,:) = 0
          endif

        endif                  ! end_if_nday

! --- radiation fluxes for other physics processes
        do i=1,im
          sfcnsw(i) = sfcfsw(i)%dnfxc - sfcfsw(i)%upfxc
          sfcdsw(i) = sfcfsw(i)%dnfxc
          sfcdswc(i)= sfcfsw(i)%dnfx0
        enddo

      endif                                ! end_if_lsswr
      
      !if lsswr is true, htrsw is recalculated, otherwise, it has the value from the previous radiation timestep
      !So, let's set the accumulated tendency at this point
      !Then, include the tendency application block here with control = 2 or 4?
      ! - if 4, htrsw (on radiation timesteps) gets passed to dcyc2t3 where it is used to calculate t_tend and applied (with control = 2)
      ! - can also use 1, 2, or 3 with radsw then
      
      ten_t = htrsw
      
      !This may belong in a separate GFS_radsw_post routine rather than here, although it would need to be created
      case_SWRAD_ten: select case (tend_opt_swrad)
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
          exit case_SWRAD_ten
        case default
          errflg = 1
          errmsg = 'A tendency application control was outside of the acceptable range (1-4)'
          return
      end select case_SWRAD_ten
      
      
      end subroutine rrtmg_sw_post_run
!> @}
      end module rrtmg_sw_post
