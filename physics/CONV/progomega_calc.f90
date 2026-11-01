!>\file progomega_calc.f90

!> This module contains the subroutine that calculates the prognostic
!! updraft velocity that is used for closure computations in
!! saSAS deep and shallow convection
!! as described in Bengtsson et al. 2026 \cite Bengtsson_2026.


module progomega

  use mo_conv_kind, only : conv_wp

  implicit none

  public progomega_calc

contains

!> This subroutine computes a prognostic updraft velocity
!! This file contains the subroutine that calculates the prognostic
!! updraft vertical velocity that is used for closure computations in
!! saSAS and C3 deep and shallow convection.
!!\section gen_progomega progomega_calc General Algorithm
    subroutine progomega_calc(first_time_step,flag_restart,im,km,kbcon1,ktcon,omegain,delt,del, &
       zi,cnvflg,omegaout,grav,buo,drag,wush,lbb1,lbb2,lbb3,conv_type)

    implicit none

    integer, intent(in) :: im,km
    integer, intent(in) :: kbcon1(im),ktcon(im)
    integer, intent(in) :: conv_type

    real(kind=conv_wp), intent(in) :: delt,grav
    real(kind=conv_wp), intent(in) :: lbb1,lbb2,lbb3
    real(kind=conv_wp), intent(in) :: omegain(im,km)
    real(kind=conv_wp), intent(in) :: del(im,km),zi(im,km)
    real(kind=conv_wp), intent(in) :: drag(im,km)
    real(kind=conv_wp), intent(in) :: buo(im,km)
    real(kind=conv_wp), intent(in) :: wush(im,km)

    real(kind=conv_wp), intent(inout) :: omegaout(im,km)
    
    logical, intent(in) :: cnvflg(im)
    logical, intent(in) :: first_time_step
    logical, intent(in) :: flag_restart

    !--------------------------------------------------------------------
    ! Local arrays
    !
    ! omega     = state at beginning of current internal substep
    ! omega_new = state at end of current internal substep
    !--------------------------------------------------------------------

    real(kind=conv_wp) :: omega(im,km)
    real(kind=conv_wp) :: omega_new(im,km)

    real(kind=conv_wp) :: termA(im,km)
    real(kind=conv_wp) :: termB(im,km)
    real(kind=conv_wp) :: termC(im,km)
    real(kind=conv_wp) :: memory_term,buoy_term,adv_term

    ! Height-dependent effective perturbation-pressure coefficients:
    ! bb2_prof = 1 - Cb for buoyancy forcing
    ! bb4_prof = 1 - Cd for vertical momentum advection
    real(kind=conv_wp) :: bb2_prof(im,km)
    real(kind=conv_wp) :: bb4_prof(im,km)
    real(kind=conv_wp) :: zkm,zbase,cb_loc,cd_loc

    integer, parameter :: nprof_deep = 7
    real(kind=conv_wp), parameter :: z_deep(nprof_deep) =       &
         (/ 0.0_conv_wp, 0.3_conv_wp, 0.8_conv_wp,         &
            1.6_conv_wp, 4.0_conv_wp, 12.8_conv_wp,        &
            16.0_conv_wp /)
    real(kind=conv_wp), parameter :: cb_deep(nprof_deep) =      &
         (/ 0.00_conv_wp, 0.30_conv_wp, 0.60_conv_wp,      &
            0.70_conv_wp, 0.55_conv_wp, 0.55_conv_wp,      &
            0.25_conv_wp /)
    real(kind=conv_wp), parameter :: cd_deep(nprof_deep) =      &
         (/ 0.00_conv_wp, 0.40_conv_wp, 0.38_conv_wp,      &
            0.35_conv_wp, 0.25_conv_wp, 0.20_conv_wp,      &
            0.05_conv_wp /)

    integer, parameter :: nprof_shal = 7
    real(kind=conv_wp), parameter :: z_shal(nprof_shal) =       &
         (/ 0.0_conv_wp, 0.2_conv_wp, 0.8_conv_wp,         &
            1.6_conv_wp, 2.4_conv_wp, 3.2_conv_wp,         &
            4.0_conv_wp /)
    real(kind=conv_wp), parameter :: cb_shal(nprof_shal) =      &
         (/ 0.00_conv_wp, 0.45_conv_wp, 0.70_conv_wp,      &
            0.85_conv_wp, 0.75_conv_wp, 0.20_conv_wp,      &
            0.00_conv_wp /)
    real(kind=conv_wp), parameter :: cd_shal(nprof_shal) =      &
         (/ 0.00_conv_wp, 0.15_conv_wp, 0.35_conv_wp,      &
            0.40_conv_wp, 0.25_conv_wp, 0.05_conv_wp,      &
            0.00_conv_wp /)
    !--------------------------------------------------------------------
    ! Scalars
    !--------------------------------------------------------------------

    real(kind=conv_wp) :: dp,dz,pi_conv,discr
    real(kind=conv_wp) :: rhs_exp

    real(kind=conv_wp) :: omega_eps
    real(kind=conv_wp) :: a_eps,b_eps,disc_eps

    real(kind=conv_wp) :: cfl_target
    real(kind=conv_wp) :: dt_sub
    real(kind=conv_wp) :: dt_remaining
    real(kind=conv_wp) :: dt_cfl
    real(kind=conv_wp) :: time_done

    integer :: i,k

    logical :: active_point_found

    !--------------------------------------------------------------------
    ! Numerical parameters
    !--------------------------------------------------------------------

    ! Target CFL for explicit pressure-coordinate vertical advection.
    cfl_target = 0.8_conv_wp

    omega_eps = 1.0e-5_conv_wp
    a_eps     = 1.0e-12_conv_wp
    b_eps     = 1.0e-12_conv_wp
    disc_eps  = 1.0e-12_conv_wp

    !--------------------------------------------------------------------
    ! Construct height-dependent perturbation-pressure coefficients.
    ! Height is measured above cloud base, consistent with Bengtsson et al.
    ! 2026
    !--------------------------------------------------------------------

    bb2_prof(:,:) = lbb2
    bb4_prof(:,:) = 1.0_conv_wp

    do k = 1,km
       do i = 1,im

          if (cnvflg(i)) then
             if (k >= kbcon1(i) .and. k < ktcon(i)) then

                ! Cloud-base layer midpoint height [m]
                zbase = 0.5_conv_wp * &
                     (zi(i,kbcon1(i)) + zi(i,kbcon1(i)+1))

                ! Layer height above cloud base [km]
                zkm = (0.5_conv_wp * (zi(i,k) + zi(i,k+1)) - zbase) &
                     * 0.001_conv_wp
                zkm = max(zkm,0.0_conv_wp)

                select case (conv_type)
                case (1)
                   ! Deep convection
                   cb_loc = interp_profile(zkm,z_deep,cb_deep,nprof_deep)
                   cd_loc = interp_profile(zkm,z_deep,cd_deep,nprof_deep)
                case (2)
                   ! Shallow convection
                   cb_loc = interp_profile(zkm,z_shal,cb_shal,nprof_shal)
                   cd_loc = interp_profile(zkm,z_shal,cd_shal,nprof_shal)
                case default
                   ! Preserve the previous UFS formulation if conv_type
                   ! is not recognized.
                   cb_loc = 1.0_conv_wp - lbb2
                   cd_loc = 0.0_conv_wp
                end select

                bb2_prof(i,k) = 1.0_conv_wp - cb_loc
                bb4_prof(i,k) = 1.0_conv_wp - cd_loc

             endif
          endif

       enddo
    enddo

    !--------------------------------------------------------------------
    ! Initialize from incoming prognostic tracer
    !--------------------------------------------------------------------

    do k = 1,km
       do i = 1,im

          termA(i,k) = 0.0_conv_wp
          termB(i,k) = 0.0_conv_wp
          termC(i,k) = 0.0_conv_wp

          omega(i,k)     = omegain(i,k)
          omega_new(i,k) = omegain(i,k)

       enddo
    enddo

    !--------------------------------------------------------------------
    ! Remove numerically negligible values for active convection
    !--------------------------------------------------------------------

    do k = 1,km
       do i = 1,im
          if (cnvflg(i)) then
             if (abs(omega(i,k)) < omega_eps) then
                omega(i,k)     = 0.0_conv_wp
                omega_new(i,k) = 0.0_conv_wp
             endif
          endif
       enddo
    enddo

    !--------------------------------------------------------------------
    !-------------------------------------------------------------------
    ! Cold-start initialization
    !--------------------------------------------------------------------

    if (first_time_step .and. .not. flag_restart) then

       do k = 1,km
          do i = 1,im

             if (cnvflg(i)) then

                if (k >= kbcon1(i) .and. k < ktcon(i)) then

                   omega(i,k)     = -1.2_conv_wp
                   omega_new(i,k) = -1.2_conv_wp
                   omegaout(i,k)  = -1.2_conv_wp

                endif

             endif

          enddo
       enddo

    endif

    !--------------------------------------------------------------------
    ! Adaptive CFL-controlled subcycling
    !--------------------------------------------------------------------

    time_done = 0.0_conv_wp

    ! Save the state entering the prognostic integration so that the
    ! host-timestep local tendency can be returned after all substeps.

    do while (time_done < delt)

       dt_remaining = delt - time_done

       !---------------------------------------------------------------
       ! Determine CFL-limited timestep from current omega profile.
       !
       ! The explicit vertical-advection term is zero at cloud base,
       ! so only levels above kbcon1 are included in the CFL estimate.
       !---------------------------------------------------------------

       dt_cfl = dt_remaining
       active_point_found = .false.

       do k = 2,km
          do i = 1,im

             if (cnvflg(i)) then
                if (k > kbcon1(i) .and. k < ktcon(i)) then
                   dp = 1000.0_conv_wp * del(i,k)
                   if (dp > 0.0_conv_wp) then
                      active_point_found = .true.
                      if (abs(omega(i,k)) > omega_eps) then
                         dt_cfl = min(dt_cfl,                       &
                              cfl_target * dp / abs(omega(i,k)))
                      endif
                   endif
                endif
             endif
          enddo
       enddo

       ! If there are no active plume points above cloud base,
       ! no prognostic subcycling is needed.
       if (.not. active_point_found) exit

       dt_sub = min(dt_remaining,dt_cfl)

       ! Numerical protection against pathological tiny timesteps.
       if (dt_sub <= 1.0e-8_conv_wp) then
          dt_sub = dt_remaining
       endif

       !---------------------------------------------------------------
       ! Start new substep from previous-substep profile.
       !
       ! This ensures every vertical level uses the same
       ! profile for the explicit vertical-advection term.
       !---------------------------------------------------------------

       omega_new(:,:) = omega(:,:)

       !---------------------------------------------------------------
       ! Solve prognostic momentum equation following
       ! Bengtsson et al. 2026
       !---------------------------------------------------------------

       do k = 2,km
          do i = 1,im

             if (cnvflg(i)) then
                if (k >= kbcon1(i) .and. k < ktcon(i)) then
                   ! Cloud-base boundary condition.
                   omega_new(i,kbcon1(i)) = 0.0_conv_wp

                   dp = 1000.0_conv_wp * del(i,k)
                   dz = zi(i,k+1) - zi(i,k)

                   if (dp <= 0.0_conv_wp) cycle
                   if (abs(dz) <= 1.0e-12_conv_wp) cycle

                   pi_conv = dp/dz

                   !----------------------------------------------------
                   ! Explicit contributions
                   !----------------------------------------------------

                   memory_term = omega(i,k)

                   buoy_term = -0.5_conv_wp * dt_sub * lbb2      &
                        * buo(i,k) * pi_conv

                   if (k == kbcon1(i)) then
                      adv_term = 0.0_conv_wp
                   else
                      adv_term = -dt_sub * omega(i,k)             &
                            * (omega(i,k-1)-omega(i,k)) / dp

                   endif
                   rhs_exp = memory_term + buoy_term + adv_term

                   !----------------------------------------------------
                   ! Quadratic coefficients
                   !
                   ! A * omega_new**2 + B * omega_new + C = 0
                   !----------------------------------------------------

                   termA(i,k) = -0.5_conv_wp * dt_sub * lbb1     &
                        * drag(i,k) / pi_conv

                   termB(i,k) = 1.0_conv_wp                       &
                        + 0.5_conv_wp * dt_sub * wush(i,k)

                   termC(i,k) = -rhs_exp

                   !----------------------------------------------------
                   ! Robust quadratic / linear solution
                   !----------------------------------------------------

                   if (abs(termA(i,k)) < a_eps) then

                      if (abs(termB(i,k)) > b_eps) then
                         omega_new(i,k) = -termC(i,k) / termB(i,k)
                      else
                         ! Degenerate case: retain previous state.
                         omega_new(i,k) = omega(i,k)
                      endif

                   else

                      discr = termB(i,k)**2                         &
                            - 4.0_conv_wp * termA(i,k) * termC(i,k)

                      if (discr >= -disc_eps) then

                         discr = max(discr,0.0_conv_wp)

                         omega_new(i,k) =                           &
                              (-termB(i,k) + sqrt(discr))          &
                              / (2.0_conv_wp * termA(i,k))

                      else

                         ! No real solution: retain previous state.
                         omega_new(i,k) = omega(i,k)

                      endif

                   endif

                   !----------------------------------------------------
                   ! Physical bounds for active updrafts
                   !----------------------------------------------------

                   omega_new(i,k) = max(                           &
                        min(omega_new(i,k),0.0_conv_wp),           &
                        -80.0_conv_wp)


                endif

             endif

          enddo
       enddo

       !---------------------------------------------------------------
       ! Advance entire vertical profile simultaneously
       !---------------------------------------------------------------

       omega(:,:) = omega_new(:,:)

       time_done = time_done + dt_sub

       ! Prevent tiny floating-point remainder from generating an
       ! unnecessary additional substep.
       if (delt-time_done < 1.0e-8_conv_wp*delt) then
          time_done = delt
       endif

    enddo

    !--------------------------------------------------------------------
    ! Return final prognostic state
    !--------------------------------------------------------------------

    do k = 1,km
       do i = 1,im

          if (cnvflg(i)) then

             omegaout(i,k) = omega(i,k)

             if (abs(omegaout(i,k)) < omega_eps) then
                omegaout(i,k) = 0.0_conv_wp
             endif

          endif

       enddo
    enddo

  end subroutine progomega_calc

  !=====================================================================
  ! Piecewise-linear interpolation of perturbation-pressure profiles.
  ! Values outside the specified height range are held at the nearest
  ! endpoint.
  !=====================================================================

  pure function interp_profile(z,zpts,cpts,npts) result(c)

    implicit none

    integer, intent(in) :: npts
    real(kind=conv_wp), intent(in) :: z
    real(kind=conv_wp), intent(in) :: zpts(npts)
    real(kind=conv_wp), intent(in) :: cpts(npts)

    real(kind=conv_wp) :: c
    real(kind=conv_wp) :: weight
    integer :: n

    if (z <= zpts(1)) then
       c = cpts(1)
       return
    endif

    if (z >= zpts(npts)) then
       c = cpts(npts)
       return
    endif

    do n = 1,npts-1
       if (z >= zpts(n) .and. z <= zpts(n+1)) then
          weight = (z-zpts(n)) / (zpts(n+1)-zpts(n))
          c = cpts(n) + weight * (cpts(n+1)-cpts(n))
          return
       endif
    enddo

    c = cpts(npts)

  end function interp_profile

end module progomega
