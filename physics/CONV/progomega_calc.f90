!>\file progomega_calc.f90

!> This module contains the subroutine that calculates the prognostic
!! updraft velocity that is used for closure computations in
!! saSAS deep and shallow convection
!! as described in Bengtsson et al. 2026 \cite Bengtsson_2026.


module progomega

  implicit none

  public progomega_calc

contains
  
!> This subroutine computes a prognostic updraft velocity
!! This file contains the subroutine that calculates the prognostic
!! updraft vertical velocity that is used for closure computations in
!! saSAS and C3 deep and shallow convection.
!!\section gen_progomega progomega_calc General Algorithm
  subroutine progomega_calc(first_time_step,flag_restart,im,km,kbcon1,ktcon,omegain,delt,del, &
       zi,cnvflg,omegaout,grav,buo,drag,wush,lbb1,lbb2,lbb3,dt_decay)

    use machine, only : kind_phys

    implicit none

    integer, intent(in) :: im,km
    integer, intent(in) :: kbcon1(im),ktcon(im)

    real(kind=kind_phys), intent(in) :: delt,grav
    real(kind=kind_phys), intent(in) :: lbb1,lbb2,lbb3,dt_decay
    real(kind=kind_phys), intent(in) :: omegain(im,km)
    real(kind=kind_phys), intent(in) :: del(im,km),zi(im,km)
    real(kind=kind_phys), intent(in) :: drag(im,km)
    real(kind=kind_phys), intent(in) :: buo(im,km)
    real(kind=kind_phys), intent(in) :: wush(im,km)

    real(kind=kind_phys), intent(inout) :: omegaout(im,km)
    
    logical, intent(in) :: cnvflg(im)
    logical, intent(in) :: first_time_step
    logical, intent(in) :: flag_restart

    !--------------------------------------------------------------------
    ! Local arrays
    !
    ! omega     = state at beginning of current internal substep
    ! omega_new = state at end of current internal substep
    !--------------------------------------------------------------------

    real(kind=kind_phys) :: omega(im,km)
    real(kind=kind_phys) :: omega_new(im,km)
    real(kind=kind_phys) :: omega_start(im,km)

    real(kind=kind_phys) :: termA(im,km)
    real(kind=kind_phys) :: termB(im,km)
    real(kind=kind_phys) :: termC(im,km)
    real(kind=kind_phys) :: memory_term,buoy_term,adv_term
    real(kind=kind_phys) :: decay_fac
    !--------------------------------------------------------------------
    ! Scalars
    !--------------------------------------------------------------------

    real(kind=kind_phys) :: dp,dz,pi_conv,discr
    real(kind=kind_phys) :: rhs_exp

    real(kind=kind_phys) :: omega_eps
    real(kind=kind_phys) :: a_eps,b_eps,disc_eps

    real(kind=kind_phys) :: cfl_target
    real(kind=kind_phys) :: dt_sub
    real(kind=kind_phys) :: dt_remaining
    real(kind=kind_phys) :: dt_cfl
    real(kind=kind_phys) :: time_done

    integer :: i,k

    logical :: active_point_found

    !--------------------------------------------------------------------
    ! Numerical parameters
    !--------------------------------------------------------------------

    ! Target CFL for explicit pressure-coordinate vertical advection.
    cfl_target = 0.8_kind_phys

    omega_eps = 1.0e-5_kind_phys
    a_eps     = 1.0e-12_kind_phys
    b_eps     = 1.0e-12_kind_phys
    disc_eps  = 1.0e-12_kind_phys

    ! Decay applied over one host physics timestep.
    decay_fac = exp(-delt/dt_decay)
    
    !--------------------------------------------------------------------
    ! Initialize from incoming prognostic tracer
    !--------------------------------------------------------------------

    do k = 1,km
       do i = 1,im

          termA(i,k) = 0.0_kind_phys
          termB(i,k) = 0.0_kind_phys
          termC(i,k) = 0.0_kind_phys

          omega(i,k)     = omegain(i,k)
          omega_new(i,k) = omegain(i,k)
          omegaout(i,k)  = omegain(i,k)

       enddo
    enddo

    !--------------------------------------------------------------------
    ! Retain memory when convection is inactive.
    !--------------------------------------------------------------------

    do k = 1,km
       do i = 1,im

          if (.not. cnvflg(i)) then

             omega(i,k)     = omegain(i,k)
             omega_new(i,k) = omega(i,k)
             omegaout(i,k)  = omega(i,k)

          endif

          ! Remove numerically negligible values.
          if (abs(omega(i,k)) < omega_eps) then

             omega(i,k)     = 0.0_kind_phys
             omega_new(i,k) = 0.0_kind_phys
             omegaout(i,k)  = 0.0_kind_phys

          endif

       enddo
    enddo

    !--------------------------------------------------------------------
    ! Decay prognostic updraft memory when convection is inactive.
    !
    ! The decay is applied once over the full host physics timestep.
    !--------------------------------------------------------------------

    do k = 1,km
       do i = 1,im

          if (.not. cnvflg(i)) then

             omega(i,k)     = omegain(i,k) * decay_fac
             omega_new(i,k) = omega(i,k)
             omegaout(i,k)  = omega(i,k)

          endif

          ! Remove numerically negligible values.
          if (abs(omega(i,k)) < omega_eps) then

             omega(i,k)     = 0.0_kind_phys
             omega_new(i,k) = 0.0_kind_phys
             omegaout(i,k)  = 0.0_kind_phys

          endif

       enddo
    enddo
    
    !--------------------------------------------------------------------
    ! Cold-start initialization
    !--------------------------------------------------------------------

    if (first_time_step .and. .not. flag_restart) then

       do k = 1,km
          do i = 1,im

             if (cnvflg(i)) then

                if (k >= kbcon1(i) .and. k < ktcon(i)) then

                   omega(i,k)     = -1.2_kind_phys
                   omega_new(i,k) = -1.2_kind_phys
                   omegaout(i,k)  = -1.2_kind_phys

                endif

             endif

          enddo
       enddo

    endif

    !--------------------------------------------------------------------
    ! Adaptive CFL-controlled subcycling
    !--------------------------------------------------------------------

    time_done = 0.0_kind_phys

    ! Save the state entering the prognostic integration so that the
    ! host-timestep local tendency can be returned after all substeps.
    omega_start(:,:) = omega(:,:)

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
                   dp = 1000.0_kind_phys * del(i,k)
                   if (dp > 0.0_kind_phys) then
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
       if (dt_sub <= 1.0e-8_kind_phys) then
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
                   omega_new(i,kbcon1(i)) = 0.0_kind_phys

                   dp = 1000.0_kind_phys * del(i,k)
                   dz = zi(i,k+1) - zi(i,k)

                   if (dp <= 0.0_kind_phys) cycle
                   if (abs(dz) <= 1.0e-12_kind_phys) cycle

                   pi_conv = dp/dz

                   !----------------------------------------------------
                   ! Explicit contributions
                   !----------------------------------------------------

                   memory_term = omega(i,k)

                   buoy_term = -0.5_kind_phys * dt_sub * lbb2      &
                        * buo(i,k) * pi_conv

                   if (k == kbcon1(i)) then
                      adv_term = 0.0_kind_phys
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

                   termA(i,k) = -0.5_kind_phys * dt_sub * lbb1     &
                        * drag(i,k) / pi_conv

                   termB(i,k) = 1.0_kind_phys                     &
                        + 0.5_kind_phys * dt_sub * wush(i,k)

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

                      discr = termB(i,k)**2                        &
                           - 4.0_kind_phys * termA(i,k) * termC(i,k)

                      if (discr >= -disc_eps) then

                         discr = max(discr,0.0_kind_phys)

                         omega_new(i,k) =                          &
                              (-termB(i,k) + sqrt(discr))          &
                              / (2.0_kind_phys * termA(i,k))

                      else

                         ! No real solution: retain previous state.
                         omega_new(i,k) = omega(i,k)

                      endif

                   endif

                   !----------------------------------------------------
                   ! Physical bounds for active updrafts
                   !----------------------------------------------------

                   omega_new(i,k) = max(                           &
                        min(omega_new(i,k),-1.2_kind_phys),         &
                        -80.0_kind_phys)

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
       if (delt-time_done < 1.0e-8_kind_phys*delt) then
          time_done = delt
       endif

    enddo

    !--------------------------------------------------------------------
    ! Return final prognostic state and host-timestep local tendency
    !--------------------------------------------------------------------

    omegaout(:,:) = omega(:,:)

  end subroutine progomega_calc

end module progomega
