      module progomega

        implicit none

        public progomega_calc

      contains

!>\file progomega_calc.f90
!! This file contains the subroutine that calculates the prognostic
!! updraft vertical velocity that is used for closure computations in 
!! saSAS and C3 deep and shallow convection. 

!>\ingroup SAMFdeep
!>\ingroup SAMF_shal
!> This subroutine computes a prognostic updraft vertical velocity
!! used in the closure computations in the samfdeepcnv.f and cu_c3_conv.f scheme
!! This subroutine computes a prognostic updraft vertical velocity
!! used in the closure computations in the samfshalcnv. and cu_c3_shal scheme
!!\section gen_progomega progomega_calc General Algorithm
       
   subroutine progomega_calc(first_time_step,flag_restart,im,km,kbcon1,ktcon,omegain,delt,del, &
        zi,cnvflg,omegaout,grav,buo,drag,wush,tentr,bb1,bb2)
     
     use machine,  only : kind_phys
     use funcphys, only : fpvs  
     implicit none

     integer, intent(in)  :: im, km
     integer, intent(in)  :: kbcon1(im),ktcon(im)
     real(kind=kind_phys), intent(in)  :: delt,grav,bb1,bb2
     real(kind=kind_phys), intent(in)  :: omegain(im,km), del(im,km),zi(im,km)
     real(kind=kind_phys), intent(in)  :: drag(im,km),buo(im,km),wush(im,km),tentr(im,km)
     real(kind=kind_phys), intent(out) :: omegaout(im,km)
     logical, intent(in)               :: cnvflg(im),first_time_step,flag_restart
     real(kind=kind_phys) :: termA(im,km),termB(im,km),termC(im,km),omega(im,km)
     real(kind=kind_phys) :: RHS(im,km),Kd(im,km)
     real(kind=kind_phys) :: dp,dz,entrn,Kdn,discr,wush_pa,lbb1,lbb2,lbb3
     integer              :: i,k

     entrn = 0.8E-4 !0.5E-4 !m^-1
     Kdn   = 0.5E-4 !2.9E-4 !m^-1
     lbb1  = 0.5 !1.0 
     lbb2  = 3.2 !3.0
     lbb3  = 0.5 !0.5
     
     
     !Initialization 2D
     do k = 1,km
        do i = 1,im
           termA(i,k)=0.
           termB(i,k)=0.
           termC(i,k)=0.
           RHS(i,k)=0.
           omega(i,k)=omegain(i,k)
        enddo
     enddo

     if(first_time_step .and. .not. flag_restart)then
        do k = 1,km
           do i = 1,im
              if(cnvflg(i))then
                 omega(i,k)=-40.0 !Pa/s 
              endif
           enddo
        enddo
     endif
     
     ! Compute RHS terms
     !Lisa Bengtsson: !  compute updraft velocity omega (Pa/s)
     !> - Expand the steady state solution of updraft velocity from Han et al.'s (2017)
     !> \cite han_et_al_2017 equation 7 to include the time-derivative, and an aerodynamic
     !> drag term from Guérémy 2016.
     !> Solve using implicit time-stepping scheme, solving the quadratic equation for omega. 
     
     do k = 2, km
        do i = 1, im
           if (cnvflg(i)) then
              if (k > kbcon1(i) .and. k < ktcon(i)) then

                 ! Aerodynamic drag parameter
                 Kd(i,k) = (tentr(i,k)/entrn)*Kdn

                 ! Scale by dp/dz to have equation in Pa/s
                 !(dp/dz > 0)
                 dp = 1000. * del(i,k)
                 dz = zi(i,k+1) - zi(i,k)
                 
                 !termA	- Ensures quadratic damping (drag).
                 !termB	- Ensures linear damping from wind shear.
                 !termC - Adds buoyancy forcing 
                 
                 !Coefficients for the quadratic equation
                 termA(i,k) = delt * ((lbb1 * drag(i,k) * (dp/dz)) + (Kd(i,k) * (dp/dz)))
                 termB(i,k) = -1.0 - delt * lbb3 * wush(i,k) * dp/dz
                 termC(i,k) = omega(i,k) - delt * lbb2 * buo(i,k) * (dp/dz) &
                      - delt * 0.5 * (omega(i,k)**2 - omega(i,k-1)**2) / dp
                 !Compute the discriminant
                 discr = termB(i,k)**2 - 4.0 * termA(i,k) * termC(i,k)

                 ! Check if discriminant is non-negative
                 if (discr >= 0.0) then
                 ! Solve quadratic equation, take the negative root
                 omegaout(i,k) = (-termB(i,k) - sqrt(discr)) / (2.0 * termA(i,k))
                 else
                 omegaout(i,k) = omega(i,k)
                 endif

                 omegaout(i,k) = MAX(MIN(omegaout(i,k), -1.2), -80.0)
                
              endif
           endif
        enddo
     enddo
     
    end subroutine progomega_calc
end module progomega
