!>\file progsigma_calc.f90

!> This module contains the subroutine that calculates the prognostic
!! updraft area fraction that is used for closure computations in
!! saSAS deep and shallow convection, based on a moisture budget
!! as described in Bengtsson et al. 2022 \cite Bengtsson_2022.
      module progsigma

        implicit none

        public progsigma_calc

      contains

!> This subroutine computes a prognostic updraft area fraction
!! used in the closure computations in the samfdeepcnv.f scheme
!! This subroutine computes a prognostic updraft area fracftion
!! used in the closure computations in the samfshalcnv. scheme
!!\section gen_progsigma progsigma_calc General Algorithm 
      subroutine progsigma_calc (im,km,flag_init,flag_restart,flag_shallow,&
           flag_mid,del,tmf,qmicro,dbyo1,zdqca,omega_u,zeta,hvap,          &
           delt,qadv,kbcon1,ktcon,cnvflg,betascu,betamcu,betadcu,          &
           sigmind,sigminm,sigmins,sigmain,sigmaout,sigmab)
!                                                           
!                                                                                                                                             
      use machine,  only : kind_phys
      use funcphys, only : fpvs

      implicit none

!     intent in
      integer, intent(in)  :: im,km,kbcon1(im),ktcon(im)
      real(kind=kind_phys), intent(in)  :: hvap,delt,betascu,betamcu,betadcu, &
                                           sigmind,sigminm,sigmins
      real(kind=kind_phys), intent(in)  :: qadv(im,km),del(im,km),    &
           qmicro(im,km),tmf(im,km),dbyo1(im,km),zdqca(im,km),           &
           omega_u(im,km),zeta(im,km)
      logical, intent(in)  :: flag_init,flag_restart,cnvflg(im),flag_shallow,flag_mid
      real(kind=kind_phys), intent(in) :: sigmain(im,km)

!     intent out
      real(kind=kind_phys), intent(out) :: sigmaout(im,km)
      real(kind=kind_phys), intent(out) :: sigmab(im)


!     Local variables
      integer              :: i,k,km1
      real(kind=kind_phys) :: termA(im),termB(im),termC(im),termD(im)
      real(kind=kind_phys) :: mcons(im),fdqa(im),form(im,km),              &
           dp(im,km),inbu(im,km)                         
                          

      real(kind=kind_phys) :: gcvalmx,epsilon,ZZ,cvg,mcon,buy2,   &
                          fdqb,dtdyn,dxlim,rmulacvg,tem,     &
                          DEN,dp1,invdelt

     !Parameters
      gcvalmx = 0.1
      rmulacvg=10.
      epsilon=1.E-11
      km1=km-1
      invdelt = 1./delt

     !Initialization 2D
      do k = 1,km
         do i = 1,im
            inbu(i,k)=0.
            form(i,k)=0.
            dp(i,k)=0.
         enddo
      enddo
     
     !Initialization 1D
      do i=1,im
         sigmab(i)=0.
         termA(i)=0.
         termB(i)=0.
         termC(i)=0.
         termD(i)=0.
         fdqa(i)=0.
         mcons(i)=0.
      enddo

      do k = 2,km1
          do i = 1,im
             if(cnvflg(i))then
                dp(i,k) = 1000. * del(i,k)
             endif
          enddo
       enddo

      !Initial computations, place maximum sigmain in sigmab
       do i=1,im
          if(cnvflg(i))then
             do k=2,km
                if(sigmain(i,k)>sigmab(i))then
                   sigmab(i)=sigmain(i,k)
                endif
             enddo
          endif
       enddo
     
      do i=1,im
         if(cnvflg(i))then
            if(sigmab(i) < 1.E-5)then !after advection
               sigmab(i)=0.                                  
            endif
         endif
      enddo
          
      !compute termD "The vertical integral of the latent heat convergence is limited to the                                        
      !buoyant layers with positive moisture convergence (accumulated from the surface).                                                       
      !Lowest level:                                                                                                               
       do i = 1,im
          dp1 = 1000. * del(i,1)
          mcons(i)=(hvap*(qadv(i,1)+tmf(i,1)+qmicro(i,1))*dp1)
       enddo
      !Levels above:
       do k = 2,km1
          do i = 1,im
             if(cnvflg(i))then
                mcon = (hvap*(qadv(i,k)+tmf(i,k)+qmicro(i,k))*dp(i,k))
                buy2 = termD(i)+mcon+mcons(i)
!               Do the integral over buoyant layers with positive mcon acc from surface
                if(dbyo1(i,k)>0 .and. buy2 > 0.)then
                   inbu(i,k)=1.
                endif
                inbu(i,k-1)=MAX(inbu(i,k-1),inbu(i,k))
                termD(i) = termD(i) + inbu(i,k-1)*mcons(i)
                mcons(i)=mcon
             endif
          enddo
       enddo

       !termA
       do k = 2,km1
          do i = 1,im
             if(cnvflg(i))then
                tem=sigmab(i)*zeta(i,k)*inbu(i,k)*dbyo1(i,k)*dp(i,k)
                termA(i)=termA(i)+tem
             endif
          enddo
       enddo

       !termB                                                                                                             
       do k = 2,km1
          do i = 1,im
             if(cnvflg(i))then
                tem=zeta(i,k)*dbyo1(i,k)*inbu(i,k)*dp(i,k)
                termB(i)=termB(i)+tem
             endif
          enddo
       enddo

      !termC
       do k = 2,km1
          do i = 1,im
             if(cnvflg(i))then
                form(i,k)=-1.0*inbu(i,k)*(omega_u(i,k)*delt)
                fdqb=0.5*((form(i,k)*zdqca(i,k)))
                termC(i)=termC(i)+inbu(i,k)*   &
                     (fdqb+fdqa(i))*hvap*zeta(i,k)
                fdqa(i)=fdqb
             endif
         enddo
      enddo

      !sigmab
      if(flag_init .and. .not. flag_restart)then
         do i = 1,im
            if(cnvflg(i))then
               sigmab(i)=0.03
            endif
         enddo
      else
         do i = 1,im
            if(cnvflg(i))then
               DEN=MIN(termC(i)+termB(i),1.E8)
               cvg=termD(i)*delt
               ZZ=MAX(0.0,SIGN(1.0,termA(i)))            &
                    *MAX(0.0,SIGN(1.0,termB(i)))         &
                    *MAX(0.0,SIGN(1.0,termC(i)-epsilon))
               cvg=MAX(0.0,cvg)
               sigmab(i)=(ZZ*(termA(i)+cvg))/(DEN+(1.0-ZZ))
               if(sigmab(i)>0.)then
                  sigmab(i)=MIN(sigmab(i),0.95)  
                  sigmab(i)=MAX(sigmab(i),0.01)
               endif
            endif!cnvflg
         enddo
      endif

      do k=1,km
         do i=1,im
            if(cnvflg(i))then
               sigmaout(i,k)=sigmab(i)
            endif
         enddo
      enddo

      !Reduce area fraction before coupling back to mass-flux computation. 
      if(flag_shallow)then
         do i= 1, im
            if(cnvflg(i)) then
               sigmab(i)=sigmab(i)/betascu
               sigmab(i)=MAX(sigmins,sigmab(i))
            endif
         enddo
      elseif(flag_mid)then
         do i= 1, im
            if(cnvflg(i)) then
               sigmab(i)=sigmab(i)/betamcu
               sigmab(i)=MAX(sigminm,sigmab(i))
            endif
         enddo
      else
         do i= 1, im
            if(cnvflg(i)) then
               sigmab(i)=sigmab(i)/betadcu
               sigmab(i)=MAX(sigmind,sigmab(i))
            endif
         enddo
      endif
      do i= 1, im
        sigmab(i) = MIN(0.95,sigmab(i))
      enddo

     end subroutine progsigma_calc

end module progsigma
