module coarsepm_settling_mod

  use machine ,      only : kind_phys
  use dust_data_mod, only : dyn_visc !den_dust, reff_dust, dyn_visc

  implicit none

CONTAINS


SUBROUTINE coarsepm_settling_driver(dt,t_phy,                         &
                                  chem,rho_phy,dz8w,p8w,p_phy,sedim,  &
                                  area,g,num_chem,                    &
                                  ids,ide, jds,jde, kds,kde,          &
                                  ims,ime, jms,jme, kms,kme,          &
                                  its,ite, jts,jte, kts,kte           )

  IMPLICIT NONE

  INTEGER,      INTENT(IN   ) ::                                           &
                                  num_chem,                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
  REAL(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),INTENT(INOUT ) :: chem
  REAL(kind_phys), DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::  t_phy,p_phy,dz8w,p8w,rho_phy
  REAL(kind_phys), DIMENSION( ims:ime ,  jms:jme ),INTENT(IN   ) ::  area
  REAL(kind_phys), INTENT(IN   ) :: dt,g

  REAL(kind_phys), DIMENSION( ims:ime, jms:jme, num_chem ), INTENT(OUT  ) :: sedim

  integer :: nv,i,j,k,kk,lmx,idust
  real(kind_phys), DIMENSION (1,1,kte-kts+1) :: tmp,airden,airmas,p_mid,delz,rh
  real(kind_phys), DIMENSION (1,1,kte-kts+1,1) :: dust
  real(kind_phys), DIMENSION (ime,jme,kme,num_chem) :: chem_before
!
! bstl is for budgets
!
! real(kind_phys), DIMENSION (5), PARAMETER :: den_dust(5)=(/2500.,2650.,2650.,2650.,2650./)
! real(kind_phys), DIMENSION (5), PARAMETER :: reff_dust(5)=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/)
  real(kind_phys), DIMENSION (1), PARAMETER :: den_dust (1)=(/2650. /)
  real(kind_phys), DIMENSION (1), PARAMETER :: reff_dust(1)=(/2.4D-6/)
  real(kind_phys), DIMENSION (1) :: bstl_dust
  real(kind_phys) conver,converi
  real(kind_phys),parameter::max_default=0.

       sedim = 0.
       conver=1.e-9
       converi=1.e9
       lmx=kte-kts+1
!
       do j=jts,jte
       do i=its,ite
! 
! initialize some met stuff
!
          kk=0
          bstl_dust(:)=0.
          do k=kts,kte 
          kk=kk+1
          p_mid(1,1,kk)=.01*p_phy(i,kte-k+kts,j) 
          delz(1,1,kk)=dz8w(i,kte-k+kts,j) 
          airmas(1,1,kk)=-(p8w(i,k+1,j)-p8w(i,k,j))*area(i,j)/g
          airden(1,1,kk)=rho_phy(i,k,j)
          tmp(1,1,kk)=t_phy(i,k,j)
          do nv = 1, num_chem
            chem_before(i,j,k,nv) =  chem(i,k,j,nv)  
          enddo
          enddo
!
! max dust in column
!
          idust=1
          kk=0
          do k=kts,kte 
            kk=kk+1
            dust(1,1,kk,1)=chem(i,k,j,1)*conver
          enddo


          call settling(1, 1, lmx, 1,g,dyn_visc, &
                    dust, tmp, p_mid, delz, airmas, &
                    den_dust, reff_dust, dt, bstl_dust, idust, airden)

           kk = 0
          do k = kts,kte
             kk = kk+1
             chem(i,k,j,1)=dust(1,1,kk,1)*converi          ! coarse dust [ug/kg]
          enddo
!
!
!
          do nv = 1, num_chem
            do k = kts,kte
              sedim(i,j,nv) = sedim(i,j,nv)+(chem_before(i,j,k,nv) - chem(i,k,j,nv))*p8w(i,k,j)/g
            enddo
            sedim(i,j,nv) = sedim(i,j,nv) / dt  !ug/m2/s
          enddo
!
!
!
       enddo
       enddo
!
!
!
END SUBROUTINE coarsepm_settling_driver


          subroutine settling(imx,jmx, lmx, nmx,g0,dyn_visc, &
                    tc, tmp, p_mid, delz, airmas, &
                    den, reff, dt, bstl, idust, airden)
! ****************************************************************************
! *                                                                          *
! *  Calculate the loss by settling, using an implicit method                *
! *                                                                          *
! *  Input variables:                                                        *
! *    SIGE(k)         - sigma coordinate of the vertical edges              *
! *    PS(i,j)         - Surface pressure (mb)                               *
! *    TMP(i,j,k)      - Air temperature  (K)                                *
! *    CT(i,j)         - Surface exchange coeff for moisture
! *                                                                          *
! **************************************************************************** 


  IMPLICIT  NONE

  INTEGER, INTENT(IN) :: imx, jmx, lmx, nmx,idust
  INTEGER :: ntdt
  REAL(kind_phys), INTENT(IN) :: dt,g0,dyn_visc
  REAL(kind_phys),    INTENT(IN) :: tmp(imx,jmx,lmx), delz(imx,jmx,lmx),  &
                         airmas(imx,jmx,lmx), &
                         den(nmx), reff(nmx),p_mid(imx,jmx,lmx),&
                         airden(imx,jmx,lmx)
  REAL(kind_phys), INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL(kind_phys), INTENT(OUT)   :: bstl(imx,jmx,nmx)

  REAL(kind_phys)    :: tc1(imx,jmx,lmx,nmx), dt_settl(nmx), rcm(nmx), rho(nmx)
  INTEGER :: ndt_settl(nmx)
  REAL(kind_phys)    :: dzmin, vsettl, dtmax, rhb, rwet(nmx), ratio_r(nmx)
  REAL(kind_phys)    :: c_stokes, free_path, c_cun, viscosity,  growth_fac
  REAL(kind_phys)    :: vd_cor(lmx),vd_wk1 
  INTEGER :: k, n, i, j, l, l2
  REAL(kind_phys)    :: transfer_to_below_level,temp_tc

  ! for OMP:
  REAL(kind_phys) :: rwet_priv(nmx), rho_priv(nmx)

  ! executable statements

  bstl = 0._kind_phys

  if(idust.ne.1.)return

!!!  WHERE (tc(:,:,:,:) < 0.0) tc(:,:,:,:) = 1.0E-32

  dzmin = MINVAL(delz(:,:,:))
  IF (idust == 1)     growth_fac = 1.0

  DO k = 1,nmx

     ! Settling velocity (m/s) for each tracer (Stokes Law)
     ! DEN         density                        (kg/m3)
     ! REFF        effective radius               (m)
     ! dyn_visc    dynamic viscosity              (kg/m/s)
     ! g0          gravity                        (m/s2)
     ! 3.0         corresponds to a growth of a factor 3 of radius with 100% RH
     ! 0.5         upper limit with temp correction

     tc1(:,:,:,k) = tc(:,:,:,k)
     vsettl = 2.0/9.0 * g0 * den(k) * (growth_fac*reff(k))**2 / &
              (0.5*dyn_visc)

     ! Determine the maximum time-step satisying the CFL condition:
     ! dt <= (dz)_min / v_settl
     ntdt=INT(dt)
     dtmax = dzmin / vsettl
     ndt_settl(k) = MAX( 1, INT( ntdt /dtmax) )
     ! limit maximum number of iterations
     IF (ndt_settl(k) > 12) ndt_settl(k) = 12
     dt_settl(k) = REAL(ntdt) / REAL(ndt_settl(k))

     ! Particles radius in centimeters
     IF (idust.eq.1)then
          rwet(k) = reff(k)
          ratio_r(k) = 1.0
          rho(k) = den(k)
      endif
  END DO

  ! Solve the bidiagonal matrix (l,l)

!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( i,   j,   l,   l2, n,   k,   rhb, rwet_priv, ratio_r, c_stokes)&
!$OMP PRIVATE( free_path, c_cun, viscosity, rho_priv, vd_cor )

  ! Loop over latitudes
  DO j = 1,jmx
 
     DO k = 1,nmx
        IF (idust.eq.1) THEN
           rwet_priv(k) = rwet(k)
           rho_priv(k)  = rho(k)
        END IF

        DO n = 1,ndt_settl(k)

           ! Solve each vertical layer successively (layer l)
        transfer_to_below_level=0
 
           DO l = lmx,1,-1
              l2 = lmx - l + 1

!           DO j = 1,jmx
              DO i = 1,imx

                 ! Dynamic viscosity
                 c_stokes = 1.458E-6 * tmp(i,j,l)**1.5/(tmp(i,j,l) + 110.4) 

                 ! Mean free path as a function of pressure (mb) and 
                 ! temperature (K)
                 ! order of p_mid is top->sfc
                 free_path = 1.1E-3/p_mid(i,j,l2)/SQRT(tmp(i,j,l))
!!!                 free_path = 1.1E-3/p_edge(i,j,l2)/SQRT(tmp(i,j,l))

                 ! Slip Correction Factor
                 c_cun = 1.0+ free_path/rwet_priv(k)* &
                      (1.257 + 0.4*EXP(-1.1*rwet_priv(k)/free_path))

                 ! Corrected dynamic viscosity (kg/m/s)
                 viscosity = c_stokes / c_cun

                 ! Settling velocity

                 vd_cor(l) = 2.0/9.0*g0*rho_priv(k)*rwet_priv(k)**2/viscosity

            ! Update mixing ratio; order of delz: top->sfc
            temp_tc=tc(i,j,l,k)      !temp_tc - for temporal storage [ug/kg]            
            vd_wk1 = dt_settl(k)*vd_cor(l)/delz(i,j,l2)   !fraction to leave level

            tc(i,j,l,k)   =  tc(i,j,l,k)*(1.- vd_wk1)+transfer_to_below_level ! [ug/kg]
            
           if (l.gt.1) transfer_to_below_level =(temp_tc*vd_wk1)*((delz(i,j,l2) &
                   *airden(i,j,l))/(delz(i,j,l2+1)*airden(i,j,l-1)))          ! [ug/kg]

              END DO   !i
!           END DO   !j
        END DO  !l

     END DO  !n
  END DO  !k

  END DO   !j
!$OMP END PARALLEL DO

  DO n = 1,nmx
     DO i = 1,imx
        DO j = 1,jmx
           bstl(i,j,n) = 0._kind_phys
           DO l = 1,lmx
              IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0D-32
              bstl(i,j,n) = bstl(i,j,n) + &
                   (tc(i,j,l,n) - tc1(i,j,l,n)) * airmas(i,j,l)
           END DO
        END DO
     END DO
  END DO
  
END SUBROUTINE settling

end module coarsepm_settling_mod
