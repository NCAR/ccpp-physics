!>\file dep_dry_mod.F90
!! This file is for the dry depostion driver.
!-------------REVISION HISTORY---------------!
! XX/XX/XXXX : original implementation (Ravan Ahmadov)
! 08/17/2023 : modified to follow Emerson et al., (2020) (Jordan Schnell)
! 08/17/2023 : gravitational settling folowing the coarse pm settling driver (Jordan Schnell)

module dep_dry_emerson_mod

  use machine ,        only : kind_phys
  use dep_data_mod     ! JLS
  use rrfs_smoke_config, only : num_chem, p_smoke, p_dust_1, p_coarse_pm, n_dbg_lines

  implicit none

  private

  public :: dry_dep_driver_emerson

contains
    subroutine dry_dep_driver_emerson(rmol,ustar,znt,ndvel,ddvel,vgrav,   &
               chem,delz,snowh,t_phy,p_phy,rho_phy,ivgtyp,g0,dt,          &
               settling_flag,drydep_flux,settling_flux,dbg_opt,           &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte, curr_secs, mpiid, xlat, xlong   )
!
! compute dry deposition velocity for aerosol particles
! Based on Emerson et al. (2020), PNAS,
! www.pnas.org/cgi/doi/10.1073/pnas.2014761117
! Code adapted from Hee-Ryu and Min, (2022):
! https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2021MS002792
!----------------------------------------------------------------------
  IMPLICIT NONE
 
       INTEGER,  INTENT(IN   ) ::  settling_flag,ndvel,                    &
                                   ids,ide, jds,jde, kds,kde,              &
                                   ims,ime, jms,jme, kms,kme,              &
                                   its,ite, jts,jte, kts,kte

       REAL(kind_phys) :: curr_secs

       REAL(kind_phys), DIMENSION( ims:ime , jms:jme )        ,            &
           INTENT(IN) :: ustar, rmol, znt, snowh
       REAL(kind_phys), DIMENSION( ims:ime , kms:kme , jms:jme ),          &
           INTENT(IN   ) :: t_phy, rho_phy, p_phy, delz                    
       INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  ivgtyp          
       REAL(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),  &
                                             INTENT(INOUT ) :: chem
       REAL(kind_phys), INTENT(IN) :: g0,dt
       LOGICAL, INTENT(IN) :: dbg_opt
 !
 ! Output arrays
       REAL(kind_phys), DIMENSION( ims:ime, jms:jme, ndvel ), INTENT(INOUT) :: ddvel
       REAL(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, ndvel), &
                                                              INTENT(OUT)   :: vgrav
       REAL(kind_phys), DIMENSION( ims:ime, jms:jme, ndvel ), INTENT(OUT)   :: settling_flux, drydep_flux
 ! Local
       REAL(kind_phys), DIMENSION( ims:ime , jms:jme )  :: aer_res
       REAL(kind_phys), DIMENSION( ndvel ) :: cblk
 ! Modpar variables, mass, density, diameter, knudsen number, mean free path
       REAL(kind_phys) :: pmasssn,pmassa,pmassc,pdensn,pdensa,pdensc,      &
                          dgnuc,dgacc,dgcor,knnuc,knacc,kncor,xlm
       real(kind_phys) :: Cc                      ! Cunningham/slip correction factor [-]
       real(kind_phys) :: DDp, Eb                 ! Brownian diffusion []
       real(kind_phys) :: Eim                     ! Impaction []
       real(kind_phys) :: Ein                     ! Interception []
       real(kind_phys) :: Sc                      ! Schmit number []
       real(kind_phys) :: St                      ! Stokes number []
       real(kind_phys) :: vg                      ! gravitational settling [cm/s]
       real(kind_phys) :: A, eps0, dumalnsg       ! land surface params [-]
       real(kind_phys) :: amu, amu_corrected      ! dynamic viscosity [g/s]
       real(kind_phys) :: airkinvisc              ! Air kinetic viscosity [cm2/s]
       real(kind_phys) :: freepath                ! Air molecular freepath [cm]
       real(kind_phys) :: dp                      ! aerosol diameter [cm]
       real(kind_phys) :: aerodens                ! aerosol density [g/cm3] 
       real(kind_phys) :: Rs                      ! Surface resistance
       real(kind_phys) :: vgpart
       real(kind_phys) :: growth_fac,vsettl,dtmax,conver,converi,dzmin
       real(kind_phys) :: rmol_local
       real(kind_phys), dimension( kts:kte) :: rho_col, delz_col
       real(kind_phys), dimension(ndvel)    :: dt_settl, chem_before, chem_after
       real(kind_phys), dimension( kts:kte, ndvel ) :: cblk_col, vg_col
       integer, dimension(ndvel) :: ndt_settl
       integer :: i, j, k, ntdt, nv
       integer :: icall=0
       integer, INTENT(IN) :: mpiid
       real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  xlat,xlong
 ! chem pointers (p_*) are not sequentially numbered, need to define so nv loops work
       integer, dimension(ndvel) :: chem_pointers
!> -- Gas constant
       real(kind_phys), parameter :: RSI = 8.314510
       chem_pointers(1) = p_smoke
       chem_pointers(2) = p_dust_1
       chem_pointers(3) = p_coarse_pm

       growth_fac = 1.0
       conver=1.e-9
       converi=1.e9
 
       if (mod(int(curr_secs),1800) .eq. 0) then
           icall = 0
       endif

       do j = jts, jte
          do i = its, ite
             aer_res(i,j) = 0.0
             rmol_local = rmol(i,j)
             do k = kts, kte
                delz_col(k) = delz(i,k,j)
                rho_col(k)  = rho_phy(i,k,j)
                do nv = 1, ndvel
                   cblk(nv) = chem(i,k,j,chem_pointers(nv))
                   if ( k == kts ) then
                      ddvel(i,j,nv) = 0.0
                      dt_settl(nv)  = 0.0
                   endif ! k==kts
                end do ! nv
                ! *** U.S. Standard Atmosphere 1962 page 14 expression
                !     for dynamic viscosity = beta * T * sqrt(T) / ( T + S)
                !     where beta = 1.458e-6 [ kg sec^-1 K**-0.5 ], s = 110.4 [ K ].
                amu = 1.458E-6 * t_phy(i,k,j) * sqrt(t_phy(i,k,j)) / ( t_phy(i,k,j) + 110.4 )
                ! Aerodynamic resistance
                call depvel( rmol_local, dep_ref_hgt, znt(i,j), ustar(i,j), vgpart, aer_res(i,j) )
                ! depvel uses meters, need to convert to s/cm
                aer_res(i,j) = max(aer_res(i,j)/100._kind_phys,0._kind_phys) 
                ! Air kinematic viscosity (cm^2/s)
                airkinvisc = ( 1.8325e-4 * ( 416.16 / ( t_phy(i,k,j) + 120.0 ) ) *   &
                             ( ( t_phy(i,k,j) / 296.16 )**1.5 ) ) / ( rho_phy(i,k,j) / 1.e3 ) ! Convert density to mol/cm^3
                ! Air molecular freepath (cm)  ! Check against XLM from above
                freepath = 7.39758e-4 * airkinvisc / sqrt( t_phy(i,k,j) )
                do nv = 1, ndvel
                   if ( chem_pointers(nv) == p_smoke ) then
                      dp = 4.E-8 !dgacc
                      aerodens = 1.4e+3 !pdensa
                   elseif ( chem_pointers(nv) == p_dust_1) then
                      dp = 1.E-6 !dgacc 
                      aerodens = 2.6e+3 !pdensa
                   elseif ( chem_pointers(nv) == p_coarse_pm ) then
                      dp = 4.5E-6 !dgcor
                      aerodens = 2.6e+3 !pdensc
                   else
                      continue
                   endif
                   ! Convert diameter to cm and aerodens to g/cm3
                   aerodens = aerodens / 1000.
                   dp = dp * 1e+2
                   ! Cunningham correction factor
                   Cc = 1. + 2. * freepath / dp * ( 1.257 + 0.4*exp( -0.55 * dp / freepath ) )
                   ! Corrected dynamic viscosity (used for settling)
                   amu_corrected = amu / Cc
                   ! Gravitational Settling
                   vg = aerodens * dp * dp * gravity * 100. * Cc / &       ! Convert gravity to cm/s^2
                          ( 18. * airkinvisc * ( rho_phy(i,k,j) / 1.e3 ) ) ! Convert density to mol/cm^3
                   ! -- Rest of loop for the surface when deposition velocity needs to be cacluated
                   if ( k == kts ) then
                      ! Brownian Diffusion
                      DDp = ( boltzmann * t_phy(i,k,j) ) * Cc / (3. * pi * airkinvisc * ( rho_phy(i,k,j) / 1.e3 )  * dp) ! Convert density to mol/cm^3
                      ! Schmit number
                      Sc = airkinvisc / DDp
                      ! Brownian Diffusion
                      Eb = Cb * Sc**(-0.666666667)
                      ! Stokes number
                      St = ( 100. * ustar(i,j) ) * ( 100.* ustar(i,j) ) * vg / airkinvisc / ( gravity * 100.) ! Convert ustar to cm/s, gravity to cm/s^2
                      ! Impaction 
                      Eim = Cim * ( St / ( alpha + St ) )**1.7
                      ! MODIS type lu, large roughness lengths (e.g., urban or forest)
                      ! -----------------------------------------------------------------------
                      ! *** TO DO -- set A and eps0 for all land surface types *** !!!
                      ! -----------------------------------------------------------------------
                      if ( ivgtyp(i,j) .eq. 13 .or. ivgtyp(i,j) .le. 5 ) then ! Forest
                         A = A_for
                         eps0 = eps0_for
                      else if ( ivgtyp(i,j) .eq. 17 ) then ! water
                         A = A_wat
                         eps0 = eps0_wat
                      else ! otherwise
                         A = A_grs
                         eps0 = eps0_grs
                      end if
                      ! Set if snow greater than 1 cm
                      ! Interception
                      Ein = Cin * ( dp / A )**vv
                      ! Surface resistance
                      Rs = 1. / ( ( ustar(i,j) * 100.) * ( Eb + Eim + Ein) * eps0 ) ! Convert ustar to cm/s
                      ! Compute final ddvel = aer_res + RS, set max at max_dep_vel in dep_data_mod.F[ m/s]
                      ! The /100. term converts from cm/s to m/s, required for MYNN.
                      if ( settling_flag == 1 ) then
                         ddvel(i,j,nv) = max(min( ( vg + 1./(aer_res(i,j)+Rs) )/100., max_dep_vel),0._kind_phys)
                      else
                         ddvel(i,j,nv) = max(min( ( 1./(aer_res(i,j)+Rs) )/100., max_dep_vel),0._kind_phys)
                      endif
                      if ( dbg_opt .and. (icall .le. n_dbg_lines) ) then
                         WRITE(1000+mpiid,*) 'dry_dep_mod_emer:xlat,xlong,curr_secs,nv',xlat(i,j),xlong(i,j),int(curr_secs),nv
                         WRITE(1000+mpiid,*) 'dry_dep_mod_emer:xlat,xlong,curr_secs,deposition velocity (m/s)',xlat(i,j),xlong(i,j),int(curr_secs),ddvel(i,j,nv)
                         icall = icall + 1
                      endif
                      drydep_flux(i,j,nv) = chem(i,kts,j,chem_pointers(nv))*rho_phy(i,k,j)*ddvel(i,j,nv)/100.0*dt
                   endif ! k == kts
                   vgrav(i,k,j,nv) = vg
                   ! Fill column variables
                   cblk_col(k,nv) = cblk(nv)
                   vg_col(k,nv)   = vg
                enddo ! nv
             enddo ! k
             ! -- Get necessary info for settling
             ! -- Determine the maximum time-step satisying the CFL condition:
             dzmin = minval(delz_col)
             ntdt=INT(dt)
             do nv = 1, ndvel
               ! -- NOTE, diameters and densities are NOT converted to cm and g/cm3 like above
               ! -- dt_settl calculations (from original coarsepm_settling)
               if ( chem_pointers(nv) == p_smoke ) then
                 dp = 4.E-8 !dgacc
                 aerodens = 1.4e+3 !pdensa
               elseif ( chem_pointers(nv) == p_dust_1) then
                 dp = 1.E-6 !dgacc 
                 aerodens = 2.6e+3 !pdensa
               elseif ( chem_pointers(nv) == p_coarse_pm ) then
                 dp = 4.5E-6 !dgcor
                 aerodens = 2.6e+3 !pdensc
               else
                 continue
               endif
               ! 1.5E-5 = dyn_visc --> dust_data_mod.F90
               vsettl = 2.0 / 9.0 * g0 * aerodens * ( growth_fac * ( 0.5 * dp ))**2.0 / ( 0.5 * 1.5E-5 )
               dtmax = dzmin / vsettl
               ndt_settl(nv) = MAX( 1, INT( ntdt /dtmax) )
               ! Limit maximum number of iterations
               IF (ndt_settl(nv) > 12) ndt_settl(nv) = 12 
               dt_settl(nv) = REAL(ntdt,kind=kind_phys) /REAL(ndt_settl(nv),kind=kind_phys)
             enddo
             ! Perform gravitational settling if desired
             if ( settling_flag == 1 ) then
                call particle_settling(cblk_col,rho_col,delz_col,vg_col,dt_settl,ndt_settl,ndvel,kts,kte)
             endif
             ! Put cblk back into chem array
             do nv= 1, ndvel
                do k = kts, kte
                   chem(i,k,j,chem_pointers(nv)) = cblk_col(k,nv)
                enddo ! k
             enddo ! nv
        end do ! j
        end do ! i
end subroutine dry_dep_driver_emerson
!
!--------------------------------------------------------------------------------
!
subroutine depvel( rmol, zr, z0, ustar, vgpart, aer_res )
!--------------------------------------------------
!     THIS FUNCTION HAS BEEN DESIGNED TO EVALUATE AN UPPER LIMIT
!     FOR THE POLLUTANT DEPOSITION VELOCITY AS A FUNCTION OF THE
!     SURFACE ROUGHNESS AND METEOROLOGICAL CONDITIONS.
!     PROGRAM WRITTEN BY GREGORY J.MCRAE (NOVEMBER 1977)
!         Modified by Darrell A. Winner  (Feb. 1991)
!                  by Winfried Seidl     (Aug. 1997)
!.....PROGRAM VARIABLES...
!     RMOL     - RECIPROCAL OF THE MONIN-OBUKHOV LENGTH
!     ZR       - REFERENCE HEIGHT
!     Z0       - SURFACE ROUGHNESS HEIGHT
!     USTAR    - FRICTION VELOCITY U*
!     AER_RES  - AERODYNAMIC RESISTANCE
!.....REFERENCES...
!     MCRAE, G.J. ET AL. (1983) MATHEMATICAL MODELING OF PHOTOCHEMICAL
!       AIR POLLUTION, ENVIRONMENTAL QUALITY LABORATORY REPORT 18,
!       CALIFORNIA INSTITUTE OF TECHNOLOGY, PASADENA, CALIFORNIA.
!.....RESTRICTIONS...
!     1. THE MODEL EDDY DIFFUSIVITIES ARE BASED ON MONIN-OBUKHOV
!        SIMILARITY THEORY AND SO ARE ONLY APPLICABLE IN THE
!        SURFACE LAYER, A HEIGHT OF O(30M).
!     2. ALL INPUT UNITS MUST BE CONSISTENT
!     3. THE PHI FUNCTIONS USED TO CALCULATE THE FRICTION
!        VELOCITY U* AND THE POLLUTANT INTEGRALS ARE BASED
!        ON THE WORK OF BUSINGER ET AL.(1971).
!     4. THE MOMENTUM AND POLLUTANT DIFFUSIVITIES ARE NOT
!        THE SAME FOR THE CASES L<0 AND L>0.
!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        REAL(kind_phys), intent(in)    :: ustar, z0, zr
        REAL(kind_phys), intent(out)   :: vgpart, aer_res
        REAL(kind_phys), intent(inout) :: rmol
!--------------------------------------------------
! .. Local Scalars ..
!--------------------------------------------------
        INTEGER :: l
        REAL(kind_phys)    :: ao, ar, polint, vk
!--------------------------------------------------
! .. Intrinsic Functions ..
!--------------------------------------------------
        INTRINSIC alog
!--------------------------------------------------
!     Set the von Karman constant
!--------------------------------------------------
        vk = karman
!--------------------------------------------------
!     DETERMINE THE STABILITY BASED ON THE CONDITIONS
!             1/L < 0 UNSTABLE
!             1/L = 0 NEUTRAL
!             1/L > 0 STABLE
!--------------------------------------------------
        if(abs(rmol) < 1.E-6 ) rmol = 0.
        IF (rmol<0) THEN
          ar = ((1.0-9.0*zr*rmol)**(0.25)+0.001)**2
          ao = ((1.0-9.0*z0*rmol)**(0.25)+0.001)**2
          polint = 0.74*(alog((ar-1.0)/(ar+1.0))-alog((ao-1.0)/(ao+1.0)))
        ELSE IF (rmol==0.) THEN
          polint = 0.74*alog(zr/z0)
        ELSE
          polint = 0.74_kind_phys*alog(zr/z0) + 4.7_kind_phys*rmol*(zr-z0)
        END IF
        vgpart = ustar*vk/polint
        aer_res = polint/(karman*max(ustar,1.0e-4))
end subroutine depvel
!
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
!
subroutine particle_settling(cblk,rho_phy,delz,vg,dt_settl,ndt_settl,ndvel,kts,kte)
     IMPLICIT NONE
     
     INTEGER, INTENT(IN ) :: kts, kte, ndvel
     REAL(kind_phys), DIMENSION(kts:kte), INTENT (IN)  :: rho_phy, delz
     REAL(kind_phys), DIMENSION(kts:kte,ndvel), INTENT(IN) :: vg
     REAL(kind_phys), DIMENSION(kts:kte,ndvel), INTENT(INOUT) :: cblk
     REAL(kind_phys), DIMENSION(ndvel), INTENT(IN) :: dt_settl
     INTEGER, DIMENSION(ndvel), INTENT(IN) :: ndt_settl
!
!--- Local------
     INTEGER :: k,nv,n,l2
     REAL(kind_phys) :: temp_tc, transfer_to_below_level, vd_wrk1
     REAL(kind_phys), DIMENSION(kts:kte) :: delz_flip 
     
     do k = kts,kte
        delz_flip(k) = delz(kte-k+kts)
     enddo

     do nv = 1, ndvel
     do n = 1,ndt_settl(nv)
     transfer_to_below_level = 0.0
     do k = kte,kts,-1
        l2 = kte - k + 1
     
        temp_tc = cblk(k,nv)
     
        vd_wrk1 = dt_settl(nv) * vg(k,nv)/100. / delz_flip(l2)               ! convert vg to m/s
     
        cblk(k,nv)= cblk(k,nv) * (1. - vd_wrk1) + transfer_to_below_level
        if (k.gt.kts) then
            transfer_to_below_level =(temp_tc*vd_wrk1)*((delz_flip(l2) &
                        *rho_phy(k))/(delz_flip(l2+1)*rho_phy(k-1)))          ! [ug/kg]
        endif
     enddo ! k
     enddo ! n
     enddo ! nv
end subroutine particle_settling

!
end module dep_dry_emerson_mod
