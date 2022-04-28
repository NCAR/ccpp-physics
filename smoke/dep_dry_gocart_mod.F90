!>\file dep_dry_gocart_mod.F90
!! This file is GOCART dry deposition module to calculate the dry deposition
!! velocities of smoke and dust.

module dep_dry_gocart_mod

  use machine ,        only : kind_phys
  use rrfs_smoke_data

  implicit none

  private

  public :: gocart_drydep_driver

CONTAINS

subroutine gocart_drydep_driver(numgas,        &
  moist,p8w,chem,rho_phy,dz8w,ddvel,xland,hfx, &
  ivgtyp,tsk,pbl,ust,znt,g,                    &
  num_moist,num_chem,                          &
  ids,ide, jds,jde, kds,kde,                   &
  ims,ime, jms,jme, kms,kme,                   &
  its,ite, jts,jte, kts,kte                    )

  IMPLICIT NONE
  INTEGER, INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,       &
                            ims,ime, jms,jme, kms,kme,       &
                            num_moist,num_chem,              &
                            its,ite, jts,jte, kts,kte,numgas
  REAL(kind_phys),    INTENT(IN   ) :: g
  REAL(kind_phys),    DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),&
           INTENT(IN   ) :: moist
  REAL(kind_phys),    DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ) ,&
           INTENT(INOUT) :: chem
  REAL(kind_phys),    DIMENSION( ims:ime , kms:kme , jms:jme )         ,&
           INTENT(IN   ) :: dz8w, p8w,rho_phy
  INTEGER, DIMENSION( ims:ime , jms:jme )                   ,&
           INTENT(IN   ) :: ivgtyp
  REAL(kind_phys),    DIMENSION( ims:ime , jms:jme )                   ,&
           INTENT(INOUT) :: tsk,                             &
                            pbl,                             &
                            ust,                             &
                            xland,znt,hfx

!! .. Local Scalars ..

  INTEGER :: iland, iprt, iseason, jce, jcs,  &
             n, nr, ipr, jpr, nvr,            &
             idrydep_onoff,imx,jmx,lmx
  integer :: ii,jj,kk,i,j,k,nv
  integer,            dimension (1,1) :: ilwi
  real(kind_phys), DIMENSION (5)   :: tc,bems
  real(kind_phys), dimension (1,1) :: z0,w10m,gwet,airden,airmas,&
                                         delz_sfc,hflux,ts,pblz,ustar,&
                                         ps,dvel,drydf
  REAL(kind_phys), DIMENSION( its:ite, jts:jte, num_chem ) :: ddvel

  do nv=1,num_chem
    do j=jts,jte
      do i=its,ite
        ddvel(i,j,nv)=0.
      enddo
    enddo
  enddo
  imx=1
  jmx=1
  lmx=1
  do j=jts,jte
  do i=its,ite
     dvel(1,1)=0.
     ilwi(1,1)=0
     if(xland(i,j).gt.1.5)ilwi=1
! for aerosols, ii=1 or ii=2
     ii=1
     if(ivgtyp(i,j).eq.19.or.ivgtyp(i,j).eq.23)ii=1
     airden(1,1)=rho_phy(i,kts,j)
     delz_sfc(1,1)=dz8w(i,kts,j)
     ustar(1,1)=ust(i,j)
     hflux(1,1)=hfx(i,j)
     pblz(1,1)=pbl(i,j)
     ps(1,1)=p8w(i,kts,j)*.01
     z0(1,1)=znt(i,j)
     ts(1,1)=tsk(i,j)

     call depvel_gocart(ii,imx,jmx,lmx,&
     airden, delz_sfc, pblz, ts, ustar, hflux, ilwi, &
     ps, z0, dvel, drydf,g)
     do nv=1,num_chem
      ddvel(i,j,nv)=dvel(1,1)
     enddo
  enddo
  enddo
end subroutine gocart_drydep_driver



SUBROUTINE depvel_gocart(      &
     ii,imx,jmx,lmx,&
     airden, delz_sfc, pblz, ts, ustar, hflux, ilwi, &
     ps, z0, dvel, drydf,g0)

! ****************************************************************************
! *                                                                          *
! *  Calculate dry deposition velocity.                                      *
! *                                                                          *
! *  Input variables:                                                        *
! *    AEROSOL(k)      - Logical, T = aerosol species, F = gas species       *
! *    IREG(i,j)       - # of landtypes in grid square                       *
! *    ILAND(i,j,ldt)  - Land type ID for element ldt =1,IREG(i,j)           *
! *    IUSE(i,j,ldt)   - Fraction of gridbox area occupied by land type      *
! *                      element ldt                                         *
! *    USTAR(i,j)      - Friction velocity (m s-1)                           *
! *    DELZ_SFC(i,j)   - Thickness of layer above surface                    *
! *    PBLZ(i,j)       - Mixing depth (m)                                    *
! *    Z0(i,j)         - Roughness height (m)                                *
! *                                                                          *
! *  Determined in this subroutine (local):                                  *
! *    OBK             - Monin-Obukhov length (m): set to 1.E5 m under       *
! *                      neutral conditions                                  *
! *    Rs(ldt)         - Bulk surface resistance(s m-1) for species k to     * 
! *                      surface ldt                                         *
! *    Ra              - Aerodynamic resistance.                             *
! *    Rb              - Sublayer resistance.                                *
! *    Rs              - Surface resistance.                                 *
! *    Rttl            - Total deposition resistance (s m-1) for species k   *
! *                      Rttl(k) = Ra + Rb + Rs.                             *
! *                                                                          *
! *  Returned:                                                               *
! *    DVEL(i,j,k)     - Deposition velocity (m s-1) of species k            *
! *    DRYDf(i,j,k)    - Deposition frequency (s-1) of species k,            *
! *                    = DVEL / DELZ_SFC                                     *
! *                                                                          *
! **************************************************************************** 


  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: imx,jmx,lmx
  REAL(kind_phys),    INTENT(IN)  :: airden(imx,jmx), delz_sfc(imx,jmx)
  REAL(kind_phys),    INTENT(IN)  :: hflux(imx,jmx), ts(imx,jmx)
  REAL(kind_phys),    INTENT(IN)  :: ustar(imx,jmx), pblz(imx,jmx)
  REAL(kind_phys),    INTENT(IN)  :: ps(imx,jmx)
  INTEGER, INTENT(IN)  :: ilwi(imx,jmx)
  REAL(kind_phys),    INTENT(IN)  :: z0(imx,jmx)
  REAL(kind=kind_phys),    INTENT(IN)  :: g0
  REAL(kind_phys),    INTENT(OUT) :: dvel(imx,jmx), drydf(imx,jmx)
  
  REAL(kind_phys)    :: obk, vds, czh, rttl, frac, logmfrac, psi_h, cz, eps
  REAL(kind_phys)    :: vd, ra, rb, rs  
  INTEGER :: i, j, k, ldt, iolson, ii
  CHARACTER(LEN=50) :: msg
  REAL(kind_phys)     :: prss, tempk, tempc, xnu, ckustr, reyno, aird, diam, xm, z
  REAL(kind_phys)     :: frpath, speed, dg, dw, rt
  REAL(kind_phys)     :: rad0, rix, gfact, gfaci, rdc, rixx, rluxx, rgsx, rclx
  REAL(kind_phys)     :: dtmp1, dtmp2, dtmp3, dtmp4
  REAL(kind_phys)     :: biofit,vk

   psi_h=0.0
  ! executable statements
  j_loop: DO j = 1,jmx               
     i_loop: DO i = 1,imx            
        vk=.4                                     
        vd    = 0.0
        ra    = 0.0
        rb    = 0.0 ! only required for gases (SO2)
        rs = 0.0

! ****************************************************************************
! *  Compute the the Monin-Obhukov length.                                   *
! *  The direct computation of the Monin-Obhukov length is:                  *
! *                                                                          *
! *           - Air density * Cp * T(surface air) * Ustar^3                  *
! *   OBK =   ----------------------------------------------                 *
! *                    vK   * g  * Sensible Heat flux                        *
! *                                                                          *
! *    Cp = 1000 J/kg/K    = specific heat at constant pressure              *
! *    vK = 0.4            = von Karman's constant                           *
! ****************************************************************************

        IF (hflux(i,j) == 0.0) THEN
           obk = 1.0E5
        ELSE
           ! MINVAL(hflux), MINVAL(airden), MINVAL(ustar) =??
           obk = -airden(i,j) * 1000.0 * ts(i,j) * (ustar(i,j))**3 &
                / (vk * g0 * hflux(i,j)) 
! -- debug:
           IF ( obk == 0.0 ) WRITE(*,211) obk, i, j
211        FORMAT(1X,'OBK=', E11.2, 1X,' i,j = ', 2I4)
           
        END IF

        cz = delz_sfc(i,j) / 2.0 ! center of the grid box above surface

! ****************************************************************************
! *  (1) Aerosodynamic resistance Ra and sublayer resistance Rb.             *
! *                                                                          *
! *  The Reynolds number REYNO diagnoses whether a surface is                *
! *  aerodynamically rough (REYNO > 10) or smooth.  Surface is               *
! *  rough in all cases except over water with low wind speeds.              *
! *                                                                          *
! *  For gas species over land and ice (REYNO >= 10) and for aerosol         *
! *  species for all surfaces:                                               *
! *                                                                          *
! *      Ra = 1./VT          (VT from GEOS Kzz at L=1, m/s).                 *
! *                                                                          *
! *  The following equations are from Walcek et al, 1986:                    *
! *                                                                          *
! *  For gas species when REYNO < 10 (smooth), Ra and Rb are combined        *
! *  as Ra:                                                                  *
! *                                                                          *
! *      Ra = { ln(ku* z1/Dg) - Sh } / ku*           eq.(13)                 *
! *                                                                          *
! *      where z1 is the altitude at the center of the lowest model layer    *
! *               (CZ);                                                      *
! *            Sh is a stability correction function;                        *
! *            k  is the von Karman constant (0.4, vK);                      *
! *            u* is the friction velocity (USTAR).                          *
! *                                                                          *
! *   Sh is computed as a function of z1 and L       eq ( 4) and (5)):       *
! *                                                                          *
! *    0 < z1/L <= 1:     Sh = -5 * z1/L                                     *
! *    z1/L < 0:          Sh = exp{ 0.598 + 0.39*ln(E) - 0.09(ln(E))^2 }     *
! *                       where E = min(1,-z1/L) (Balkanski, thesis).        *
! *                                                                          *
! *   For gas species when REYNO >= 10,                                      *
! *                                                                          *
! *      Rb = 2/ku* (Dair/Dg)**(2/3)                 eq.(12)                 *
! *      where Dg is the gas diffusivity, and                                *
! *            Dair is the air diffusivity.                                  *
! *                                                                          *
! *  For aerosol species, Rb is combined with surface resistance as Rs.      *
! *                                                                          *
! ****************************************************************************

        frac = cz / obk
        IF (frac > 1.0) frac = 1.0
        IF (frac > 0.0 .AND. frac <= 1.0) THEN 
           psi_h = -5.0*frac
        ELSE IF (frac < 0.0) THEN
           eps = MIN(1.0D0, -frac)
           logmfrac = LOG(eps)
           psi_h = EXP( 0.598 + 0.39 * logmfrac - 0.09 * (logmfrac)**2 )
        END IF
           !--------------------------------------------------------------
           !  Aerosol species, Rs here is the combination of Rb and Rs.

              ra = (LOG(cz/z0(i,j)) - psi_h) / (vk*ustar(i,j))

           vds = 0.002*ustar(i,j)
           IF (obk < 0.0) &
                vds = vds * (1.0+(-300.0/obk)**0.6667)

           czh  = pblz(i,j)/obk
           IF (czh < -30.0) vds = 0.0009*ustar(i,j)*(-czh)**0.6667

           ! --Set Vds to be less than VDSMAX (entry in input file divided --
           !   by 1.E4). VDSMAX is taken from Table 2 of Walcek et al. [1986].
           !   Invert to get corresponding R
           if(ii.eq.1)then
              rs=1.0/MIN(vds,2.0D-2)
           else
              rs=1.0/MIN(vds,2.0D-3)
           endif
           

        ! ------ Set max and min values for bulk surface resistances ------

           rs= MAX(1.0D0, MIN(rs, 9.9990D+3))

! ****************************************************************************
! *                                                                          *
! *  Compute dry deposition velocity.                                        *
! *                                                                          *
! *  IUSE is the fraction of the grid square occupied by surface ldt in      *
! *  units of per mil (IUSE=500 -> 50% of the grid square).  Add the         *
! *  contribution of surface type ldt to the deposition velocity; this is    *
! *  a loop over all surface types in the gridbox.                           *
! *                                                                          *
! *  Total resistance = Ra + Rb + Rs.
! *                                                                          *
! ****************************************************************************

           rttl = ra + rb + rs
           vd   = vd + 1./rttl

        ! ------ Load array DVEL ------
           dvel(i,j) = vd * 1.2

        ! -- Set a minimum value for DVEL
        !    MIN(VdSO2)      = 2.0e-3 m/s  over ice
        !                    = 3.0e-3 m/s  over land
        !    MIN(vd_aerosol) = 1.0e-4 m/s

           IF (dvel(i,j) < 1.0E-4) dvel(i,j) = 1.0E-4
        drydf(i,j) = dvel(i,j) / delz_sfc(i,j)

     END DO i_loop
  END DO j_loop

END SUBROUTINE depvel_gocart

end module dep_dry_gocart_mod
