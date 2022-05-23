!>\file  dep_simple_mod.F90
!! This file contains the Wesely dry deposition module.

module dep_simple_mod

  use rrfs_smoke_data
  use rrfs_smoke_config, GOCART_SIMPLE => CHEM_OPT_GOCART, chem_opt=>chem_opt
!  use chem_tracers_mod, config_flags => chem_config

! USE module_data_sorgam

  implicit none

!--------------------------------------------------
! .. Default Accessibility ..
!--------------------------------------------------
    PUBLIC


    CONTAINS

SUBROUTINE wesely_driver( data, ktau, dtstep, current_month,  &
                          gmt, julday, t_phy,moist, p8w, t8w, raincv,     &
                          p_phy, chem, rho_phy, dz8w, ddvel, aer_res_def, &
                          aer_res_zcen, ivgtyp, tsk, gsw, vegfra, pbl,    &
                          rmol, ust, znt, xlat, xlong,                    &
                          z, z_at_w, snowh, numgas,                       &
                          ids,ide, jds,jde, kds,kde,                      &
                          ims,ime, jms,jme, kms,kme,                      &
                          its,ite, jts,jte, kts,kte                       )
  implicit none
!--------------------------------------------------
!  Wesely dry dposition driver
!--------------------------------------------------

!  USE module_model_constants 
!  USE module_wrf_control,only:num_moist,num_chem
!  USE module_state_description                       
! USE module_initial_chem_namelists
! USE module_data_sorgam
!  USE module_state_description, only:  param_first_scalar        
  type(smoke_data), intent(inout), pointer :: data
   INTEGER,      INTENT(IN   ) :: julday,                              &
                                  numgas, current_month,                  &
                                  ids,ide, jds,jde, kds,kde,              &
                                  ims,ime, jms,jme, kms,kme,              &
                                  its,ite, jts,jte, kts,kte     
   INTEGER,      INTENT(IN   ) :: ktau            
      REAL(kind_phys),      INTENT(IN   ) :: dtstep,gmt

!--------------------------------------------------
! advected moisture variables
!--------------------------------------------------
   REAL(KIND_PHYS), DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ), INTENT(IN ) :: &
                                                      moist  
!--------------------------------------------------
! advected chemical species
!--------------------------------------------------
   REAL(KIND_PHYS), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), INTENT(INOUT ) :: &
                                                      chem
!--------------------------------------------------
! deposition velocities
!--------------------------------------------------
   REAL(KIND_PHYS), DIMENSION( its:ite, jts:jte, num_chem ), INTENT(INOUT ) ::      &
                                                      ddvel                     
!--------------------------------------------------
! input from met model
!--------------------------------------------------
   REAL(KIND_PHYS),  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) ::      &
                                                      t_phy,              &
                                                      p_phy,              &
                                                      dz8w,               &
                                                      z,                  &
                                                      t8w,                &
                                                      p8w,                &
                                                      z_at_w,             &
                                                      rho_phy
   INTEGER,DIMENSION( ims:ime , jms:jme ), INTENT(IN   ) ::               &
                                                     ivgtyp
   REAL(KIND_PHYS),  DIMENSION( ims:ime , jms:jme ), INTENT(INOUT   ) ::             &
                                                     tsk,                 &
                                                     gsw,                 &
                                                     vegfra,              &
                                                     pbl,                 &
                                                     rmol,                &
                                                     ust,                 &
                                                     xlat,                &
                                                     xlong,               &
                                                     raincv,              &
                                                     znt
   REAL(KIND_PHYS), intent(inout) ::                            aer_res_def(its:ite,jts:jte)
   REAL(KIND_PHYS), intent(inout) ::                            aer_res_zcen(its:ite,jts:jte)
   REAL(KIND_PHYS), INTENT(IN)    ::                            snowh(ims:ime,jms:jme)

!--------------------------------------------------
! .. Local Scalars
!--------------------------------------------------
      REAL(kind_phys)    ::  clwchem, dvfog, dvpart, pa, rad, dep_vap
      REAL(KIND_PHYS) ::  rhchem, ta, ustar, vegfrac, z1, zntt
      INTEGER :: i, iland, iprt, iseason, j, jce, jcs, n, nr, ipr,jpr,nvr
      LOGICAL :: highnh3, rainflag, vegflag, wetflag
!--------------------------------------------------
! .. Local Arrays
!--------------------------------------------------
      REAL(KIND_PHYS) :: p(kts:kte)
      REAL(KIND_PHYS) :: srfres(numgas)
      REAL(KIND_PHYS) :: ddvel0d(numgas)

!-----------------------------------------------------------
! necessary for aerosols (module dependent)         
!-----------------------------------------------------------
      real(kind_phys) :: rcx(numgas)


!-----------------------------------------------------------
! .. Intrinsic Functions
!-----------------------------------------------------------
!      integer :: chem_opt 
 
      INTRINSIC max, min                             

      data => get_thread_smoke_data()

!      chem_opt = chem_opt

      dep_vap  = depo_fact
      !print*,'hli simple chem_opt',chem_opt

!      CALL wrf_debug(15,'in dry_dep_wesely')

         if( julday < 90 .or. julday > 270 ) then
            iseason = 2
!            CALL wrf_debug(15,'setting iseason to 2')
         else
            iseason = 1
         endif


tile_lat_loop : &
      do j = jts,jte 
tile_lon_loop : &
         do i = its,ite 
            iprt  = 0

            iland = luse2usgs( ivgtyp(i,j) )
!--

            ta    = tsk(i,j)
            rad   = gsw(i,j)
            vegfrac = vegfra(i,j)
            pa      = .01*p_phy(i,kts,j)
            clwchem = moist(i,kts,j,p_qc)
            ustar = ust(i,j)
            zntt  = znt(i,j)
            z1    = z_at_w(i,kts+1,j) - z_at_w(i,kts,j)
!-----------------------------------------------------------
!     Set logical default values
!-----------------------------------------------------------
            rainflag = .FALSE.
            wetflag  = .FALSE.
            highnh3  = .FALSE.
!            if(p_qr > 1) then
!               if(moist(i,kts,j,p_qr) > 1.e-18 .or. raincv(i,j) > 0.) then
!                  rainflag = .true.
!               endif
!            endif
            rhchem = MIN( 100.,100. * moist(i,kts,j,p_qv) / &
                     (3.80*exp(17.27*(t_phy(i,kts,j)-273.)/(t_phy(i,kts,j)-36.))/pa))
            rhchem = MAX(5.,RHCHEM)
            if (rhchem >= 95.) wetflag = .true.

!-----------------------------------------------------------
!--- deposition
!-----------------------------------------------------------
!     if(snowc(i,j).gt.0.)iseason=4
            CALL rc( data, rcx, ta, rad, rhchem, iland, &
                     iseason, numgas, wetflag, rainflag, highnh3, &
                     iprt, moist(i,kts,j,p_qv), p8w(i,kts,j) )
               srfres(1:numgas-2) = rcx(1:numgas-2)
               srfres(numgas-1:numgas) = 0.
            CALL deppart( data, rmol(i,j), ustar, rhchem, clwchem, iland, dvpart, dvfog )
            ddvel0d(1:numgas) = 0.
            aer_res_def(i,j)  = 0.
            aer_res_zcen(i,j) = 0.
            CALL landusevg( data, ddvel0d, ustar, rmol(i,j), zntt, z1, dvpart, iland,  &
                            numgas, srfres, aer_res_def(i,j), aer_res_zcen(i,j), p_sulf )

!-----------------------------------------------------------
!wig: CBMZ does not have HO and HO2 last so need to copy all species
!      ddvel(i,j,1:numgas-2)=ddvel0d(1:numgas-2)
!-----------------------------------------------------------
            ddvel(i,j,1:numgas) = ddvel0d(1:numgas)
         end do tile_lon_loop
      end do tile_lat_loop
     
!-----------------------------------------------------------
! For the additional CBMZ species, assign similar RADM counter parts for
! now. Short lived species get a zero velocity since dry dep should be
! unimportant.  **ALSO**, treat p_sulf as h2so4 vapor, not aerosol sulfate
!-----------------------------------------------------------
!

!-----------------------------------------------------------
! For gocartsimple : need msa. On the other hand sulf comes from aerosol routine
!-----------------------------------------------------------
      if  (chem_opt == GOCART_SIMPLE          )   then
         do j=jts,jte
            do i=its,ite
               ddvel(i,j,p_msa)         = ddvel(i,j,p_sulf)
               ddvel(i,j,p_sulf)        = 0.
               ddvel(i,j,p_dms)         = 0.
            end do
         end do
      end if
      
END SUBROUTINE wesely_driver

      SUBROUTINE rc( data, rcx, t, rad, rh, iland, &
                     iseason, numgas, wetflag, rainflag, highnh3, &
                     iprt, spec_hum, p_srf )
!----------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES SURFACE RESISTENCES ACCORDING
!     TO THE MODEL OF
!     M. L. WESELY,
!     ATMOSPHERIC ENVIRONMENT 23 (1989), 1293-1304
!     WITH SOME ADDITIONS ACCORDING TO
!     J. W. ERISMAN, A. VAN PUL, AND P. WYERS,
!     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607
!     WRITTEN BY  WINFRIED SEIDL, APRIL 1997
!     MODYFIED BY WINFRIED SEIDL, MARCH 2000
!                    FOR MM5 VERSION 3
!----------------------------------------------------------------------

!  USE module_state_description                       
!  USE module_initial_chem_namelists
        implicit none
        type(smoke_data), pointer, intent(inout) :: data
!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
        INTEGER, intent(in) :: iland, iseason, numgas
        INTEGER, intent(in) :: iprt
        REAL(KIND_PHYS), intent(in)    :: rad, rh
        REAL(KIND_PHYS), intent(in)    :: t                            ! surface temp (K)
        REAL(KIND_PHYS), intent(in)    :: p_srf                        ! surface pressure (Pa)
        REAL(KIND_PHYS), intent(in)    :: spec_hum                     ! surface specific humidity (kg/kg)
        real(kind_phys), intent(out)   :: rcx(numgas)
        LOGICAL, intent(in) :: highnh3, rainflag, wetflag

!----------------------------------------------------------------------
! .. Local Scalars ..
!----------------------------------------------------------------------
        REAL(KIND_PHYS), parameter :: t0    = 298.
        REAL(KIND_PHYS), parameter :: tmelt = 273.16
        INTEGER :: lt, n
        INTEGER :: chem_opt
        REAL(KIND_PHYS) :: rclx, rdc, resice, rgsx, rluo1, rluo2
        REAL(KIND_PHYS) :: rlux, rmx, rs, rsmx, rdtheta, z, wrk
        REAL(KIND_PHYS) :: qs, es, ws, dewm, dv_pan, drat
        REAL(KIND_PHYS) :: crs, tc
        REAL(KIND_PHYS) :: rs_pan, tc_pan
        LOGICAL :: has_dew
!----------------------------------------------------------------------
! .. Local Arrays ..
!----------------------------------------------------------------------
        REAL(KIND_PHYS) :: hstary(numgas)

!----------------------------------------------------------------------
! .. Intrinsic Functions ..
!----------------------------------------------------------------------
        INTRINSIC exp

        chem_opt = chem_opt

        rcx(1:numgas) = 1.

        tc = t - 273.15
        rdtheta = 0.

        z = 200./(rad+0.1)

!!!  HARDWIRE VALUES FOR TESTING
!       z=0.4727409
!       tc=22.76083
!       t=tc+273.15
!       rad = 412.8426
!       rainflag=.false.
!       wetflag=.false.

        IF ( tc<=0. .OR. tc>=40. ) THEN
          rs = 9999.
        ELSE
          rs = data%ri(iland,iseason)*(1+z*z)*(400./(tc*(40.-tc)))
        END IF
        rdc   = 100.*(1. + 1000./(rad + 10.))/(1. + 1000.*rdtheta)
        rluo1 = 1./(1./3000. + 3./data%rlu(iland,iseason))
        rluo2 = 1./(1./1000. + 3./data%rlu(iland,iseason))
        resice = 1000.*exp( -(tc + 4.) )
        wrk    = (t0 - t)/(t0*t)


        DO n = 1, numgas
          IF( data%hstar(n) /= 0. ) then
             hstary(n) = data%hstar(n)*exp( data%dhr(n)*wrk )
!----------------------------------------------------------------------
!     SPECIAL TREATMENT FOR HNO3, HNO4, H2O2, PAA
!----------------------------------------------------------------------
             rmx = 1./(hstary(n)/3000. + 100.*data%f0(n))
             rsmx = rs*data%dratio(n) + rmx
             rclx = 1./(1.e-5*hstary(n)/data%rcls(iland,iseason) &
                        + data%f0(n)/data%rclo(iland,iseason)) + resice
             rgsx = 1./(1.e-5*hstary(n)/data%rgss(iland,iseason) &
                        + data%f0(n)/data%rgso(iland,iseason)) + resice
             rlux = data%rlu(iland,iseason)/(1.e-5*hstary(n) + data%f0(n)) + resice
             IF( wetflag ) THEN
               rlux = 1./(1./(3.*data%rlu(iland,iseason)) + 1.e-7*hstary(n) + data%f0(n)/rluo1)
             END IF
             IF( rainflag ) THEN
               rlux = 1./(1./(3.*data%rlu(iland,iseason)) + 1.e-7*hstary(n) + data%f0(n)/rluo2)
             END IF
             rcx(n) = 1./(1./rsmx + 1./rlux + 1./(rdc + rclx) + 1./(data%rac(iland,iseason) + rgsx))
             rcx(n) = max( 1.,rcx(n) )
          end IF
        END DO

!--------------------------------------------------
!     SPECIAL TREATMENT FOR OZONE
!--------------------------------------------------
!     SPECIAL TREATMENT FOR SO2 (Wesely)
!       HSTARY(P_SO2)=DATA%HSTAR(P_SO2)*EXP(DATA%DHR(P_SO2)*(1./T-1./298.))
!       RMX=1./(HSTARY(P_SO2)/3000.+100.*DATA%F0(P_SO2))
!       RSMX=RS*DATA%DRATIO(P_SO2)+RMX
!       RLUX=DATA%RLU(ILAND,ISEASON)/(1.E-5*HSTARY(P_SO2)+DATA%F0(P_SO2))
!    &       +RESICE
!       RCLX=DATA%RCLS(ILAND,ISEASON)+RESICE
!       RGSX=DATA%RGSS(ILAND,ISEASON)+RESICE
!       IF ((wetflag).OR.(RAINFLAG)) THEN
!         IF (ILAND.EQ.1) THEN
!           RLUX=50.
!         ELSE
!           RLUX=100.
!         END IF
!       END IF
!       RCX(P_SO2)=1./(1./RSMX+1./RLUX+1./(RDC+RCLX)
!    &                +1./(DATA%RAC(ILAND,ISEASON)+RGSX))
!       IF (RCX(P_SO2).LT.1.) RCX(P_SO2)=1.

!--------------------------------------------------
!     SO2 according to Erisman et al. 1994
!       R_STOM
!--------------------------------------------------
is_so2 : &
     if( p_so2 > 1 ) then
        rsmx = rs*data%dratio(p_so2)
!--------------------------------------------------
!       R_EXT
!--------------------------------------------------
        IF (tc> -1. ) THEN
          IF (rh<81.3) THEN
            rlux = 25000.*exp(-0.0693*rh)
          ELSE
            rlux = 0.58E12*exp(-0.278*rh)
          END IF
        END IF
        IF (((wetflag) .OR. (rainflag)) .AND. (tc> -1. )) THEN
          rlux = 1.
        END IF
        IF ((tc>= -5. ) .AND. (tc<= -1. )) THEN
          rlux = 200.
        END IF
        IF (tc< -5. ) THEN
          rlux = 500.
        END IF
!--------------------------------------------------
!       INSTEAD OF R_INC R_CL and R_DC of Wesely are used
!--------------------------------------------------
        rclx = data%rcls(iland,iseason)
!--------------------------------------------------
!       DRY SURFACE
!--------------------------------------------------
        rgsx = 1000.
!--------------------------------------------------
!       WET SURFACE
!--------------------------------------------------
        IF ((wetflag) .OR. (rainflag)) THEN
          IF (highnh3) THEN
            rgsx = 0.
          ELSE
            rgsx = 500.
          END IF
        END IF
!--------------------------------------------------
!       WATER
!--------------------------------------------------
        IF (iland==iswater_temp) THEN
          rgsx = 0.
        END IF
!--------------------------------------------------
!       SNOW
!--------------------------------------------------
        IF( iseason==4 .OR. iland==isice_temp ) THEN
          IF( tc > 2. ) THEN
            rgsx = 0.
          else IF ( tc >= -1. .AND. tc <= 2. ) THEN
            rgsx = 70.*(2. - tc)
          else IF ( tc < -1. ) THEN
            rgsx = 500.
          END IF
        END IF
!--------------------------------------------------
!       TOTAL SURFACE RESISTENCE
!--------------------------------------------------
        IF ((iseason/=4) .AND. (data%ixxxlu(iland)/=1) .AND. (iland/=iswater_temp) .AND. &
            (iland/=isice_temp)) THEN
          rcx(p_so2) = 1./(1./rsmx+1./rlux+1./(rclx+rdc+rgsx))
        ELSE
          rcx(p_so2) = rgsx
        END IF
        rcx(p_so2) = max( 1.,rcx(p_so2) )
     end if is_so2
!--------------------------------------------------
!     NH3 according to Erisman et al. 1994
!       R_STOM
!--------------------------------------------------
      END SUBROUTINE rc

      SUBROUTINE deppart( data, rmol, ustar, rh, clw, iland, &
                          dvpart, dvfog )
!--------------------------------------------------
!     THIS SUBROUTINE CALCULATES SURFACE DEPOSITION VELOCITIES
!     FOR FINE AEROSOL PARTICLES ACCORDING TO THE MODEL OF
!     J. W. ERISMAN, A. VAN PUL, AND P. WYERS,
!     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607
!     WRITTEN BY WINFRIED SEIDL, APRIL 1997
!     MODIFIED BY WINFRIED SEIDL, MARCH 2000
!            FOR MM5 VERSION 3
!--------------------------------------------------
        implicit none
        type(smoke_data), pointer, intent(inout) :: data

!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        INTEGER, intent(in) :: iland
        REAL(KIND_PHYS), intent(in)    :: clw, rh, rmol, ustar
        REAL(KIND_PHYS), intent(out)   :: dvfog, dvpart

!--------------------------------------------------
! .. Intrinsic Functions ..
!--------------------------------------------------
        INTRINSIC exp

        dvpart = ustar/data%kpart(iland)
        IF (rmol<0.) THEN
!--------------------------------------------------
!         UNSTABLE LAYERING CORRECTION
!--------------------------------------------------
          dvpart = dvpart*(1.+(-300.*rmol)**0.66667)
        END IF
        IF (rh>80.) THEN
!--------------------------------------------------
!         HIGH RELATIVE HUMIDITY CORRECTION
!         ACCORDING TO J. W. ERISMAN ET AL.
!         ATMOSPHERIC ENVIRONMENT 31 (1997), 321-332
!--------------------------------------------------
          dvpart = dvpart*(1.+0.37*exp((rh-80.)/20.))
        END IF

!--------------------------------------------------
!       SEDIMENTATION VELOCITY OF FOG WATER ACCORDING TO
!       R. FORKEL, W. SEIDL, R. DLUGI AND E. DEIGELE
!       J. GEOPHYS. RES. 95D (1990), 18501-18515
!--------------------------------------------------
        dvfog = 0.06*clw
        IF (data%ixxxlu(iland)==5) THEN
!--------------------------------------------------
!         TURBULENT DEPOSITION OF FOG WATER IN CONIFEROUS FOREST ACCORDI
!         A. T. VERMEULEN ET AL.
!         ATMOSPHERIC ENVIRONMENT 31 (1997), 375-386
!--------------------------------------------------
          dvfog = dvfog + 0.195*ustar*ustar
        END IF

      END SUBROUTINE deppart

      SUBROUTINE landusevg( data, vgs, ustar, rmol, z0, zz, &
                            dvparx, iland, numgas, srfres, aer_res_def, &
                            aer_res_zcen, p_sulf )
!--------------------------------------------------
!     This subroutine calculates the species specific deposition velocit
!     as a function of the local meteorology and land use.  The depositi
!     Velocity is also landuse specific.
!     Reference: Hsieh, C.M., Wesely, M.L. and Walcek, C.J. (1986)
!                A Dry Deposition Module for Regional Acid Deposition
!                EPA report under agreement DW89930060-01
!     Revised version by Darrell Winner (January 1991)
!        Environmental Engineering Science 138-78
!           California Institute of Technology
!              Pasadena, CA  91125
!     Modified by Winfried Seidl (August 1997)
!       Fraunhofer-Institut fuer Atmosphaerische Umweltforschung
!                    Garmisch-Partenkirchen, D-82467
!          for use of Wesely and Erisman surface resistances
!     Inputs:
!        Ustar  : The grid average friction velocity (m/s)
!        Rmol   : Reciprocal of the Monin-Obukhov length (1/m)
!        Z0     : Surface roughness height for the grid square (m)
!        SrfRes : Array of landuse/atmospheric/species resistances (s/m)
!        Slist  : Array of chemical species codes
!        Dvparx : Array of surface deposition velocity of fine aerosol p
!     Outputs:
!        Vgs    : Array of species and landuse specific deposition
!                 velocities (m/s)
!        Vg     : Cell-average deposition velocity by species (m/s)
!     Variables used:
!        SCPR23  : (Schmidt #/Prandtl #)**(2/3) Diffusion correction fac
!        Zr      : Reference Height (m)
!        Iatmo   : Parameter specifying the stabilty class (Function of
!        Z0      : Surface roughness height (m)
!        karman  : Von Karman constant (from module_model_constants)
!--------------------------------------------------

!        USE module_model_constants, only: karman
        implicit none

        type(smoke_data), pointer, intent(inout) :: data

!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        INTEGER, intent(in) :: iland, numgas, p_sulf
        REAL(KIND_PHYS), intent(in)    :: dvparx, ustar, z0, zz
        REAL(KIND_PHYS), intent(inout) :: rmol
        REAL(KIND_PHYS), intent(inout) :: aer_res_def
        REAL(KIND_PHYS), intent(inout) :: aer_res_zcen
!--------------------------------------------------
! .. Array Arguments ..
!--------------------------------------------------
        REAL(KIND_PHYS), intent(in)  :: srfres(numgas)
        REAL(KIND_PHYS), intent(out) :: vgs(numgas)

!--------------------------------------------------
! .. Local Scalars ..
!--------------------------------------------------
        INTEGER :: jspec
        REAL(KIND_PHYS) :: vgp, vgpart, zr
        REAL(KIND_PHYS) :: rmol_tmp
!--------------------------------------------------
! .. Local Arrays ..
!--------------------------------------------------
        REAL(KIND_PHYS) :: vgspec(numgas)

!--------------------------------------------------
!   Calculate aerodynamic resistance for reference
!   height = layer center
!--------------------------------------------------
        zr = zz*.5
        rmol_tmp = rmol
        CALL depvel( data, numgas, rmol_tmp, zr, z0, ustar, &
                     vgspec, vgpart, aer_res_zcen )
!--------------------------------------------------
!   Set the reference height (2.0 m)
!--------------------------------------------------
!       zr = 10.0
        zr = 2.0

!--------------------------------------------------
!   CALCULATE THE DEPOSITION VELOCITY without any surface
!   resistance term, i.e. 1 / (ra + rb)
!--------------------------------------------------
        CALL depvel( data, numgas, rmol, zr, z0, ustar, &
                     vgspec, vgpart, aer_res_def )

!--------------------------------------------------
!   Calculate the deposition velocity for each species
!   and grid cell by looping through all the possibile combinations
!   of the two
!--------------------------------------------------
        vgp = 1.0/((1.0/vgpart)+(1.0/dvparx))
!--------------------------------------------------
!   Loop through the various species
!--------------------------------------------------
        DO jspec = 1, numgas
!--------------------------------------------------
!   Add in the surface resistance term, rc (SrfRes)
!--------------------------------------------------
          vgs(jspec) = 1.0/(1.0/vgspec(jspec) + srfres(jspec))
        END DO
        vgs(p_sulf) = vgp

        CALL cellvg( data, vgs, ustar, zz, zr, rmol, numgas )

      END SUBROUTINE landusevg

      SUBROUTINE cellvg( data, vgtemp, ustar, dz, zr, rmol, nspec )
!--------------------------------------------------
!     THIS PROGRAM HAS BEEN DESIGNED TO CALCULATE THE CELL AVERAGE
!     DEPOSITION VELOCITY GIVEN THE VALUE OF VG AT SOME REFERENCE
!     HEIGHT ZR WHICH IS MUCH SMALLER THAN THE CELL HEIGHT DZ.
!       PROGRAM WRITTEN BY GREGORY J.MCRAE (NOVEMBER 1977)
!         Modified by Darrell A. Winner    (February 1991)
!.....PROGRAM VARIABLES...
!     VgTemp   - DEPOSITION VELOCITY AT THE REFERENCE HEIGHT
!     USTAR    - FRICTION VELOCITY
!     RMOL     - RECIPROCAL OF THE MONIN-OBUKHOV LENGTH
!     ZR       - REFERENCE HEIGHT
!     DZ       - CELL HEIGHT
!     CELLVG   - CELL AVERAGE DEPOSITION VELOCITY
!     VK       - VON KARMAN CONSTANT
!--------------------------------------------------

!        USE module_model_constants, only: karman
        implicit none
        type(smoke_data), pointer, intent(inout) :: data

!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        INTEGER, intent(in) :: nspec
        REAL(KIND_PHYS), intent(in)    :: dz, rmol, ustar, zr
!--------------------------------------------------
! .. Array Arguments ..
!--------------------------------------------------
        REAL(KIND_PHYS), intent(out) :: vgtemp(nspec)
!--------------------------------------------------
! .. Local Scalars ..
!--------------------------------------------------
        INTEGER :: nss
        REAL(KIND_PHYS) :: a, fac, pdz, pzr, vk
!--------------------------------------------------
! .. Intrinsic Functions ..
!--------------------------------------------------
        INTRINSIC alog, sqrt

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
        DO nss = 1, nspec
          IF (rmol < 0.) THEN
            pdz = sqrt(1.0 - 9.0*dz*rmol)
            pzr = sqrt(1.0 - 9.0*zr*rmol)
            fac = ((pdz - 1.0)/(pzr - 1.0))*((pzr + 1.0)/(pdz + 1.0))
            a   = 0.74*dz*alog(fac) + (0.164/rmol)*(pdz-pzr)
          ELSE IF (rmol == 0.) THEN
            a = 0.74*(dz*alog(dz/zr) - dz + zr)
          ELSE
            a = 0.74*(dz*alog(dz/zr) - dz + zr) + (2.35*rmol)*(dz - zr)**2
          END IF
!--------------------------------------------------
!     CALCULATE THE DEPOSITION VELOCITIY
!--------------------------------------------------
          vgtemp(nss) = vgtemp(nss)/(1.0 + vgtemp(nss)*a/(vk*ustar*(dz - zr)))
        END DO

      END SUBROUTINE cellvg

      SUBROUTINE depvel( data, numgas, rmol, zr, z0, ustar, &
                         depv, vgpart, aer_res )
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
!     SCPR23   - (Schmidt #/Prandtl #)**(2/3) Diffusion correction fact
!     UBAR     - ABSOLUTE VALUE OF SURFACE WIND SPEED
!     DEPVEL   - POLLUTANT DEPOSITION VELOCITY
!     Vk       - VON KARMAN CONSTANT
!     USTAR    - FRICTION VELOCITY U*
!     POLINT   - POLLUTANT INTEGRAL
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

!        USE module_model_constants, only: karman
        implicit none
        type(smoke_data), pointer, intent(inout) :: data

!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        INTEGER, intent(in) :: numgas
        REAL(KIND_PHYS), intent(in)    :: ustar, z0, zr
        REAL(KIND_PHYS), intent(out)   :: vgpart, aer_res
        REAL(KIND_PHYS), intent(inout) :: rmol
!--------------------------------------------------
! .. Array Arguments ..
!--------------------------------------------------
        REAL(KIND_PHYS), intent(out) :: depv(numgas)
!--------------------------------------------------
! .. Local Scalars ..
!--------------------------------------------------
        INTEGER :: l
        REAL(KIND_PHYS) :: ao, ar, polint, vk
!--------------------------------------------------
! .. Intrinsic Functions ..
!--------------------------------------------------
        INTRINSIC alog
!--------------------------------------------------
!     Set the von Karman constant
!--------------------------------------------------
        vk = karman

!--------------------------------------------------
!     Calculate the diffusion correction factor
!     SCPR23 is calculated as (Sc/Pr)**(2/3) using Sc= 1.15 and Pr= 1.0
!     DATA%SCPR23 = 1.10
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
          polint = 0.74*alog(zr/z0) + 4.7*rmol*(zr-z0)
        END IF

!--------------------------------------------------
!     CALCULATE THE Maximum DEPOSITION VELOCITY
!--------------------------------------------------
        DO l = 1, numgas
          depv(l) = ustar*vk/(2.0*data%scpr23(l)+polint)
        END DO
        vgpart = ustar*vk/polint
        aer_res = polint/(karman*max(ustar,1.0e-4))

      END SUBROUTINE depvel

      ! NOTE: dep_init is now in rrfs_smoke_data

end module dep_simple_mod
