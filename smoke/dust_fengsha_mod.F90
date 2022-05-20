!>\file  dust_fengsha_mod.F90
!! This file contains the FENGSHA dust scheme.

module dust_fengsha_mod
!
!  This module developed by Barry Baker (NOAA ARL)
!  For serious questions contact barry.baker@noaa.gov
!
!  07/16/2019 - Adapted for NUOPC/GOCART, R. Montuoro
!  02/01/2020 - Adapted for FV3/CCPP, Haiqin Li

  use rrfs_smoke_data
  use machine ,        only : kind_phys
  use dust_data_mod

  implicit none

  private

  public :: gocart_dust_fengsha_driver

contains

  subroutine gocart_dust_fengsha_driver(data, dt,           &
       chem,rho_phy,smois,p8w,ssm,                       &
       isltyp,vegfra,snowh,xland,area,g,emis_dust,       &
       ust,znt,clay,sand,rdrag,uthr,                     &
       num_emis_dust,num_moist,num_chem,num_soil_layers, &
       ids,ide, jds,jde, kds,kde,                        &
       ims,ime, jms,jme, kms,kme,                        &
       its,ite, jts,jte, kts,kte)
    IMPLICIT NONE
    type(smoke_data), intent(inout) :: data
    INTEGER,      INTENT(IN   ) ::                       &
         ids,ide, jds,jde, kds,kde,                      &
         ims,ime, jms,jme, kms,kme,                      &
         its,ite, jts,jte, kts,kte,                      &
         num_emis_dust,num_moist,num_chem,num_soil_layers
    INTEGER,DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: isltyp
    REAL(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), INTENT(INOUT) :: chem
    REAL(kind_phys), DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL, INTENT(INOUT) :: emis_dust
    REAL(kind_phys), DIMENSION( ims:ime, num_soil_layers, jms:jme ), INTENT(IN) :: smois
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: ssm
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: vegfra,     &
                                                        snowh,      &
                                                        xland,      &
                                                        area,       &
                                                        ust,        &
                                                        znt,        &
                                                        clay,       &
                                                        sand,       &
                                                        rdrag,      &
                                                        uthr
    REAL(kind_phys), DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: &
         p8w,             &
         rho_phy
    REAL(kind_phys), INTENT(IN) :: dt,g

    ! Local variables

    integer :: nmx,smx,i,j,k,imx,jmx,lmx
    integer,dimension (1,1) :: ilwi
    real(kind_phys), DIMENSION (1,1) :: erodtot
    REAL(kind_phys), DIMENSION (1,1) :: gravsm
    REAL(kind_phys), DIMENSION (1,1) :: drylimit
    real(kind_phys), DIMENSION (5)   :: tc,bems
    real(kind_phys), dimension (1,1) :: airden,airmas,ustar
    real(kind_phys), dimension (1) :: dxy
    real(kind_phys), dimension (3) :: massfrac
    real(kind_phys) :: conver,converi
    real(kind_phys) :: R

    ! threshold values
    conver=1.e-9
    converi=1.e9

    ! Number of dust bins

    imx=1
    jmx=1
    lmx=1
    nmx=ndust
    smx=nsalt

    k=kts
    do j=jts,jte
       do i=its,ite

          ! Don't do dust over water!!!

          ilwi(1,1)=0
          if(xland(i,j).lt.1.5)then
             ilwi(1,1)=1

             ! Total concentration at lowest model level. This is still hardcoded for 5 bins.

             !    if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
             !       tc(:)=1.e-16*conver
             !    else
             tc(1)=chem(i,kts,j,p_dust_1)*conver
             tc(2)=chem(i,kts,j,p_dust_2)*conver
             tc(3)=chem(i,kts,j,p_dust_3)*conver
             tc(4)=chem(i,kts,j,p_dust_4)*conver
             tc(5)=chem(i,kts,j,p_dust_5)*conver
             !    endif

             ! Air mass and density at lowest model level.

             airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g
             airden(1,1)=rho_phy(i,kts,j)
             ustar(1,1)=ust(i,j)
             dxy(1)=area(i,j)

             ! Mass fractions of clay, silt, and sand.
             massfrac(1)=clay(i,j)
             massfrac(2)=1-(clay(i,j)+sand(i,j))
             massfrac(3)=sand(i,j)


             ! Total erodibility.
             
             erodtot(1,1) = ssm(i,j) ! SUM(erod(i,j,:))
             
             ! Don't allow roughness lengths greater than 20 cm to be lofted.
             ! This kludge accounts for land use types like urban areas and
             ! forests which would otherwise show up as high dust emitters.
             ! This is a placeholder for a more widely accepted kludge
             ! factor in the literature, which reduces lofting for rough areas.
             ! Forthcoming...

             IF (znt(i,j) .gt. 0.2) then
                ilwi(1,1)=0
             endif

             ! limit where there is lots of vegetation
             if (vegfra(i,j) .gt. .17) then
                ilwi(1,1) = 0
             endif

             ! limit where there is snow on the ground
             if (snowh(i,j) .gt. 0) then
                ilwi(1,1) = 0
             endif

             ! Do not allow areas with bedrock, lava, or land-ice to loft

             IF (isltyp(i,j) .eq. 15 .or. isltyp(i,j) .eq. 16. .or. &
                  isltyp(i,j) .eq. 18) then
                ilwi(1,1)=0
             ENDIF
             IF (isltyp(i,j) .eq. 0)then
                ilwi(1,1)=0
             endif
             if(ilwi(1,1) == 0 ) cycle

             ! Calculate gravimetric soil moisture and drylimit.
             gravsm(1,1)=100.*smois(i,1,j)/((1.-maxsmc(isltyp(i,j)))*(2.65*(1.-clay(i,j))+2.50*clay(i,j)))
             drylimit(1,1)=14.0*clay(i,j)*clay(i,j)+17.0*clay(i,j)

             ! get drag partition
             ! FENGSHA uses the drag partition correction of MacKinnon et al 2004
             !     doi:10.1016/j.geomorph.2004.03.009
             if (dust_calcdrag .ne. 1) then
                call fengsha_drag(data,znt(i,j),R)
             else
                ! use the precalculated version derived from ASCAT; Prigent et al. (2012,2015)
                ! doi:10.1109/TGRS.2014.2338913 & doi:10.5194/amt-5-2703-2012
                ! pick only valid values
                if (rdrag(i,j) > 0.) then
                  R = real(rdrag(i,j), kind=kind_phys)
                else
                  cycle
                endif
             endif

             ! Call dust emission routine.
             call source_dust(data, imx, jmx, lmx, nmx, smx, dt, tc, ustar, massfrac, & 
                  erodtot, dxy, gravsm, airden, airmas, &
                  bems, g, drylimit, dust_alpha, dust_gamma, R, uthr(i,j))

             !    if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
             !     dustin(i,j,1:5)=tc(1:5)*converi
             !    else
             chem(i,kts,j,p_dust_1)=tc(1)*converi
             chem(i,kts,j,p_dust_2)=tc(2)*converi
             chem(i,kts,j,p_dust_3)=tc(3)*converi
             chem(i,kts,j,p_dust_4)=tc(4)*converi
             chem(i,kts,j,p_dust_5)=tc(5)*converi
             !    endif

             !     chem(i,kts,j,p_dust_1)=tc(1)*converi
             !     chem(i,kts,j,p_dust_2)=tc(2)*converi
             !     chem(i,kts,j,p_dust_3)=tc(3)*converi
             !     chem(i,kts,j,p_dust_4)=tc(4)*converi
             !     chem(i,kts,j,p_dust_5)=tc(5)*converi

             ! For output diagnostics

             emis_dust(i,1,j,p_edust1)=bems(1)
             emis_dust(i,1,j,p_edust2)=bems(2)
             emis_dust(i,1,j,p_edust3)=bems(3)
             emis_dust(i,1,j,p_edust4)=bems(4)
             emis_dust(i,1,j,p_edust5)=bems(5)
          endif
       enddo
    enddo
    !

  end subroutine gocart_dust_fengsha_driver


  SUBROUTINE source_dust(data, imx, jmx, lmx, nmx, smx, dt1, tc, ustar, massfrac, &
       erod, dxy, gravsm, airden, airmas, bems, g0, drylimit, alpha,  &
       gamma, R, uthres)

    ! ****************************************************************************
    ! *  Evaluate the source of each dust particles size bin by soil emission
    ! *
    ! *  Input:
    ! *         EROD      Fraction of erodible grid cell                (-)
    ! *         GRAVSM    Gravimetric soil moisture                     (g/g)
    ! *         DRYLIMIT  Upper GRAVSM limit for air-dry soil           (g/g)
    ! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
    ! *         GAMMA     Tuning constant for erodibility               (-)
    ! *         DXY       Surface of each grid cell                     (m2)
    ! *         AIRMAS    Mass of air for each grid box                 (kg)
    ! *         AIRDEN    Density of air for each grid box              (kg/m3)
    ! *         USTAR     Friction velocity                             (m/s)
    ! *         DT1       Time step                                     (s)
    ! *         NMX       Number of dust bins                           (-)
    ! *         SMX       Number of saltation bins                      (-)
    ! *         IMX       Number of I points                            (-)
    ! *         JMX       Number of J points                            (-)
    ! *         LMX       Number of L points                            (-)
    ! *         R         Drag Partition                                (-)
    ! *         UTHRES    FENGSHA Dry Threshold Velocities              (m/s)
    ! *
    ! *  Data:
    ! *         MASSFRAC  Fraction of mass in each of 3 soil classes    (-)
    ! *         SPOINT    Pointer to 3 soil classes                     (-)
    ! *         DEN_DUST  Dust density                                  (kg/m3)
    ! *         DEN_SALT  Saltation particle density                    (kg/m3)
    ! *         REFF_SALT Reference saltation particle diameter         (m)
    ! *         REFF_DUST Reference dust particle diameter              (m)
    ! *         LO_DUST   Lower diameter limits for dust bins           (m)
    ! *         UP_DUST   Upper diameter limits for dust bins           (m)
    ! *         FRAC_SALT Soil class mass fraction for saltation bins   (-)
    ! *
    ! *  Parameters:
    ! *         CMB       Constant of proportionality                   (-)
    ! *         MMD_DUST  Mass median diameter of dust                  (m)
    ! *         GSD_DUST  Geometric standard deviation of dust          (-)
    ! *         LAMBDA    Side crack propagation length                 (m)
    ! *         CV        Normalization constant                        (-)
    ! *         G0        Gravitational acceleration                    (m/s2)
    ! *         G         Gravitational acceleration in cgs             (cm/s2)
    ! *
    ! *  Working:
    ! *         U_TS0     "Dry" threshold friction velocity             (m/s)
    ! *         U_TS      Moisture-adjusted threshold friction velocity (m/s)
    ! *         RHOA      Density of air in cgs                         (g/cm3)
    ! *         DEN       Dust density in cgs                           (g/cm3)
    ! *         DIAM      Dust diameter in cgs                          (cm)
    ! *         DMASS     Saltation mass distribution                   (-)
    ! *         DSURFACE  Saltation surface area per unit mass          (m2/kg)
    ! *         DS_REL    Saltation surface area distribution           (-)
    ! *         SALT      Saltation flux                                (kg/m/s)
    ! *         DLNDP     Dust bin width                                (-)
    ! *         EMIT      Total vertical mass flux                      (kg/m2/s)
    ! *         EMIT_VOL  Total vertical volume flux                    (m/s)
    ! *         DSRC      Mass of emitted dust               (kg/timestep/cell)
    ! *
    ! *  Output:
    ! *         TC        Total concentration of dust        (kg/kg/timestep/cell)
    ! *         BEMS      Source of each dust type           (kg/timestep/cell)
    ! *
    ! ****************************************************************************
    implicit none
    type(smoke_data), intent(inout) :: data

    INTEGER,            INTENT(IN)    :: imx,jmx,lmx,nmx,smx
    REAL(kind_phys), INTENT(IN)    :: dt1
    REAL(kind_phys), INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
    REAL(kind_phys), INTENT(IN)    :: ustar(imx,jmx)
    REAL(kind_phys), INTENT(IN)    :: massfrac(3)
    REAL(kind_phys), INTENT(IN)    :: erod(imx,jmx)
    REAL(kind_phys), INTENT(IN)    :: dxy(jmx)
    REAL(kind_phys), INTENT(IN)    :: gravsm(imx,jmx)
    REAL(kind_phys), INTENT(IN)    :: airden(imx,jmx,lmx)
    REAL(kind_phys), INTENT(IN)    :: airmas(imx,jmx,lmx)
    REAL(kind_phys), INTENT(OUT)   :: bems(imx,jmx,nmx)
    REAL(kind_phys), INTENT(IN)    :: g0
    REAL(kind_phys), INTENT(IN)    :: drylimit(imx,jmx)
    !! Sandblasting mass efficiency, aka "fudge factor" (based on Tegen et al,
    !! 2006 and Hemold et al, 2007)
    !
    !  REAL, PARAMETER :: alpha=1.8E-8  ! (m^-1)
    REAL(kind_phys), INTENT(IN)    :: alpha
    ! Experimental optional exponential tuning constant for erodibility.
    ! 0 < gamma < 1 -> more relative impact by low erodibility regions.
    REAL(kind_phys), INTENT(IN)    :: gamma
    REAL(kind_phys), INTENT(IN)    :: R
    REAL(kind_phys), INTENT(IN)    :: uthres

    REAL(kind_phys)    :: den(smx), diam(smx)
    REAL(kind_phys)    :: dvol(nmx), distr_dust(nmx), dlndp(nmx)
    REAL(kind_phys)    :: dsurface(smx), ds_rel(smx)
    REAL(kind_phys)    :: u_ts0, u_ts, dsrc, dmass, dvol_tot
    REAL(kind_phys)    :: salt,emit, emit_vol, stotal
    REAL(kind_phys)    :: rhoa, g
    INTEGER   :: i, j, n

    ! Sandblasting mass efficiency, beta.
    ! Beta maxes out for clay fractions above 0.2 = betamax.

    REAL(kind_phys), PARAMETER :: betamax=5.25E-4
    REAL(kind_phys) :: beta
    integer :: styp

    ! Constant of proportionality from Marticorena et al, 1997 (unitless)
    ! Arguably more ~consistent~ fudge than alpha, which has many walnuts
    ! sprinkled throughout the literature. - GC

    REAL(kind_phys), PARAMETER :: cmb=1.0
    ! REAL, PARAMETER :: cmb=2.61   ! from White,1979

    ! Parameters used in Kok distribution function. Advise not to play with
    ! these without the expressed written consent of someone who knows what
    ! they're doing. - GC

    REAL(kind_phys), PARAMETER :: mmd_dust=3.4D-6  ! median mass diameter (m)
    REAL(kind_phys), PARAMETER :: gsd_dust=3.0     ! geom. std deviation
    REAL(kind_phys), PARAMETER :: lambda=12.0D-6   ! crack propagation length (m)
    REAL(kind_phys), PARAMETER :: cv=12.62D-6      ! normalization constant

    ! Calculate saltation surface area distribution from sand, silt, and clay
    ! mass fractions and saltation bin fraction. This will later become a
    ! modifier to the total saltation flux.  The reasoning here is that the
    ! size and availability of saltators affects saltation efficiency. Based
    ! on Eqn. (32) in Marticorena & Bergametti, 1995 (hereon, MB95).

    DO n=1,smx
       dmass=massfrac(spoint(n))*frac_salt(n)
       dsurface(n)=0.75*dmass/(den_salt(n)*reff_salt(n))
    ENDDO

    ! The following equation yields relative surface area fraction.  It will only
    ! work if you are representing the "full range" of all three soil classes.
    ! For this reason alone, we have incorporated particle sizes that encompass
    ! the clay class, to account for the its relative area over the basal
    ! surface, even though these smaller bins would be unlikely to play any large
    ! role in the actual saltation process. - GC

    stotal=SUM(dsurface(:))
    DO n=1,smx
       ds_rel(n)=dsurface(n)/stotal
    ENDDO

    ! Calculate total dust emission due to saltation of sand sized particles.
    ! Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
    ! u_ts0 for moisture to get threshold friction velocity (u_ts). Then
    ! calculate saltation flux (salt) where ustar has exceeded u_ts.  Finally,
    ! calculate total dust emission (tot_emit), taking into account erodibility.

    ! Set DRY threshold friction velocity to input value
    u_ts0 = uthres

    g = g0*1.0E2
    emit=0.0

    DO n = 1, smx
       den(n) = den_salt(n)*1.0D-3         ! (g cm^-3)
       diam(n) = 2.0*reff_salt(n)*1.0D2    ! (cm)
       DO i = 1,imx
          DO j = 1,jmx
             rhoa = airden(i,j,1)*1.0D-3       ! (g cm^-3)

             ! FENGSHA uses the 13 category soil type from the USDA
             ! call calc_fengsha_styp(massfrac(1),massfrac(3),massfrac(2),styp)
             ! Fengsha uses threshold velocities based on dale gilletes data
             ! call fengsha_utst(styp,uthres,u_ts0)

             ! Friction velocity threshold correction function based on physical
             ! properties related to moisture tension. Soil moisture greater than
             ! dry limit serves to increase threshold friction velocity (making
             ! it more difficult to loft dust). When soil moisture has not reached
             ! dry limit, treat as dry

             IF (gravsm(i,j) > drylimit(i,j)) THEN
                u_ts = MAX(0.0D+0,u_ts0*(sqrt(1.0+1.21*(gravsm(i,j)-drylimit(i,j))**0.68)) / R)
             ELSE
                u_ts = u_ts0 / R
             END IF

             ! Calculate total vertical mass flux (note beta has units of m^-1)
             ! Beta acts to tone down dust in areas with so few dust-sized particles that the
             ! lofting efficiency decreases.  Otherwise, super sandy zones would be huge dust
             ! producers, which is generally not the case.  Equation derived from wind-tunnel
             ! experiments (see MB95).

             beta=10**(13.6*massfrac(1)-6.0)  ! (unitless)
             if (massfrac(1) <= 0.2) then
                beta=10**(13.4*massfrac(1)-6.0)
             else
                beta = 2.E-4
             endif

             !---------------------------------------------------------------------
             ! formula of Draxler & Gillette (2001) Atmos. Environ.
             ! F   =  K A (r/g) U* ( U*^2 - Ut*^2 )
             !
             ! where:
             !     F   = vertical emission flux  [g/m**2-s]
             !     K   = constant 2.0E-04                      [1/m]
             !     A   = 0~3.5  mean = 2.8  (fudge factor)
             !     U*  = friction velocity                     [m/s]
             !     Ut* = threshold friction velocity           [m/s]
             !
             !--------------------------------------------------------------------

             IF (ustar(i,j) .gt. u_ts) then
                call fengsha_hflux(data,ustar(i,j),u_ts,beta, salt)
                salt = alpha * cmb * ds_rel(n) * airden(i,j,1) / g0 * salt * (erod(i,j)**gamma) * beta
             else
                salt = 0.
             endif
             ! EROD is taken into account above
             emit = emit + salt 
          END DO
       END DO
    END DO

    ! Now that we have the total dust emission, distribute into dust bins using
    ! lognormal distribution (Dr. Jasper Kok, in press), and
    ! calculate total mass emitted over the grid box over the timestep.
    !
    ! In calculating the Kok distribution, we assume upper and lower limits to each bin.
    ! For reff_dust=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/) (default),
    ! lower limits were ASSUMED at lo_dust=(/0.1D-6,1.0D-6,1.8D-6,3.0D-6,6.0D-6/)
    ! upper limits were ASSUMED at up_dust=(/1.0D-6,1.8D-6,3.0D-6,6.0D-6,10.0D-6/)
    ! These may be changed within module_data_gocart_dust.F, but make sure it is
    ! consistent with reff_dust values.  These values were taken from the original
    ! GOCART bin configuration. We use them here to calculate dust bin width, dlndp.
    ! dVol is the volume distribution. You know...if you were wondering. GC

    dvol_tot=0.
    DO n=1,nmx
       dlndp(n)=LOG(up_dust(n)/lo_dust(n))
       dvol(n)=(2.0*reff_dust(n)/cv)*(1.+ERF(LOG(2.0*reff_dust(n)/mmd_dust)/(SQRT(2.)*LOG(gsd_dust))))*&
            EXP(-(2.0*reff_dust(n)/lambda)**3.0)*dlndp(n)
       dvol_tot=dvol_tot+dvol(n)
       ! Convert mass flux to volume flux
       !emit_vol=emit/den_dust(n) ! (m s^-1)
    END DO
    DO n=1,nmx
       distr_dust(n)=dvol(n)/dvol_tot
       !print *,"distr_dust(",n,")=",distr_dust(n)
    END DO

    ! Now distribute total vertical emission into dust bins and update concentration.

    DO n=1,nmx
       DO i=1,imx
          DO j=1,jmx
             ! Calculate total mass emitted
             dsrc = emit*distr_dust(n)*dxy(j)*dt1  ! (kg)
             IF (dsrc < 0.0) dsrc = 0.0

             ! Update dust mixing ratio at first model level.
             tc(i,j,1,n) = tc(i,j,1,n) + dsrc / airmas(i,j,1) ! (kg/kg)
             !   bems(i,j,n) = dsrc  ! diagnostic
             !bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) ! diagnostic (g/m2/s)
             bems(i,j,n) = 1.e+9*dsrc/(dxy(j)*dt1) ! diagnostic (ug/m2/s) !lzhang
          END DO
       END DO
    END DO

  END SUBROUTINE source_dust

  subroutine fengsha_utst(data,styp,uth, ut)
    implicit none
    type(smoke_data), intent(inout) :: data

    integer,                            intent(in)  :: styp
    real(kind_phys), dimension(fengsha_maxstypes), intent(in)  :: uth
    real(kind_phys),                 intent(out) :: ut
    ut = uth(styp)
!     real (kind_phys) :: uth(13) = &
!          (/ 0.08,   & ! Sand          - 1
!          0.20,    & ! Loamy Sand      - 2
!          0.30,    & ! Sandy Loam      - 3
!          0.30,    & ! Silt Loam       - 4
!          0.35,    & ! Silt            - 5
!          0.60,    & ! Loam            - 6
!          0.30,    & ! Sandy Clay Loam - 7
!          0.35,    & ! Silty Clay Loam - 8
!          0.45,    & ! Clay Loam       - 9
!          0.45,    & ! Sandy Clay      - 10
!          0.45,    & ! Silty Clay      - 11
!          0.60,    & ! Clay            - 12
!          9.999 /)   ! Other           - 13
    return
  end subroutine fengsha_utst

  subroutine calc_fengsha_styp(data, clay, sand, silt, type)
    implicit none
    type(smoke_data), intent(inout) :: data

    !---------------------------------------------------------------
    ! Function: calculate soil type based on USDA definition.
    ! Source: USDA soil texture calculator
    !
    ! Defintion of soil types:
    !
    !
    ! NOAH 1      2             3           4           5      6      7                 8                9           10           11           12
    ! PX   1      2             3           4           -      5      6                 7                8           9            10           11
    ! Soil "Sand" "Loamy Sand" "Sandy Loam" "Silt Loam" "Silt" "Loam" "Sandy Clay Loam" "Silt Clay Loam" "Clay Loam" "Sandy Clay" "Silty Clay" "Clay"
    !---------------------------------------------------------------
    REAL(kind_phys), intent(in) ::  clay, sand, silt
    integer, intent(out) ::  type
    real(kind_phys) :: cly, snd, slt

    type = 0

    snd = sand * 100.
    cly = clay * 100.
    slt = silt * 100.
    if (slt+1.5*cly .lt. 15)                                                                type = 1      ! snd
    if (slt+1.5*cly .ge. 15 .and.slt+1.5*cly .lt. 30)                                       type = 2      ! loamy snd
    if (cly .ge. 7 .and. cly .lt. 20 .and. snd .gt. 52 .and. slt+2*cly .ge. 30)             type = 3      ! sndy loam (cond 1)
    if (cly .lt. 7 .and. slt .lt. 50 .and. slt+2*cly .ge. 30)                               type = 3      ! sndy loam (cond 2)
    if (slt .ge. 50 .and. cly .ge. 12 .and.cly .lt. 27 )                                    type = 4      ! slt loam (cond 1)
    if (slt .ge. 50 .and. slt .lt. 80 .and.cly .lt. 12)                                     type = 4      ! slt loam (cond 2)
    if (slt .ge. 80 .and. cly .lt. 12)                                                      type = 5      ! slt
    if (cly .ge. 7  .and. cly .lt. 27 .and.slt .ge. 28 .and. slt .lt. 50 .and.snd .le. 52)  type = 6      ! loam
    if (cly .ge. 20 .and. cly .lt. 35 .and.slt .lt. 28 .and. snd .gt. 45)                   type = 7      ! sndy cly loam
    if (cly .ge. 27 .and. cly .lt. 40 .and.snd .lt. 20)                                     type = 8      ! slt cly loam
    if (cly .ge. 27 .and. cly .lt. 40 .and.snd .ge. 20 .and. snd .le. 45)                   type = 9      ! cly loam
    if (cly .ge. 35 .and. snd .gt. 45)                                                      type = 10     ! sndy cly
    if (cly .ge. 40 .and. slt .ge. 40)                                                      type = 11     ! slty cly
    if (cly .ge. 40 .and. snd .le. 45 .and.slt .lt. 40)                                     type = 12     ! clay
    return
  end subroutine calc_fengsha_styp

  subroutine fengsha_drag(data,z0,R)
    implicit none
    type(smoke_data), intent(inout) :: data

    real(kind_phys), intent(in) :: z0
    real(kind_phys), intent(out) :: R
    real(kind_phys), parameter :: z0s = 1.0e-04 !Surface roughness for ideal bare surface [m]
    ! ------------------------------------------------------------------------
    ! Function: Calculates the MacKinnon et al. 2004 Drag Partition Correction
    !
    !   R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)
    !
    !--------------------------------------------------------------------------
    ! Drag partition correction. See MacKinnon et al. (2004),
    !     doi:10.1016/j.geomorph.2004.03.009
    R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)

    ! Drag partition correction. See Marticorena et al. (1997),
    !     doi:10.1029/96JD02964
    !R = 1.0 - log(z0 / z0s) / log( 0.7 * (10./z0s) ** 0.8)

    return
  end subroutine fengsha_drag

  subroutine fengsha_hflux(data,ust,utst, kvh, salt)
    !---------------------------------------------------------------------
    ! Function: Calculates the Horizontal Saltation Flux, Q, and then
    !           calculates the vertical flux.
    !
    ! formula of Draxler & Gillette (2001) Atmos. Environ.
    ! F   =  K A (r/g) U* ( U*^2 - Ut*^2 )
    !
    ! where:
    !     F   = vertical emission flux  [g/m**2-s]
    !     K   = constant 2.0E-04                      [1/m]
    !     A   = 0~3.5  mean = 2.8  (fudge factor)
    !     U*  = friction velocity                     [m/s]
    !     Ut* = threshold friction velocity           [m/s]
    !
    !--------------------------------------------------------------------
    implicit none
    type(smoke_data), intent(inout) :: data
    real(kind_phys), intent(in) :: ust, & ! friction velocity
                                     utst, & ! threshold friction velocity
                                      kvh    ! vertical to horizontal mass flux ratio

    real(kind_phys), intent(out) :: salt
    real(kind_phys) :: Q
    Q = ust * (ust * ust - utst * utst)
    salt = Q ! sdep * kvh * Q

    return
  end subroutine fengsha_hflux


end module dust_fengsha_mod
