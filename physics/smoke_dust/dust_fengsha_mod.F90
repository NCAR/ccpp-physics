!>\file  dust_fengsha_mod.F90
!! This file contains the FENGSHA dust scheme.

module dust_fengsha_mod
!
!  This module developed by Barry Baker (NOAA ARL)
!  For serious questions contact barry.baker@noaa.gov
!
!  07/16/2019 - Adapted for NUOPC/GOCART, R. Montuoro
!  02/01/2020 - Adapted for FV3/CCPP, Haiqin Li

  use machine ,        only : kind_phys
  use dust_data_mod

  implicit none

  private

  public :: gocart_dust_fengsha_driver

contains

  subroutine gocart_dust_fengsha_driver(dt,              &
       chem,rho_phy,smois,stemp,p8w,ssm,                 &
       isltyp,snowh,xland,area,g,emis_dust,              &
       ust,znt,clay,sand,rdrag,uthr,                     &
       num_emis_dust,num_chem,num_soil_layers,           &
       ids,ide, jds,jde, kds,kde,                        &
       ims,ime, jms,jme, kms,kme,                        &
       its,ite, jts,jte, kts,kte)
    IMPLICIT NONE
    INTEGER,      INTENT(IN   ) ::                       &
         ids,ide, jds,jde, kds,kde,                      &
         ims,ime, jms,jme, kms,kme,                      &
         its,ite, jts,jte, kts,kte,                      &
         num_emis_dust,num_chem,num_soil_layers

    ! 2d input variables
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: ssm     ! Sediment supply map
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: snowh   ! snow height (m)
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: xland   ! dominant land use type
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: area    ! area of grid cell
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: ust     ! friction velocity
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: znt     ! Surface Roughness length (m)
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: clay    ! Clay Fraction (-)
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: sand    ! Sand Fraction (-)
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: rdrag   ! Drag Partition (-)
    REAL(kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: uthr    ! Dry Threshold Velocity (m/s)

    INTEGER,         DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: isltyp  ! soil type

    ! 3d input variables
    REAL(kind_phys), DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN) :: p8w
    REAL(kind_phys), DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN) :: rho_phy
    REAL(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), INTENT(INOUT) :: chem
    REAL(kind_phys), DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL, INTENT(INOUT) :: emis_dust
    REAL(kind_phys), DIMENSION( ims:ime, num_soil_layers, jms:jme ), INTENT(IN) :: smois, stemp

    !0d input variables 
    REAL(kind_phys), INTENT(IN) :: dt ! time step
    REAL(kind_phys), INTENT(IN) :: g  ! gravity (m/s**2)



    ! Local variables
    integer :: nmx,i,j,k,imx,jmx,lmx
    integer :: ilwi
    real(kind_phys) :: airden ! air density
    REAL(kind_phys) :: airmas ! dry air mass
    real(kind_phys) :: dxy
    real(kind_phys) :: conver,converi ! conversion values 
    real(kind_phys) :: R ! local drag partition
    real(kind_phys) :: ustar
    real(kind_phys), DIMENSION (num_emis_dust) :: tc
    real(kind_phys), DIMENSION (num_emis_dust) :: bems
    real(kind_phys), DIMENSION (num_emis_dust) :: distribution
    real(kind_phys), dimension (3) :: massfrac
    real(kind_phys) :: erodtot
    real(kind_phys) :: moist_volumetric

    ! conversion values
    conver=1.e-9
    converi=1.e9

    ! Number of dust bins

    imx=1
    jmx=1
    lmx=1
    nmx=ndust

    k=kts
    do j=jts,jte
       do i=its,ite

          ! Don't do dust over water!!!

          ilwi=0
          if(xland(i,j).lt.1.5)then
             ilwi=1

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

             airmas=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g
             airden=rho_phy(i,kts,j)
             ustar=ust(i,j)
             dxy=area(i,j)

             ! Mass fractions of clay, silt, and sand.
             massfrac(1)=clay(i,j)
             massfrac(2)=1-(clay(i,j)+sand(i,j))
             massfrac(3)=sand(i,j)


             ! Total erodibility.
             
             erodtot = ssm(i,j) ! SUM(erod(i,j,:))
             
             ! Don't allow roughness lengths greater than 20 cm to be lofted.
             ! This kludge accounts for land use types like urban areas and
             ! forests which would otherwise show up as high dust emitters.
             ! This is a placeholder for a more widely accepted kludge
             ! factor in the literature, which reduces lofting for rough areas.
             ! Forthcoming...

             IF (znt(i,j) .gt. 0.2) then
                ilwi=0
             endif

             ! limit where there is lots of vegetation

             ! limit where there is snow on the ground
             if (snowh(i,j) .gt. 0) then
                ilwi = 0
             endif

             ! Don't emit over frozen soil
             if (stemp(i,1,j) < 268.0) then ! -5C
                ilwi = 0
             endif

             ! Do not allow areas with bedrock, lava, or land-ice to loft

             IF (isltyp(i,j) .eq. 15 .or. isltyp(i,j) .eq. 16. .or. &
                  isltyp(i,j) .eq. 18) then
                ilwi=0
             ENDIF
             IF (isltyp(i,j) .eq. 0)then
                ilwi=0
             endif
             if(ilwi == 0 ) cycle

             ! get drag partition
             ! FENGSHA uses the drag partition correction of MacKinnon et al 2004
             !     doi:10.1016/j.geomorph.2004.03.009
             if (dust_calcdrag .ne. 1) then
                call fengsha_drag(znt(i,j),R)
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

             ! soil moisture correction factor 
             moist_volumetric = dust_moist_correction * smois(i,2,j) 

             ! Call dust emission routine.
             
             call source_dust(imx,jmx, lmx, nmx, dt, tc, ustar, massfrac, & 
                  erodtot, dxy, moist_volumetric, airden, airmas, bems, g, dust_alpha, dust_gamma, &
                  R, uthr(i,j))

             ! convert back to concentration

             chem(i,kts,j,p_dust_1)=tc(1)*converi
             chem(i,kts,j,p_dust_2)=tc(2)*converi
             chem(i,kts,j,p_dust_3)=tc(3)*converi
             chem(i,kts,j,p_dust_4)=tc(4)*converi
             chem(i,kts,j,p_dust_5)=tc(5)*converi

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


  subroutine source_dust(imx, jmx, lmx, nmx, dt1, tc, ustar, massfrac, &
                  erod, dxy, smois, airden, airmas, bems, g0, alpha, gamma, &
                  R, uthres)

    ! ****************************************************************************
    ! *  Evaluate the source of each dust particles size bin by soil emission
    ! *
    ! *  Input:
    ! *         EROD      Fraction of erodible grid cell                (-)
    ! *         smois     Volumetric  soil moisture                     (m3/m3)
    ! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
    ! *         GAMMA     Tuning constant for erodibility               (-)
    ! *         DXY       Surface of each grid cell                     (m2)
    ! *         AIRMAS    Mass of air for each grid box                 (kg)
    ! *         AIRDEN    Density of air for each grid box              (kg/m3)
    ! *         USTAR     Friction velocity                             (m/s)
    ! *         DT1       Time step                                     (s)
    ! *         NMX       Number of dust bins                           (-)
    ! *         IMX       Number of I points                            (-)
    ! *         JMX       Number of J points                            (-)
    ! *         LMX       Number of L points                            (-)
    ! *         R         Drag Partition                                (-)
    ! *         UTHRES    FENGSHA Dry Threshold Velocities              (m/s)
    ! *
    ! *  Data:
    ! *         MASSFRAC  Fraction of mass in each of 3 soil classes    (-) (clay silt sand) 
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
    ! *
    ! *  Working:
    ! *         RHOA      Density of air in cgs                         (g/cm3)
    ! *         DS_REL    Saltation surface area distribution           (-)
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

    ! Input
    INTEGER,            INTENT(IN)    :: imx,jmx,lmx,nmx
    REAL(kind_phys), INTENT(IN)    :: dt1
    REAL(kind_phys), INTENT(IN)    :: ustar
    REAL(kind_phys), INTENT(IN)    :: massfrac(3)
    REAL(kind_phys), INTENT(IN)    :: erod
    REAL(kind_phys), INTENT(IN)    :: dxy
    REAL(kind_phys), INTENT(IN)    :: smois
    REAL(kind_phys), INTENT(IN)    :: airden
    REAL(kind_phys), INTENT(IN)    :: airmas
    REAL(kind_phys), INTENT(IN)    :: g0
    REAL(kind_phys), INTENT(IN)    :: alpha
    REAL(kind_phys), INTENT(IN)    :: gamma
    REAL(kind_phys), INTENT(IN)    :: R
    REAL(kind_phys), INTENT(IN)    :: uthres

    ! Output
    REAL(kind_phys), INTENT(INOUT) :: tc(nmx)

    ! Local Variables
    REAL(kind_phys), INTENT(OUT)   :: bems(nmx)
    
    REAL(kind_phys) :: dvol(nmx)
    REAL(kind_phys) :: distr_dust(nmx)
    REAL(kind_phys) :: dlndp(nmx)
    REAL(kind_phys) :: dsrc
    REAL(kind_phys) :: dvol_tot
    REAL(kind_phys) :: emit
    REAL(kind_phys) :: emit_vol
    REAL(kind_phys) :: rhoa
    INTEGER   :: i, j, n

    ! Constant of proportionality from Marticorena et al, 1997 (unitless)
    ! Arguably more ~consistent~ fudge than alpha, which has many walnuts
    ! sprinkled throughout the literature. - GC

    REAL(kind_phys), PARAMETER :: cmb=1.0
    REAL(kind_phys), PARAMETER :: kvhmax=2.0e-4

    ! Parameters used in Kok distribution function. Advise not to play with
    ! these without the expressed written consent of someone who knows what
    ! they're doing. - GC

    REAL(kind_phys), PARAMETER :: mmd_dust=3.4D-6  ! median mass diameter (m)
    REAL(kind_phys), PARAMETER :: gsd_dust=3.0     ! geom. std deviation
    REAL(kind_phys), PARAMETER :: lambda=12.0D-6   ! crack propagation length (m)
    REAL(kind_phys), PARAMETER :: cv=12.62D-6      ! normalization constant
    REAL(kind_phys), PARAMETER :: RHOSOIL=2650.


    ! calculate the total vertical dust flux 

    emit = 0.0

    call DustEmissionFENGSHA(smois,massfrac(1),massfrac(3), massfrac(2), &
                                erod, R, airden, ustar, uthres, alpha, gamma, kvhmax, &
                                g0, RHOSOIL, emit)

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
       ! Calculate total mass emitted
       dsrc = emit*distr_dust(n)*dxy*dt1  ! (kg)
       IF (dsrc < 0.0) dsrc = 0.0
       
       ! Update dust mixing ratio at first model level.
       tc(n) = tc(n) + dsrc / airmas ! (kg/kg)
       !   bems(i,j,n) = dsrc  ! diagnostic
       !bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) ! diagnostic (g/m2/s)
       bems(n) = 1.e+9*dsrc/(dxy*dt1) ! diagnostic (ug/m2/s) !lzhang
       
    END DO
    tc(1)=tc(1)+0.286*tc(2)       ! This is just for RRFS-SD. DO NOT use in other models!!!
    tc(5)=0.714*tc(2)+tc(3)+tc(4) ! This is just for RRFS-SD. DO NOT use in other models!!!

  END SUBROUTINE source_dust


  subroutine fengsha_drag(z0,R)
    implicit none

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

  subroutine DustEmissionFENGSHA(slc, clay, sand, silt,  &
                                  ssm, rdrag, airdens, ustar, uthrs, alpha, gamma, &
                                  kvhmax, grav, rhop, emissions)
    
    ! !USES:
    implicit NONE
    
! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
    REAL(kind_phys), intent(in) :: clay     ! fractional clay content [1]
    REAL(kind_phys), intent(in) :: sand     ! fractional sand content [1]
    REAL(kind_phys), intent(in) :: silt     ! fractional silt content [1]
    REAL(kind_phys), intent(in) :: ssm      ! erosion map [1]
    REAL(kind_phys), intent(in) :: rdrag    ! drag partition [1/m]
    REAL(kind_phys), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
    REAL(kind_phys), intent(in) :: ustar    ! friction velocity [m/sec]
    REAL(kind_phys), intent(in) :: uthrs    ! threshold velocity [m/2]
    REAL(kind_phys), intent(in) :: alpha    ! scaling factor [1]
    REAL(kind_phys), intent(in) :: gamma    ! scaling factor [1]
    REAL(kind_phys), intent(in) :: kvhmax   ! max. vertical to horizontal mass flux ratio [1]
    REAL(kind_phys), intent(in) :: grav     ! gravity [m/sec^2]
    REAL(kind_phys), intent(in) :: rhop     ! soil class density [kg/m^3]
    
    ! !OUTPUT PARAMETERS:
    REAL(kind_phys), intent(inout) :: emissions ! binned surface emissions [kg/(m^2 sec)]
    
    ! !DESCRIPTION: Compute dust emissions using NOAA/ARL FENGSHA model
    !
    ! !REVISION HISTORY:
    !
    ! 22Feb2020 B.Baker/NOAA    - Original implementation
    ! 29Mar2021 R.Montuoro/NOAA - Refactored for process library
    ! 09Aug2022 B.Baker/NOAA    - Adapted for CCPP-Physics
    
    ! !Local Variables
    real(kind_phys)                  :: alpha_grav
    real(kind_phys)                  :: h
    real(kind_phys)                  :: kvh
    real(kind_phys)                  :: q
    real(kind_phys)                  :: rustar
    real(kind_phys)                  :: total_emissions
    real(kind_phys)                  :: u_sum, u_thresh
    
!EOP
!-------------------------------------------------------------------------
!  Begin

!  Initialize emissions
!  --------------------
   emissions = 0.

!  Prepare scaling factor
!  ----------------------
   alpha_grav = alpha / grav

   ! Compute vertical-to-horizontal mass flux ratio
   ! ----------------------------------------------
   kvh = DustFluxV2HRatioMB95(clay, kvhmax)

   ! Compute total emissions
   ! -----------------------
   emissions = alpha_grav * (ssm ** gamma) * airdens * kvh

   !  Compute threshold wind friction velocity using drag partition
   !  -------------------------------------------------------------
   rustar = rdrag * ustar

   !  Now compute size-dependent total emission flux
   !  ----------------------------------------------

   if (dust_moist_opt .eq. 1) then

      ! Fecan moisture correction
      ! -------------------------
      h = moistureCorrectionFecan(slc, sand, clay)
   else
      ! shao soil moisture correction 
      h = moistureCorrectionShao(slc)
   end if
   ! Adjust threshold
   ! ----------------
   u_thresh = uthrs * h
   
   u_sum = rustar + u_thresh
   
   ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
   ! ---------------------------------------------------------------------------
   q = max(0., rustar - u_thresh) * u_sum * u_sum
   
   ! Distribute emissions to bins and convert to mass flux (kg s-1)
   ! --------------------------------------------------------------
   emissions = emissions * q


 end subroutine DustEmissionFENGSHA
!-----------------------------------------------------------------
  real function soilMoistureConvertVol2Grav(vsoil, sandfrac)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: vsoil       ! volumetric soil moisture fraction [1]
    REAL(kind_phys), intent(in) :: sandfrac    ! fractional sand content [1]

! !DESCRIPTION: Convert soil moisture fraction from volumetric to gravimetric.
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: vsat

!  !CONSTANTS:
    REAL(kind_phys), parameter :: rhow = 1000.    ! density of water [kg m-3]
    REAL(kind_phys), parameter :: rhop = 1700.    ! density of dry soil 
!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Saturated volumetric water content (sand-dependent) ! [m3 m-3]
    vsat = 0.489 - 0.126 * sandfrac 
    

!  Gravimetric soil content
    soilMoistureConvertVol2Grav = 100.0 * (vsoil * rhow / rhop / ( 1. - vsat))

  end function soilMoistureConvertVol2Grav
!----------------------------------------------------------------
  real function moistureCorrectionFecan(slc, sand, clay)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: slc     ! liquid water content of top soil layer, volumetric fraction [1]
    REAL(kind_phys), intent(in) :: sand    ! fractional sand content [1]
    REAL(kind_phys), intent(in) :: clay    ! fractional clay content [1]

! !DESCRIPTION: Compute correction factor to account for Fecal soil moisture
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: grvsoilm
    real :: drylimit

!EOP
!---------------------------------------------------------------
!  Begin...

!  Convert soil moisture from volumetric to gravimetric
    grvsoilm = soilMoistureConvertVol2Grav(slc, sand)

!  Compute fecan dry limit
    drylimit = dust_drylimit_factor * clay * (14.0 * clay + 17.0)

!  Compute soil moisture correction
    moistureCorrectionFecan = sqrt(1.0 + 1.21 * max(0., grvsoilm - drylimit)**0.68)

  end function moistureCorrectionFecan
!----------------------------------------------------------------
  real function moistureCorrectionShao(slc)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: slc     ! liquid water content of top soil layer, volumetric fraction [1]

! !DESCRIPTION: Compute correction factor to account for Fecal soil moisture
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: grvsoilm
    real :: drylimit

!EOP
!---------------------------------------------------------------
!  Begin...

    if (slc < 0.03) then
       moistureCorrectionShao = exp(22.7 * slc) 
    else
       moistureCorrectionShao = exp(95.3 * slc - 2.029)
    end if

  end function moistureCorrectionShao
!---------------------------------------------------------------
  real function DustFluxV2HRatioMB95(clay, kvhmax)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: clay      ! fractional clay content [1]
    REAL(kind_phys), intent(in) :: kvhmax    ! maximum flux ratio [1]

!  !CONSTANTS:
    REAL(kind_phys), parameter :: clay_thresh = 0.2    ! clay fraction above which the maximum flux ratio is returned

! !DESCRIPTION: Computes the vertical-to-horizontal dust flux ratio according to
!               B.Marticorena, G.Bergametti, J.Geophys.Res., 100(D8), 164!               doi:10.1029/95JD00690
!
! !REVISION HISTORY:
!
! 22Feb2020 B.Baker/NOAA    - Original implementation
! 01Apr2021 R.Montuoro/NOAA - Adapted for GOCART process library
!
!EOP
!-------------------------------------------------------------------------
!  Begin...

    if (clay > clay_thresh) then
       DustFluxV2HRatioMB95 = kvhmax
    else
       DustFluxV2HRatioMB95 = 10.0**(13.4*clay-6.0)
    end if

  end function DustFluxV2HRatioMB95
  
end module dust_fengsha_mod
