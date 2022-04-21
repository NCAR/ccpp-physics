!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!         Adapted by NOAA/GSD/ESRL                                       !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  seas_ngac_mod.F90 --- Calculate the Seasalt Emissions
!
! !INTERFACE:
!

   module seas_ngac_mod

! !USES:

!  use chem_comm_mod,  only : chem_comm_isroot
   use machine ,       only : kind_phys
   use physcons,       only : pi=>con_pi

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  SeasaltEmission


! !CONSTANTS
   real(kind=kind_phys), parameter    :: r80fac = 1.65     ! ratio of radius(RH=0.8)/radius(RH=0.) [Gerber]
   real(kind=kind_phys), parameter    :: rhop = 2200.      ! dry seasalt density [kg m-3]

!
! !DESCRIPTION:
!
!  This module implements the sea salt aerosol emission parameterizations.
!  For all variants, emissions are some function of wind speed (and possibly
!  other dynamical parameters) and the sea salt particle radius.  Here,
!  we assume the model passes in dry radius (or dry radius of size bin edges).
!  Output is the mass emission flux (kg m-2 s-1) into that radius bin.
!
! !REVISION HISTORY:
!
!  30Mar2010 Colarco    First crack!
!
!EOP
!-------------------------------------------------------------------------
CONTAINS
!
! !IROUTINE:  SeasaltEmission - Master driver to compute the sea salt emissions
!
! !INTERFACE:
!
   subroutine SeasaltEmission ( rLow, rUp, method, w10m, ustar, &
                                memissions, nemissions, rc )

! !DESCRIPTION: Calculates the seasalt mass emission flux every timestep.
!  The particular method (algorithm) used for the calculation is based
!  on the value of "method" passed on input.  Mostly these algorithms are
!  a function of wind speed and particle size (nominally at 80% RH).
!  Routine is called once for each size bin, passing in the edge radii
!  "rLow" and "rUp" (in dry radius, units of um).  Returned in the emission
!  mass flux [kg m-2 s-1].  A sub-bin assumption is made to break (possibly)
!  large size bins into a smaller space.
!
! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   real(kind=kind_phys),    intent(in)           :: rLow, rUp   ! Dry particle bin edge radii [um]
   real(kind=kind_phys),    intent(in)           :: w10m        ! 10-m wind speed [m s-1]
   real(kind=kind_phys),    intent(in)           :: ustar       ! friction velocity [m s-1]
   integer, intent(in)           :: method      ! Algorithm to use

! !OUTPUT PARAMETERS:

   real(kind=kind_phys),    intent(inout)        :: memissions      ! Mass Emissions Flux [kg m-2 s-1]
   real(kind=kind_phys),    intent(inout)        :: nemissions      ! Number Emissions Flux [# m-2 s-1]
   integer, intent(out)          :: rc              ! Error return code:
                                                    !  0 - all is well
                                                    !  1 - 
! !Local Variables
   integer       :: ir
   real(kind=kind_phys)          :: w                               ! Intermediary wind speed [m s-1]
   real(kind=kind_phys)          :: r, dr                           ! sub-bin radius spacing (dry, um)
   real(kind=kind_phys)          :: rwet, drwet                     ! sub-bin radius spacing (rh=80%, um)
   real(kind=kind_phys)          :: aFac, bFac, scalefac, rpow, exppow, wpow

   integer, parameter :: nr = 10                    ! Number of (linear) sub-size bins

   character(len=*), parameter :: myname = 'SeasaltEmission'

!  Define the sub-bins (still in dry radius)
   dr = (rUp - rLow)/nr
   r  = rLow + 0.5*dr

!  Loop over size bins
   nemissions = 0.
   memissions = 0.

   do ir = 1, nr

    rwet  = r80fac * r
    drwet = r80fac * dr

    select case(method)

     case(1)  ! Gong 2003
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 1.
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41
      w        =  w10m

     case(2)  ! Gong 1997
      aFac     = 3.
      bFac     = (0.380-log10(rwet))/0.650
      scalefac = 1.
      rpow     = 1.05
      exppow   = 1.19
      wpow     = 3.41
      w        =  w10m

     case(3)  ! GEOS5 2012
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 33.0e3
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41 - 1.
      w        =  ustar

     case default
!     if(chem_comm_isroot()) print *, 'SeasaltEmission missing algorithm method'
      rc = 1
      return

    end select


!   Number emissions flux (# m-2 s-1)
    nemissions = nemissions + SeasaltEmissionGong( rwet, drwet, w, scalefac, aFac, bFac, rpow, exppow, wpow )
!   Mass emissions flux (kg m-2 s-1)
    scalefac = scalefac * 4./3.*pi*rhop*r**3.*1.e-18 
    memissions = memissions + SeasaltEmissionGong( rwet, drwet, w, scalefac, aFac, bFac, rpow, exppow, wpow )

    r = r + dr

   end do

   rc = 0

  end subroutine SeasaltEmission


! Function to compute sea salt emissions following the Gong style
! parameterization.  Functional form is from Gong 2003:
!  dN/dr = scalefac * 1.373 * (w^wpow) * (r^-aFac) * (1+0.057*r^rpow) * 10^(exppow*exp(-bFac^2))
! where r is the particle radius at 80% RH, dr is the size bin width at 80% RH, and w is the wind speed

  function SeasaltEmissionGong ( r, dr, w, scalefac, aFac, bFac, rpow, exppow, wpow )

   real(kind=kind_phys), intent(in)    :: r, dr     ! Wet particle radius, bin width [um]
   real(kind=kind_phys), intent(in)    :: w         ! Grid box mean wind speed [m s-1] (10-m or ustar wind)
   real(kind=kind_phys), intent(in)    :: scalefac, aFac, bFac, rpow, exppow, wpow
   real(kind=kind_phys)                :: SeasaltEmissionGong

!  Initialize
   SeasaltEmissionGong = 0.

!  Particle size distribution function
   SeasaltEmissionGong = scalefac * 1.373*r**(-aFac)*(1.+0.057*r**rpow) &
                         *10**(exppow*exp(-bFac**2.))*dr
!  Apply wind speed function
   SeasaltEmissionGong = w**wpow * SeasaltEmissionGong

  end function SeasaltEmissionGong


  end module seas_ngac_mod
