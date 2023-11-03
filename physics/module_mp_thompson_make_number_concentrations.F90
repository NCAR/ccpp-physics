!>\file module_mp_thompson_make_number_concentrations.F90
!! This file contains

!>\ingroup aathompson
module module_mp_thompson_make_number_concentrations

      use physcons, only: PI => con_pi

      implicit none

      private

      public make_IceNumber, make_DropletNumber, make_RainNumber

!      Q_ice              is cloud ice mixing ratio, units of kg/m3
!      Q_cloud            is cloud water mixing ratio, units of kg/m3
!      Q_rain             is rain mixing ratio, units of kg/m3
!      temp               is air temperature in Kelvin
!      make_IceNumber     is cloud droplet number mixing ratio, units of number per m3
!      make_DropletNumber is rain number mixing ratio, units of number per kg of m3
!      make_RainNumber    is rain number mixing ratio, units of number per kg of m3
!      qnwfa              is number of water-friendly aerosols in number per kg

!+---+-----------------------------------------------------------------+ 
!+---+-----------------------------------------------------------------+ 

   contains
!>\ingroup aathompson
!!Table of lookup values of radiative effective radius of ice crystals
!! as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
!! radiation code where it is attributed to Jon Egill Kristjansson
!! and coauthors.
      elemental real function make_IceNumber (Q_ice, temp)

      !IMPLICIT NONE
      REAL, PARAMETER:: Ice_density = 890.0
      !REAL, PARAMETER:: PI = 3.1415926536
      real, intent(in):: Q_ice, temp
      integer idx_rei
      real corr, reice, deice
      double precision lambda

!+---+-----------------------------------------------------------------+ 
!..Table of lookup values of radiative effective radius of ice crystals
!.. as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
!.. radiation code where it is attributed to Jon Egill Kristjansson
!.. and coauthors.
!+---+-----------------------------------------------------------------+ 

      !real retab(95)
      !data retab /                                                      &
      !   5.92779, 6.26422, 6.61973, 6.99539, 7.39234,                   &
      !   7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,          &
      !   10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,          &
      !   15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,          &
      !   20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,          &
      !   27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,          &
      !   31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,          &
      !   34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,          &
      !   38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,          &
      !   42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,          &
      !   50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,          &
      !   65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,          &
      !   93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,          &
      !   124.954, 130.630, 136.457, 142.446, 148.608, 154.956,          &
      !   161.503, 168.262, 175.248, 182.473, 189.952, 197.699,          &
      !   205.728, 214.055, 222.694, 231.661, 240.971, 250.639/
      real, dimension(95), parameter:: retab = (/                       &
         5.92779, 6.26422, 6.61973, 6.99539, 7.39234,                   &
         7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,          &
         10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,          &
         15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,          &
         20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,          &
         27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,          &
         31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,          &
         34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,          &
         38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,          &
         42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,          &
         50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,          &
         65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,          &
         93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,          &
         124.954, 130.630, 136.457, 142.446, 148.608, 154.956,          &
         161.503, 168.262, 175.248, 182.473, 189.952, 197.699,          &
         205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)

      if (Q_ice == 0) then
         make_IceNumber = 0
         return
      end if

!+---+-----------------------------------------------------------------+ 
!..From the model 3D temperature field, subtract 179K for which
!.. index value of retab as a start.  Value of corr is for
!.. interpolating between neighboring values in the table.
!+---+-----------------------------------------------------------------+ 

      idx_rei = int(temp-179.)
      idx_rei = min(max(idx_rei,1),94)
      corr = temp - int(temp)
      reice = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
      deice = 2.*reice * 1.E-6

!+---+-----------------------------------------------------------------+ 
!..Now we have the final radiative effective size of ice (as function
!.. of temperature only).  This size represents 3rd moment divided by
!.. second moment of the ice size distribution, so we can compute a
!.. number concentration from the mean size and mass mixing ratio.
!.. The mean (radiative effective) diameter is 3./Slope for an inverse
!.. exponential size distribution.  So, starting with slope, work
!.. backwords to get number concentration.
!+---+-----------------------------------------------------------------+ 

      lambda = 3.0 / deice
      make_IceNumber = Q_ice * lambda*lambda*lambda / (PI*Ice_density)

!+---+-----------------------------------------------------------------+ 
!..Example1: Common ice size coming from Thompson scheme is about 30 microns.
!.. An example ice mixing ratio could be 0.001 g/kg for a temperature of -50C.
!.. Remember to convert both into MKS units.  This gives N_ice=357652 per kg.
!..Example2: Lower in atmosphere at T=-10C matching ~162 microns in retab,
!.. and assuming we have 0.1 g/kg mixing ratio, then N_ice=28122 per kg, 
!.. which is 28 crystals per liter of air if the air density is 1.0.
!+---+-----------------------------------------------------------------+ 

      return
      end function make_IceNumber

!+---+-----------------------------------------------------------------+ 
!+---+-----------------------------------------------------------------+ 

!>\ingroup aathompson
!!
      elemental real function make_DropletNumber (Q_cloud, qnwfa)

      !IMPLICIT NONE

      real, intent(in):: Q_cloud, qnwfa

      !real, parameter:: PI = 3.1415926536
      real, parameter:: am_r = PI*1000./6.
      real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336,   &
     &                504,720,990,1320,1716,2184,2730,3360,4080,4896/)
      double precision:: lambda, qnc
      real:: q_nwfa, x1, xDc
      integer:: nu_c

      if (Q_cloud == 0) then
         make_DropletNumber = 0
         return
      end if

!+---+

      q_nwfa = MAX(99.E6, MIN(qnwfa,5.E10))
      nu_c = MAX(2, MIN(NINT(2.5E10/q_nwfa), 15))

      x1 = MAX(1., MIN(q_nwfa*1.E-9, 10.)) - 1.
      xDc = (30. - x1*20./9.) * 1.E-6

      lambda = (4.0D0 + nu_c) / xDc
      qnc = Q_cloud / g_ratio(nu_c) * lambda*lambda*lambda / am_r
      make_DropletNumber = SNGL(qnc)

      return
      end function make_DropletNumber

!+---+-----------------------------------------------------------------+ 
!+---+-----------------------------------------------------------------+ 

!>\ingroup aathompson
!!
      elemental real function make_RainNumber (Q_rain, temp)

      IMPLICIT NONE

      real, intent(in):: Q_rain, temp
      double precision:: lambda, N0, qnr
      !real, parameter:: PI = 3.1415926536
      real, parameter:: am_r = PI*1000./6.

      if (Q_rain == 0) then
         make_RainNumber = 0
         return
      end if

      !+---+-----------------------------------------------------------------+ 
      !.. Not thrilled with it, but set Y-intercept parameter to Marshal-Palmer value
      !.. that basically assumes melting snow becomes typical rain. However, for
      !.. -2C < T < 0C, make linear increase in exponent to attempt to keep
      !.. supercooled collision-coalescence (warm-rain) similar to drizzle rather
      !.. than bigger rain drops.  While this could also exist at T>0C, it is
      !.. more difficult to assume it directly from having mass and not number.
      !+---+-----------------------------------------------------------------+ 

      N0 = 8.E6

      if (temp .le. 271.15) then
         N0 = 8.E8
      elseif (temp .gt. 271.15 .and. temp.lt.273.15) then
         N0 = 8. * 10**(279.15-temp)
      endif

      lambda = SQRT(SQRT(N0*am_r*6.0/Q_rain))
      qnr = Q_rain / 6.0 * lambda*lambda*lambda / am_r
      make_RainNumber = SNGL(qnr)

      return
      end function make_RainNumber

!+---+-----------------------------------------------------------------+ 
!+---+-----------------------------------------------------------------+ 

end module module_mp_thompson_make_number_concentrations
