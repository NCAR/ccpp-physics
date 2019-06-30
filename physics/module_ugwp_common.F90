    module ugwp_common
!
     implicit none
     
      real, parameter ::  grav =9.80665, cpd = 1004.6, grcp = grav/cpd 
      real, parameter ::  rd = 287.05 , rv =461.5      
      real, parameter ::  rgrav = 1.0/grav
           
      real, parameter ::  fv   = rv/rd - 1.0  
      real, parameter ::  rdi  = 1.0 / rd
      real, parameter ::  gor  = grav/rd
      real, parameter ::  gr2  = grav*gor
      real, parameter ::  gocp = grav/cpd
      real, parameter ::  pi = 4.*atan(1.0), pi2 = 2.*pi
!   
      real, parameter ::  rad_to_deg=180.0/pi, deg_to_rad=pi/180.0

      real, parameter ::  arad = 6370.e3
      real, parameter ::  rcpd2 = 0.5/cpd,  rcpd = 1./cpd                  
      real, parameter ::  dw2min=1.0
      real, parameter ::  bnv2min=1.e-6
      real, parameter ::  velmin=sqrt(dw2min)     
      real, parameter ::  omega1 = pi2/86400.
      real, parameter ::  omega2 = 2.*omega1   
     end module ugwp_common
!
