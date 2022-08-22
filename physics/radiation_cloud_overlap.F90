!>\file radiation_cloud_overlap.F90
!!

!>\defgroup rad_cld_ovr_mod Radiation Cloud Overlap Module
!! This module contains the calculation of cloud overlap parameters for both RRTMG and RRTMGP. 
module module_radiation_cloud_overlap
  use physparam,        only : kind_phys
  implicit none
  public :: cmp_dcorr_lgth, get_alpha_exper, get_beta_exper
  
  interface cmp_dcorr_lgth
     module procedure cmp_dcorr_lgth_hogan
     module procedure cmp_dcorr_lgth_oreopoulos
  end interface
  
contains

!>\defgroup rad_cld_ovr_mod Radiation Cloud Overlap Module
!! This module contains the calculation of cloud overlap parameters for both RRTMG and RRTMGP.
!>@{
  ! ######################################################################################
  ! Hogan et al. (2010)
  ! "Effect of improving representation of horizontal and vertical cloud structure on the 
  ! Earth's global radiation budget. Part I: Review and parametrization"
  ! https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.647
  ! ######################################################################################
!> see Shonk et al.(2010) \cite shonk_et_al_2010
  subroutine cmp_dcorr_lgth_hogan(nCol, lat, con_pi, dcorr_lgth)
    ! Inputs
    integer, intent(in) :: &
         nCol       ! Number of horizontal grid-points
    real(kind_phys), intent(in) :: &
         con_pi     ! Physical constant: Pi
    real(kind_phys), dimension(:), intent(in) :: &
         lat        ! Latitude  
    ! Outputs
    real(kind_phys), dimension(:),intent(out) :: &
         dcorr_lgth ! Decorrelation length  
    
    ! Local variables
    integer :: iCol
    
    ! Parameters
    real(kind_phys),parameter :: min_dcorr = 0.6 ! (see section 2.3)
    
    do iCol =1,nCol
       dcorr_lgth(iCol) = max(min_dcorr, 2.78-4.6*abs(lat(iCol)/con_pi))        ! (eq. 13)
    enddo
    
  end subroutine cmp_dcorr_lgth_hogan
  ! ######################################################################################
  ! Oreopoulos et al. (2012)
  ! "Radiative impacts of cloud heterogeneity and overlap in an
  ! atmospheric General Circulation Model"
  ! 10.5194/acp-12-9097-2012
  ! ######################################################################################
!>\see Oreopoulos et al.(2012) \cite Oreopoulos_2012
  subroutine cmp_dcorr_lgth_oreopoulos(nCol, icond, lat, juldat, yearlength, &
                                       dcorr_lgth, dcorr_lgth_cond)

    ! Inputs
    integer, intent(in) :: &
         nCol,        & ! Number of horizontal grid-points
         icond,       & ! Cloud condensate overlap flag
                        ! 0 = Condensate vertical overlap off
                        ! 1 = Condensate vertical overlap on
         yearlength     ! Number of days in year

    real(kind_phys), intent(in) :: &
         juldat         ! Julian date
    real(kind_phys), dimension(:), intent(in) :: &
         lat            ! Latitude  
    
    ! Outputs
    real(kind_phys), dimension(:),intent(out) :: &
         dcorr_lgth,  &  ! Cloud fraction decorrelation length (km)
         dcorr_lgth_cond ! Cloud condensate decorrelation length (km)
    
    ! Parameters for the Gaussian fits per Eqs. (10) and (11) (See Table 1)
    ! Cloud fraction (Note: amr is incorrect in Table 1)
    real(kind_phys), parameter ::  & !  
         am1 = 1.4315_kind_phys,   & ! 
         am2 = 2.1219_kind_phys,   & !
         am4 = -25.584_kind_phys,  & !
         amr = 7.0_kind_phys         !
    ! Parameters for the Gaussian fits per Eqs. (10) and (11) (See Table 1)
    ! Cloud condensate (Note: bmr is incorrect in Table 1)
    real(kind_phys), parameter ::  & !  
         bm1 = 0.72192_kind_phys,  & ! 
         bm2 = 0.78996_kind_phys,  & !
         bm4 = 40.404_kind_phys,   & !
         bmr = 8.5_kind_phys         !
    
    ! Local variables
    integer :: iCol     
    real(kind_phys) :: am3, bm3
    
    do iCol = 1, nCol
       if (juldat .gt. 181._kind_phys) then
          am3 = -4._kind_phys * amr * (juldat - 272._kind_phys) / yearlength  ! (eq. 11a)
       else
          am3 =  4._kind_phys * amr * (juldat - 91._kind_phys)  / yearlength  ! (eq. 11b)
       endif
       dcorr_lgth(iCol) = am1 + am2 * exp( -(lat(iCol) - am3)**2 / am4**2)    ! (eq. 10)
    enddo
    
    dcorr_lgth_cond(:nCol) = 0.0_kind_phys
    if (icond == 1) then    
       do iCol = 1, nCol
          if (juldat .gt. 181._kind_phys) then
             bm3 = -4._kind_phys * bmr * (juldat - 272._kind_phys) / yearlength    ! (eq. 11a)
          else
             bm3 =  4._kind_phys * bmr * (juldat - 91._kind_phys)  / yearlength    ! (eq. 11b)
          endif
          dcorr_lgth_cond(iCol) = bm1 + bm2 * exp( -(lat(iCol) - bm3)**2 / bm4**2) ! (eq. 10)
       enddo
    endif
    
  end subroutine cmp_dcorr_lgth_oreopoulos
  
  ! ######################################################################################
  !
  ! ######################################################################################
!>This subroutine provides the alpha cloud overlap parameter for both RRTMG and RRTMGP
  subroutine get_alpha_exper(nCol, nLay, iovr, iovr_exprand, dzlay,    &
                             dcorr_lgth, cld_frac, alpha)
    
    ! Inputs
    integer, intent(in) :: &
         nCol,     & ! Number of horizontal grid points
         nLay        ! Number of vertical grid points
    integer, intent(in) :: &
         iovr,     &
         iovr_exprand
    real(kind_phys), dimension(:), intent(in) :: &
         dcorr_lgth  ! Cloud fraction decorrelation length (km)
    real(kind_phys), dimension(:,:), intent(in) :: &
         dzlay       !
    real(kind_phys), dimension(:,:), intent(in) ::  &
         cld_frac
    
    ! Outputs
    real(kind_phys), dimension(:,:) :: &
         alpha       ! Cloud fraction overlap parameter
    
    ! Local variables
    integer :: iCol,iLay
    
    do iCol = 1, nCol
       alpha(iCol,1) = 0.0d0
       do iLay = 2, nLay
          alpha(iCol,iLay) = exp( -(dzlay(iCol,iLay)) / dcorr_lgth(iCol))
       enddo
    enddo
   
    ! Revise alpha for exponential-random cloud overlap
    ! Decorrelate layers when a clear layer follows a cloudy layer to enforce
    ! random correlation between non-adjacent blocks of cloudy layers
    if (iovr == iovr_exprand) then
      do iLay = 2, nLay
        do iCol = 1, nCol
          if (cld_frac(iCol,iLay) == 0.0 .and. cld_frac(iCol,iLay-1) > 0.0) then
            alpha(iCol,iLay) = 0.0
          endif
        enddo
      enddo
    endif
 
    return
    
  end subroutine get_alpha_exper
  
  ! ######################################################################################
  !
  ! ######################################################################################
  subroutine get_beta_exper(nCol, nLay, iovr, iovr_exprand, dzlay,    &
                            dcorr_lgth_cond, cld_frac, beta)
    
    ! Inputs
    integer, intent(in) :: &
         nCol,     & ! Number of horizontal grid points
         nLay        ! Number of vertical grid points
    integer, intent(in) :: &
         iovr,     &
         iovr_exprand
    real(kind_phys), dimension(:), intent(in) :: &
         dcorr_lgth_cond  ! Cloud condensate decorrelation length (km)
    real(kind_phys), dimension(:,:), intent(in) :: &
         dzlay       !
    real(kind_phys), dimension(:,:), intent(in) ::  &
         cld_frac
    
    ! Outputs
    real(kind_phys), dimension(:,:), intent(out) :: &
         beta       ! Cloud condensate overlap parameter
    
    ! Local variables
    integer :: iCol,iLay

    do iCol = 1, nCol
       beta(iCol,1) = 0.0
       do iLay = 2, nLay
          beta(iCol,iLay) = exp( -(dzlay(iCol,iLay)) / dcorr_lgth_cond(iCol))
       enddo
    enddo
   
    ! Revise beta for exponential-random cloud overlap
    ! Decorrelate layers when a clear layer follows a cloudy layer to enforce
    ! random correlation between non-adjacent blocks of cloudy layers
    if (iovr == iovr_exprand) then
      do iLay = 2, nLay
        do iCol = 1, nCol
          if (cld_frac(iCol,iLay) == 0.0 .and. cld_frac(iCol,iLay-1) > 0.0) then
            beta(iCol,iLay) = 0.0
          endif
        enddo
      enddo
    endif
 
    return
    
  end subroutine get_beta_exper
!>@}  
end module module_radiation_cloud_overlap
