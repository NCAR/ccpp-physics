!>\file GFS_cloud_diagnostics.F90
!!

module GFS_cloud_diagnostics
  use machine,                 only: kind_phys
  use module_radiation_clouds, only: gethml

  ! Module parameters (imported directly from radiation_cloud.f)
  integer, parameter :: &
       NF_CLDS = 9,        &  ! Number of fields in cloud array
       NK_CLDS = 3            ! Number of cloud vertical domains
  real(kind_phys), parameter :: &
       climit = 0.001,      & ! Lowest allowable cloud-fraction
       ovcst  = 1.0 - 1.0e-8  ! Overcast cloud-fraction 0.999999999
  real(kind_phys), parameter, dimension(NK_CLDS+1,2) :: &
       ptopc = reshape(source=(/ 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /),  &
                       shape=(/NK_CLDS+1,2/)) 
                      
  ! Version tag and last revision date
  character(40), parameter :: VTAGCLD='UFS-cloud-diagnostics    vX.x May 2020 '
     
  public GFS_cloud_diagnostics_run

contains

!>\defgroup gfs_cloud_diagnostics_mod GFS Cloud Diagnostics Module
!> This module contains code to produce the UFS High/Mid/Low cloud-diagnostics.
!! This was bundled together with the prognostic cloud modules within the RRTMG implementation.
!! For the RRTMGP implementation we propose to keep these diagnostics independent.
!> @{
!> \section arg_table_GFS_cloud_diagnostics_run
!! \htmlinclude GFS_cloud_diagnostics_run.html
!!  
  subroutine GFS_cloud_diagnostics_run(nCol, nLev, iovr, iovr_rand, iovr_maxrand,        &
       iovr_max, iovr_dcorr, iovr_exp, iovr_exprand, lsswr, lslwr, lat, de_lgth, p_lay,  &
       cld_frac, p_lev, deltaZ, cloud_overlap_param, precip_overlap_param, con_pi,       &
       top_at_1, si, mtopa, mbota, cldsa, errmsg, errflg)
    implicit none
     
    ! Inputs 
    integer, intent(in) ::     &
         nCol,                 & ! Number of horizontal grid-points
         nLev                    ! Number of vertical-layers
    integer, intent(in) ::     &
         iovr,                 & ! Choice of cloud-overlap method
         iovr_rand,            & ! Flag for random cloud overlap method
         iovr_maxrand,         & ! Flag for maximum-random cloud overlap method
         iovr_max,             & ! Flag for maximum cloud overlap method
         iovr_dcorr,           & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,             & ! Flag for exponential cloud overlap method
         iovr_exprand            ! Flag for exponential-random cloud overlap method
    logical, intent(in) ::     &
    	 lsswr,                & ! Call SW radiation?
    	 lslwr,                & ! Call LW radiation?
         top_at_1                ! Vertical ordering flag
    real(kind_phys), intent(in) :: &
         con_pi                  ! Physical constant: pi
    real(kind_phys), dimension(:), intent(in) ::   &
         lat,                  & ! Latitude
         de_lgth,              & ! Decorrelation length
         si                      ! Vertical sigma coordinate
    real(kind_phys), dimension(:,:), intent(in) :: &
         p_lay,                & ! Pressure at model-layer
         cld_frac                ! Total cloud fraction
    real(kind_phys), dimension(:,:), intent(in) :: &
         p_lev                   ! Pressure at model interfaces         
    real(kind_phys), dimension(:,:), intent(in) :: &
    	 deltaZ,               & ! Layer thickness (m)
         cloud_overlap_param,  & ! Cloud-overlap parameter
         precip_overlap_param    ! Precipitation overlap parameter
    
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                  ! Error message
    integer, intent(out) :: &  
         errflg                  ! Error flag
    integer,dimension(:,:),intent(out) :: &
         mbota,                & ! Vertical indices for cloud tops
         mtopa                   ! Vertical indices for cloud bases
    real(kind_phys),dimension(:,:), intent(out) :: &
         cldsa                   ! Fraction of clouds for low, middle, high, total and BL 
    
    ! Local variables
    integer i,id,iCol,iLay,icld
    real(kind_phys) :: tem1
    real(kind_phys),dimension(nCol,NK_CLDS+1) :: ptop1
    real(kind_phys),dimension(nCol) :: rlat
    real(kind_phys),dimension(nCol,nLev) :: cldcnv
	
    if (.not. (lsswr .or. lslwr)) return
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0	
    
    ! This is set to zero in all of the progcld() routines and passed to gethml().
    cldcnv(:,:) = 0._kind_phys
    
    do icld = 1, NK_CLDS+1
       tem1 = ptopc(icld,2) - ptopc(icld,1)
       do i=1,nCol
          rlat(i) = abs(lat(i) / con_pi )
          ptop1(i,icld) = ptopc(icld,1) + tem1*max( 0.0, 4.0*rlat(i)-1.0 )
       enddo
    enddo
	
    ! Compute low, mid, high, total, and boundary layer cloud fractions and clouds top/bottom 
    ! layer indices for low, mid, and high clouds. The three cloud domain boundaries are 
    ! defined by ptopc. The cloud overlapping method is defined by control flag 'iovr', which may
    ! be different for lw and sw radiation programs.
    call gethml(p_lay*0.01, ptop1, cld_frac, cldcnv, deltaZ, de_lgth, cloud_overlap_param,&
         nCol, nLev, iovr, iovr_rand, iovr_maxrand, iovr_max, iovr_dcorr, iovr_exp,       &
         iovr_exprand, top_at_1, si, cldsa, mtopa, mbota)
    
  end subroutine GFS_cloud_diagnostics_run
!> @}
end module GFS_cloud_diagnostics
