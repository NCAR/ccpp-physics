! ########################################################################################
! This module contains code to produce the UFS High/Mid/Low cloud-diagnostics. 
! This was bundled together with the prognostic cloud modules within the RRTMG implementation.
! For the RRTMGP implementation we propose to keep these diagnostics independent.
! ########################################################################################
module GFS_cloud_diagnostics
  use machine,                 only: kind_phys
  use physparam,               only: icldflg
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
  
  ! Module variables
  integer :: &
       llyr    = 2            ! Upper limit of boundary layer clouds
     
  public GFS_cloud_diagnostics_run, GFS_cloud_diagnostics_init,&
       GFS_cloud_diagnostics_finalize, hml_cloud_diagnostics_init
contains
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_cloud_diagnostics_init()
  end subroutine GFS_cloud_diagnostics_init
  
  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_GFS_cloud_diagnostics_run
!! \htmlinclude GFS_cloud_diagnostics_run.html
!!  
  subroutine GFS_cloud_diagnostics_run(nCol, nLev, lsswr, lslwr, lat, de_lgth, p_lay,    &
       cld_frac, p_lev, deltaZ, cloud_overlap_param, precip_overlap_param, con_pi,       &
       mbota, mtopa, cldsa, errmsg, errflg)
    implicit none
     
    ! Inputs 
    integer, intent(in) :: &
         nCol,              & ! Number of horizontal grid-points
         nLev                 ! Number of vertical-layers
    logical, intent(in) :: &
    	 lsswr,             & ! Call SW radiation?
    	 lslwr                ! Call LW radiation 
    real(kind_phys), intent(in) :: &
         con_pi               ! Physical constant: pi  
    real(kind_phys), dimension(nCol), intent(in) :: &
         lat,               & ! Latitude       
         de_lgth              ! Decorrelation length     
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &
         p_lay,             & ! Pressure at model-layer
         cld_frac             ! Total cloud fraction
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &
         p_lev                ! Pressure at model interfaces         
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &
    	 deltaZ,              & ! Layer thickness (km)
         cloud_overlap_param, & ! Cloud-overlap parameter
         precip_overlap_param   ! Precipitation overlap parameter
    
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                 ! Error message
    integer, intent(out) :: &  
         errflg                 ! Error flag
    integer,dimension(ncol,3),intent(out) :: &
         mbota,               & ! Vertical indices for cloud tops
         mtopa                  ! Vertical indices for cloud bases
    real(kind_phys), dimension(ncol,5), intent(out) :: &
         cldsa                  ! Fraction of clouds for low, middle, high, total and BL 
    
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
    call gethml(p_lay/100., ptop1, cld_frac, cldcnv, deltaZ, de_lgth, cloud_overlap_param,&
         nCol, nLev, cldsa, mtopa, mbota)	
    
  end subroutine GFS_cloud_diagnostics_run
  
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_cloud_diagnostics_finalize()
  end subroutine GFS_cloud_diagnostics_finalize
  
  ! ######################################################################################
  ! Initialization routine for High/Mid/Low cloud diagnostics. 
  ! ######################################################################################
  subroutine hml_cloud_diagnostics_initialize(imp_physics, imp_physics_fer_hires,        &
          imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,                      &
          imp_physics_zhao_carr, imp_physics_zhao_carr_pdf, imp_physics_mg, nLev,        &
          mpi_rank, sigmainit, errflg)
    implicit none
    ! Inputs
    integer, intent(in) :: &
          imp_physics,               & ! Flag for MP scheme
          imp_physics_fer_hires,     & ! Flag for fer-hires scheme
          imp_physics_gfdl,          & ! Flag for gfdl scheme
          imp_physics_thompson,      & ! Flag for thompsonscheme
          imp_physics_wsm6,          & ! Flag for wsm6 scheme
          imp_physics_zhao_carr,     & ! Flag for zhao-carr scheme
          imp_physics_zhao_carr_pdf, & ! Flag for zhao-carr+PDF scheme
          imp_physics_mg               ! Flag for MG scheme
    integer, intent(in) :: &
         nLev,        & ! Number of vertical-layers
         mpi_rank 
    real(kind_phys), dimension(nLev+1), intent(in) :: &
         sigmainit
    ! Outputs
    integer, intent(out) :: &
    	errflg
    
    ! Local variables
    integer :: iLay, kl
 
    ! Initialize error flag
    errflg = 0
    
    if (mpi_rank == 0) print *, VTAGCLD      !print out version tag
    
    if ( icldflg == 0 ) then
       print *,' - Diagnostic Cloud Method has been discontinued'
       errflg = 1
    else
       if (mpi_rank == 0) then
          print *,' - Using Prognostic Cloud Method'
          if (imp_physics == imp_physics_zhao_carr) then
             print *,'   --- Zhao/Carr/Sundqvist microphysics'
          elseif (imp_physics == imp_physics_zhao_carr_pdf) then
             print *,'   --- zhao/carr/sundqvist + pdf cloud'
          elseif (imp_physics == imp_physics_gfdl) then
             print *,'   --- GFDL Lin cloud microphysics'
          elseif (imp_physics == imp_physics_thompson) then
             print *,'   --- Thompson cloud microphysics'
          elseif (imp_physics == imp_physics_wsm6) then
             print *,'   --- WSM6 cloud microphysics'
          elseif (imp_physics == imp_physics_mg) then
             print *,'   --- MG cloud microphysics'
          elseif (imp_physics == imp_physics_fer_hires) then
             print *,'   --- Ferrier-Aligo cloud microphysics'
          else
             print *,'  !!! ERROR in cloud microphysc specification!!!', &
                  '  imp_physics (NP3D) =',imp_physics
             errflg = 1
          endif
       endif
    endif
    
    ! Compute the top of BL cld (llyr), which is the topmost non cld(low) layer for 
    ! stratiform (at or above lowest 0.1 of the atmosphere).
    lab_do_k0 : do iLay = nLev, 2, -1
       kl = iLay
       if (sigmainit(iLay) < 0.9e0) exit lab_do_k0
    enddo  lab_do_k0
    llyr = kl      
    
    return
  end subroutine hml_cloud_diagnostics_initialize
end module GFS_cloud_diagnostics
