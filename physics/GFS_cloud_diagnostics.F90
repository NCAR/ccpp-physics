! ########################################################################################
! This module contains code to produce the UFS High/Mid/Low cloud-diagnostics. 
! This was bundled together with the prognostic cloud modules within the RRTMG implementation.
! For the RRTMGP implementation we propose to keep these diagnostics independent.
! ########################################################################################
module GFS_cloud_diagnostics
  use machine,                 only: kind_phys
  use physparam,               only: iovrlw, iovrsw, ivflip, icldflg, idcor
  
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
       iovr    = 1,        &  ! Cloud overlap used for diagnostic HML cloud outputs
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
    
    ! Cloud overlap used for diagnostic HML cloud outputs      
    iovr = max(iovrsw,iovrlw)  
    
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
  
  ! #########################################################################################
  ! #########################################################################################
  subroutine gethml(plyr, ptop1, cldtot, cldcnv, dz, de_lgth, alpha, IX, NLAY, clds, mtop, mbot)
    !  ===================================================================  !
    !                                                                       !
    ! abstract: compute high, mid, low, total, and boundary cloud fractions !
    !   and cloud top/bottom layer indices for model diagnostic output.     !
    !   the three cloud domain boundaries are defined by ptopc.  the cloud  !
    !   overlapping method is defined by control flag 'iovr', which is also !
    !   used by lw and sw radiation programs.                               !
    !                                                                       !
    ! usage:         call gethml                                            !
    !                                                                       !
    ! subprograms called:  none                                             !
    !                                                                       !
    ! attributes:                                                           !
    !   language:   fortran 90                                              !
    !   machine:    ibm-sp, sgi                                             !
    !                                                                       !
    !                                                                       !
    !  ====================  definition of variables  ====================  !
    !                                                                       !
    ! input variables:                                                      !
    !   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
    !   ptop1 (IX,4)    : pressure limits of cloud domain interfaces        !
    !                     (sfc,low,mid,high) in mb (100Pa)                  !
    !   cldtot(IX,NLAY) : total or straiform cloud profile in fraction      !
    !   cldcnv(IX,NLAY) : convective cloud (for diagnostic scheme only)     !
    !   dz    (ix,nlay) : layer thickness (km)                              !
    !   de_lgth(ix)     : clouds vertical de-correlation length (km)        !
    !   alpha(ix,nlay)  : alpha decorrelation parameter                     !
    !   IX              : horizontal dimention                              !
    !   NLAY            : vertical layer dimensions                         !
    !                                                                       !
    ! output variables:                                                     !
    !   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
    !   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
    !   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
    !                                                                       !
    ! external module variables:  (in physparam)                            !
    !   ivflip          : control flag of vertical index direction          !
    !                     =0: index from toa to surface                     !
    !                     =1: index from surface to toa                     !
    !                                                                       !
    ! internal module variables:                                            !
    !   iovr            : control flag for cloud overlap                    !
    !                     =0 random overlapping clouds                      !
    !                     =1 max/ran overlapping clouds                     !
    !                     =2 maximum overlapping  ( for mcica only )        !
    !                     =3 decorr-length ovlp   ( for mcica only )        !
    !                     =4 exponential cloud overlap  (AER; mcica only)  !
    !                     =5 exponential-random overlap (AER; mcica only)  !
    !                                                                       !
    !  ====================    end of description    =====================  !
    !
    implicit none!
    
    !  ---  inputs:
    integer, intent(in) :: IX, NLAY
    
    real (kind=kind_phys), dimension(:,:), intent(in) :: plyr, ptop1, &
         cldtot, cldcnv, dz
    real (kind=kind_phys), dimension(:),   intent(in) :: de_lgth
    real (kind=kind_phys), dimension(:,:), intent(in) :: alpha
    
    !  ---  outputs
    real (kind=kind_phys), dimension(:,:), intent(out) :: clds
    
    integer,               dimension(:,:), intent(out) :: mtop, mbot
    
    !  ---  local variables:
    real (kind=kind_phys) :: cl1(IX), cl2(IX), dz1(ix)
    real (kind=kind_phys) :: pcur, pnxt, ccur, cnxt, alfa
    
    integer, dimension(IX):: idom, kbt1, kth1, kbt2, kth2
    integer :: i, k, id, id1, kstr, kend, kinc
    
    !
    !===> ... begin here
    !
    clds(:,:) = 0.0
    
    do i = 1, IX
       cl1(i) = 1.0
       cl2(i) = 1.0
    enddo
    
    !  ---  total and bl clouds, where cl1, cl2 are fractions of clear-sky view
    !       layer processed from surface and up
    
    !> - Calculate total and BL cloud fractions (maximum-random cloud
    !!    overlapping is operational).
    
    if ( ivflip == 0 ) then                   ! input data from toa to sfc
       kstr = NLAY
       kend = 1
       kinc = -1
    else                                      ! input data from sfc to toa
       kstr = 1
       kend = NLAY
       kinc = 1
    endif                                     ! end_if_ivflip
    
    if ( iovr == 0 ) then                     ! random overlap
       
       do k = kstr, kend, kinc
          do i = 1, IX
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) cl1(i) = cl1(i) * (1.0 - ccur)
          enddo
          
          if (k == llyr) then
             do i = 1, IX
                clds(i,5) = 1.0 - cl1(i)          ! save bl cloud
             enddo
          endif
       enddo
       
       do i = 1, IX
          clds(i,4) = 1.0 - cl1(i)              ! save total cloud
       enddo
       
    elseif ( iovr == 1 ) then                 ! max/ran overlap
       
       do k = kstr, kend, kinc
          do i = 1, IX
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) then             ! cloudy layer
                cl2(i) = min( cl2(i), (1.0 - ccur) )
             else                                ! clear layer
                cl1(i) = cl1(i) * cl2(i)
                cl2(i) = 1.0
             endif
          enddo
          
          if (k == llyr) then
             do i = 1, IX
                clds(i,5) = 1.0 - cl1(i) * cl2(i) ! save bl cloud
             enddo
          endif
       enddo
       
       do i = 1, IX
          clds(i,4) = 1.0 - cl1(i) * cl2(i)     ! save total cloud
       enddo
       
    elseif ( iovr == 2 ) then                 ! maximum overlap all levels
       
       cl1(:) = 0.0
       
       do k = kstr, kend, kinc
          do i = 1, IX
             ccur = min( ovcst,  max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) cl1(i) = max( cl1(i), ccur )
          enddo
          
          if (k == llyr) then
             do i = 1, IX
                clds(i,5) = cl1(i)    ! save bl cloud
             enddo
          endif
       enddo
       
       do i = 1, IX
          clds(i,4) = cl1(i)        ! save total cloud
       enddo
       
    elseif ( iovr == 3 ) then                 ! random if clear-layer divided,
                                              ! otherwise de-corrlength method
       do i = 1, ix
          dz1(i) = - dz(i,kstr)
       enddo
       
       do k = kstr, kend, kinc
          do i = 1, ix
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) then                           ! cloudy layer
                alfa = exp( -0.5*((dz1(i)+dz(i,k)))/de_lgth(i) )
                dz1(i) = dz(i,k)
                cl2(i) =    alfa      * min(cl2(i), (1.0 - ccur))         & ! maximum part
                     + (1.0 - alfa) * (cl2(i) * (1.0 - ccur))             ! random part
             else                                               ! clear layer
                cl1(i) = cl1(i) * cl2(i)
                cl2(i) = 1.0
                if (k /= kend) dz1(i) = -dz(i,k+kinc)
             endif
          enddo
          
          if (k == llyr) then
             do i = 1, ix
                clds(i,5) = 1.0 - cl1(i) * cl2(i) ! save bl cloud
             enddo
          endif
       enddo
       
       do i = 1, ix
          clds(i,4) = 1.0 - cl1(i) * cl2(i)     ! save total cloud
       enddo
       
    elseif ( iovr == 4 .or. iovr == 5 ) then  ! exponential overlap (iovr=4), or
                                              ! exponential-random  (iovr=5);
                                              ! distinction defined by alpha
       do k = kstr, kend, kinc
          do i = 1, ix
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) then                                   ! cloudy layer
                cl2(i) =   alpha(i,k) * min(cl2(i), (1.0 - ccur))     & ! maximum part
                     + (1.0 - alpha(i,k)) * (cl2(i) * (1.0 - ccur))     ! random part
             else                                                       ! clear layer
                cl1(i) = cl1(i) * cl2(i)
                cl2(i) = 1.0
             endif
          enddo
          if (k == llyr) then
             do i = 1, ix
                clds(i,5) = 1.0 - cl1(i) * cl2(i) ! save bl cloud
             enddo
          endif
       enddo
       do i = 1, ix
          clds(i,4) = 1.0 - cl1(i) * cl2(i)     ! save total cloud
       enddo
    endif                                     ! end_if_iovr
    
    !  ---  high, mid, low clouds, where cl1, cl2 are cloud fractions
    !       layer processed from one layer below llyr and up
    !  ---  change! layer processed from surface to top, so low clouds will
    !       contains both bl and low clouds.
    
    !> - Calculte high, mid, low cloud fractions and vertical indices of
    !!    cloud tops/bases.
    if ( ivflip == 0 ) then                   ! input data from toa to sfc
       
       do i = 1, IX
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = NLAY
          kbt2(i) = NLAY
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = NLAY
          mtop(i,1) = NLAY
          mbot(i,2) = NLAY - 1
          mtop(i,2) = NLAY - 1
          mbot(i,3) = NLAY - 1
          mtop(i,3) = NLAY - 1
       enddo
       
       !org    do k = llyr-1, 1, -1
       do k = NLAY, 1, -1
          do i = 1, IX
             id = idom(i)
             id1= id + 1
             
             pcur = plyr(i,k)
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             
             if (k > 1) then
                pnxt = plyr(i,k-1)
                cnxt = min( ovcst, max( cldtot(i,k-1), cldcnv(i,k-1) ))
             else
                pnxt = -1.0
                cnxt = 0.0
             endif
             
             if (pcur < ptop1(i,id1)) then
                id = id + 1
                id1= id1 + 1
                idom(i) = id
             endif
             
             if (ccur >= climit) then
                if (kth2(i) == 0) kbt2(i) = k
                kth2(i) = kth2(i) + 1
                
                if ( iovr == 0 ) then
                   cl2(i) = cl2(i) + ccur - cl2(i)*ccur
                else
                   cl2(i) = max( cl2(i), ccur )
                endif
                
                if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                   kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i) )      &
                                     / (cl1(i) + cl2(i)) )
                   kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i) )      &
                                     / (cl1(i) + cl2(i)) )
                   cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)
                   
                   kbt2(i) = k - 1
                   kth2(i) = 0
                   cl2 (i) = 0.0
                endif   ! end_if_cnxt_or_pnxt
             endif     ! end_if_ccur
             
             if (pnxt < ptop1(i,id1)) then
                clds(i,id) = cl1(i)
                mtop(i,id) = min( kbt1(i), kbt1(i)-kth1(i)+1 )
                mbot(i,id) = kbt1(i)
                
                cl1 (i) = 0.0
                kbt1(i) = k - 1
                kth1(i) = 0
                
                if (id1 <= NK_CLDS) then
                   mbot(i,id1) = kbt1(i)
                   mtop(i,id1) = kbt1(i)
                endif
             endif     ! end_if_pnxt
             
          enddo       ! end_do_i_loop
       enddo         ! end_do_k_loop
       
    else                                      ! input data from sfc to toa
       
       do i = 1, IX
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = 1
          kbt2(i) = 1
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = 1
          mtop(i,1) = 1
          mbot(i,2) = 2
          mtop(i,2) = 2
          mbot(i,3) = 2
          mtop(i,3) = 2
       enddo
       
       !org    do k = llyr+1, NLAY
       do k = 1, NLAY
          do i = 1, IX
             id = idom(i)
             id1= id + 1
             
             pcur = plyr(i,k)
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             
             if (k < NLAY) then
                pnxt = plyr(i,k+1)
                cnxt = min( ovcst, max( cldtot(i,k+1), cldcnv(i,k+1) ))
             else
                pnxt = -1.0
                cnxt = 0.0
             endif
             
             if (pcur < ptop1(i,id1)) then
                id = id + 1
                id1= id1 + 1
                idom(i) = id
             endif
             
             if (ccur >= climit) then
                if (kth2(i) == 0) kbt2(i) = k
                kth2(i) = kth2(i) + 1
                
                if ( iovr == 0 ) then
                   cl2(i) = cl2(i) + ccur - cl2(i)*ccur
                else
                   cl2(i) = max( cl2(i), ccur )
                endif
                
                if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                   kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))       &
                                           / (cl1(i) + cl2(i)) )
                   kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))       &
                                           / (cl1(i) + cl2(i)) )
                   cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)
                   
                   kbt2(i) = k + 1
                   kth2(i) = 0
                   cl2 (i) = 0.0
                endif     ! end_if_cnxt_or_pnxt
             endif       ! end_if_ccur
             
             if (pnxt < ptop1(i,id1)) then
                clds(i,id) = cl1(i)
                mtop(i,id) = max( kbt1(i), kbt1(i)+kth1(i)-1 )
                mbot(i,id) = kbt1(i)
                
                cl1 (i) = 0.0
                kbt1(i) = min(k+1, nlay)
                kth1(i) = 0
                
                if (id1 <= NK_CLDS) then
                   mbot(i,id1) = kbt1(i)
                   mtop(i,id1) = kbt1(i)
                endif
             endif     ! end_if_pnxt
             
          enddo       ! end_do_i_loop
       enddo         ! end_do_k_loop
       
    endif                                     ! end_if_ivflip
    
    !
    return
    !...................................
  end subroutine gethml
end module GFS_cloud_diagnostics
