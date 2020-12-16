! ########################################################################################
!
! ########################################################################################
module GFS_rrtmgp_cloud_overlap_pre
  use machine,      only: kind_phys
  use rrtmgp_aux,   only: check_error_msg
  use module_radiation_cloud_overlap, only: cmp_dcorr_lgth, get_alpha_exp  

  public GFS_rrtmgp_cloud_overlap_pre_init, GFS_rrtmgp_cloud_overlap_pre_run, GFS_rrtmgp_cloud_overlap_pre_finalize

contains  
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_cloud_overlap_pre_init()
  end subroutine GFS_rrtmgp_cloud_overlap_pre_init

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_GFS_rrtmgp_cloud_overlap_pre_run
!! \htmlinclude GFS_rrtmgp_cloud_overlap_pre_run.html
!!  
  subroutine GFS_rrtmgp_cloud_overlap_pre_run(nCol, nLev, yearlen, doSWrad, doLWrad,     &
       julian, lat, p_lev, p_lay, tv_lay, con_pi, con_g, con_rd, con_epsq, dcorr_con,    &
       idcor, iovr, iovr_dcorr, iovr_exprand, iovr_exp, idcor_con, idcor_hogan,          &
       idcor_oreopoulos, cld_frac,                                                       &
       cloud_overlap_param, precip_overlap_param, de_lgth, deltaZc, errmsg, errflg)
    implicit none
    
    ! Inputs   
    integer, intent(in)    :: &
         nCol,                 & ! Number of horizontal grid points
         nLev,                 & ! Number of vertical layers
         yearlen,              & ! Length of current year (365/366) WTF?
         iovr,                 & ! Choice of cloud-overlap method
         iovr_dcorr,           & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,             & ! Flag for exponential cloud overlap method
         iovr_exprand,         & ! Flag for exponential-random cloud overlap method
         idcor,                & ! Choice of method for decorrelation length computation
         idcor_con,            & ! Flag for decorrelation-length. Use constant value
         idcor_hogan,          & ! Flag for decorrelation-length. (https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.647)
         idcor_oreopoulos        ! Flag for decorrelation-length. (10.5194/acp-12-9097-2012) 
    logical, intent(in) :: &
    	 doSWrad,              & ! Call SW radiation?
    	 doLWrad                 ! Call LW radiation
    real(kind_phys), intent(in) :: &
         julian,               & ! Julian day 
         con_pi,               & ! Physical constant: pi
         con_g,                & ! Physical constant: gravitational constant
         con_rd,               & ! Physical constant: gas-constant for dry air
         con_epsq,             & ! Physical constant: Minimum value for specific humidity
         dcorr_con               ! Decorrelation-length (used if idcor = idcor_con)
    real(kind_phys), dimension(nCol), intent(in) :: &
         lat                     ! Latitude             
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &         
         tv_lay,               & ! Virtual temperature (K)
         p_lay,                & ! Pressure at model-layers (Pa)
         cld_frac                ! Total cloud fraction
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &         
         p_lev                   ! Pressure at model-level interfaces (Pa) 
    
    ! Outputs     
    real(kind_phys), dimension(nCol),intent(out) :: &
         de_lgth                 ! Decorrelation length     
    real(kind_phys), dimension(nCol,nLev),intent(out) :: &
         cloud_overlap_param,  & ! Cloud-overlap parameter
         precip_overlap_param, & ! Precipitation overlap parameter  
         deltaZc                 ! Layer thickness (from layer-centers)(km)          
    character(len=*), intent(out) :: &
         errmsg                  ! Error message
    integer, intent(out) :: &  
         errflg                  ! Error flag
    
    ! Local variables
    real(kind_phys) :: tem1,pfac
    real(kind_phys), dimension(nLev+1) :: hgtb
    real(kind_phys), dimension(nLev)   :: hgtc
    integer :: iCol,iLay,l,iSFC,iTOA
    real(kind_phys), dimension(nCol,nLev) :: deltaZ
    logical :: top_at_1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. (doSWrad .or. doLWrad)) return

    ! What is vertical ordering?                                                                                                                                                              
    top_at_1 = (p_lev(1,1) .lt.  p_lev(1, nLev))
    if (top_at_1) then
       iSFC = nLev
       iTOA = 1
    else
       iSFC = 1
       iTOA = nLev
    endif

    !                                                                                                                                                                                        
    ! Compute layer-thickness between layer boundaries (deltaZ) and layer centers (deltaZc)                                                                                                   
    !                                                                                                                                                                                         
    do iCol=1,nCol 
       if (top_at_1) then
          ! Layer thickness (km)                                                                                                                                                              
          do iLay=1,nLev
             deltaZ(iCol,iLay) = ((con_rd/con_g)*0.001) * abs(log(p_lev(iCol,iLay+1)) - log(p_lev(iCol,iLay))) * tv_lay(iCol,iLay)
          enddo
          ! Height at layer boundaries                                                                                                                                                        
          hgtb(nLev+1) = 0._kind_phys
          do iLay=nLev,1,-1
             hgtb(iLay)= hgtb(iLay+1) + deltaZ(iCol,iLay)
          enddo
          ! Height at layer centers                                                                                                                                                           
          do iLay = nLev, 1, -1
             pfac = abs(log(p_lev(iCol,iLay+1)) - log(p_lay(iCol,iLay))) /  &
                  abs(log(p_lev(iCol,iLay+1)) - log(p_lev(iCol,iLay)))
             hgtc(iLay) = hgtb(iLay+1) + pfac * (hgtb(iLay) - hgtb(iLay+1))
          enddo
          ! Layer thickness between centers                                                                                                                                                   
          do iLay = nLev-1, 1, -1
             deltaZc(iCol,iLay) = hgtc(iLay) - hgtc(iLay+1)
          enddo
          deltaZc(iCol,nLev) = hgtc(nLev) - hgtb(nLev+1)
       else
          do iLay=nLev,1,-1
             deltaZ(iCol,iLay) = ((con_rd/con_g)*0.001) * abs(log(p_lev(iCol,iLay))  - log(p_lev(iCol,iLay+1))) * tv_lay(iCol,iLay)
          enddo
          ! Height at layer boundaries                                                                                                                                                        
          hgtb(1) = 0._kind_phys
          do iLay=1,nLev
             hgtb(iLay+1)= hgtb(iLay) + deltaZ(iCol,iLay)
          enddo
          ! Height at layer centers                                                                                                                                                           
          do iLay = 1, nLev
             pfac = abs(log(p_lev(iCol,iLay)) - log(p_lay(iCol,iLay)  )) /  &
                  abs(log(p_lev(iCol,iLay)) - log(p_lev(iCol,iLay+1)))
             hgtc(iLay) = hgtb(iLay) + pfac * (hgtb(iLay+1) - hgtb(iLay))
          enddo
          ! Layer thickness between centers                                                                                                                                                   
          do iLay = 2, nLev
             deltaZc(iCol,iLay) = hgtc(iLay) - hgtc(iLay-1)
          enddo
          deltaZc(iCol,1) = hgtc(1) - hgtb(1)
       endif
    enddo

    !                                                                                                                                                                                         
    ! Cloud decorrelation length                                                                                                                                                              
    !                                                                                                                                                                                         
    if (idcor == idcor_hogan) then
       call cmp_dcorr_lgth(nCol, lat, con_pi, de_lgth)
    endif
    if (idcor == idcor_oreopoulos) then
       call cmp_dcorr_lgth(nCol, lat*(180._kind_phys/con_pi), julian, yearlen, de_lgth)
    endif
    if (idcor == idcor_con) then
       de_lgth(:) = dcorr_con
    endif

    !
    ! Cloud overlap parameter
    !
    if (iovr == iovr_dcorr .or. iovr == iovr_exp .or. iovr == iovr_exprand) then
       call get_alpha_exp(nCol, nLev, deltaZc, de_lgth, cloud_overlap_param)
    else
       de_lgth(:)               = 0.
       cloud_overlap_param(:,:) = 0.
    endif

    ! For exponential random overlap...                                                                                                                                                      
    ! Decorrelate layers when a clear layer follows a cloudy layer to enforce                                                                                                                 
    ! random correlation between non-adjacent blocks of cloudy layers                                                                                                                         
    if (iovr == iovr_exprand) then
       do iLay = 1, nLev
          do iCol = 1, nCol
             if (cld_frac(iCol,iLay) .eq. 0. .and. cld_frac(iCol,iLay-1) .gt. 0.) then
                cloud_overlap_param(iCol,iLay) = 0._kind_phys
             endif
          enddo
       enddo
    endif

    ! 
    ! Compute precipitation overlap parameter (Hack. Using same as cloud for now)
    !
    precip_overlap_param = cloud_overlap_param    
    
  end subroutine GFS_rrtmgp_cloud_overlap_pre_run
  
  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_cloud_overlap_pre_finalize()
  end subroutine GFS_rrtmgp_cloud_overlap_pre_finalize
end module GFS_rrtmgp_cloud_overlap_pre
