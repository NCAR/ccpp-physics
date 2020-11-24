! ########################################################################################
! This module contains the interface between the GFDL macrophysics and the RRTMGP radiation
! schemes. Only compatable with Model%imp_physics = Model%imp_physics_gfdl
! ########################################################################################
module GFS_rrtmgp_gfdlmp_pre
  use machine,      only: kind_phys
  use rrtmgp_aux,   only: check_error_msg
  use module_radiation_cloud_overlap, only: cmp_dcorr_lgth, get_alpha_exp  

  ! Parameters
  real(kind_phys), parameter :: &
       reliq_def  = 10.0 ,      & ! Default liq radius to 10 micron (used when effr_in=F)
       reice_def  = 50.0,       & ! Default ice radius to 50 micron (used when effr_in=F)
       rerain_def = 1000.0,     & ! Default rain radius to 1000 micron (used when effr_in=F)
       resnow_def = 250.0,      & ! Default snow radius to 250 micron (used when effr_in=F)    
       reice_min  = 10.0,       & ! Minimum ice size allowed by scheme
       reice_max  = 150.0,      & ! Maximum ice size allowed by scheme
       cllimit    = 0.001         ! Lowest cloud fraction in GFDL MP scheme
         
   public GFS_rrtmgp_gfdlmp_pre_init, GFS_rrtmgp_gfdlmp_pre_run, GFS_rrtmgp_gfdlmp_pre_finalize

contains  
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_gfdlmp_pre_init()
  end subroutine GFS_rrtmgp_gfdlmp_pre_init

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_GFS_rrtmgp_gfdlmp_pre_run
!! \htmlinclude GFS_rrtmgp_gfdlmp_pre_run.html
!!  
  subroutine GFS_rrtmgp_gfdlmp_pre_run(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice,   &
       i_cldrain, i_cldsnow, i_cldgrpl, i_cldtot, yearlen, doSWrad, doLWrad, effr_in,    &
       julian, lat, p_lev, p_lay, tv_lay, effrin_cldliq, effrin_cldice, effrin_cldrain,  &
       effrin_cldsnow, tracer, con_pi, con_g, con_rd, con_epsq, dcorr_con, idcor, iovr,  &
       iovr_dcorr, iovr_exprand, iovr_exp, idcor_con, idcor_hogan, idcor_oreopoulos,     &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,   &
       cld_rerain, precip_frac, cloud_overlap_param, precip_overlap_param, de_lgth,      &
       deltaZb, errmsg, errflg)
    implicit none
    
    ! Inputs   
    integer, intent(in)    :: &
         nCol,                 & ! Number of horizontal grid points
         nLev,                 & ! Number of vertical layers
         ncnd,                 & ! Number of cloud condensation types.
         nTracers,             & ! Number of tracers from model. 
         i_cldliq,             & ! Index into tracer array for cloud liquid. 
         i_cldice,             & ! Index into tracer array for cloud ice.
         i_cldrain,            & ! Index into tracer array for cloud rain.
         i_cldsnow,            & ! Index into tracer array for cloud snow.
         i_cldgrpl,            & ! Index into tracer array for cloud groupel.
         i_cldtot,             & ! Index into tracer array for cloud total amount.
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
    	 doLWrad,              & ! Call LW radiation
    	 effr_in                 ! Provide hydrometeor radii from macrophysics?
    real(kind_phys), intent(in) :: &
         julian,               & ! Julian day 
         con_pi,               & ! Physical constant: pi
         con_g,                & ! Physical constant: gravitational constant
         con_rd,               & ! Physical constant: gas-constant for dry air
         con_epsq,             & ! Physical constant(?): Minimum value for specific humidity
         dcorr_con               ! Decorrelation-length (used if idcor = idcor_con)
    real(kind_phys), dimension(nCol), intent(in) :: &
         lat                     ! Latitude             
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &         
         tv_lay,               & ! Virtual temperature (K)
         p_lay,                & ! Pressure at model-layers (Pa)
         effrin_cldliq,        & ! Effective radius for liquid cloud-particles (microns)
         effrin_cldice,        & ! Effective radius for ice cloud-particles (microns)
         effrin_cldrain,       & ! Effective radius for rain cloud-particles (microns)
         effrin_cldsnow          ! Effective radius for snow cloud-particles (microns)
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &         
         p_lev                   ! Pressure at model-level interfaces (Pa)
    real(kind_phys), dimension(nCol, nLev, nTracers),intent(in) :: &
         tracer                  ! Cloud condensate amount in layer by type ()         
    
    ! Outputs     
    real(kind_phys), dimension(nCol),intent(out) :: &
         de_lgth                 ! Decorrelation length     
    real(kind_phys), dimension(nCol,nLev),intent(out) :: &
         cld_frac,             & ! Total cloud fraction
         cld_lwp,              & ! Cloud liquid water path
         cld_reliq,            & ! Cloud liquid effective radius
         cld_iwp,              & ! Cloud ice water path
         cld_reice,            & ! Cloud ice effecive radius
         cld_swp,              & ! Cloud snow water path
         cld_resnow,           & ! Cloud snow effective radius
         cld_rwp,              & ! Cloud rain water path
         cld_rerain,           & ! Cloud rain effective radius       
         precip_frac,          & ! Precipitation fraction
         cloud_overlap_param,  & ! Cloud-overlap parameter
         precip_overlap_param, & ! Precipitation overlap parameter  
         deltaZb                 ! Layer thickness (km)          
    character(len=*), intent(out) :: &
         errmsg                  ! Error message
    integer, intent(out) :: &  
         errflg                  ! Error flag
    
    ! Local variables
    real(kind_phys) :: tem1,pfac
    real(kind_phys), dimension(nLev+1) :: hgtb
    real(kind_phys), dimension(nLev)   :: hgtc
    real(kind_phys), dimension(nCol, nLev, min(4,ncnd)) :: cld_condensate
    integer :: iCol,iLay,l,ncndl,iSFC,iTOA
    real(kind_phys), dimension(nCol,nLev) :: deltaP,deltaZ
    logical :: top_at_1

    if (.not. (doSWrad .or. doLWrad)) return
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    ! Test inputs
    if (ncnd .ne. 5) then
       errmsg = 'Incorrect number of cloud condensates provided'
       errflg = 1
       call check_error_msg('GFS_rrtmgp_gfdlmp_pre_run',errmsg)
       return
    endif

    ! What is vertical ordering?                                                                                                                                                              
    top_at_1 = (p_lev(1,1) .lt.  p_lev(1, nLev))
    if (top_at_1) then
       iSFC = nLev
       iTOA = 1
    else
       iSFC = 1
       iTOA = nLev
    endif

    ! Initialize outputs                                                                                                                                                                      
    cld_lwp(:,:)    = 0.0
    cld_reliq(:,:)  = reliq_def
    cld_iwp(:,:)    = 0.0
    cld_reice(:,:)  = reice_def
    cld_rwp(:,:)    = 0.0
    cld_rerain(:,:) = rerain_def
    cld_swp(:,:)    = 0.0
    cld_resnow(:,:) = resnow_def
    
    ! ####################################################################################
    ! Pull out cloud information for GFDL MP scheme.
    ! ####################################################################################
    ! Condensate
    cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,i_cldliq)     ! -liquid water
    cld_condensate(1:nCol,1:nLev,2) = tracer(1:nCol,1:nLev,i_cldice)     ! -ice water
    cld_condensate(1:nCol,1:nLev,3) = tracer(1:nCol,1:nLev,i_cldrain)    ! -rain water
    cld_condensate(1:nCol,1:nLev,4) = tracer(1:nCol,1:nLev,i_cldsnow) + &! -snow + grapuel
                                      tracer(1:nCol,1:nLev,i_cldgrpl) 
                           
    ! Since we combine the snow and grapuel, define local variable for number of condensate types.
    ncndl = min(4,ncnd)

    ! Set really tiny suspended particle amounts to clear
    do l=1,ncndl
       do iLay=1,nLev
          do iCol=1,nCol   
             if (cld_condensate(iCol,iLay,l) < con_epsq) cld_condensate(iCol,iLay,l) = 0.0
          enddo
       enddo
    enddo
    
    ! Cloud-fraction
    cld_frac(1:nCol,1:nLev)   = tracer(1:nCol,1:nLev,i_cldtot)
    
    ! Precipitation fraction (Hack. For now use cloud-fraction)
    precip_frac(1:nCol,1:nLev) = tracer(1:nCol,1:nLev,i_cldtot)
    
    ! Condensate and effective size
    deltaP = abs(p_lev(:,2:nLev+1)-p_lev(:,1:nLev))/100.  
    do iLay = 1, nLev
       do iCol = 1, nCol
          ! Compute liquid/ice condensate path from mixing ratios (kg/kg)->(g/m2)   
          if (cld_frac(iCol,iLay) .ge. cllimit) then         
             tem1          = (1.0e5/con_g) * deltaP(iCol,iLay)
             cld_lwp(iCol,iLay)  = cld_condensate(iCol,iLay,1) * tem1
             cld_iwp(iCol,iLay)  = cld_condensate(iCol,iLay,2) * tem1
             cld_rwp(iCol,iLay)  = cld_condensate(iCol,iLay,3) * tem1
             cld_swp(iCol,iLay)  = cld_condensate(iCol,iLay,4) * tem1  
          endif
          ! Use radii provided from the macrophysics        
          if (effr_in) then
             cld_reliq(iCol,iLay)  = effrin_cldliq(iCol,iLay)
             cld_reice(iCol,iLay)  = max(reice_min, min(reice_max,effrin_cldice(iCol,iLay)))
             cld_rerain(iCol,iLay) = effrin_cldrain(iCol,iLay)
             cld_resnow(iCol,iLay) = effrin_cldsnow(iCol,iLay)
          else
             cld_reliq(iCol,iLay)  = reliq_def
             cld_reice(iCol,iLay)  = reice_def
             cld_rerain(iCol,iLay) = rerain_def
             cld_resnow(iCol,iLay) = resnow_def
          endif
       enddo
    enddo
    
    ! ####################################################################################
    ! Cloud (and precipitation) overlap
    ! ####################################################################################
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
             deltaZb(iCol,iLay) = hgtc(iLay) - hgtc(iLay+1)
          enddo
          deltaZb(iCol,nLev) = hgtc(nLev) - hgtb(nLev+1)
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
             deltaZb(iCol,iLay) = hgtc(iLay) - hgtc(iLay-1)
          enddo
          deltaZb(iCol,1) = hgtc(1) - hgtb(1)
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
       call get_alpha_exp(nCol, nLev, deltaZb, de_lgth, cloud_overlap_param)
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
    
  end subroutine GFS_rrtmgp_gfdlmp_pre_run
  
  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_gfdlmp_pre_finalize()
  end subroutine GFS_rrtmgp_gfdlmp_pre_finalize
end module GFS_rrtmgp_gfdlmp_pre
