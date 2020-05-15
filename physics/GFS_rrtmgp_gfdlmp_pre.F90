module GFS_rrtmgp_gfdlmp_pre
  use machine,      only: kind_phys
  use GFS_typedefs, only: GFS_control_type, GFS_tbd_type 
  use physcons,     only: con_ttp, & ! Temperature at h2o 3pt (K)
  						  con_rd,  & ! Gas constant for dry air (J/KgK)
  						  con_pi,  & ! PI
  						  con_g      ! Gravity (m/s2)
  use physparam,    only: lcnorm,lcrick
  ! Parameters
  real(kind_phys), parameter :: &
  	 reliq_def = 10.0,     & ! Default liq radius to 10 micron
     reice_def = 50.0,     & ! Default ice radius to 50 micron
     rrain_def = 1000.0,   & ! Default rain radius to 1000 micron
     rsnow_def = 250.0,    & ! Default snow radius to 250 micron  
     epsq      = 1.0e-12,  & ! Tiny value
     cllimit   = 0.001,    & ! Lowest cloud fraction in GFDL MP scheme
     gfac=1.0e5/con_g        !
       
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
  subroutine GFS_rrtmgp_gfdlmp_pre_run(Model, Tbd, nCol, nLev, slmsk, lat, p_lay, p_lev, &
      t_lay, tv_lay, tracer,                                                             &
      cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,    &
      cld_rerain, errmsg, errflg)
    implicit none
    
    ! Inputs  
    type(GFS_control_type), intent(in) :: &
         Model                ! DDT: FV3-GFS model control parameters      
    type(GFS_tbd_type), intent(in) :: &
         Tbd                  ! DDT: FV3-GFS data not yet assigned to a defined container         
    integer, intent(in) :: &
         nCol,              & ! Number of horizontal grid-points
         nLev                 ! Number of vertical-layers
    real(kind_phys), dimension(nCol), intent(in) :: &
         slmsk,             & ! Land/sea/sea-ice mask
         lat                  ! Latitude
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &
         p_lay,             & ! Pressure at model-layer (Pa)
         t_lay,             & ! Temperature at model layer  (K)
         tv_lay               ! Virtual temperature at model-layers (K)
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &         
         p_lev                ! Pressure at model-level interfaces (Pa)
    real(kind_phys), dimension(nCol, nLev, Model%ntrac),intent(in) :: &
         tracer               ! Cloud condensate amount in layer by type ()         
    
    ! Outputs     
    real(kind_phys), dimension(nCol,nLev),intent(out) :: &
         cld_frac,          & ! Total cloud fraction
         cld_lwp,           & ! Cloud liquid water path
         cld_reliq,         & ! Cloud liquid effective radius
         cld_iwp,           & ! Cloud ice water path
         cld_reice,         & ! Cloud ice effecive radius
         cld_swp,           & ! Cloud snow water path
         cld_resnow,        & ! Cloud snow effective radius
         cld_rwp,           & ! Cloud rain water path
         cld_rerain           ! Cloud rain effective radius             
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag
                      
    ! Local variables
    real(kind_phys) :: tem1, tem2, tem3, clwt
    real(kind_phys), dimension(nCol) :: rlat
    real(kind_phys), dimension(nCol, nLev, min(4,Model%ncnd)) :: cld_condensate, clwf
    integer :: i,k,l,ncndl,icnd
    real(kind_phys), dimension(nCol,nLev) :: deltaP, cldcov
    real(kind_phys), dimension(nCol,nLev,9) :: clouds         
    
    if (.not. (Model%lsswr .or. Model%lslwr)) return
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Initialize outputs
    cld_lwp(:,:)    = 0.0
    cld_reliq(:,:)  = 0.0
    cld_iwp(:,:)    = 0.0
    cld_reice(:,:)  = 0.0
    cld_rwp(:,:)    = 0.0
    cld_rerain(:,:) = 0.0
    cld_swp(:,:)    = 0.0
    cld_resnow(:,:) = 0.0

    ! Compute layer pressure thickness (hPa)
    deltaP = abs(p_lev(:,2:nLev+1)-p_lev(:,1:nLev))/100.  

	! ####################################################################################
    ! Pull out cloud information for GFDL MP scheme.
    ! ####################################################################################
    ! Cloud hydrometeors
    cld_condensate(:,:,:) = 0._kind_phys
    if (Model%ncnd .eq. 2) then
       cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,Model%ntcw)     ! -liquid water
       cld_condensate(1:nCol,1:nLev,2) = tracer(1:nCol,1:nLev,Model%ntiw)     ! -ice water
       ncndl = Model%ncnd
    endif   
    if (Model%ncnd .eq. 5) then
       cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,Model%ntcw)     ! -liquid water
       cld_condensate(1:nCol,1:nLev,2) = tracer(1:nCol,1:nLev,Model%ntiw)     ! -ice water
       cld_condensate(1:nCol,1:nLev,3) = tracer(1:nCol,1:nLev,Model%ntrw)     ! -rain water
       cld_condensate(1:nCol,1:nLev,4) = tracer(1:nCol,1:nLev,Model%ntsw) + & ! -snow + grapuel
                                         tracer(1:nCol,1:nLev,Model%ntgl) 
                           
       ! Since we combine the snow and grapuel, define local variable for number of condensate types.
       ncndl = min(4,Model%ncnd)       
    endif
    
	! Cloud-fraction
    cld_frac(1:nCol,1:nLev) = tracer(1:nCol,1:nLev,Model%ntclamt)
                           			    
	! Set really tiny suspended particle amounts to clear
    do l=1,ncndl
       do k=1,nLev
          do i=1,nCol   
             if (cld_condensate(i,k,l) < epsq) cld_condensate(i,k,l) = 0.0
	      enddo
       enddo
    enddo
            
    ! DJS asks. Do we need lcrick? If not replace clwf with cld_condensate(:,:,1)
    if ( lcrick ) then
       do icnd=1,ncndl
          do i = 1, nCol
             clwf(i,1,icnd)    = 0.75*cld_condensate(i,1,icnd)    + 0.25*cld_condensate(i,2,icnd)
             clwf(i,nlev,icnd) = 0.75*cld_condensate(i,nLev,icnd) + 0.25*cld_condensate(i,nLev-1,icnd)
          enddo
          do k = 2, nLev-1
             do i = 1, nCol
               clwf(i,k,icnd) = 0.25*cld_condensate(i,k-1,icnd) + 0.5*cld_condensate(i,k,icnd) + &
            	   		        0.25*cld_condensate(i,k+1,icnd)
             enddo
          enddo
        enddo
    else
       do icnd=1,ncndl
          do k = 1, nLev
             do i = 1, nCol
                clwf(i,k,icnd) = cld_condensate(i,k,icnd)
             enddo
           enddo
        enddo
    endif    

    ! ####################################################################################
	! A) Compute Liquid/Ice/Rain/Snow(+groupel) cloud condensate paths
    ! ####################################################################################
    
    ! ####################################################################################
    ! i) This option uses the mixing-ratios and effective radii for 5 cloud hydrometeor types,
    !    Liquid, Ice, Rain, and Snow(+groupel), to determine cloud properties. 
    !    Formerly progclduni()
    ! ####################################################################################
    if (Model%lgfdlmprad) then    	
      ! Compute liquid/ice condensate path from mixing ratios (kg/kg)->(g/m2)   
      do k = 1, nLev
         do i = 1, nCol
         	if (cld_frac(i,k) .ge. cllimit) then         
               tem1          = gfac * deltaP(i,k)
               cld_lwp(i,k)  = clwf(i,k,1) * tem1
               cld_iwp(i,k)  = clwf(i,k,2) * tem1
               ! Also Rain and Snow(+groupel) if provided
               if (ncndl .eq. 4) then
                  cld_rwp(i,k)  = clwf(i,k,3) * tem1
                  cld_swp(i,k)  = clwf(i,k,4) * tem1  
               endif                
            endif
          enddo
       enddo            
    ! ####################################################################################
    ! ii) This option uses only a single mixing-ratio and partitions into liquid/ice cloud
    !     properties by phase.
    !     Formerly progcld4()
    ! ####################################################################################
    else
	   ! Compute total-cloud suspended water.
       clwf(:,:,1) = sum(clwf,dim=3)

      ! Compute liquid/ice condensate path (g/m2)
      do k = 1, nLev
         do i = 1, nCol
         	if (cld_frac(i,k) .ge. cllimit) then
               clwt         = max(0.0,clwf(i,k,1)) * gfac * deltaP(i,k)
               tem2         = min( 1.0, max( 0.0, (con_ttp-t_lay(i,k))*0.05 ) )
               cld_iwp(i,k) = clwt * tem2
               cld_lwp(i,k) = clwt - cld_iwp(i,k)
            endif
         enddo
      enddo
   endif

    ! ####################################################################################
    ! B) Particle sizes
    ! ####################################################################################

    ! ####################################################################################
    ! i) Use radii provided from the macrophysics        
    ! ####################################################################################
    if (Model%effr_in) then
       do k=1,nLev
          do i=1,nCol
            cld_reliq(i,k)  = Tbd%phy_f3d(i,k,1)
            cld_reice(i,k)  = max(10.0, min(150.0,Tbd%phy_f3d(i,k,2)))
            cld_rerain(i,k) = Tbd%phy_f3d(i,k,3)
            cld_resnow(i,k) = Tbd%phy_f3d(i,k,4)
          enddo
        enddo
    ! ####################################################################################
    ! ii) Start with default values. Modify liquid sizes over land. Adjust ice sizes following
    !     Hemsfield and McFarquhar (1996) https://doi.org/10.1175/1520-0469
    ! ####################################################################################
    else
       cld_reliq(:,:)  = reliq_def
       cld_reice(:,:)  = reice_def
       cld_rerain(:,:) = rrain_def     
       cld_resnow(:,:) = rsnow_def     
        
	   ! Compute effective liquid cloud droplet radius over land.
       do i = 1, nCol
         if (nint(slmsk(i)) == 1) then
           do k = 1, nLev
             tem2           = min( 1.0, max( 0.0, (con_ttp-t_lay(i,k))*0.05 ) )
             cld_reliq(i,k) = 5.0 + 5.0 * tem2
           enddo
         endif
       enddo

       ! Compute effective ice cloud droplet radius.
       do k = 1, nLev
          do i = 1, nCol
             tem2 = t_lay(i,k) - con_ttp
             if (cld_iwp(i,k) > 0.0) then
               tem3 = (con_g/con_rd)* cld_iwp(i,k) * (p_lay(i,k)/100.) / (deltaP(i,k)*tv_lay(i,k))
               if (tem2 < -50.0) then
                 cld_reice(i,k) = (1250.0/9.917) * tem3 ** 0.109
               elseif (tem2 < -40.0) then
                 cld_reice(i,k) = (1250.0/9.337) * tem3 ** 0.08
               elseif (tem2 < -30.0) then
                 cld_reice(i,k) = (1250.0/9.208) * tem3 ** 0.055
               else
                 cld_reice(i,k) = (1250.0/9.387) * tem3 ** 0.031
               endif
               cld_reice(i,k)   = max(10.0, min(cld_reice(i,k), 150.0))
             endif
           enddo
       enddo        
    endif  
    
   ! Normalize cloud-condensate by cloud-cover?
   if ( lcnorm ) then
      do k = 1, nLev
         do i = 1, nCol
            if (cld_frac(i,k) >= cllimit) then
               tem1 = 1.0 / max(0.05, cld_frac(i,k))
               cld_lwp(i,k) = cld_lwp(i,k) * tem1
               cld_iwp(i,k) = cld_iwp(i,k) * tem1
               cld_rwp(i,k) = cld_rwp(i,k) * tem1
               cld_swp(i,k) = cld_swp(i,k) * tem1
            endif
         enddo
      enddo
    endif  
    
  end subroutine GFS_rrtmgp_gfdlmp_pre_run

  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_gfdlmp_pre_finalize()
  end subroutine GFS_rrtmgp_gfdlmp_pre_finalize

end module GFS_rrtmgp_gfdlmp_pre