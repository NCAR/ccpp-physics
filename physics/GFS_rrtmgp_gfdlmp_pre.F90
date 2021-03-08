! ########################################################################################
! This module contains the interface between the GFDL macrophysics and the RRTMGP radiation
! schemes. Only compatable with Model%imp_physics = Model%imp_physics_gfdl
! ########################################################################################
module GFS_rrtmgp_gfdlmp_pre
  use machine,      only: kind_phys
  use rrtmgp_aux,   only: check_error_msg
  use module_radiation_cloud_overlap, only: cmp_dcorr_lgth, get_alpha_exp  
  use rrtmgp_lw_cloud_optics, only: radliq_lwr => radliq_lwrLW, radliq_upr => radliq_uprLW,&
                                    radice_lwr => radice_lwrLW, radice_upr => radice_uprLW

  ! Parameters
  real(kind_phys), parameter :: &
       reliq_def  = 10.0 ,      & ! Default liq radius to 10 micron (used when effr_in=F)
       reice_def  = 50.0,       & ! Default ice radius to 50 micron (used when effr_in=F)
       rerain_def = 1000.0,     & ! Default rain radius to 1000 micron (used when effr_in=F)
       resnow_def = 250.0,      & ! Default snow radius to 250 micron (used when effr_in=F)
       reice_min  = 10.0,       & ! Minimum ice size allowed by GFDL MP scheme
       reice_max  = 150.0         ! Maximum ice size allowed by GFDL MP scheme
  ! NOTE: When using RRTMGP cloud-optics, the min/max particle size allowed are imported 
  !       from initialization.

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
       i_cldrain, i_cldsnow, i_cldgrpl, i_cldtot, doSWrad, doLWrad, effr_in, kdt,        &
       do_mynnedmf, p_lev, p_lay, tv_lay, effrin_cldliq, effrin_cldice, effrin_cldrain,  &
       effrin_cldsnow, tracer, con_g, con_rd, doGP_cldoptics_PADE, doGP_cldoptics_LUT,   &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,   &
       cld_rerain, precip_frac, errmsg, errflg)
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
         kdt                     ! Current forecast iteration
    logical, intent(in) :: &
    	 doSWrad,              & ! Call SW radiation?
    	 doLWrad,              & ! Call LW radiation
    	 effr_in,              & ! Provide hydrometeor radii from macrophysics?
         do_mynnedmf,          & ! Flag to activate MYNN-EDMF 
         doGP_cldoptics_LUT,   & ! Flag to do GP cloud-optics (LUTs)
         doGP_cldoptics_PADE     !                            (PADE approximation)
    real(kind_phys), intent(in) :: &
         con_g,                & ! Physical constant: gravitational constant
         con_rd                  ! Physical constant: gas-constant for dry air
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
    real(kind_phys), dimension(nCol,nLev),intent(inout) :: &
         cld_frac,             & ! Total cloud fraction
         cld_lwp,              & ! Cloud liquid water path
         cld_reliq,            & ! Cloud liquid effective radius
         cld_iwp,              & ! Cloud ice water path
         cld_reice,            & ! Cloud ice effecive radius
         cld_swp,              & ! Cloud snow water path
         cld_resnow,           & ! Cloud snow effective radius
         cld_rwp,              & ! Cloud rain water path
         cld_rerain,           & ! Cloud rain effective radius       
         precip_frac             ! Precipitation fraction         
    character(len=*), intent(out) :: &
         errmsg                  ! Error message
    integer, intent(out) :: &  
         errflg                  ! Error flag
    
    ! Local variables
    real(kind_phys) :: tem1,pfac
    real(kind_phys), dimension(nCol, nLev, min(4,ncnd)) :: cld_condensate
    integer :: iCol,iLay,l,ncndl
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

    ! Initialize outputs
    cld_reliq(:,:)  = reliq_def
    cld_reice(:,:)  = reice_def
    cld_rerain(:,:) = rerain_def
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
                           
    ! Cloud water path (g/m2)
    deltaP = abs(p_lev(:,2:nLev+1)-p_lev(:,1:nLev))/100.  
    do iLay = 1, nLev
       do iCol = 1, nCol
          ! Compute liquid/ice condensate path from mixing ratios (kg/kg)->(g/m2)   
          tem1                = (1.0e5/con_g) * deltaP(iCol,iLay)
          cld_lwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,1) * tem1)
          cld_iwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,2) * tem1)
          cld_rwp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,3) * tem1)
          cld_swp(iCol,iLay)  = max(0., cld_condensate(iCol,iLay,4) * tem1) 
       enddo
    enddo

    ! Particle size
    do iLay = 1, nLev
       do iCol = 1, nCol
          ! Use radii provided from the macrophysics        
          if (effr_in) then
             cld_reliq(iCol,iLay)  = effrin_cldliq(iCol,iLay)
             cld_reice(iCol,iLay)  = max(reice_min, min(reice_max,effrin_cldice(iCol,iLay)))
             cld_rerain(iCol,iLay) = effrin_cldrain(iCol,iLay)
             cld_resnow(iCol,iLay) = effrin_cldsnow(iCol,iLay)
          endif
       enddo
    enddo

    ! Bound effective radii for RRTMGP, LUT's for cloud-optics go from 
    !   2.5 - 21.5 microns for liquid clouds,
    !   10  - 180  microns for ice-clouds    
    if (doGP_cldoptics_PADE .or. doGP_cldoptics_LUT) then
       where(cld_reliq .lt. radliq_lwr) cld_reliq = radliq_lwr
       where(cld_reliq .gt. radliq_upr) cld_reliq = radliq_upr
       where(cld_reice .lt. radice_lwr) cld_reice = radice_lwr
       where(cld_reice .gt. radice_upr) cld_reice = radice_upr
    endif
    
    ! Cloud-fraction
    if (do_mynnedmf .and. kdt .gt. 1) then
       do iLay = 1, nLev
          do iCol = 1, nCol
             if (tracer(iCol,iLay,i_cldrain) > 1.0e-7 .OR. tracer(iCol,iLay,i_cldsnow)>1.0e-7) then
                cld_frac(iCol,iLay) = tracer(iCol,iLay,i_cldtot)
             endif
          enddo
       enddo
    else
       cld_frac(1:nCol,1:nLev) = tracer(1:nCol,1:nLev,i_cldtot)
    endif
    
    ! Precipitation fraction (Hack. For now use cloud-fraction)
    precip_frac(1:nCol,1:nLev) = cld_frac(1:nCol,1:nLev)

  end subroutine GFS_rrtmgp_gfdlmp_pre_run
  
  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_gfdlmp_pre_finalize()
  end subroutine GFS_rrtmgp_gfdlmp_pre_finalize
end module GFS_rrtmgp_gfdlmp_pre
