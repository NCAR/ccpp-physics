! ########################################################################################
! This module contains the interface between the Zhao-Carr macrophysics and the RRTMGP 
! radiation schemes. Only compatable with imp_physics = imp_physics_zhaocarr
! ########################################################################################
module GFS_rrtmgp_zhaocarr_pre
  use machine,      only: kind_phys
  use rrtmgp_aux,   only: check_error_msg
  use funcphys,     only: fpvs
  use module_radiation_clouds, only: get_alpha_dcorr

  ! Zhao-Carr MP parameters.
  real(kind_phys), parameter :: &
       reliq_def = 10.0 ,       & ! Default liq radius to 10 micron
       reice_def = 50.0,        & ! Default ice radius to 50 micron
       rerain_def = 1000.0,     & ! Default rain radius to 1000 micron
       resnow_def = 250.0         ! Default snow radius to 250 micron
    
   public GFS_rrtmgp_zhaocarr_pre_init, GFS_rrtmgp_zhaocarr_pre_run, GFS_rrtmgp_zhaocarr_pre_finalize

contains  
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_zhaocarr_pre_init()
  end subroutine GFS_rrtmgp_zhaocarr_pre_init

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_GFS_rrtmgp_zhaocarr_pre_run
!! \htmlinclude GFS_rrtmgp_zhaocarr_pre_run.html
!!  
  subroutine GFS_rrtmgp_zhaocarr_pre_run(nCol, nLev, nCnd, nTracers, i_cldliq, lsswr,    &
       lslwr, effr_in, uni_cld, lmfshal, lat, lsmask, p_lev, p_lay, t_lay, relhum,       &
       tv_lay, effrin_cldliq, effrin_cldice, effrin_cldrain, effrin_cldsnow,             &
       shoc_sgs_cldfrac, cncvw, tracer,                                                  &
       con_eps, con_epsq, con_epsqs, con_epsm1, con_g, con_ttp, con_rd, con_pi,          &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,   &
       cld_rerain, de_lgth, deltaZ, cloud_overlap_param, errmsg, errflg)
    implicit none

    ! Inputs   
    integer, intent(in)    :: &
         nCol,              & ! Number of horizontal grid points
         nLev,              & ! Number of vertical layers
         nCnd,              & ! Number of cloud condensation types.
         nTracers,          & ! Number of tracers from model. 
         i_cldliq             ! Index into tracer array for cloud liquid. 
    logical, intent(in) :: &
    	 lsswr,             & ! Call SW radiation?
    	 lslwr,             & ! Call LW radiation
    	 effr_in,           & ! Provide hydrometeor radii from macrophysics?         
         uni_cld,           & !
         lmfshal            
    real(kind_phys), intent(in) :: &
         con_eps,           & ! rd/rv
         con_epsm1,         & ! (rd/rv) - 1
         con_epsq,          & ! Floor value for specific humidity
         con_epsqs,         & ! Floor value for saturation mixing ratio
         con_g,             & ! Gravitational acceleration (m/s2)
         con_ttp,           & ! Triple point temperature of water (K)
         con_rd,            & ! Ideal gas constant for dry air (J/kg/K)
         con_pi               ! Pi
    real(kind_phys), dimension(nCol), intent(in) :: &
         lsmask,            & ! Land/Sea mask
         lat                  ! Latitude             
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &         
         tv_lay,            & ! Virtual temperature (K)
         p_lay,             & ! Pressure at model-layers (Pa)
         t_lay,             & ! Temperature at model-layers (K)
         relhum,            & ! Relative humidity at model-layers ()
         effrin_cldliq,     & ! Effective radius for liquid cloud-particles (microns)
         effrin_cldice,     & ! Effective radius for ice cloud-particles (microns)
         effrin_cldrain,    & ! Effective radius for rain cloud-particles (microns)
         effrin_cldsnow,    & ! Effective radius for snow cloud-particles (microns)         
         shoc_sgs_cldfrac,  & ! Subgrid-scale cloud fraction from the SHOC scheme
         cncvw                ! Convective cloud water mixing ratio (kg/kg)
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &         
         p_lev                ! Pressure at model-level interfaces (Pa)
    real(kind_phys), dimension(nCol, nLev, nTracers),intent(in) :: &
         tracer               ! Cloud condensate amount in layer by type ()         
    
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
         deltaZ,               & ! Layer thickness (km)
         cloud_overlap_param     ! Cloud-overlap parameter
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag
    
    ! Local variables
    real(kind_phys) :: tem1,tem2,tem3,clwt,onemrh,clwm,clwmin,es,qs,value
    real(kind_phys), dimension(nCol, nLev, min(4,nCnd)) :: cld_condensate
    integer :: iCol,iLay
    real(kind_phys), dimension(nCol,nLev) :: deltaP
    
    if (.not. (lsswr .or. lslwr)) return
    
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
    
    ! ####################################################################################
    ! Pull out cloud information for Zhao-Carr MP scheme.
    ! ####################################################################################
    ! Condensate
    cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,i_cldliq)     ! Liquid water

    ! Set really tiny suspended particle amounts to clear
    do iLay=1,nLev
       do iCol=1,nCol   
          if (cld_condensate(iCol,iLay,1) < con_epsq) cld_condensate(iCol,iLay,1) = 0.0
       enddo
    enddo
    
    ! Use radii provided from the macrophysics?  
    if (effr_in) then
       cld_reliq(1:nCol,1:nLev)  = effrin_cldliq(1:nCol,1:nLev)
       cld_reice(1:nCol,1:nLev)  = effrin_cldice(1:nCol,1:nLev)
       cld_rerain(1:nCol,1:nLev) = effrin_cldrain(1:nCol,1:nLev)
       cld_resnow(1:nCol,1:nLev) = effrin_cldsnow(1:nCol,1:nLev)
    endif
    
    ! Use cloud-fraction from SHOC?
    if (uni_cld) then
       cld_frac(1:nCol,1:nLev) = shoc_sgs_cldfrac(1:nCol,1:nLev)
    ! Compute cloud-fraction?
    else
       clwmin = 0.0e-6
       if (.not. lmfshal) then
          do iLay = 1,nLev
             do iCol = 1, nCol
                es = min( p_lay(iCol,iLay),  fpvs( t_lay(iCol,iLay) ) )  ! fpvs and prsl in pa
                qs = max( con_epsqs, con_eps * es / (p_lay(iCol,iLay) + con_epsm1*es) )
                clwt = 1.0e-6 * (p_lay(iCol,iLay)*0.00001)
                if (cld_condensate(iCol,iLay,1) > clwt) then
                   onemrh= max( 1.e-10, 1.0-relhum(iCol,iLay) )
                   clwm  = clwmin / max( 0.01, p_lay(iCol,iLay)*0.00001 )
                   tem1  = min(max(sqrt(sqrt(onemrh*qs)),0.0001),1.0)
                   tem1  = 2000.0 / tem1
                   value = max( min( tem1*(cld_condensate(iCol,iLay,1)-clwm), 50.0 ), 0.0 )
                   tem2  = sqrt( sqrt(relhum(iCol,iLay)) )
                   cld_frac(iCol,iLay) = max( tem2*(1.0-exp(-value)), 0.0 )
                endif
             enddo
          enddo
       else
          do iLay=1,nLev
             do iCol = 1, nCol
                es = min( p_lay(iCol,iLay),  fpvs( t_lay(iCol,iLay) ) )  ! fpvs and prsl in pa
                qs = max( con_epsqs, con_eps * es / (p_lay(iCol,iLay) + con_epsm1*es) )
                clwt = 1.0e-6 * (p_lay(iCol,iLay)*0.00001)
                if (cld_condensate(iCol,iLay,1) > clwt) then
                   onemrh= max( 1.e-10, 1.0-relhum(iCol,iLay) )
                   clwm  = clwmin / max( 0.01, p_lay(iCol,iLay)*0.00001 )
                   tem1  = min(max((onemrh*qs)**0.49,0.0001),1.0)  !jhan
                   tem1  = 100.0 / tem1
                   value = max( min( tem1*(cld_condensate(iCol,iLay,1)-clwm), 50.0 ), 0.0 )
                   tem2  = sqrt( sqrt(relhum(iCol,iLay)) )                   
                   cld_frac(iCol,iLay) = max( tem2*(1.0-exp(-value)), 0.0 )
                endif
             enddo
          enddo
       endif
    endif
    
    ! Add suspended convective cloud water to grid-scale cloud water only for cloud 
    ! fraction & radiation computation it is to enhance cloudiness due to suspended convec 
    ! cloud water for zhao/moorthi's (imp_phys=99)    
    cld_condensate(1:nCol,1:nLev,1) = cld_condensate(1:nCol,1:nLev,1) + cncvw(1:nCol,1:nLev)
    
    ! Compute cloud liquid/ice condensate path.
    deltaP = abs(p_lev(:,2:nLev+1)-p_lev(:,1:nLev))/100.
    do iLay=1,nLev
       do iCol=1,nCol
          tem1               = max(0.0, cld_condensate(iCol,iLay,1)) * (1.0e5/con_g) * deltaP(iCol,iLay)
          cld_iwp(iCol,iLay) = tem1*(t_lay(iCol,iLay) - 273.16)
          cld_lwp(iCol,iLay) = tem1 - cld_iwp(iCol,iLay)
       enddo
    enddo

    ! Compute effective liquid cloud droplet radius over land.
    if(.not. effr_in) then
        do iCol = 1, nCol
           if (nint(lsmask(iCol)) == 1) then
              do iLay = 1, nLev
                 cld_reliq(iCol,iLay) = 5.0 + 5.0 * (t_lay(iCol,iLay) - 273.16)
              enddo
           endif
        enddo

       ! Compute effective ice cloud droplet radius following Heymsfield 
       ! and McFarquhar (1996) \cite heymsfield_and_mcfarquhar_1996.
       do iLay=1,nLev
          do iCol=1,nCol
             tem2 = t_lay(iCol,iLay) - con_ttp
             if (cld_iwp(iCol,iLay) > 0.0) then
                tem3 = (con_g/con_rd ) * cld_iwp(iCol,iLay) * (0.01*p_lay(iCol,iLay)) / (deltaP(iCol,iLay)*tv_lay(iCol,iLay))
                if (tem2 < -50.0) then
                   cld_reice(iCol,iLay) = (1250.0/9.917) * tem3 ** 0.109
                elseif (tem2 < -40.0) then
                   cld_reice(iCol,iLay) = (1250.0/9.337) * tem3 ** 0.08
                elseif (tem2 < -30.0) then
                   cld_reice(iCol,iLay) = (1250.0/9.208) * tem3 ** 0.055
                else
                   cld_reice(iCol,iLay) = (1250.0/9.387) * tem3 ** 0.031
                endif
                cld_reice(iCol,iLay)   = max(10.0, min(cld_reice(iCol,iLay), 150.0))
             endif
          enddo
       enddo
    endif

    ! #################################################################################### 
    ! Cloud (and precipitation) overlap                                                                                                                                                                    ! ####################################################################################
    ! Compute layer-thickness    
    do iCol=1,nCol
       do iLay=1,nLev
          deltaZ(iCol,iLay) = ((con_rd/con_g)*0.001) * abs(log(p_lev(iCol,iLay)) - log(p_lev(iCol,iLay+1))) * tv_lay(iCol,iLay)
       enddo
    enddo

    ! Cloud overlap parameter
    call get_alpha_dcorr(nCol, nLev, lat, con_pi, deltaZ, de_lgth, cloud_overlap_param)
        
  end subroutine GFS_rrtmgp_zhaocarr_pre_run
  
  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_zhaocarr_pre_finalize()
  end subroutine GFS_rrtmgp_zhaocarr_pre_finalize
 
end module GFS_rrtmgp_zhaocarr_pre
