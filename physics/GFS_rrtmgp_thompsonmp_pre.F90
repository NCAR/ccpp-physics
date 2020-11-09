! ########################################################################################
! This module contains the interface between the THOMPSON macrophysics and the RRTMGP radiation
! schemes. Only compatable with Model%imp_physics = Model%imp_physics_thompson
! ########################################################################################
module GFS_rrtmgp_thompsonmp_pre
  use machine, only: &
  	  kind_phys
  use rrtmgp_aux, only: &
  	 check_error_msg
  use module_radiation_cloud_overlap, only: &
  	 cmp_dcorr_lgth, &
     get_alpha_exp  
  use module_mp_thompson, only: &
  	 calc_effectRad, &
  	 Nt_c
  use module_mp_thompson_make_number_concentrations, only: &
  	 make_IceNumber,      &
     make_DropletNumber, &
     make_RainNumber
  implicit none
  
  ! Parameters specific to THOMPSONMP scheme.
  real(kind_phys), parameter :: &
       reliq_def  = 10.0 ,      & ! Default liq radius to 10 micron (used when effr_in=F)
       reice_def  = 50.0,       & ! Default ice radius to 50 micron (used when effr_in=F)
       rerain_def = 1000.0,     & ! Default rain radius to 1000 micron (used when effr_in=F)
       resnow_def = 250.0,      & ! Default snow radius to 250 micron (used when effr_in=F)  
       cllimit    = 0.001         ! Lowest cloud fraction in GFDL MP scheme
         
   public GFS_rrtmgp_thompsonmp_pre_init, GFS_rrtmgp_thompsonmp_pre_run, GFS_rrtmgp_thompsonmp_pre_finalize

contains  
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_rrtmgp_thompsonmp_pre_init()
  end subroutine GFS_rrtmgp_thompsonmp_pre_init

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_GFS_rrtmgp_thompsonmp_pre_run
!! \htmlinclude GFS_rrtmgp_thompsonmp_pre_run.html
!!  
  subroutine GFS_rrtmgp_thompsonmp_pre_run(nCol, nLev, nTracers, ncnd, i_cldliq, i_cldice,&
       i_cldrain, i_cldsnow, i_cldgrpl, i_cldtot, i_cldliq_nc, i_cldice_nc, i_twa, yearlen, doSWrad, doLWrad, effr_in, julian,          &
       lat, p_lev, p_lay, tv_lay, t_lay, effrin_cldliq, effrin_cldice, effrin_cldsnow, tracer,   &
       qs_lay, q_lay, relhum, cld_frac_mg, con_pi, con_g, con_rd, con_epsq, iovr, iovr_dcorr, uni_cld, lmfshal, lmfdeep2, ltaerosol,      &
       iovr_exp,iovr_exprand, idcor, dcorr_con, idcor_con, idcor_hogan, idcor_oreopoulos, &
       do_mynnedmf, imfdeepcnv, imfdeepcnv_gf, &
       cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp,    &
       cld_rerain, precip_frac, cloud_overlap_param, precip_overlap_param, de_lgth,       &
       deltaZb, errmsg, errflg)
    implicit none
    
    ! Inputs   
    integer, intent(in)    :: &
         nCol,              & ! Number of horizontal grid points
         nLev,              & ! Number of vertical layers
         ncnd,              & ! Number of cloud condensation types.
         nTracers,          & ! Number of tracers from model. 
         i_cldliq,          & ! Index into tracer array for cloud liquid amount. 
         i_cldice,          & !                             cloud ice amount.
         i_cldrain,         & !                             cloud rain amount.
         i_cldsnow,         & !                             cloud snow amount.
         i_cldgrpl,         & !                             cloud groupel amount.
         i_cldtot,          & !                             cloud total amount.
         i_cldliq_nc,       & !                             cloud liquid number concentration.
         i_cldice_nc,       & !                             cloud ice number concentration.
         i_twa,             & !                             water friendly aerosol.
         yearlen,           & ! Length of current year (365/366) WTF?             
         iovr,              & ! Choice of cloud-overlap method
         iovr_dcorr,        & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,          & ! Flag for exponential cloud overlap method
         iovr_exprand,      & ! Flag for exponential-random cloud overlap method
         idcor,             & ! Choice of method for decorrelation length computation
         idcor_con,         & ! Flag for decorrelation-length. Use constant value
         idcor_hogan,       & ! Flag for decorrelation-length. (https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.647)
         idcor_oreopoulos,  & ! Flag for decorrelation-length. (10.5194/acp-12-9097-2012)
         imfdeepcnv,        & ! Choice of mass-flux deep convection scheme
         imfdeepcnv_gf        ! Flag for Grell-Freitas deep convection scheme
    logical, intent(in) :: &
    	 doSWrad,           & ! Call SW radiation?
    	 doLWrad,           & ! Call LW radiation
    	 effr_in,           & ! Use cloud effective radii provided by model?
         uni_cld,           & !
         lmfshal,           & !
         lmfdeep2,          & !
         ltaerosol,         & !
         do_mynnedmf          ! Flag to activate MYNN-EDMF
    real(kind_phys), intent(in) :: &
         julian,            & ! Julian day 
         con_pi,            & ! Physical constant: pi
         con_g,             & ! Physical constant: gravitational constant
         con_rd,            & ! Physical constant: gas-constant for dry air
         con_epsq,          & ! Physical constant(?): Minimum value for specific humidity
         dcorr_con            ! Decorrelation-length (used if idcor = 0, default is idcor = 1)
    real(kind_phys), dimension(nCol), intent(in) :: &
         lat                  ! Latitude (radians)      
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &         
         tv_lay,            & ! Virtual temperature (K)
         t_lay,             & ! Temperature (K)
         qs_lay,            & ! Saturation vapor pressure (Pa)
         q_lay,             & ! water-vapor mixing ratio (kg/kg)
         relhum,            & ! Relative humidity
         p_lay,             & ! Pressure at model-layers (Pa)
         cld_frac_mg          ! Cloud-fraction from MG scheme. WTF?????
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &         
         p_lev                ! Pressure at model-level interfaces (Pa)
    real(kind_phys), dimension(nCol, nLev, nTracers),intent(in) :: &
         tracer               ! Cloud condensate amount in layer by type ()         

	! In/Outs
    real(kind_phys), dimension(nCol,nLev), intent(inout) :: &  
         cld_frac,          & ! Total cloud fraction
         cld_lwp,           & ! Cloud liquid water path
         cld_reliq,         & ! Cloud liquid effective radius
         cld_iwp,           & ! Cloud ice water path
         cld_reice,         & ! Cloud ice effecive radius
         effrin_cldliq,     & ! Effective radius for liquid cloud-particles (microns)
         effrin_cldice,     & ! Effective radius for ice cloud-particles (microns)
         effrin_cldsnow       ! Effective radius for snow cloud-particles (microns)	         
    
    ! Outputs     
    real(kind_phys), dimension(nCol),intent(out) :: &
         de_lgth                 ! Decorrelation length     
    real(kind_phys), dimension(nCol,nLev),intent(out) :: &
         cld_swp,              & ! Cloud snow water path
         cld_resnow,           & ! Cloud snow effective radius
         cld_rwp,              & ! Cloud rain water path
         cld_rerain,           & ! Cloud rain effective radius       
         precip_frac,          & ! Precipitation fraction
         cloud_overlap_param,  & ! Cloud-overlap parameter
         precip_overlap_param, & ! Precipitation overlap parameter  
         deltaZb                 ! Layer thickness (km)          
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag
    
    ! Local variables
    real(kind_phys) :: tem0, tem1, tem2, pfac, clwt, clwm, onemrh, clwmin, clwf
    real(kind_phys), dimension(nLev+1) :: hgtb
    real(kind_phys), dimension(nLev)   :: hgtc
    real(kind_phys), dimension(nCol, nLev, min(4,ncnd)) :: cld_condensate
    integer :: iCol,iLay,l,iSFC,iTOA
    real(kind_phys), dimension(nCol,nLev) :: deltaP, deltaZ, rho, orho, re_cloud, re_ice,&
       re_snow, qv_mp, qc_mp, qi_mp, qs_mp, nc_mp, ni_mp, nwfa
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
    
    ! ####################################################################################
    ! Pull out cloud information for THOMPSON MP scheme.
    ! ####################################################################################

    ! Cloud condensate
    cld_condensate(1:nCol,1:nLev,1) = tracer(1:nCol,1:nLev,i_cldliq)     ! -liquid water
    cld_condensate(1:nCol,1:nLev,2) = tracer(1:nCol,1:nLev,i_cldice)     ! -ice water
    cld_condensate(1:nCol,1:nLev,3) = tracer(1:nCol,1:nLev,i_cldrain)    ! -rain water
    cld_condensate(1:nCol,1:nLev,4) = tracer(1:nCol,1:nLev,i_cldsnow) + &! -snow + grapuel
                                      tracer(1:nCol,1:nLev,i_cldgrpl) 
                                      
    ! Cloud particle size
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
        	
    ! First, prepare cloud mixing-ratios and number concentrations for Calc_Re
    rho  = p_lay(1:nCol,1:nLev)/(con_rd*t_lay(1:nCol,1:nLev))    
    orho = 1./rho    
    do iLay = 1, nLev
       do iCol = 1, nCol    
         qv_mp(iCol,iLay) = q_lay(iCol,iLay)/(1.-q_lay(iCol,iLay))
         qc_mp(iCol,iLay) = tracer(iCol,iLay,i_cldliq)       / (1.-q_lay(iCol,iLay))
         qi_mp(iCol,iLay) = tracer(iCol,iLay,i_cldice)       / (1.-q_lay(iCol,iLay))
         qs_mp(iCol,iLay) = tracer(iCol,iLay,i_cldsnow)      / (1.-q_lay(iCol,iLay))
         nc_mp(iCol,iLay) = tracer(iCol,iLay,i_cldliq_nc)    / (1.-q_lay(iCol,iLay))
         if (ltaerosol) then
            ni_mp(iCol,iLay) = tracer(iCol,iLay,i_cldice_nc) / (1.-q_lay(iCol,iLay))
            nwfa(iCol,iLay)  = tracer(iCol,iLay,i_twa)
         else
            nc_mp(iCol,iLay) = nt_c*orho(iCol,iLay)
            ni_mp(iCol,iLay) = tracer(iCol,iLay,i_cldice_nc) / (1.-q_lay(iCol,iLay)) 
         endif        
       enddo
    enddo

    ! Update number concentration, consistent with sub-grid clouds
    do iLay = 1, nLev
       do iCol = 1, nCol  
          if (ltaerosol .and. qc_mp(iCol,iLay) > 1.e-12 .and. nc_mp(iCol,iLay) < 100.) then
             nc_mp(iCol,iLay) = make_DropletNumber(qc_mp(iCol,iLay)*rho(iCol,iLay), nwfa(iCol,iLay)) * orho(iCol,iLay)
          endif
          if (qi_mp(iCol,iLay) > 1.e-12 .and. ni_mp(iCol,iLay) < 100.) then
             ni_mp(iCol,iLay) = make_IceNumber(qi_mp(iCol,iLay)*rho(iCol,iLay), t_lay(iCol,iLay)) * orho(iCol,iLay)
          endif
       enddo
    enddo

    ! Compute effective radii for liquid/ice/snow using subgrid scale clouds
    ! Call Thompson's subroutine to compute effective radii
    do iCol=1,nCol
       call calc_effectRad (t_lay(iCol,:), p_lay(iCol,:), qv_mp(iCol,:), qc_mp(iCol,:),  &
                           nc_mp(iCol,:), qi_mp(iCol,:), ni_mp(iCol,:), qs_mp(iCol,:),   &
                           re_cloud(iCol,:), re_ice(iCol,:), re_snow(iCol,:), 1, nLev )
    enddo
    
    ! Scale Thompson's effective radii from meter to micron and update global effective radii.
    effrin_cldliq(1:nCol,1:nLev)  = re_cloud(1:nCol,1:nLev)*1.e6
    effrin_cldice(1:nCol,1:nLev)  = re_ice(1:nCol,1:nLev)*1.e6
    effrin_cldsnow(1:nCol,1:nLev) = re_snow(1:nCol,1:nLev)*1.e6
    cld_reliq(1:nCol,1:nLev)      = effrin_cldliq(1:nCol,1:nLev)
    cld_reice(1:nCol,1:nLev)      = effrin_cldice(1:nCol,1:nLev)
    cld_resnow(1:nCol,1:nLev)     = effrin_cldsnow(1:nCol,1:nLev)
    cld_rerain(1:nCol,1:nLev)     = rerain_def

	! Compute cloud-fraction. The logic is a mess here. I don't have any idea where these
	! magic numbers are coming from.
    if(.not. do_mynnedmf .or. imfdeepcnv .ne. imfdeepcnv_gf ) then ! MYNN PBL or GF conv   	
	   ! Cloud-fraction
	   if (uni_cld) then
          cld_frac(1:nCol,1:nLev) = cld_frac_mg(1:nCol,1:nLev)    
       else
         clwmin = 0.0
         if (.not. lmfshal) then
            do iLay = 1, nLev
               do iCol = 1, nCol         
                  clwf = tracer(iCol,iLay,i_cldliq) + tracer(iCol,iLay,i_cldice) + &
                         tracer(iCol,iLay,i_cldsnow)
                  clwt = 1.0e-6 * (p_lay(iCol,iLay)*0.001)
                  if (clwf > clwt) then
                     onemrh= max( 1.e-10, 1.0-relhum(iCol,iLay) )
                     clwm  = clwmin / max( 0.01, p_lay(iCol,iLay)*0.001 )
                     tem1  = 2000.0 / min(max(sqrt(sqrt(onemrh*qs_lay(iCol,iLay))),0.0001),1.0)
                     tem1  = max( min( tem1*(clwf-clwm), 50.0 ), 0.0 )
                     tem2  = sqrt( sqrt(relhum(iCol,iLay)) )
                     !
                     cld_frac(iCol,iLay) = max( tem2*(1.0-exp(-tem1)), 0.0 )
                   endif
               enddo
          enddo
        else
           do iLay = 1, nLev
              do iCol = 1, nCol              
                 clwf = tracer(iCol,iLay,i_cldliq) + tracer(iCol,iLay,i_cldice) + &
                         tracer(iCol,iLay,i_cldsnow)              	 
                 clwt = 1.0e-6 * (p_lay(iCol,iLay)*0.001)

                 if (clwf > clwt) then
                    onemrh= max( 1.e-10, 1.0-relhum(iCol,iLay) )
                    clwm  = clwmin / max( 0.01, p_lay(iCol,iLay)*0.001 )
                    tem1  = 100.0 / min(max((onemrh*qs_lay(iCol,iLay))**0.49,0.0001),1.0)  !jhan
                    tem1  = max( min( tem1*(clwf-clwm), 50.0 ), 0.0 )
                    tem2  = sqrt( sqrt(relhum(iCol,iLay)) )
                    !
                    cld_frac(iCol,iLay) = max( tem2*(1.0-exp(-tem1)), 0.0 )
                  endif
               enddo
             enddo
          endif     
       endif  
    endif     
        
    ! Precipitation fraction (Hack. For now use cloud-fraction)
    precip_frac(1:nCol,1:nLev) = cld_frac(1:nCol,1:nLev)
        
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
       call cmp_dcorr_lgth(nCol, abs(lat/con_pi), con_pi, de_lgth)      
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
    call get_alpha_exp(nCol, nLev, deltaZb, de_lgth, cloud_overlap_param)
    
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
    
  end subroutine GFS_rrtmgp_thompsonmp_pre_run
  
  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_thompsonmp_pre_finalize()
  end subroutine GFS_rrtmgp_thompsonmp_pre_finalize
end module GFS_rrtmgp_thompsonmp_pre