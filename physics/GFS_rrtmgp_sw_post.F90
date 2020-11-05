module GFS_rrtmgp_sw_post 
  use machine,                   only: kind_phys
  use module_radiation_aerosols, only: NSPC1
  use module_radsw_parameters,   only: topfsw_type, sfcfsw_type, profsw_type, cmpfsw_type
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_fluxes_byband,          only: ty_fluxes_byband
  use mo_heating_rates,          only: compute_heating_rate
  use rrtmgp_aux,                only: check_error_msg
  implicit none
  
  public GFS_rrtmgp_sw_post_init,GFS_rrtmgp_sw_post_run,GFS_rrtmgp_sw_post_finalize

contains

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_post_init
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_post_init()
  end subroutine GFS_rrtmgp_sw_post_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_post_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_sw_post_run
!! \htmlinclude GFS_rrtmgp_sw_post_run.html
!!
  subroutine GFS_rrtmgp_sw_post_run (nCol, nLev, nDay, idxday, lsswr, do_sw_clrsky_hr,      &
       save_diag, fhswr,  coszen, coszdg, t_lay, p_lev, sfc_alb_nir_dir, sfc_alb_nir_dif,   &
       sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, sw_gas_props, fluxswUP_allsky,                 &
       fluxswDOWN_allsky, fluxswUP_clrsky, fluxswDOWN_clrsky, raddt, aerodp, cldsa, mbota,  &
       mtopa, cld_frac, cldtausw, fluxr,                                                    &
       nirbmdi, nirdfdi, visbmdi, visdfdi, nirbmui, nirdfui, visbmui, visdfui, sfcnsw,      &
       sfcdsw, htrsw, sfcfsw, topfsw, htrswc, flxprf_sw, scmpsw, errmsg, errflg)

    ! Inputs      
    integer, intent(in) :: &
         nCol,              & ! Horizontal loop extent 
         nLev,              & ! Number of vertical layers
         nDay                 ! Number of daylit columns
    integer, intent(in), dimension(nday) :: &
         idxday               ! Index array for daytime points
    logical, intent(in) :: &
    	 lsswr,             & ! Call SW radiation?
    	 do_sw_clrsky_hr,   & ! Output clear-sky SW heating-rate?         
    	 save_diag            ! Output radiation diagnostics?
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props         ! DDT containing SW spectral information
    real(kind_phys), intent(in) :: &
         fhswr                ! Frequency for SW radiation
    real(kind_phys), dimension(nCol), intent(in) :: &
         t_lay,             & ! Temperature at model layer centers (K)
         coszen,            & ! Cosine(SZA)     
         coszdg               ! Cosine(SZA), daytime     
    real(kind_phys), dimension(nCol, nLev+1), intent(in) :: &
         p_lev                ! Pressure @ model layer-interfaces    (Pa)
    real(kind_phys), dimension(sw_gas_props%get_nband(),ncol), intent(in) :: &
         sfc_alb_nir_dir,   & ! Surface albedo (direct) 
         sfc_alb_nir_dif,   & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir, & ! Surface albedo (direct)
         sfc_alb_uvvis_dif    ! Surface albedo (diffuse)
    real(kind_phys), dimension(nCol, nLev+1), intent(in) :: &
         fluxswUP_allsky,   & ! SW All-sky flux                    (W/m2)
         fluxswDOWN_allsky, & ! SW All-sky flux                    (W/m2)
         fluxswUP_clrsky,   & ! SW Clear-sky flux                  (W/m2)
         fluxswDOWN_clrsky    ! SW All-sky flux                    (W/m2)
    real(kind_phys), intent(in) :: &
         raddt                ! Radiation time step
    real(kind_phys), dimension(nCol,NSPC1), intent(in) :: &
         aerodp               ! Vertical integrated optical depth for various aerosol species  
    real(kind_phys), dimension(nCol,5), intent(in) :: &
         cldsa                ! Fraction of clouds for low, middle, high, total and BL 
    integer,         dimension(nCol,3), intent(in) ::&
         mbota,             & ! vertical indices for low, middle and high cloud tops 
         mtopa                ! vertical indices for low, middle and high cloud bases
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &
         cld_frac,          & ! Total cloud fraction in each layer
         cldtausw             ! approx .55mu band layer cloud optical depth
    
    ! Inputs (optional)     
    type(cmpfsw_type), dimension(nCol), intent(in), optional :: &
         scmpsw           ! 2D surface fluxes, components:
                          ! uvbfc - total sky downward uv-b flux at  (W/m2)
                          ! uvbf0 - clear sky downward uv-b flux at  (W/m2)
                          ! nirbm - downward nir direct beam flux    (W/m2)
                          ! nirdf - downward nir diffused flux       (W/m2)
                          ! visbm - downward uv+vis direct beam flux (W/m2)
                          ! visdf - downward uv+vis diffused flux    (W/m2)           
    
    real(kind=kind_phys), dimension(:,:), intent(inout) :: fluxr
    
    ! Outputs (mandatory)
    real(kind_phys), dimension(nCol), intent(out) :: &
         nirbmdi,           & ! sfc nir beam sw downward flux    (W/m2)
         nirdfdi,           & ! sfc nir diff sw downward flux    (W/m2)
         visbmdi,           & ! sfc uv+vis beam sw downward flux (W/m2)
         visdfdi,           & ! sfc uv+vis diff sw downward flux (W/m2)
         nirbmui,           & ! sfc nir beam sw upward flux      (W/m2)
         nirdfui,           & ! sfc nir diff sw upward flux      (W/m2)
         visbmui,           & ! sfc uv+vis beam sw upward flux   (W/m2)
         visdfui,           & ! sfc uv+vis diff sw upward flux   (W/m2)    
         sfcnsw,            & ! total sky sfc netsw flx into ground
         sfcdsw               !
    real(kind_phys), dimension(nCol,nLev), intent(out) :: &
         htrsw                ! SW all-sky heating rate
    type(sfcfsw_type), dimension(nCol), intent(out) :: &
         sfcfsw               ! sw radiation fluxes at sfc
    type(topfsw_type), dimension(nCol), intent(out) :: &
         topfsw               ! sw_fluxes_top_atmosphere
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg

    ! Outputs (optional)
    type(profsw_type), dimension(nCol, nLev), intent(out), optional :: &
         flxprf_sw        ! 2D radiative fluxes, components:
                          ! upfxc - total sky upward flux            (W/m2)
                          ! dnfxc - total sky dnward flux            (W/m2)
                          ! upfx0 - clear sky upward flux            (W/m2)
                          ! dnfx0 - clear sky dnward flux            (W/m2)
    real(kind_phys),dimension(nCol, nLev),intent(out),optional :: &
         htrswc           ! Clear-sky heating rate (K/s)
	
    ! Local variables
    integer :: i, j, k, iSFC, iTOA, itop, ibtc
    real(kind_phys) :: tem0d, tem1, tem2
    real(kind_phys), dimension(nDay, nLev) :: thetaTendClrSky, thetaTendAllSky
    logical :: l_fluxessw2d, top_at_1, l_scmpsw

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. lsswr) return
    if (nDay .gt. 0) then

       ! Are any optional outputs requested?
       l_fluxessw2d = present(flxprf_sw)
       
       ! Are the components of the surface fluxes provided?
       l_scmpsw = present(scmpsw)

       ! #######################################################################################
       ! What is vertical ordering?
       ! #######################################################################################
       top_at_1 = (p_lev(1,1) .lt. p_lev(1, nLev))
       if (top_at_1) then 
          iSFC = nLev+1
          iTOA = 1
       else
          iSFC = 1
          iTOA = nLev+1
       endif
       
       ! #######################################################################################
       ! Compute SW heating-rates
       ! #######################################################################################
       ! Clear-sky heating-rate (optional)
       if (do_sw_clrsky_hr) then
          htrswc(:,:) = 0._kind_phys
          call check_error_msg('GFS_rrtmgp_post',compute_heating_rate( &
               fluxswUP_clrsky(idxday(1:nDay),:),   & ! IN  - Shortwave upward clear-sky flux profiles (W/m2)
               fluxswDOWN_clrsky(idxday(1:nDay),:), & ! IN  - Shortwave downward clear-sky flux profiles (W/m2)
               p_lev(idxday(1:nDay),:),             & ! IN  - Pressure at model-interface (Pa)
               thetaTendClrSky))                      ! OUT - Clear-sky heating-rate (K/sec)
          htrswc(idxday(1:nDay),:)=thetaTendClrSky !**NOTE** GP doesn't use radiation levels, it uses the model fields. Not sure if this is necessary
       endif

       ! All-sky heating-rate (mandatory)
       htrsw(:,:) = 0._kind_phys
       call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(    &
            fluxswUP_allsky(idxday(1:nDay),:),      & ! IN  - Shortwave upward all-sky flux profiles (W/m2)
            fluxswDOWN_allsky(idxday(1:nDay),:),    & ! IN  - Shortwave downward all-sky flux profiles (W/m2)
            p_lev(idxday(1:nDay),:),                & ! IN  - Pressure at model-interface (Pa)
            thetaTendAllSky))                         ! OUT - All-sky heating-rate (K/sec)
       htrsw(idxday(1:nDay),:) = thetaTendAllSky

       ! #######################################################################################
       ! Save SW outputs
       ! #######################################################################################
       ! Copy fluxes from RRTGMP types into model radiation types.
       ! Mandatory outputs
       topfsw(:)%upfxc = fluxswUP_allsky(:,iTOA)
       topfsw(:)%upfx0 = fluxswUP_clrsky(:,iTOA)
       topfsw(:)%dnfxc = fluxswDOWN_allsky(:,iTOA)
       sfcfsw(:)%upfxc = fluxswUP_allsky(:,iSFC)
       sfcfsw(:)%upfx0 = fluxswUP_clrsky(:,iSFC)
       sfcfsw(:)%dnfxc = fluxswDOWN_allsky(:,iSFC)
       sfcfsw(:)%dnfx0 = fluxswDOWN_clrsky(:,iSFC)

       ! Optional output
       if(l_fluxessw2D) then
          flxprf_sw(:,:)%upfxc = fluxswUP_allsky(:,:)
          flxprf_sw(:,:)%dnfxc = fluxswDOWN_allsky(:,:)
          flxprf_sw(:,:)%upfx0 = fluxswUP_clrsky(:,:)
          flxprf_sw(:,:)%dnfx0 = fluxswDOWN_clrsky(:,:)
       endif
       
       ! Surface down and up spectral component fluxes
       ! - Save two spectral bands' surface downward and upward fluxes for output.
	   if (l_scmpsw) then
           do i=1,nCol
              nirbmdi(i) = scmpsw(i)%nirbm
              nirdfdi(i) = scmpsw(i)%nirdf
              visbmdi(i) = scmpsw(i)%visbm
              visdfdi(i) = scmpsw(i)%visdf
              nirbmui(i) = scmpsw(i)%nirbm * sfc_alb_nir_dir(1,i)
              nirdfui(i) = scmpsw(i)%nirdf * sfc_alb_nir_dif(1,i)
              visbmui(i) = scmpsw(i)%visbm * sfc_alb_uvvis_dir(1,i)
              visdfui(i) = scmpsw(i)%visdf * sfc_alb_uvvis_dif(1,i)
           enddo
        else
          nirbmdi(:) = 0.0
          nirdfdi(:) = 0.0
          visbmdi(:) = 0.0
          visdfdi(:) = 0.0
          nirbmui(:) = 0.0
          nirdfui(:) = 0.0
          visbmui(:) = 0.0
          visdfui(:) = 0.0        
        endif	
    else                   ! if_nday_block
       ! #######################################################################################
       ! Dark everywhere
       ! #######################################################################################
       htrsw(:,:) = 0.0
       sfcfsw     = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
       topfsw     = topfsw_type( 0.0, 0.0, 0.0 )
       nirbmdi(:) = 0.0
       nirdfdi(:) = 0.0
       visbmdi(:) = 0.0
       visdfdi(:) = 0.0
       nirbmui(:) = 0.0
       nirdfui(:) = 0.0
       visbmui(:) = 0.0
       visdfui(:) = 0.0
       
       if (do_sw_clrsky_hr) then
          htrswc(:,:) = 0
       endif
    endif                  ! end_if_nday

    ! Radiation fluxes for other physics processes
    do i=1,nCol
       sfcnsw(i) = sfcfsw(i)%dnfxc - sfcfsw(i)%upfxc
       sfcdsw(i) = sfcfsw(i)%dnfxc
    enddo

    ! #######################################################################################
    ! Save SW diagnostics
    ! - For time averaged output quantities (including total-sky and clear-sky SW and LW 
    !   fluxes at TOA and surface; conventional 3-domain cloud amount, cloud top and base 
    !   pressure, and cloud top temperature; aerosols AOD, etc.), store computed results in
    !   corresponding slots of array fluxr with appropriate time weights.
    ! - Collect the fluxr data for wrtsfc
    ! #######################################################################################
    if (save_diag) then
       do i=1,nCol
          fluxr(i,34) = fluxr(i,34) + fhswr*aerodp(i,1)  ! total aod at 550nm
          fluxr(i,35) = fluxr(i,35) + fhswr*aerodp(i,2)  ! DU aod at 550nm
          fluxr(i,36) = fluxr(i,36) + fhswr*aerodp(i,3)  ! BC aod at 550nm
          fluxr(i,37) = fluxr(i,37) + fhswr*aerodp(i,4)  ! OC aod at 550nm
          fluxr(i,38) = fluxr(i,38) + fhswr*aerodp(i,5)  ! SU aod at 550nm
          fluxr(i,39) = fluxr(i,39) + fhswr*aerodp(i,6)  ! SS aod at 550nm
          if (coszen(i) > 0.) then
             ! SW all-sky fluxes
             tem0d = fhswr * coszdg(i) / coszen(i)
             fluxr(i,2 ) = fluxr(i,2)  + topfsw(i)%upfxc * tem0d  ! total sky top sw up
             fluxr(i,3 ) = fluxr(i,3)  + sfcfsw(i)%upfxc * tem0d  
             fluxr(i,4 ) = fluxr(i,4)  + sfcfsw(i)%dnfxc * tem0d  ! total sky sfc sw dn
             ! SW uv-b fluxes
             fluxr(i,21) = fluxr(i,21) + scmpsw(i)%uvbfc * tem0d          ! total sky uv-b sw dn
             fluxr(i,22) = fluxr(i,22) + scmpsw(i)%uvbf0 * tem0d          ! clear sky uv-b sw dn
             ! SW TOA incoming fluxes
             fluxr(i,23) = fluxr(i,23) + topfsw(i)%dnfxc * tem0d     ! top sw dn 
             ! SW SFC flux components
             fluxr(i,24) = fluxr(i,24) + visbmdi(i) * tem0d          ! uv/vis beam sw dn
             fluxr(i,25) = fluxr(i,25) + visdfdi(i) * tem0d          ! uv/vis diff sw dn
             fluxr(i,26) = fluxr(i,26) + nirbmdi(i) * tem0d          ! nir beam sw dn
             fluxr(i,27) = fluxr(i,27) + nirdfdi(i) * tem0d          ! nir diff sw dn
             ! SW clear-sky fluxes
             fluxr(i,29) = fluxr(i,29) + topfsw(i)%upfx0 * tem0d
             fluxr(i,31) = fluxr(i,31) + sfcfsw(i)%upfx0 * tem0d 
             fluxr(i,32) = fluxr(i,32) + sfcfsw(i)%dnfx0 * tem0d
          endif
       enddo

       ! Save total and boundary-layer clouds
       do i=1,nCol
          fluxr(i,17) = fluxr(i,17) + raddt * cldsa(i,4)
          fluxr(i,18) = fluxr(i,18) + raddt * cldsa(i,5)
       enddo

       ! Save cld frac,toplyr,botlyr and top temp, note that the order of h,m,l cloud 
       ! is reversed for the fluxr output. save interface pressure (pa) of top/bot
       do j = 1, 3
          do i = 1, nCol
             tem0d = raddt * cldsa(i,j)
             itop  = mtopa(i,j)
             ibtc  = mbota(i,j)
             fluxr(i, 8-j) = fluxr(i, 8-j) + tem0d
             fluxr(i,11-j) = fluxr(i,11-j) + tem0d * p_lev(i,itop)
             fluxr(i,14-j) = fluxr(i,14-j) + tem0d * p_lev(i,ibtc)
             fluxr(i,17-j) = fluxr(i,17-j) + tem0d * p_lev(i,itop)

             ! Add optical depth and emissivity output
             tem1 = 0.
             do k=ibtc,itop
                tem1 = tem1 + cldtausw(i,k)      ! approx .55 mu channel
             enddo
             fluxr(i,43-j) = fluxr(i,43-j) + tem0d * tem1
          enddo
       enddo
    endif
  end subroutine GFS_rrtmgp_sw_post_run

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_post_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_post_finalize ()
  end subroutine GFS_rrtmgp_sw_post_finalize

end module GFS_rrtmgp_sw_post
