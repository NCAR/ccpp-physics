! ############################################################################################# 
!> \file GFS_radiation_post.F90
!!
!! Radiation post-processing routine.
!!
!! This module has two purposes:
!! 1*) Perform coupling from the radiation scheme(s) to other physical parameterizations.
!! 2) Compute diagnostics
!!
!! *For RRTMG,  this coupling is handled in the SCHEME.
!! *For RRTMGP, this coupling is handled HERE (more on this below).
!!
! ############################################################################################# 
module GFS_radiation_post
  use machine,                   only: kind_phys
  use module_radlw_parameters,   only: topflw_type, sfcflw_type
  use module_radsw_parameters,   only: topfsw_type, sfcfsw_type, cmpfsw_type
  use mo_heating_rates,          only: compute_heating_rate
  use radiation_tools,           only: check_error_msg
  implicit none

  public GFS_radiation_post_run

contains
! #############################################################################################
!> \section arg_table_GFS_radiation_post_run Argument Table
!! \htmlinclude GFS_radiation_post_run.html
!!
!! This routine needs to be called AFTER the RRTMG (radlw_main.F90 and radsw_main.F90)
!! or the RRTMGP (rrtmgp_lw_main.F90 and rrtmgp_sw_main.F90) radiaiton schemes in the
!! CCPP enabled UFS.
!!
!! For RRTMG, not much is done here, since the scheme outputs the fields needed by the
!! UFS. For example, RRTMG provides the heating-rate profiles and has been modified to use 
!! UFS native DDTs for storing the fluxes.
!!
!! For RRTMGP*:
!! - The all-sky radiation tendency is computed. The clear-sky tendency is computed, if
!!   requested.
!! - Surface and TOA fluxes are copied to UFS native DDTs that persist between radiation/physics
!!   calls.
!!
!! *Note on RTE-RRTMGP implementation in CCPP
!!  This is done in an attempt to make the CCPP enabled RRTMGP LW/SW drivers more host agnostic.
!!  The drivers are outputting the same fields as RTE, flux profiles, maintaining the same scheme
!!  interface at the lowest CCPP entrypoint, the "Scheme-level" interstitial. Any host specific
!!  coupling to the scheme happens here, a layer above, within the "Suite-level" interstitial.
!!
!! For ALL Radiaiton Schemes:
!! - Compute SW total cloud albedo
!! - Compute diagnostics
!!
! ############################################################################################# 

  ! ###########################################################################################
  ! GFS_radiation_post_run
  ! ###########################################################################################
  subroutine GFS_radiation_post_run(doLWrad, doSWrad, lssav, total_albedo, topfsw, fhlwr, fhswr,&
      coszen, coszdg, raddt, aerodp, cldsa, mtopa, mbota, cldtausw, cldtaulw, p_lev, tgrs, kb,  &
      kd, kt, sfcflw, sfcfsw, topflw, scmpsw, nCol, nLev, lmk, nDay, nfxr, nspc1, fluxr,        &
      do_RRTMGP, do_lw_clrsky_hr, fluxlwUP_clrsky, fluxlwDOWN_clrsky, htrlwc, fluxlwUP_allsky,  &
      fluxlwDOWN_allsky, htrlw, do_sw_clrsky_hr, htrswc, fluxswUP_clrsky, idxday,               &
      fluxswDOWN_clrsky, htrsw, fluxswUP_allsky, fluxswDOWN_allsky, iSFC, iTOA, tsflw, tsfa,    &
      sfcdlw, sfculw, htrlwu, nirbmdi, nirdfdi, visbmdi, visdfdi, nirbmui, nirdfui, visbmui,    &
      visdfui, sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, sfcnsw,  &
      sfcdsw, errmsg, errflg) 

    ! Inputs
    integer, intent(in) :: &
         nCol,              & !< Horizontal loop extent 
         nLev,              & !< Number of vertical layers
         lmk,               & !< Number of vertical layers for radiation (adjusted)
         nDay,              & !< Number of daylit columns
         nfxr,              & !< Number of variables stored in the fluxr array
         nspc1,             & !< Number of species for output aerosol optical depth
         kb,                & !< Vertical index difference between layer and lower bound (H/M/L diag)
         kd,                & !< Vertical index difference between in/out and local  (H/M/L diag)
         kt,                & !< Vertical index difference between layer and upper bound (H/M/L diag)
         iSFC,              & !< Vertical index for surface level
         iTOA                 !< Vertical index for TOA level
    integer, intent(in), dimension(:) :: &
         idxday               !< Index array for daytime points
    logical, intent(in) :: & 
         doLWrad,           & !< Logical flags for lw radiation calls
         doSWrad,           & !< Logical flags for sw radiation calls
         lssav,             & !< Flag for radiation diagnostics
         do_RRTMGP,         & !< Flag for using RRTMGP scheme
         do_lw_clrsky_hr,   & !< Output clear-sky LW heating-rate?
         do_sw_clrsky_hr      !< Output clear-sky SW heating-rate? 
    real(kind_phys), intent(in) ::  &
         fhlwr,             & !< Frequency for longwave radiation  (sec)
         fhswr,             & !< Frequency for shortwave radiation (sec)
         raddt                !< Radiation time step               (sec)
    real(kind_phys), dimension(:), intent(in) ::  &
         coszen,            & !< Mean cos of zenith angle over rad call period
         coszdg               !< Daytime mean cosz over rad call period
    real(kind_phys), dimension(:), intent(in) ::  &
         tsfa,              & !< Lowest model layer air temperature for radiation (K)
         sfc_alb_nir_dir,   & !< Surface albedo (direct) 
         sfc_alb_nir_dif,   & !< Surface albedo (diffuse)
         sfc_alb_uvvis_dir, & !< Surface albedo (direct)
         sfc_alb_uvvis_dif    !< Surface albedo (diffuse)
    real(kind_phys), dimension(:,:), intent(in) :: &
         p_lev                !< Pressure @ model layer-interfaces (Pa)
    real(kind_phys), dimension(:,:), intent(in) :: & 
         tgrs                 !< Temperature @ model layer-centers (K)
    real(kind_phys), dimension(:,:), intent(in), optional :: &
         fluxlwUP_clrsky,   & !< RRTMGP longwave clear-sky flux    (W/m2)
         fluxlwDOWN_clrsky, & !< RRTMGP longwave clear-sky flux    (W/m2)
         fluxlwUP_allsky,   & !< RRTMGP longwave all-sky flux      (W/m2)
         fluxlwDOWN_allsky, & !< RRTMGP longwave all-sky flux      (W/m2)
         fluxswUP_clrsky,   & !< RRTMGP shortwave clear-sky flux   (W/m2)
         fluxswDOWN_clrsky, & !< RRTMGP shortwave clear-sky flux   (W/m2)
         fluxswUP_allsky,   & !< RRTMGP shortwave all-sky flux     (W/m2)
         fluxswDOWN_allsky    !< RRTMGP shortwave all-sky flux     (W/m2)
    real(kind_phys), dimension(:,:), intent(in) ::  &
         aerodp               !< Vertical integrated optical depth for <nspc1> aerosol species
    real(kind_phys), dimension(:,:), intent(in) ::  & 
         cldtausw,          & !< .55mu band layer cloud optical depth (SW)
         cldtaulw             !< 10mu  band layer cloud optical depth (LW)
    real(kind_phys), dimension(:,:), intent(in) ::  &
         cldsa                !< Fraction of clouds for High/Mid/Low diagnostics:
                              !<   low(1), middle(2), high(3), total(4) and BL(5)
    integer, intent(in), dimension(:,:) :: &
         mtopa,             & !< Vertical indices for low, middle and high cloud tops  (H/M/L diag)
         mbota                !< Vertical indices for low, middle and high cloud bases (H/M/L diag)
    type(cmpfsw_type), dimension(:), intent(in) :: &
         scmpsw               !< 2D surface fluxes, components:
                              !!\n uvbfc - total sky downward uv-b flux at  (W/m2)
                              !!\n uvbf0 - clear sky downward uv-b flux at  (W/m2)
                              !!\n nirbm - downward nir direct beam flux    (W/m2)
                              !!\n nirdf - downward nir diffused flux       (W/m2)
                              !!\n visbm - downward uv+vis direct beam flux (W/m2)
                              !!\n visdf - downward uv+vis diffused flux    (W/m2)
    ! Outputs (mandatory)
    real(kind_phys), dimension(:), intent(inout) :: &
         tsflw,             & !< LW sfc air temp during calculation (K)
         sfcdlw,            & !< LW sfc all-sky     downward flux   (W/m2)
         sfculw,            & !< LW sfc all-sky     upward   flux   (W/m2)
         nirbmdi,           & !< SW sfc nir    beam downward flux   (W/m2)
         nirdfdi,           & !< SW sfc nir    diff downward flux   (W/m2)
         visbmdi,           & !< SW sfc uv+vis beam downward flux   (W/m2)
         visdfdi,           & !< SW sfc uv+vis diff downward flux   (W/m2)
         nirbmui,           & !< SW sfc nir    beam upward   flux   (W/m2)
         nirdfui,           & !< SW sfc nir    diff upward   flux   (W/m2)
         visbmui,           & !< SW sfc uv+vis beam upward   flux   (W/m2)
         visdfui,           & !< SW sfc uv+vis diff upward   flux   (W/m2)
         sfcnsw,            & !< SW sfc all-sky     net      flux   (W/m2) flux into ground
         sfcdsw               !< SW sfc all-sky     downward flux   (W/m2)
    real(kind_phys), dimension(:,:), intent(inout) :: &
          htrlw,             & !< LW all-sky heating rate (K/s)
          htrsw                !< SW all-sky heating rate (K/s)
    real(kind_phys), dimension(:), intent(inout) :: &
         total_albedo         !< Total sky albedo at TOA (W/m2)
    real(kind_phys), dimension(:,:), intent(inout), optional :: &
         htrlwu               !< LW all-sky heating-rate updated in-between radiation calls.
    type(sfcflw_type), dimension(:), intent(inout) :: &
         sfcflw               !< LW radiation fluxes at sfc
    type(sfcfsw_type), dimension(:), intent(inout) :: &
         sfcfsw               !< SW radiation fluxes at sfc
    type(topfsw_type), dimension(:), intent(inout) :: &
         topfsw               !< SW fluxes at top atmosphere
    type(topflw_type), dimension(:), intent(inout) :: &
         topflw               !< LW  fluxes at top atmosphere
    real(kind_phys), dimension(:,:), intent(inout) :: &
         fluxr                !< LW/SW diagnostics
    character(len=*), intent(out) :: &
         errmsg               !< CCPP error message
    integer, intent(out) :: &
         errflg               !< CCPP error code
    ! Outputs (optional)
    real(kind_phys),dimension(:,:),intent(inout),optional  :: &
         htrlwc,            & !< LW clear-sky heating-rate (K/s)
         htrswc               !< SW clear-sky heating rate (K/s)

    ! Local variables
    integer :: i
    real(kind_phys), dimension(nDay, nLev) :: thetaTendClrSky, thetaTendAllSky

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Only proceed if radiation is being called.
    if (.not. (doLWrad .or. doSWrad)) return

    ! #######################################################################################
    ! Longwave Radiation
    ! ####################################################################################### 
    if (doLWRad) then
      if (do_RRTMGP) then
        ! Clear-sky heating-rate (optional)
        if (do_lw_clrsky_hr) then
           call check_error_msg('GFS_radiation_post',compute_heating_rate(  &
                fluxlwUP_clrsky,   & ! IN  - RRTMGP upward longwave clear-sky flux profiles (W/m2)
                fluxlwDOWN_clrsky, & ! IN  - RRTMGP downward longwave clear-sky flux profiles (W/m2)
                p_lev,             & ! IN  - Pressure @ layer-interfaces (Pa)
                htrlwc))             ! OUT - Longwave clear-sky heating rate (K/sec)
        endif
    
          ! All-sky heating-rate (mandatory)
        call check_error_msg('GFS_radiation_post',compute_heating_rate(     &
             fluxlwUP_allsky,      & ! IN  - RRTMGP upward longwave all-sky flux profiles (W/m2)
             fluxlwDOWN_allsky,    & ! IN  - RRTMGP downward longwave all-sky flux profiles (W/m2)
             p_lev,                & ! IN  - Pressure @ layer-interfaces (Pa)
             htrlw))                 ! OUT - Longwave all-sky heating rate (K/sec)
    
        ! (Copy fluxes from RRTMGP types into model radiation types.) 
        ! TOA fluxes
        topflw(:)%upfxc = fluxlwUP_allsky(:,iTOA)
        topflw(:)%upfx0 = fluxlwUP_clrsky(:,iTOA)

        ! Surface fluxes
        sfcflw(:)%upfxc = fluxlwUP_allsky(:,iSFC)
        sfcflw(:)%upfx0 = fluxlwUP_clrsky(:,iSFC)
        sfcflw(:)%dnfxc = fluxlwDOWN_allsky(:,iSFC)
        sfcflw(:)%dnfx0 = fluxlwDOWN_clrsky(:,iSFC)

        ! Save surface air temp for diurnal adjustment at model t-steps
        tsflw(:) = tsfa(:)

        ! Radiation fluxes for other physics processes
        sfcdlw(:) = sfcflw(:)%dnfxc
        sfculw(:) = sfcflw(:)%upfxc

        ! Heating-rate at radiation timestep, used for adjustment between radiation calls.
        htrlwu = htrlw
      endif ! RRTMGP Longwave Radiaiton
    endif    ! ALL Longwave Radiation

    ! #######################################################################################
    ! Shortwave Radiation
    ! #######################################################################################
    if (doSWRad) then
      if (do_RRTMGP) then
        if (nDay .gt. 0) then
          ! Clear-sky heating-rate (optional)
          if (do_sw_clrsky_hr) then
             htrswc(:,:) = 0._kind_phys
             call check_error_msg('GFS_radiation_post',compute_heating_rate( &
                  fluxswUP_clrsky(idxday(1:nDay),:),   & ! IN  - Shortwave upward clear-sky flux profiles (W/m2)
                  fluxswDOWN_clrsky(idxday(1:nDay),:), & ! IN  - Shortwave downward clear-sky flux profiles (W/m2)
                  p_lev(idxday(1:nDay),:),             & ! IN  - Pressure at model-interface (Pa)
                  thetaTendClrSky))                      ! OUT - Clear-sky heating-rate (K/sec)
             htrswc(idxday(1:nDay),:)=thetaTendClrSky !**NOTE** GP doesn't use radiation levels, it uses the model fields. Not sure if this is necessary
          endif

          ! All-sky heating-rate (mandatory)
          htrsw(:,:) = 0._kind_phys
          call check_error_msg('GFS_radiation_post',compute_heating_rate(    &
               fluxswUP_allsky(idxday(1:nDay),:),      & ! IN  - Shortwave upward all-sky flux profiles (W/m2)
               fluxswDOWN_allsky(idxday(1:nDay),:),    & ! IN  - Shortwave downward all-sky flux profiles (W/m2)
               p_lev(idxday(1:nDay),:),                & ! IN  - Pressure at model-interface (Pa)
               thetaTendAllSky))                         ! OUT - All-sky heating-rate (K/sec)
          htrsw(idxday(1:nDay),:) = thetaTendAllSky

          ! (Copy fluxes from RRTMGP types into model radiation types.)

          ! TOA fluxes
          topfsw(:)%upfxc = fluxswUP_allsky(:,iTOA)
          topfsw(:)%upfx0 = fluxswUP_clrsky(:,iTOA)
          topfsw(:)%dnfxc = fluxswDOWN_allsky(:,iTOA)

          ! Surface fluxes
          sfcfsw(:)%upfxc = fluxswUP_allsky(:,iSFC)
          sfcfsw(:)%upfx0 = fluxswUP_clrsky(:,iSFC)
          sfcfsw(:)%dnfxc = fluxswDOWN_allsky(:,iSFC)
          sfcfsw(:)%dnfx0 = fluxswDOWN_clrsky(:,iSFC)

          ! Surface down and up spectral component fluxes
          ! - Save two spectral bands' surface downward and upward fluxes for output.
          do i=1,nCol
             nirbmdi(i) = scmpsw(i)%nirbm
             nirdfdi(i) = scmpsw(i)%nirdf
             visbmdi(i) = scmpsw(i)%visbm
             visdfdi(i) = scmpsw(i)%visdf
             nirbmui(i) = scmpsw(i)%nirbm * sfc_alb_nir_dir(i)
             nirdfui(i) = scmpsw(i)%nirdf * sfc_alb_nir_dif(i)
             visbmui(i) = scmpsw(i)%visbm * sfc_alb_uvvis_dir(i)
             visdfui(i) = scmpsw(i)%visdf * sfc_alb_uvvis_dif(i)
          enddo
        else  ! if_nday_block
          ! Dark everywhere
          htrsw(:,:) = 0.0
          sfcfsw     = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          topfsw     = topfsw_type( 0.0, 0.0, 0.0 )
          do i=1,nCol
             nirbmdi(i) = 0.0
             nirdfdi(i) = 0.0
             visbmdi(i) = 0.0
             visdfdi(i) = 0.0
             nirbmui(i) = 0.0
             nirdfui(i) = 0.0
             visbmui(i) = 0.0
             visdfui(i) = 0.0
          enddo

          if (do_sw_clrsky_hr) then
             htrswc(:,:) = 0
          endif
        endif                  ! end_if_nday

        ! Radiation fluxes for other physics processes
        ! *NOTE* For RRTMG, sfcnsw and sfcdsw are provided.
        !        For RRTMGP, we compute them here.
        do i=1,nCol
           sfcnsw(i) = sfcfsw(i)%dnfxc - sfcfsw(i)%upfxc
           sfcdsw(i) = sfcfsw(i)%dnfxc
        enddo
      endif ! RRTMGP Shortwave Radiaiton
    endif ! ALL Shortwave Radiation  

    ! The total sky (with clouds) shortwave albedo
    total_albedo = 0.0
    where(topfsw(:)%dnfxc>0) total_albedo(:) = topfsw(:)%upfxc/topfsw(:)%dnfxc

    ! #########################################################################################
    ! Compute radiation diagnostics
    ! #########################################################################################
    if (lssav) then
       call GFS_radiation_diagnostics(doLWrad, doSWrad, fhlwr, fhswr, coszen, coszdg, raddt,  &
            aerodp, cldsa, mtopa, mbota, cldtausw, cldtaulw, p_lev, tgrs, kb, kd, kt, sfcflw, &
            sfcfsw, topfsw, topflw, scmpsw, nCol, nDay, nLev, lmk, nfxr, nspc1, fluxr)
    endif

  end subroutine GFS_radiation_post_run
  
  ! ###########################################################################################
  ! GFS_radiation_diagnostics
  !
  ! For time averaged output quantities (including total-sky and clear-sky SW and LW fluxes at
  ! TOA and surface; conventional 3-domain cloud amount, cloud top and base pressure, and cloud
  ! top temperature; aerosols AOD, etc.), store computed results in corresponding slots of
  ! array <fluxr> with appropriate time weights.
  !
  ! ###########################################################################################
  subroutine GFS_radiation_diagnostics(doLWrad, doSWrad, fhlwr, fhswr, coszen, coszdg, raddt, &
       aerodp, cldsa, mtopa, mbota, cldtausw, cldtaulw, p_lev, tgrs, kb, kd, kt, sfcflw,      &
       sfcfsw, topfsw, topflw, scmpsw, nCol, nDay, nLev, lmk, nfxr, nspc1, fluxr)
    ! Inputs
    logical,           intent(in) :: doLWrad, doSWrad
    integer,           intent(in) :: nCol, nLev, lmk, nfxr, nspc1, nDay
    real(kind_phys),   intent(in) :: fhlwr, fhswr, coszen(nCol), coszdg(nCol), raddt
    real(kind_phys),   intent(in) :: aerodp(nCol,nspc1)
    real(kind_phys),   intent(in) :: cldtausw(nCol,lmk), cldtaulw(nCol,lmk)
    real(kind_phys),   intent(in) :: p_lev(nCol,nLev+1), tgrs(nCol,nLev)
    type(cmpfsw_type), intent(in) :: scmpsw(nCol)
    type(sfcflw_type), intent(in) :: sfcflw(nCol)
    type(sfcfsw_type), intent(in) :: sfcfsw(nCol)
    type(topfsw_type), intent(in) :: topfsw(nCol)
    type(topflw_type), intent(in) :: topflw(nCol)
    ! For High/Mid/Low cloud flux diagnsotics
    integer,           intent(in) :: kb, kd, kt
    integer,           intent(in) :: mtopa(nCol,3), mbota(nCol,3)
    real(kind_phys),   intent(in) :: cldsa(nCol,5)
    
    ! Outputs
    real(kind_phys), intent(inout) :: fluxr(nCol,nfxr)
    ! Locals
    integer :: i, j, k, itop, ibtc
    real(kind_phys) :: tem0d, tem1, tem2

    ! Save LW toa and sfc fluxes
    if (doLWrad) then
       do i=1,nCol
          ! LW total-sky fluxes
          fluxr(i,1 ) = fluxr(i,1 ) + fhlwr * topflw(i)%upfxc      ! total sky TOA LW up
          fluxr(i,19) = fluxr(i,19) + fhlwr * sfcflw(i)%dnfxc      ! total sky SFC LW down
          fluxr(i,20) = fluxr(i,20) + fhlwr * sfcflw(i)%upfxc      ! total sky SFC LW up
          ! LW clear-sky fluxes
          fluxr(i,28) = fluxr(i,28) + fhlwr * topflw(i)%upfx0      ! clear sky TOA LW up
          fluxr(i,30) = fluxr(i,30) + fhlwr * sfcflw(i)%dnfx0      ! clear sky SFC LW down
          fluxr(i,33) = fluxr(i,33) + fhlwr * sfcflw(i)%upfx0      ! clear sky SFC LW up
       enddo
    endif  ! END DOLWRAD

    ! Save SW toa and sfc fluxes with proper diurnal sw wgt. coszen=mean cosz over daylight
    ! part of sw calling interval, while coszdg= mean cosz over entire interval
    if (doSWrad) then
       do i=1,nCol
          ! Aerosol optical-depths
          fluxr(i,34) = aerodp(i,1)  ! Total aod at 550nm
          fluxr(i,35) = aerodp(i,2)  ! Dust aod at 550nm
          fluxr(i,36) = aerodp(i,3)  ! Soot aod at 550nm
          fluxr(i,37) = aerodp(i,4)  ! Waso aod at 550nm
          fluxr(i,38) = aerodp(i,5)  ! Suso aod at 550nm
          fluxr(i,39) = aerodp(i,6)  ! Salt aod at 550nm
          
          if (coszen(i) > 0.) then
             ! SW total-sky fluxes
             tem0d = fhswr * coszdg(i) / coszen(i)
             fluxr(i,2 ) = fluxr(i,2)  + topfsw(i)%upfxc * tem0d   ! total sky TOA SW up
             fluxr(i,3 ) = fluxr(i,3)  + sfcfsw(i)%upfxc * tem0d   ! total sky SFC SW up
             fluxr(i,4 ) = fluxr(i,4)  + sfcfsw(i)%dnfxc * tem0d   ! total sky SFC SW down
             ! SW uv-b fluxes
             fluxr(i,21) = fluxr(i,21) + scmpsw(i)%uvbfc * tem0d   ! total sky uv-b SW down
             fluxr(i,22) = fluxr(i,22) + scmpsw(i)%uvbf0 * tem0d   ! clear sky uv-b SW down
             ! SW toa incoming fluxes
             fluxr(i,23) = fluxr(i,23) + topfsw(i)%dnfxc * tem0d   ! TOA SW down
             ! SW sfc flux components
             fluxr(i,24) = fluxr(i,24) + scmpsw(i)%visbm * tem0d   ! uv/vis beam SW down
             fluxr(i,25) = fluxr(i,25) + scmpsw(i)%visdf * tem0d   ! uv/vis diff SW down
             fluxr(i,26) = fluxr(i,26) + scmpsw(i)%nirbm * tem0d   ! nir beam SW down
             fluxr(i,27) = fluxr(i,27) + scmpsw(i)%nirdf * tem0d   ! nir diff SW down
             ! SW clear-sky fluxes
             fluxr(i,29) = fluxr(i,29) + topfsw(i)%upfx0 * tem0d   ! clear sky TOA SW up
             fluxr(i,31) = fluxr(i,31) + sfcfsw(i)%upfx0 * tem0d   ! clear sky SFC SW up
             fluxr(i,32) = fluxr(i,32) + sfcfsw(i)%dnfx0 * tem0d   ! clear sky SFC SW down
          endif
       enddo
    endif  ! END DOSWRAD

    !
    ! High/Mid/Low diagnostics
    !
    if (doLWrad .or. doSWrad) then
       ! Save total and boundary layer clouds
       do i=1,nCol
          fluxr(i,17) = fluxr(i,17) + raddt * cldsa(i,4)
          fluxr(i,18) = fluxr(i,18) + raddt * cldsa(i,5)
       enddo

       ! Save cld frac,toplyr,botlyr and top temp, note that the order
       ! of h,m,l cloud is reversed for the fluxr output.
       ! save interface pressure (pa) of top/bot
       do j = 1, 3
          do i = 1, nCol
             tem0d = raddt * cldsa(i,j)
             itop  = mtopa(i,j) - kd
             ibtc  = mbota(i,j) - kd
             fluxr(i, 8-j) = fluxr(i, 8-j) + tem0d
             fluxr(i,11-j) = fluxr(i,11-j) + tem0d * p_lev(i,itop+kt)
             fluxr(i,14-j) = fluxr(i,14-j) + tem0d * p_lev(i,ibtc+kb)
             fluxr(i,17-j) = fluxr(i,17-j) + tem0d * tgrs(i,itop)
          enddo
       enddo

       ! In-cloud (shortwave) optical depth at approx .55 um channel
       if (doSWrad .and. (nDay > 0)) then
          do j = 1, 3
             do i = 1, nCol
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j) - kd
                ibtc  = mbota(i,j) - kd
                tem1 = 0.
                do k=ibtc,itop
                   tem1 = tem1 + cldtausw(i,k)
                enddo
                fluxr(i,43-j) = fluxr(i,43-j) + tem0d * tem1
             enddo
          enddo
       endif  ! END DOSWRAD

       ! In-cloud (longwave) optical depth at approx 10. um channel
       if (doLWrad) then
          do j = 1, 3
             do i = 1, nCol
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j) - kd
                ibtc  = mbota(i,j) - kd
                tem2 = 0.
                do k=ibtc,itop
                   tem2 = tem2 + cldtaulw(i,k)
                enddo
                fluxr(i,46-j) = fluxr(i,46-j) + tem0d * (1.0-exp(-tem2))
             enddo
          enddo
       endif  ! END DOLWRAD
    endif     ! END DOSWRAD OR DOLWRAD

  end subroutine GFS_radiation_diagnostics

end module GFS_radiation_post
