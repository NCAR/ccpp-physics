!> \file GFS_rrtmgp_post.F90
!!
!> \defgroup GFS_rrtmgp_post GFS_rrtmgp_post.F90
!!
!! \brief RRTMGP post-processing routine.
!!
module GFS_rrtmgp_post
  use machine,                   only: kind_phys
  use module_radlw_parameters,   only: topflw_type, sfcflw_type
  use module_radsw_parameters,   only: topfsw_type, sfcfsw_type, cmpfsw_type
  use mo_heating_rates,          only: compute_heating_rate
  use radiation_tools,           only: check_error_msg
  implicit none

  public GFS_rrtmgp_post_run

contains
  ! ######################################################################################## 
!>\defgroup gfs_rrtmgp_post_mod GFS RRTMGP Post Module
!> \section arg_table_GFS_rrtmgp_post_run
!! \htmlinclude GFS_rrtmgp_post.html
!!
!! \ingroup GFS_rrtmgp_post
!!
!! \brief The all-sky radiation tendency is computed, the clear-sky tendency is computed 
!! if requested.
!!
!! RRTMGP surface and TOA fluxes are copied to fields that persist between radiation/physics
!! calls.
!! 
!! (optional) Save additional diagnostics.
!!
!! \section GFS_rrtmgp_post_run
  ! ######################################################################################## 
  subroutine GFS_rrtmgp_post_run (nCol, nLev, nDay, iSFC, iTOA, idxday, doLWrad, doSWrad,  &
       do_lw_clrsky_hr, do_sw_clrsky_hr, save_diag, fhlwr, fhswr, sfc_alb_nir_dir,         &
       sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, p_lev, tsfa, coszen, coszdg, &
       fluxlwDOWN_clrsky, fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky,             &
       fluxswDOWN_clrsky, fluxswUP_allsky, fluxswDOWN_allsky, fluxswUP_clrsky,             &
       raddt, aerodp, cldsa, mtopa, mbota, cld_frac, cldtaulw, cldtausw, scmpsw, fluxr,    &
       sfcdlw, sfculw, sfcflw, tsflw, htrlw, htrlwu, topflw, nirbmdi, nirdfdi, visbmdi,    &
       visdfdi, nirbmui, nirdfui, visbmui, visdfui, sfcnsw, sfcdsw, htrsw, sfcfsw, topfsw, &
       htrswc, htrlwc, errmsg, errflg)

    ! Inputs
    integer, intent(in) ::  &
         nCol,              & ! Horizontal loop extent 
         nLev,              & ! Number of vertical layers
         nDay,              & ! Number of daylit columns
         iSFC,              & ! Vertical index for surface level
         iTOA                 ! Vertical index for TOA level
    integer, intent(in), dimension(:) :: &
         idxday               ! Index array for daytime points
    integer, intent(in), dimension(:,:) :: &
         mbota,             & ! Vertical indices for low, middle and high cloud tops
         mtopa                ! ertical indices for low, middle and high cloud bases
    logical, intent(in) :: & 
         doLWrad,           & ! Logical flags for lw radiation calls
         doSWrad,           & ! Logical flags for sw radiation calls
         do_lw_clrsky_hr,   & ! Output clear-sky LW heating-rate?
         do_sw_clrsky_hr,   & ! Output clear-sky SW heating-rate? 
         save_diag            ! Output radiation diagnostics?
    real(kind_phys), intent(in) :: &
         fhlwr,             & ! Frequency for LW radiation calls
         fhswr                ! Frequency for SW radiation calls
    real(kind_phys), dimension(:), intent(in) ::  &
         tsfa,              & ! Lowest model layer air temperature for radiation (K)
         coszen,            & ! Cosine(SZA)
         coszdg,            & ! Cosine(SZA), daytime
         sfc_alb_nir_dir,   & ! Surface albedo (direct) 
         sfc_alb_nir_dif,   & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir, & ! Surface albedo (direct)
         sfc_alb_uvvis_dif    ! Surface albedo (diffuse)
    real(kind_phys), dimension(:,:), intent(in) :: &
         p_lev,             & ! Pressure @ model layer-interfaces (Pa)
         fluxlwUP_allsky,   & ! RRTMGP longwave all-sky flux      (W/m2)
         fluxlwDOWN_allsky, & ! RRTMGP longwave all-sky flux      (W/m2)
         fluxlwUP_clrsky,   & ! RRTMGP longwave clear-sky flux    (W/m2)
         fluxlwDOWN_clrsky, & ! RRTMGP longwave clear-sky flux    (W/m2)
         fluxswUP_allsky,   & ! RRTMGP shortwave all-sky flux     (W/m2)
         fluxswDOWN_allsky, & ! RRTMGP shortwave all-sky flux     (W/m2)
         fluxswUP_clrsky,   & ! RRTMGP shortwave clear-sky flux   (W/m2)
         fluxswDOWN_clrsky    ! RRTMGP shortwave clear-sky flux   (W/m2)
    real(kind_phys), intent(in) :: &
         raddt                ! Radiation time step
    real(kind_phys), dimension(:,:), intent(in) :: &
         aerodp,            & ! Vertical integrated optical depth for various aerosol species
         cldsa,             & ! Fraction of clouds for low, middle, high, total and BL
         cld_frac,          & ! Total cloud fraction in each layer
         cldtaulw,          & ! approx 10.mu band layer cloud optical depth
         cldtausw             ! approx .55mu band layer cloud optical depth
    type(cmpfsw_type), dimension(:), intent(in) :: &
         scmpsw               ! 2D surface fluxes, components:
                              ! uvbfc - total sky downward uv-b flux at  (W/m2)
                              ! uvbf0 - clear sky downward uv-b flux at  (W/m2)
                              ! nirbm - downward nir direct beam flux    (W/m2)
                              ! nirdf - downward nir diffused flux       (W/m2)
                              ! visbm - downward uv+vis direct beam flux (W/m2)
                              ! visdf - downward uv+vis diffused flux    (W/m2)


    real(kind=kind_phys), dimension(:,:), intent(inout) :: fluxr

    ! Outputs (mandatory)
    real(kind_phys), dimension(:), intent(inout) :: &
         tsflw,             & ! LW sfc air temp during calculation (K)
         sfcdlw,            & ! LW sfc all-sky     downward flux   (W/m2)
         sfculw,            & ! LW sfc all-sky     upward   flux   (W/m2)
         nirbmdi,           & ! SW sfc nir    beam downward flux   (W/m2)
         nirdfdi,           & ! SW sfc nir    diff downward flux   (W/m2)
         visbmdi,           & ! SW sfc uv+vis beam downward flux   (W/m2)
         visdfdi,           & ! SW sfc uv+vis diff downward flux   (W/m2)
         nirbmui,           & ! SW sfc nir    beam upward   flux   (W/m2)
         nirdfui,           & ! SW sfc nir    diff upward   flux   (W/m2)
         visbmui,           & ! SW sfc uv+vis beam upward   flux   (W/m2)
         visdfui,           & ! SW sfc uv+vis diff upward   flux   (W/m2)
         sfcnsw,            & ! SW sfc all-sky     net      flux   (W/m2) flux into ground
         sfcdsw               ! SW sfc all-sky     downward flux   (W/m2)
    real(kind_phys), dimension(:,:), intent(inout) :: &
         htrlw,             & ! LW all-sky heating rate (K/s)
         htrsw,             & ! SW all-sky heating rate (K/s)
         htrlwu               ! LW all-sky heating-rate updated in-between radiation calls.
    type(sfcflw_type), dimension(:), intent(inout) :: &
         sfcflw               ! LW radiation fluxes at sfc
    type(sfcfsw_type), dimension(:), intent(inout) :: &
         sfcfsw               ! SW radiation fluxes at sfc
    type(topfsw_type), dimension(:), intent(inout) :: &
         topfsw               ! SW fluxes at top atmosphere
    type(topflw_type), dimension(:), intent(inout) :: &
         topflw               ! LW  fluxes at top atmosphere
    character(len=*), intent(out) :: &
         errmsg               ! CCPP error message
    integer, intent(out) :: &
         errflg               ! CCPP error code

    ! Outputs (optional)
    real(kind_phys),dimension(:,:),intent(inout),optional  :: &
         htrlwc,            & ! LW clear-sky heating-rate (K/s)
         htrswc               ! SW clear-sky heating rate (K/s)

    ! Local variables
    integer :: i, j, k, itop, ibtc
    real(kind_phys) :: tem0d, tem1, tem2
    real(kind_phys), dimension(nDay, nLev) :: thetaTendClrSky, thetaTendAllSky

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. (doLWrad .or. doSWrad)) return

    if (doLWRad) then
       ! #######################################################################################
       ! Compute LW heating-rates.
       ! #######################################################################################

       ! Clear-sky heating-rate (optional)
       if (do_lw_clrsky_hr) then
          call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(  &
               fluxlwUP_clrsky,   & ! IN  - RRTMGP upward longwave clear-sky flux profiles (W/m2)
               fluxlwDOWN_clrsky, & ! IN  - RRTMGP downward longwave clear-sky flux profiles (W/m2)
               p_lev,             & ! IN  - Pressure @ layer-interfaces (Pa)
               htrlwc))             ! OUT - Longwave clear-sky heating rate (K/sec)
       endif

       ! All-sky heating-rate (mandatory)
       call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(     &
            fluxlwUP_allsky,      & ! IN  - RRTMGP upward longwave all-sky flux profiles (W/m2)
            fluxlwDOWN_allsky,    & ! IN  - RRTMGP downward longwave all-sky flux profiles (W/m2)
            p_lev,                & ! IN  - Pressure @ layer-interfaces (Pa)
            htrlw))                 ! OUT - Longwave all-sky heating rate (K/sec)

       ! #######################################################################################
       ! Save LW outputs.
       ! (Copy fluxes from RRTMGP types into model radiation types.)
       ! #######################################################################################
       ! TOA fluxes

       topflw(:)%upfxc = fluxlwUP_allsky(:,iTOA)
       topflw(:)%upfx0 = fluxlwUP_clrsky(:,iTOA)
       
       ! Surface fluxes
       sfcflw(:)%upfxc = fluxlwUP_allsky(:,iSFC)
       sfcflw(:)%upfx0 = fluxlwUP_clrsky(:,iSFC)
       sfcflw(:)%dnfxc = fluxlwDOWN_allsky(:,iSFC)
       sfcflw(:)%dnfx0 = fluxlwDOWN_clrsky(:,iSFC)
       
       ! Save surface air temp for diurnal adjustment at model t-steps
       tsflw (:) = tsfa(:)
       
       ! Radiation fluxes for other physics processes
       sfcdlw(:) = sfcflw(:)%dnfxc
       sfculw(:) = sfcflw(:)%upfxc
       
       ! Heating-rate at radiation timestep, used for adjustment between radiation calls.
       htrlwu = htrlw
       
       ! #######################################################################################
       ! Save LW diagnostics
       ! - For time averaged output quantities (including total-sky and clear-sky SW and LW
       !   fluxes at TOA and surface; conventional 3-domain cloud amount, cloud top and base
       !   pressure, and cloud top temperature; aerosols AOD, etc.), store computed results in
       !   corresponding slots of array fluxr with appropriate time weights.
       ! - Collect the fluxr data for wrtsfc
       ! #######################################################################################
       if (save_diag) then
          do i=1,nCol
             ! LW all-sky fluxes
             fluxr(i,1 ) = fluxr(i,1 ) + fhlwr * fluxlwUP_allsky(  i,iTOA)   ! total sky top lw up
             fluxr(i,19) = fluxr(i,19) + fhlwr * fluxlwDOWN_allsky(i,iSFC)   ! total sky sfc lw dn
             fluxr(i,20) = fluxr(i,20) + fhlwr * fluxlwUP_allsky(  i,iSFC)   ! total sky sfc lw up
             ! LW clear-sky fluxes
             fluxr(i,28) = fluxr(i,28) + fhlwr * fluxlwUP_clrsky(  i,iTOA)   ! clear sky top lw up
             fluxr(i,30) = fluxr(i,30) + fhlwr * fluxlwDOWN_clrsky(i,iSFC)   ! clear sky sfc lw dn
             fluxr(i,33) = fluxr(i,33) + fhlwr * fluxlwUP_clrsky(  i,iSFC)   ! clear sky sfc lw up
          enddo
          
          ! Save cld frac,toplyr,botlyr and top temp, note that the order of h,m,l cloud is reversed for 
          ! the fluxr output. save interface pressure (pa) of top/bot
          do j = 1, 3
             do i = 1, nCol
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j)
                ibtc  = mbota(i,j)
                
                ! Add optical depth and emissivity output
                tem2 = 0.
                do k=ibtc,itop
                   tem2 = tem2 + cldtaulw(i,k)      ! approx 10. mu channel
                enddo
                fluxr(i,46-j) = fluxr(i,46-j) + tem0d * (1.0-exp(-tem2))
             enddo
          enddo
       endif
    endif
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    if (doSWRad) then
       if (nDay .gt. 0) then
          ! #################################################################################
          ! Compute SW heating-rates
          ! #################################################################################

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
          
          ! #################################################################################
          ! Save SW outputs
          ! (Copy fluxes from RRTMGP types into model radiation types.)
          ! #################################################################################
          
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
       else                   ! if_nday_block
          ! #################################################################################
          ! Dark everywhere
          ! #################################################################################
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
       do i=1,nCol
          sfcnsw(i) = sfcfsw(i)%dnfxc - sfcfsw(i)%upfxc
          sfcdsw(i) = sfcfsw(i)%dnfxc
       enddo
       
       ! #################################################################################
       ! Save SW diagnostics
       ! - For time averaged output quantities (including total-sky and clear-sky SW and LW 
       !   fluxes at TOA and surface; conventional 3-domain cloud amount, cloud top and base 
       !   pressure, and cloud top temperature; aerosols AOD, etc.), store computed results in
       !   corresponding slots of array fluxr with appropriate time weights.
       ! - Collect the fluxr data for wrtsfc
       ! #################################################################################
       if (save_diag) then
          do i=1,nCol
             fluxr(i,34) = aerodp(i,1)  ! total aod at 550nm
             fluxr(i,35) = aerodp(i,2)  ! DU aod at 550nm
             fluxr(i,36) = aerodp(i,3)  ! BC aod at 550nm
             fluxr(i,37) = aerodp(i,4)  ! OC aod at 550nm
             fluxr(i,38) = aerodp(i,5)  ! SU aod at 550nm
             fluxr(i,39) = aerodp(i,6)  ! SS aod at 550nm
             if (coszen(i) > 0.) then
                ! SW all-sky fluxes
                tem0d = fhswr * coszdg(i) / coszen(i)
                fluxr(i,2 ) = fluxr(i,2)  + topfsw(i)%upfxc * tem0d  ! total sky top sw up
                fluxr(i,3 ) = fluxr(i,3)  + sfcfsw(i)%upfxc * tem0d  
                fluxr(i,4 ) = fluxr(i,4)  + sfcfsw(i)%dnfxc * tem0d  ! total sky sfc sw dn
                ! SW uv-b fluxes
                fluxr(i,21) = fluxr(i,21) + scmpsw(i)%uvbfc * tem0d  ! total sky uv-b sw dn
                fluxr(i,22) = fluxr(i,22) + scmpsw(i)%uvbf0 * tem0d  ! clear sky uv-b sw dn
                ! SW TOA incoming fluxes
                fluxr(i,23) = fluxr(i,23) + topfsw(i)%dnfxc * tem0d  ! top sw dn 
                ! SW SFC flux components
                fluxr(i,24) = fluxr(i,24) + visbmdi(i) * tem0d       ! uv/vis beam sw dn
                fluxr(i,25) = fluxr(i,25) + visdfdi(i) * tem0d       ! uv/vis diff sw dn
                fluxr(i,26) = fluxr(i,26) + nirbmdi(i) * tem0d       ! nir beam sw dn
                fluxr(i,27) = fluxr(i,27) + nirdfdi(i) * tem0d       ! nir diff sw dn
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
    endif

  end subroutine GFS_rrtmgp_post_run
end module GFS_rrtmgp_post
