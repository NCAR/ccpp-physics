!> \file GFS_radiation_post.F90
!! RRTMGP post-processing routine.
!!
module GFS_radiation_post
  use machine,                   only: kind_phys
  use module_radlw_parameters,   only: topflw_type, sfcflw_type
  use module_radsw_parameters,   only: topfsw_type, sfcfsw_type, cmpfsw_type
  use mo_heating_rates,          only: compute_heating_rate
  use radiation_tools,           only: check_error_msg
  implicit none

  public GFS_radiation_post_run

contains

!> \section arg_table_GFS_radiation_post_run Argument Table
!! \htmlinclude GFS_radiation_post_run.html
!!
!!The all-sky radiation tendency is computed, the clear-sky tendency is computed 
!! if requested.
!!
!! RRTMGP surface and TOA fluxes are copied to fields that persist between radiation/physics
!! calls.
!!
!! (optional) Save additional diagnostics.
  subroutine GFS_radiation_post_run (nCol, nLev, nDay, iSFC, iTOA, idxday, doLWrad, doSWrad,  &
       do_lw_clrsky_hr, do_sw_clrsky_hr, do_RRTMGP, sfc_alb_nir_dir,                       &
       sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, p_lev, tsfa,                 &
       fluxlwDOWN_clrsky, fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky,             &
       fluxswDOWN_clrsky, fluxswUP_allsky, fluxswDOWN_allsky, fluxswUP_clrsky, scmpsw,     &
       sfcdlw, sfculw, sfcflw, tsflw, htrlw, htrlwu, topflw, nirbmdi, nirdfdi, visbmdi,    &
       visdfdi, nirbmui, nirdfui, visbmui, visdfui, sfcnsw, sfcdsw, htrsw, sfcfsw, topfsw, &
       htrswc, htrlwc, total_albedo, errmsg, errflg)

    ! Inputs
    integer, intent(in) ::  &
         nCol,              & !< Horizontal loop extent 
         nLev,              & !< Number of vertical layers
         nDay,              & !< Number of daylit columns
         iSFC,              & !< Vertical index for surface level
         iTOA                 !< Vertical index for TOA level
    integer, intent(in), dimension(:) :: &
         idxday               ! Index array for daytime points
    logical, intent(in) :: & 
         doLWrad,           & ! Logical flags for lw radiation calls
         doSWrad,           & ! Logical flags for sw radiation calls
         do_lw_clrsky_hr,   & ! Output clear-sky LW heating-rate?
         do_sw_clrsky_hr,   & ! Output clear-sky SW heating-rate? 
         do_RRTMGP            ! Flag for using RRTMGP scheme
    real(kind_phys), dimension(:), intent(in) ::  &
         tsfa,              & ! Lowest model layer air temperature for radiation (K)
         sfc_alb_nir_dir,   & ! Surface albedo (direct) 
         sfc_alb_nir_dif,   & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir, & ! Surface albedo (direct)
         sfc_alb_uvvis_dif    ! Surface albedo (diffuse)
    real(kind_phys), dimension(:,:), intent(in), optional :: &
         p_lev,             & ! Pressure @ model layer-interfaces (Pa)
         fluxlwUP_allsky,   & ! RRTMGP longwave all-sky flux      (W/m2)
         fluxlwDOWN_allsky, & ! RRTMGP longwave all-sky flux      (W/m2)
         fluxlwUP_clrsky,   & ! RRTMGP longwave clear-sky flux    (W/m2)
         fluxlwDOWN_clrsky, & ! RRTMGP longwave clear-sky flux    (W/m2)
         fluxswUP_allsky,   & ! RRTMGP shortwave all-sky flux     (W/m2)
         fluxswDOWN_allsky, & ! RRTMGP shortwave all-sky flux     (W/m2)
         fluxswUP_clrsky,   & ! RRTMGP shortwave clear-sky flux   (W/m2)
         fluxswDOWN_clrsky    ! RRTMGP shortwave clear-sky flux   (W/m2)
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
         htrlw,             & ! LW all-sky heating rate (K/s)
         htrsw                ! SW all-sky heating rate (K/s)
    real(kind_phys), dimension(nCol), intent(inout) :: &
         total_albedo         ! Total sky albedo at TOA (W/m2)
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
    character(len=*), intent(out) :: &
         errmsg               !< CCPP error message
    integer, intent(out) :: &
         errflg               !< CCPP error code

    ! Outputs (optional)
    real(kind_phys),dimension(:,:),intent(inout),optional  :: &
         htrlwc,            & !< LW clear-sky heating-rate (K/s)
         htrswc               !< SW clear-sky heating rate (K/s)

    ! Local variables
    integer :: i, j, k, itop, ibtc
    real(kind_phys) :: tem0d, tem1, tem2
    real(kind_phys), dimension(nDay, nLev) :: thetaTendClrSky, thetaTendAllSky

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. (doLWrad .or. doSWrad)) return

    if (doLWRad) then
       if (do_RRTMGP) then
          ! #######################################################################################
          ! Compute LW heating-rates.
          ! #######################################################################################

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
       endif
    
!  ---  The total sky (with clouds) shortwave albedo
       total_albedo = 0.0
       where(topfsw(:)%dnfxc>0) total_albedo(:) = topfsw(:)%upfxc/topfsw(:)%dnfxc
    endif
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    if (doSWRad .and. do_RRTMGP) then
       if (nDay .gt. 0) then
          ! #################################################################################
          ! Compute SW heating-rates
          ! #################################################################################

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
       
    endif

  end subroutine GFS_radiation_post_run
end module GFS_radiation_post
