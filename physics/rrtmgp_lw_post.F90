! ###########################################################################################
! ###########################################################################################
module rrtmgp_lw_post
  use machine,                 only: kind_phys
  use mo_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp
  use mo_fluxes_byband,        only: ty_fluxes_byband
  use mo_heating_rates,        only: compute_heating_rate
  use rrtmgp_lw,               only: check_error_msg
  use module_radlw_parameters, only: topflw_type, sfcflw_type, proflw_type

  implicit none

  ! Logical flags for optional output fields in rrtmgp_lw_post_run(), default=.false.
  logical :: &
       l_AllSky_HR_byband  = .false., & ! 2D [ncol,nlay] all-sky heating rates, in each band [ncol,nlay,nBandsLW]?
       l_ClrSky_HR         = .false., & ! 2D [ncol,nlay] clear-sky heating rate?
       l_fluxes2D          = .false.    ! 2D [ncol,nlay] radiative fluxes? *Note* fluxes is a DDT w/ 4 fields.

  public rrtmgp_lw_post_init, rrtmgp_lw_post_run, rrtmgp_lw_post_finalize
contains

  subroutine rrtmgp_lw_post_init()
  end subroutine rrtmgp_lw_post_init
  
  ! #########################################################################################
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_post_run Argument Table
!! | local_name            | standard_name                                                                                 | long_name                                                          | units | rank | type                  |    kind   | intent | optional |
!! |-----------------------|-----------------------------------------------------------------------------------------------|--------------------------------------------------------------------|-------|------|-----------------------|-----------|--------|----------|
!! | ncol                  | horizontal_loop_extent                                                                        | horizontal dimension                                               | count |    0 | integer               |           | in     | F        |
!! | nlay                  | adjusted_vertical_layer_dimension_for_radiation                                               | number of vertical layers for radiation                            | count |    0 | integer               |           | in     | F        |
!! | p_lev                 | air_pressure_at_interface_for_radiation_in_hPa                                                | air pressure level                                                 | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | fluxLW_allsky         | lw_flux_profiles_byband_allsky                                                                | Fortran DDT containing RRTMGP 3D fluxes                            | DDT   |    0 | ty_fluxes_byband      |           | in     | F        |
!! | fluxLW_clrsky         | lw_flux_profiles_byband_clrsky                                                                | Fortran DDT containing RRTMGP 3D fluxes                            | DDT   |    0 | ty_fluxes_byband      |           | in     | F        |
!! | kdist_lw              | K_distribution_file_for_RRTMGP_LW_scheme                                                      | DDT containing spectral information for RRTMGP LW radiation scheme | DDT   |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | errmsg                | ccpp_error_message                                                                            | error message for error handling in CCPP                           | none  |    0 | character             | len=*     | out    | F        |
!! | errflg                | ccpp_error_flag                                                                               | error flag for error handling in CCPP                              | flag  |    0 | integer               |           | out    | F        |
!! | hlwc                  | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step                    | longwave total sky heating rate                                    | K s-1 |    2 | real                  | kind_phys | out    | F        |
!! | topflx                | lw_fluxes_top_atmosphere                                                                      | longwave total sky fluxes at the top of the atm                    | W m-2 |    1 | topflw_type           |           | inout  | F        |
!! | sfcflx                | lw_fluxes_sfc                                                                                 | longwave total sky fluxes at the Earth surface                     | W m-2 |    1 | sfcflw_type           |           | inout  | F        |
!! | hlw0                  | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step | longwave clear sky heating rate                                    | K s-1 |    2 | real                  | kind_phys | inout  | T        |
!! | hlwb                  | lw_heating_rate_spectral                                                                      | longwave total sky heating rate (spectral)                         | K s-1 |    3 | real                  | kind_phys | inout  | T        |
!! | flxprf                | lw_fluxes                                                                                     | lw fluxes total sky / csk and up / down at levels                  | W m-2 |    2 | proflw_type           |           | inout  | T        |
!!
  subroutine rrtmgp_lw_post_run(ncol, nlay, p_lev, kdist_lw, fluxLW_allsky, fluxLW_clrsky, &
              hlwc, topflx, sfcflx,  hlw0, hlwb, flxprf, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: &
         ncol,          & ! Number of horizontal gridpoints
         nlay             ! Number of vertical layers
    real(kind_phys), dimension(ncol,nlay+1), intent(in) :: &
         p_lev            ! Pressure @ model layer-interfaces (hPa)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         kdist_lw         ! DDT containing LW spectral information
    type(ty_fluxes_byband),intent(in) :: &
         fluxLW_allsky, & ! All-sky flux                      (W/m2)
         fluxLW_clrsky    ! Clear-sky flux                    (W/m2)

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg
    real(kind_phys),dimension(ncol,nlay),intent(inout) :: &
         hlwc             ! All-sky heating-rate              (K/sec)
    type(topflw_type), dimension(ncol), intent(inout) :: &
         topflx           ! radiation fluxes at top, components:
                          ! upfxc - total sky upward flux at top   (w/m2)
                          ! upfx0 - clear sky upward flux at top   (w/m2)
    type(sfcflw_type), dimension(ncol), intent(inout) :: &
         sfcflx           ! radiation fluxes at sfc, components:
                          ! upfxc - total sky upward flux at sfc   (w/m2)  
                          ! upfx0 - clear sky upward flux at sfc   (w/m2)
                          ! dnfxc - total sky downward flux at sfc (w/m2)
                          ! dnfx0 - clear sky downward flux at sfc (w/m2)

    ! Outputs (optional)
    real(kind_phys), dimension(ncol,nlay,kdist_lw%get_nband()), optional, intent(inout) :: &
         hlwb             ! All-sky heating rate, by band     (K/sec)
    real(kind_phys), dimension(ncol,nlay), optional, intent(inout) :: &
         hlw0             ! Clear-sky heating rate            (K/sec)
    type(proflw_type), dimension(ncol,nlay+1), optional, intent(inout) :: &
         flxprf           ! 2D radiative fluxes, components:
                          ! upfxc - total sky upward flux     (W/m2)
                          ! dnfxc - total sky dnward flux     (W/m2)
                          ! upfx0 - clear sky upward flux     (W/m2)
                          ! dnfx0 - clear sky dnward flux     (W/m2)

    ! Local variables
    integer :: iBand, iTOA, iSFC
    logical :: top_at_1
    real(kind_phys), dimension(ncol,nlay) :: thetaTendClrSky, thetaTendAllSky
    real(kind_phys), dimension(ncol,nlay,kdist_lw%get_nband()) :: thetaTendByBandAllSky
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! What is vertical ordering?
    top_at_1 = (p_lev(1,1) .lt. p_lev(1,nlay))
    if (top_at_1) then 
       iSFC = nlay+1
       iTOA = 1
    else
       iSFC = 1
       iTOA = nlay+1
    endif

    ! Are any optional outputs requested?
    l_ClrSky_HR        = present(hlw0)
    l_AllSky_HR_byband = present(hlwb)
    l_fluxes2D         = present(flxprf)

    ! #######################################################################################
    ! Compute heating rates
    ! #######################################################################################
    if (l_ClrSky_HR) then
       call check_error_msg(compute_heating_rate(     &
            fluxLW_clrsky%flux_up,                    &
            fluxLW_clrsky%flux_dn,                    &
            p_lev(1:ncol,1:nlay+1),                   &
            thetaTendClrSky))
    endif
    if (l_AllSky_HR_byband) then
       do iBand=1,kdist_lw%get_nband()
          call check_error_msg(compute_heating_rate(  &
               fluxLW_allsky%bnd_flux_up(:,:,iBand),  &
               fluxLW_allsky%bnd_flux_dn(:,:,iBand),  &
               p_lev(1:ncol,1:nlay+1),                &
               thetaTendByBandAllSky(:,:,iBand)))
       enddo
    else
       call check_error_msg(compute_heating_rate(     &
            fluxLW_allsky%flux_up,                    &
            fluxLW_allsky%flux_dn,                    &
            p_lev(1:ncol,1:nlay+1),                   &
            thetaTendAllSky))
    endif    

    ! #######################################################################################
    ! Copy fluxes from RRTGMP types into model radiation types.
    ! #######################################################################################
    ! Mandatory outputs
    topflx%upfxc = fluxLW_allsky%flux_up(:,iTOA)
    topflx%upfx0 = fluxLW_clrsky%flux_up(:,iTOA)
    sfcflx%upfxc = fluxLW_allsky%flux_up(:,iSFC)
    sfcflx%upfx0 = fluxLW_clrsky%flux_up(:,iSFC)
    sfcflx%dnfxc = fluxLW_allsky%flux_dn(:,iSFC)
    sfcflx%dnfx0 = fluxLW_clrsky%flux_dn(:,iSFC)
    !cldtau       = optical_props_cldy%tau(:,:,7)
    hlwc         = thetaTendAllSky

    ! Optional output
    if(l_fluxes2D) then
       flxprf%upfxc = fluxLW_allsky%flux_up
       flxprf%dnfxc = fluxLW_allsky%flux_dn
       flxprf%upfx0 = fluxLW_clrsky%flux_up
       flxprf%dnfx0 = fluxLW_clrsky%flux_dn
    endif
    if (l_AllSky_HR_byband) then
       hlwb = thetaTendByBandAllSky
    endif
    if (l_ClrSky_HR) then
       hlw0 = thetaTendClrSky
    endif


  end subroutine rrtmgp_lw_post_run
  
  subroutine rrtmgp_lw_post_finalize()
  end subroutine rrtmgp_lw_post_finalize
  


end module rrtmgp_lw_post
