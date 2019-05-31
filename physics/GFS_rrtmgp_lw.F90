module GFS_rrtmgp_lw
 use GFS_typedefs,               only: GFS_control_type
  use machine,                   only: kind_phys
  use physparam,                 only: isubclw, iovrlw
  use rrtmgp_lw,                 only: nrghice_lw => nrghice, ipsdlw0
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,           only: ty_cloud_optics
  use mo_optical_props,          only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_cloud_sampling,         only: sampled_mask_max_ran, sampled_mask_exp_ran, draw_samples
  use mo_gas_concentrations,     only: ty_gas_concs
  use mersenne_twister,          only: random_setseed, random_number, random_stat
  use mo_rrtmgp_lw_cloud_optics, only: rrtmgp_lw_cloud_optics
  public GFS_rrtmgp_lw_run,GFS_rrtmgp_lw_init,GFS_rrtmgp_lw_finalize
  
contains

  subroutine GFS_rrtmgp_lw_init()
  end subroutine GFS_rrtmgp_lw_init
  ! #########################################################################################
  ! #########################################################################################
!! \section arg_table_GFS_rrtmgp_lw_run Argument Table
!! | local_name            | standard_name                                       | long_name                                                                    | units   | rank | type                  |    kind   | intent | optional |
!! |-----------------------|-----------------------------------------------------|------------------------------------------------------------------------------|---------|------|-----------------------|-----------|--------|----------|
!! | Model                 | GFS_control_type_instance                           | Fortran DDT containing FV3-GFS model control parameters                      | DDT     |    0 | GFS_control_type      |           | in     | F        |
!! | ncol                  | horizontal_loop_extent                              | horizontal dimension                                                         | count   |    0 | integer               |           | in     | F        |
!! | p_lay                 | air_pressure_at_layer_for_RRTMGP_in_hPa             | air pressure layer                                                           | hPa     |    2 | real                  | kind_phys | in     | F        |
!! | t_lay                 | air_temperature_at_layer_for_RRTMGP                 | air temperature layer                                                        | K       |    2 | real                  | kind_phys | in     | F        |
!! | p_lev                 | air_pressure_at_interface_for_RRTMGP_in_hPa         | air pressure level                                                           | hPa     |    2 | real                  | kind_phys | in     | F        |
!! | cld_frac              | total_cloud_fraction                                | layer total cloud fraction                                                   | frac    |    2 | real                  | kind_phys | in     | F        |
!! | cld_lwp               | cloud_liquid_water_path                             | layer cloud liquid water path                                                | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_reliq             | mean_effective_radius_for_liquid_cloud              | mean effective radius for liquid cloud                                       | micron  |    2 | real                  | kind_phys | in     | F        |
!! | cld_iwp               | cloud_ice_water_path                                | layer cloud ice water path                                                   | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_reice             | mean_effective_radius_for_ice_cloud                 | mean effective radius for ice cloud                                          | micron  |    2 | real                  | kind_phys | in     | F        |
!! | cld_swp               | cloud_snow_water_path                               | layer cloud snow water path                                                  | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_resnow            | mean_effective_radius_for_snow_flake                | mean effective radius for snow cloud                                         | micron  |    2 | real                  | kind_phys | in     | F        |
!! | cld_rwp               | cloud_rain_water_path                               | layer cloud rain water path                                                  | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_rerain            | mean_effective_radius_for_rain_drop                 | mean effective radius for rain cloud                                         | micron  |    2 | real                  | kind_phys | in     | F        |
!! | gas_concentrations    | Gas_concentrations_for_RRTMGP_suite                 | DDT containing gas concentrations for RRTMGP radiation scheme                | DDT     |    0 | ty_gas_concs          |           | in     | F        |
!! | icseed_lw             | seed_random_numbers_sw                              | seed for random number generation for shortwave radiation                    | none    |    1 | integer               |           | in     | F        |
!! | kdist_lw              | K_distribution_file_for_RRTMGP_LW_scheme            | DDT containing spectral information for RRTMGP LW radiation scheme           | DDT     |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | aerosols              | aerosol_optical_properties_for_longwave_bands_01-16 | aerosol optical properties for longwave bands 01-16                          | various |    4 | real                  | kind_phys | in     | F        |
!! | kdist_cldy_lw         | K_distribution_file_for_cloudy_RRTMGP_LW_scheme     | DDT containing spectral information for cloudy RRTMGP LW radiation scheme    | DDT     |    0 | ty_cloud_optics       |           | in     | F        |
!! | optical_props_clouds  | longwave_optical_properties_for_cloudy_atmosphere   | Fortran DDT containing RRTMGP optical properties                             | DDT     |    0 | ty_optical_props_1scl |           | out    | F        |
!! | optical_props_aerosol | longwave_optical_properties_for_aerosols            | Fortran DDT containing RRTMGP optical properties                             | DDT     |    0 | ty_optical_props_1scl |           | out    | F        |
!! | cldtaulw              | cloud_optical_depth_layers_at_10mu_band             | approx 10mu band layer cloud optical depth                                   | none    |    2 | real                  | kind_phys | out    | F        |
!! | errmsg                | ccpp_error_message                                  | error message for error handling in CCPP                                     | none    |    0 | character             | len=*     | out    | F        |
!! | errflg                | ccpp_error_flag                                     | error flag for error handling in CCPP                                        | flag    |    0 | integer               |           | out    | F        |
!!
  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_lw_run(Model, ncol, icseed_lw, p_lay, t_lay, p_lev, cld_frac,       &
       cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp, cld_rerain,    &
       gas_concentrations, kdist_lw, aerosols, kdist_cldy_lw,                               &
       optical_props_clouds, optical_props_aerosol, cldtaulw, errmsg, errflg)
    
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model
    integer, intent(in) :: &
         ncol                ! Number of horizontal gridpoints
    integer,intent(in),dimension(ncol) :: &
         icseed_lw           ! auxiliary special cloud related array when module 
                             ! variable isubclw=2, it provides permutation seed 
                             ! for each column profile that are used for generating 
                             ! random numbers. when isubclw /=2, it will not be used.
    real(kind_phys), dimension(ncol,model%levs), intent(in) :: &
         p_lay,            & ! Pressure @ model layer-centers         (hPa)
         t_lay               ! Temperature                            (K)
    real(kind_phys), dimension(ncol,model%levs+1), intent(in) :: &
         p_lev               ! Pressure @ model layer-interfaces      (hPa)
    real(kind_phys), dimension(ncol,model%levs),intent(in) :: &
         cld_frac,         & ! Total cloud fraction by layer
         cld_lwp,          & ! Cloud liquid water path
         cld_reliq,        & ! Cloud liquid effective radius
         cld_iwp,          & ! Cloud ice water path
         cld_reice,        & ! Cloud ice effective radius
         cld_swp,          & ! Cloud snow water path
         cld_resnow,       & ! Cloud snow effective radius
         cld_rwp,          & ! Cloud rain water path
         cld_rerain          ! Cloud rain effective radius
    type(ty_gas_concs),intent(in) :: &
         gas_concentrations  !
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         kdist_lw            ! RRTMGP DDT containing spectral information for LW calculation
    type(ty_cloud_optics),intent(in) :: &
         kdist_cldy_lw       !
    real(kind_phys), intent(in),dimension(ncol, model%levs, kdist_lw%get_nband(),3) :: &
         aerosols            !
    real(kind_phys), dimension(ncol,Model%levs), intent(out) :: &
         cldtaulw            ! approx 10.mu band layer cloud optical depth  

    ! Outputs
    type(ty_optical_props_1scl),intent(out) :: &
         optical_props_clouds, &
         optical_props_aerosol
    integer, intent(out) :: errflg
    character(len=*), intent(out) :: errmsg

    ! Local variables
    integer :: iCol
    integer,dimension(ncol) :: ipseed_lw
    logical,dimension(ncol,model%levs) :: liqmask, icemask
    type(ty_optical_props_1scl) :: optical_props_cloudsByBand
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(kdist_lw%get_ngpt(),model%levs,ncol) :: rng3D
    real(kind_phys), dimension(kdist_lw%get_ngpt()*model%levs) :: rng1D
    logical, dimension(ncol,model%levs,kdist_lw%get_ngpt()) :: cldfracMCICA
    real(kind_phys), dimension(ncol,model%levs,kdist_lw%get_nband()) :: &
         tau_cld

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. Model%lslwr) return

    ! #######################################################################################
    ! Change random number seed value for each radiation invocation (isubclw =1 or 2).
    ! #######################################################################################
    if(isubclw == 1) then      ! advance prescribed permutation seed
       do iCol = 1, nCol
          ipseed_lw(iCol) = ipsdlw0 + iCol
       enddo
    elseif (isubclw == 2) then ! use input array of permutaion seeds
       do iCol = 1, nCol
          ipseed_lw(iCol) = icseed_lw(iCol)
       enddo
    endif
    
    ! #######################################################################################
    ! Compute ice/liquid cloud masks, needed by rrtmgp_cloud_optics
    ! #######################################################################################
    liqmask = (cld_frac .gt. 0 .and. cld_lwp .gt. 0)
    icemask = (cld_frac .gt. 0 .and. cld_iwp .gt. 0)

    ! #######################################################################################
    ! Allocate space for RRTMGP DDTs containing cloud and aerosol radiative properties
    ! #######################################################################################
    ! Cloud optics [nCol,model%levs,nBands]
    call check_error_msg('GFS_rrtmgp_lw_run',optical_props_cloudsByBand%alloc_1scl(ncol, model%levs, kdist_lw%get_band_lims_wavenumber()))
    ! Aerosol optics [Ccol,model%levs,nBands]
    call check_error_msg('GFS_rrtmgp_lw_run',optical_props_aerosol%alloc_1scl(ncol, model%levs, kdist_lw%get_band_lims_wavenumber()))
    ! Cloud optics [nCol,model%levs,nGpts]
    call check_error_msg('GFS_rrtmgp_lw_run',optical_props_clouds%alloc_1scl(ncol, model%levs, kdist_lw))

    ! #######################################################################################
    ! Copy aerosol optical information to RRTMGP DDT
    ! #######################################################################################
    optical_props_aerosol%tau = aerosols(:,:,:,1) * (1. - aerosols(:,:,:,2))

    ! #######################################################################################
    ! Compute cloud-optics for RTE.
    ! #######################################################################################
    if (Model%rrtmgp_cld_phys .gt. 0) then
       ! i) RRTMGP cloud-optics.
       call check_error_msg('GFS_rrtmgp_lw_run',kdist_cldy_lw%cloud_optics(&
            ncol,                       & ! IN  - Number of horizontal gridpoints 
            model%levs,                 & ! IN  - Number of vertical layers
            kdist_lw%get_nband(),       & ! IN  - Number of LW bands
            nrghice_lw,                 & ! IN  - Number of ice-roughness categories
            liqmask,                    & ! IN  - Liquid-cloud mask
            icemask,                    & ! IN  - Ice-cloud mask
            cld_lwp,                    & ! IN  - Cloud liquid water path
            cld_iwp,                    & ! IN  - Cloud ice water path
            cld_reliq,                  & ! IN  - Cloud liquid effective radius
            cld_reice,                  & ! IN  - Cloud ice effective radius
            optical_props_cloudsByBand))  ! OUT - RRTMGP DDT containing cloud radiative properties
                                          !       in each band
    else
       ! ii) RRTMG cloud-optics.
       if (any(cld_frac .gt. 0)) then
          call rrtmgp_lw_cloud_optics(ncol, model%levs, kdist_lw%get_nband(), cld_lwp,     &
               cld_reliq, cld_iwp, cld_reice, cld_rwp, cld_rerain, cld_swp, cld_resnow,    &
               cld_frac, tau_cld)
          optical_props_cloudsByBand%tau = tau_cld
       endif
    endif

    ! #######################################################################################
    ! Call McICA to generate subcolumns.
    ! #######################################################################################
    ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
    ! and layers. ([nGpts,model%levs,nColumn]-> [nGpts*model%levs]*nColumn)
    do iCol=1,ncol
       call random_setseed(ipseed_lw(icol),rng_stat)
       call random_number(rng1D,rng_stat)
       rng3D(:,:,iCol) = reshape(source = rng1D,shape=[kdist_lw%get_ngpt(),model%levs])
    enddo
    
    ! Call McICA
    select case ( iovrlw )
       ! Maximumn-random 
    case(1)
       call check_error_msg('GFS_rrtmgp_lw_run',sampled_mask_max_ran(rng3D,cld_frac,cldfracMCICA))       
    end select
    
    ! Map band optical depth to each g-point using McICA
    call check_error_msg('GFS_rrtmgp_lw_run',draw_samples(cldfracMCICA,optical_props_cloudsByBand,optical_props_clouds))
    
    ! GFS_RRTMGP_POST_RUN() requires the LW optical depth ~10microns
    cldtaulw = optical_props_cloudsByBand%tau(:,:,7)

  end subroutine GFS_rrtmgp_lw_run
  
  subroutine GFS_rrtmgp_lw_finalize()
  end subroutine GFS_rrtmgp_lw_finalize

  subroutine check_error_msg(routine_name, error_msg)
    character(len=*), intent(in) :: &
         error_msg, routine_name
    
    if(error_msg /= "") then
       print*,"ERROR("//trim(routine_name)//"): "
       print*,trim(error_msg)
       return
    end if
  end subroutine check_error_msg    
end module GFS_rrtmgp_lw
