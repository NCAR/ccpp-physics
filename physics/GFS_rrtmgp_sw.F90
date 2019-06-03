module GFS_rrtmgp_sw
 use GFS_typedefs,               only: GFS_control_type
  use machine,                   only: kind_phys
  use physparam,                 only: isubcsw, iovrsw
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,           only: ty_cloud_optics
  use mo_optical_props,          only: ty_optical_props_2str
  use mo_cloud_sampling,         only: sampled_mask_max_ran, sampled_mask_exp_ran, draw_samples
  use mersenne_twister,          only: random_setseed, random_number, random_stat
  use mo_rrtmgp_sw_cloud_optics, only: rrtmgp_sw_cloud_optics
  use rrtmgp_aux,                only: check_error_msg, nrghice_sw, ipsdsw0

  public GFS_rrtmgp_sw_run,GFS_rrtmgp_sw_init,GFS_rrtmgp_sw_finalize
  
contains

  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_init()
  end subroutine GFS_rrtmgp_sw_init

!! \section arg_table_GFS_rrtmgp_sw_run Argument Table
!! | local_name            | standard_name                                        | long_name                                                                    | units   | rank | type                  |    kind   | intent | optional |
!! |-----------------------|------------------------------------------------------|------------------------------------------------------------------------------|---------|------|-----------------------|-----------|--------|----------|
!! | Model                 | GFS_control_type_instance                            | Fortran DDT containing FV3-GFS model control parameters                      | DDT     |    0 | GFS_control_type      |           | in     | F        |
!! | ncol                  | horizontal_loop_extent                               | horizontal dimension                                                         | count   |    0 | integer               |           | in     | F        |
!! | p_lay                 | air_pressure_at_layer_for_RRTMGP_in_hPa              | air pressure layer                                                           | hPa     |    2 | real                  | kind_phys | in     | F        |
!! | t_lay                 | air_temperature_at_layer_for_RRTMGP                  | air temperature layer                                                        | K       |    2 | real                  | kind_phys | in     | F        |
!! | p_lev                 | air_pressure_at_interface_for_RRTMGP_in_hPa          | air pressure level                                                           | hPa     |    2 | real                  | kind_phys | in     | F        |
!! | cld_frac              | total_cloud_fraction                                 | layer total cloud fraction                                                   | frac    |    2 | real                  | kind_phys | in     | F        |
!! | cld_lwp               | cloud_liquid_water_path                              | layer cloud liquid water path                                                | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_reliq             | mean_effective_radius_for_liquid_cloud               | mean effective radius for liquid cloud                                       | micron  |    2 | real                  | kind_phys | in     | F        |
!! | cld_iwp               | cloud_ice_water_path                                 | layer cloud ice water path                                                   | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_reice             | mean_effective_radius_for_ice_cloud                  | mean effective radius for ice cloud                                          | micron  |    2 | real                  | kind_phys | in     | F        |
!! | cld_swp               | cloud_snow_water_path                                | layer cloud snow water path                                                  | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_resnow            | mean_effective_radius_for_snow_flake                 | mean effective radius for snow cloud                                         | micron  |    2 | real                  | kind_phys | in     | F        |
!! | cld_rwp               | cloud_rain_water_path                                | layer cloud rain water path                                                  | g m-2   |    2 | real                  | kind_phys | in     | F        |
!! | cld_rerain            | mean_effective_radius_for_rain_drop                  | mean effective radius for rain cloud                                         | micron  |    2 | real                  | kind_phys | in     | F        |
!! | icseed_sw             | seed_random_numbers_sw                               | seed for random number generation for shortwave radiation                    | none    |    1 | integer               |           | in     | F        |
!! | sw_gas_props          | coefficients_for_sw_gas_optics                       | DDT containing spectral information for RRTMGP SW radiation scheme           | DDT     |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | aerosols              | aerosol_optical_properties_for_shortwave_bands_01-16 | aerosol optical properties for shortwave bands 01-16                         | various |    4 | real                  | kind_phys | in     | F        |
!! | sw_cloud_props        | coefficients_for_sw_cloud_optics                     | DDT containing spectral information for cloudy RRTMGP SW radiation scheme    | DDT     |    0 | ty_cloud_optics       |           | in     | F        |
!! | nday                  | daytime_points_dimension                             | daytime points dimension                                                     | count   |    0 | integer               |           | in     | F        |
!! | idxday                | daytime_points                                       | daytime points                                                               | index   |    1 | integer               |           | in     | F        |
!! | optical_props_clouds  | shortwave_optical_properties_for_cloudy_atmosphere   | Fortran DDT containing RRTMGP optical properties                             | DDT     |    0 | ty_optical_props_2str |           | out    | F        |
!! | optical_props_aerosol | shortwave_optical_properties_for_aerosols            | Fortran DDT containing RRTMGP optical properties                             | DDT     |    0 | ty_optical_props_2str |           | out    | F        |
!! | cldtausw              | cloud_optical_depth_layers_at_0.55mu_band            | approx .55mu band layer cloud optical depth                                  | none    |    2 | real                  | kind_phys | out    | F        |
!! | errmsg                | ccpp_error_message                                   | error message for error handling in CCPP                                     | none    |    0 | character             | len=*     | out    | F        |
!! | errflg                | ccpp_error_flag                                      | error flag for error handling in CCPP                                        | flag    |    0 | integer               |           | out    | F        |
!!
  ! #########################################################################################
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_run(Model, ncol, icseed_sw, p_lay, t_lay, p_lev, cld_frac,       & ! IN
       cld_lwp, cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp, cld_rerain,    & ! IN
       sw_gas_props, aerosols, sw_cloud_props, nday, idxday,                                     & ! IN
       optical_props_clouds, optical_props_aerosol, cldtausw, errmsg, errflg)                 ! OUT
    
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model
    integer, intent(in) :: &
         ncol,             & ! Number of horizontal gridpoints
         nday                ! Number of daylit points.
    integer,intent(in),dimension(nday) :: &
         idxday              ! Indices for daylit points.
    integer,intent(in),dimension(ncol) :: &
         icseed_sw           ! auxiliary special cloud related array when module 
                             ! variable isubcsw=2, it provides permutation seed 
                             ! for each column profile that are used for generating 
                             ! random numbers. when isubcsw /=2, it will not be used.
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
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props            ! RRTMGP DDT containing spectral information for SW calculation
    type(ty_cloud_optics),intent(in) :: &
         sw_cloud_props       ! 
    real(kind_phys), intent(in),dimension(ncol, model%levs, sw_gas_props%get_nband(),3) :: &
         aerosols            !

    ! Outputs
    type(ty_optical_props_2str),intent(out) :: &
         optical_props_clouds, &
         optical_props_aerosol
    real(kind_phys), dimension(ncol,Model%levs), intent(out) :: &
         cldtausw            ! approx 10.mu band layer cloud optical depth  
    integer, intent(out) :: errflg
    character(len=*), intent(out) :: errmsg

    ! Local variables
    integer :: iCol
    integer,dimension(ncol) :: ipseed_sw
    logical,dimension(ncol,model%levs) :: liqmask, icemask
    type(ty_optical_props_2str) :: optical_props_cloudsByBand
    type(random_stat) :: rng_stat
    real(kind_phys), dimension(sw_gas_props%get_ngpt(),model%levs,ncol) :: rng3D
    real(kind_phys), dimension(sw_gas_props%get_ngpt()*model%levs) :: rng1D
    logical, dimension(ncol,model%levs,sw_gas_props%get_ngpt()) :: cldfracMCICA
    real(kind_phys), dimension(nday,model%levs,sw_gas_props%get_nband()) :: &
         tau_cld, ssa_cld, asy_cld

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. Model%lsswr) return

    ! #######################################################################################
    ! Change random number seed value for each radiation invocation (isubcsw =1 or 2).
    ! #######################################################################################
    if(isubcsw == 1) then      ! advance prescribed permutation seed
       do iCol = 1, ncol
          ipseed_sw(iCol) = ipsdsw0 + iCol
       enddo
    elseif (isubcsw == 2) then ! use input array of permutaion seeds
       do iCol = 1, ncol
          ipseed_sw(iCol) = icseed_sw(iCol)
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
    ! Cloud optics [ncol,model%levs,nBands]
    call check_error_msg('GFS_rrtmgp_sw_run',optical_props_cloudsByBand%alloc_2str(ncol, model%levs, sw_gas_props%get_band_lims_wavenumber()))
    ! Aerosol optics [ncol,model%levs,nBands]
    call check_error_msg('GFS_rrtmgp_sw_run',optical_props_aerosol%alloc_2str(ncol, model%levs, sw_gas_props%get_band_lims_wavenumber()))
    ! Cloud optics [ncol,model%levs,nGpts]
    call check_error_msg('GFS_rrtmgp_sw_run',optical_props_clouds%alloc_2str(ncol, model%levs, sw_gas_props))

    ! #######################################################################################
    ! Copy aerosol optical information to RRTMGP DDT
    ! #######################################################################################
    optical_props_aerosol%tau = aerosols(:,:,:,1)
    optical_props_aerosol%ssa = aerosols(:,:,:,2)
    optical_props_aerosol%g   = aerosols(:,:,:,3)

    ! #######################################################################################
    ! Compute cloud-optics for RTE.
    ! #######################################################################################
    if (Model%rrtmgp_cld_optics .gt. 0) then
       ! RRTMGP cloud-optics.
       call check_error_msg('GFS_rrtmgp_sw_run',sw_cloud_props%cloud_optics(&
            ncol,                       & ! IN  - Number of daylit gridpoints
            model%levs,                 & ! IN  - Number of vertical layers
            sw_gas_props%get_nband(),       & ! IN  - Number of SW bands
            nrghice_sw,                 & ! IN  - Number of ice-roughness categories
            liqmask,                    & ! IN  - Liquid-cloud mask
            icemask,                    & ! IN  - Ice-cloud mask
            cld_lwp,                    & ! IN  - Cloud liquid water path
            cld_iwp,                    & ! IN  - Cloud ice water path
            cld_reliq,                  & ! IN  - Cloud liquid effective radius
            cld_reice,                  & ! IN  - Cloud ice effective radius
            optical_props_cloudsByBand))  ! OUT - RRTMGP DDT containing cloud radiative properties
                                          !       in each band
    else
       ! RRTMG cloud-optics
       if (any(cld_frac .gt. 0)) then
          optical_props_cloudsByBand%tau(:,:,:) = 0._kind_phys
          optical_props_cloudsByBand%ssa(:,:,:) = 0._kind_phys
          optical_props_cloudsByBand%g(:,:,:)   = 0._kind_phys
          call rrtmgp_sw_cloud_optics(nday, model%levs, sw_gas_props%get_nband(), cld_lwp(idxday,:), &
               cld_reliq(idxday,:), cld_iwp(idxday,:), cld_reice(idxday,:), cld_rwp(idxday,:),   &
               cld_rerain(idxday,:), cld_swp(idxday,:), cld_resnow(idxday,:), cld_frac(idxday,:),&
               tau_cld, ssa_cld, asy_cld)
          optical_props_cloudsByBand%tau(idxday,:,:) = tau_cld
          optical_props_cloudsByBand%ssa(idxday,:,:) = ssa_cld
          optical_props_cloudsByBand%g(idxday,:,:)   = asy_cld
       endif
    endif
    ! #######################################################################################
    ! Call McICA to generate subcolumns.
    ! #######################################################################################
    ! Call RNG. Mersennse Twister accepts 1D array, so loop over columns and collapse along G-points 
    ! and layers. ([nGpts,model%levs,nColumn]-> [nGpts*model%levs]*nColumn)
    do iCol=1,ncol
       call random_setseed(ipseed_sw(icol),rng_stat)
       call random_number(rng1D,rng_stat)
       rng3D(:,:,iCol) = reshape(source = rng1D,shape=[sw_gas_props%get_ngpt(),model%levs])
    enddo
    
    ! Call McICA
    select case ( iovrsw )
       ! Maximumn-random 
    case(1)
       call check_error_msg('GFS_rrtmgp_sw_run',sampled_mask_max_ran(rng3D,cld_frac,cldfracMCICA))       
    end select
    
    ! Map band optical depth to each g-point using McICA
    call check_error_msg('GFS_rrtmgp_sw_run',draw_samples(cldfracMCICA,optical_props_cloudsByBand,optical_props_clouds))

    ! GFS_RRTMGP_POST_RUN() requires the SW optical depth ~0.55microns
    cldtausw = optical_props_cloudsByBand%tau(:,:,11)    

  end subroutine GFS_rrtmgp_sw_run
  
  subroutine GFS_rrtmgp_sw_finalize()
  end subroutine GFS_rrtmgp_sw_finalize
 
end module GFS_rrtmgp_sw
