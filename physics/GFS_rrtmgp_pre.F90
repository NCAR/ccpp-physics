!> \file GFS_rrtmgp_pre.f90
!! This file contains
      module GFS_rrtmgp_pre

      public GFS_rrtmgp_pre_run

      contains

!> \defgroup GFS_rrtmgp_pre GFS RRTMGP Scheme Pre
!! @{
!! \section arg_table_GFS_rrtmgp_pre_init Argument Table
!!
      subroutine GFS_rrtmgp_pre_init ()
        open(58,file='GFS_rrtmgp_aux_dump.txt',status='unknown')
      end subroutine GFS_rrtmgp_pre_init

!> \section arg_table_GFS_rrtmgp_pre_run Argument Table
!! | local_name        | standard_name                                                 | long_name                                                                     | units    | rank |  type            |   kind    | intent | optional |
!! |-------------------|---------------------------------------------------------------|-------------------------------------------------------------------------------|----------|------|------------------|-----------|--------|----------|
!! | Model             | GFS_control_type_instance                                     | Fortran DDT containing FV3-GFS model control parameters                       | DDT      |    0 | GFS_control_type |           | in     | F        |
!! | Grid              | GFS_grid_type_instance                                        | Fortran DDT containing FV3-GFS grid and interpolation related data            | DDT      |    0 | GFS_grid_type    |           | in     | F        |
!! | Sfcprop           | GFS_sfcprop_type_instance                                     | Fortran DDT containing FV3-GFS surface fields                                 | DDT      |    0 | GFS_sfcprop_type |           | in     | F        |
!! | Statein           | GFS_statein_type_instance                                     | Fortran DDT containing FV3-GFS prognostic state data in from dycore           | DDT      |    0 | GFS_statein_type |           | in     | F        |
!! | Tbd               | GFS_tbd_type_instance                                         | Fortran DDT containing FV3-GFS data not yet assigned to a defined container   | DDT      |    0 | GFS_tbd_type     |           | in     | F        |
!! | Cldprop           | GFS_cldprop_type_instance                                     | Fortran DDT containing FV3-GFS cloud fields needed by radiation from physics  | DDT      |    0 | GFS_cldprop_type |           | in     | F        |
!! | Coupling          | GFS_coupling_type_instance                                    | Fortran DDT containing FV3-GFS fields needed for coupling                     | DDT      |    0 | GFS_coupling_type|           | in     | F        |
!! | Radtend           | GFS_radtend_type_instance                                     | Fortran DDT containing FV3-GFS radiation tendencies                           | DDT      |    0 | GFS_radtend_type |           | inout  | F        |
!! | lm                | vertical_layer_dimension_for_radiation                        | number of vertical layers for radiation calculation                           | count    |    0 | integer          |           | in     | F        |
!! | im                | horizontal_loop_extent                                        | horizontal loop extent                                                        | count    |    0 | integer          |           | in     | F        |
!! | lmk               | adjusted_vertical_layer_dimension_for_radiation               | number of vertical layers for radiation                                       | count    |    0 | integer          |           | in     | F        |
!! | lmp               | adjusted_vertical_level_dimension_for_radiation               | number of vertical levels for radiation                                       | count    |    0 | integer          |           | in     | F        |
!! | kd                | vertical_index_difference_between_inout_and_local             | vertical index difference between in/out and local                            | index    |    0 | integer          |           | out    | F        |
!! | kt                | vertical_index_difference_between_layer_and_upper_bound       | vertical index difference between layer and upper bound                       | index    |    0 | integer          |           | out    | F        |
!! | kb                | vertical_index_difference_between_layer_and_lower_bound       | vertical index difference between layer and lower bound                       | index    |    0 | integer          |           | out    | F        |
!! | raddt             | time_step_for_radiation                                       | radiation time step                                                           | s        |    0 | real             | kind_phys | out    | F        |
!! | delp              | layer_pressure_thickness_for_radiation                        | layer pressure thickness on radiation levels                                  | hPa      |    2 | real             | kind_phys | out    | F        |
!! | dz                | layer_thickness_for_radiation                                 | layer thickness on radiation levels                                           | km       |    2 | real             | kind_phys | out    | F        |
!! | plvl              | air_pressure_at_interface_for_radiation_in_hPa                | air pressure at vertical interface for radiation calculation                  | hPa      |    2 | real             | kind_phys | out    | F        |
!! | plyr              | air_pressure_at_layer_for_radiation_in_hPa                    | air pressure at vertical layer for radiation calculation                      | hPa      |    2 | real             | kind_phys | out    | F        |
!! | tlvl              | air_temperature_at_interface_for_radiation                    | air temperature at vertical interface for radiation calculation               | K        |    2 | real             | kind_phys | out    | F        |
!! | tlyr              | air_temperature_at_layer_for_radiation                        | air temperature at vertical layer for radiation calculation                   | K        |    2 | real             | kind_phys | out    | F        |
!! | tsfg              | surface_ground_temperature_for_radiation                      | surface ground temperature for radiation                                      | K        |    1 | real             | kind_phys | out    | F        |
!! | tsfa              | surface_air_temperature_for_radiation                         | lowest model layer air temperature for radiation                              | K        |    1 | real             | kind_phys | out    | F        |
!! | qlyr              | water_vapor_specific_humidity_at_layer_for_radiation          | water vapor specific humidity at vertical layer for radiation calculation     | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | olyr              | ozone_concentration_at_layer_for_radiation                    | ozone concentration                                                           | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_co2        | volume_mixing_ratio_co2                                       | CO2 volume mixing ratio                                                       | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_n2o        | volume_mixing_ratio_n2o                                       | N2O volume mixing ratio                                                       | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_ch4        | volume_mixing_ratio_ch4                                       | CH4 volume mixing ratio                                                       | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_o2         | volume_mixing_ratio_o2                                        | O2 volume mixing ratio                                                        | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_co         | volume_mixing_ratio_co                                        | CO volume mixing ratio                                                        | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_cfc11      | volume_mixing_ratio_cfc11                                     | CFC11 volume mixing ratio                                                     | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_cfc12      | volume_mixing_ratio_cfc12                                     | CFC12 volume mixing ratio                                                     | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_cfc22      | volume_mixing_ratio_cfc22                                     | CFC22 volume mixing ratio                                                     | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_ccl4       | volume_mixing_ratio_ccl4                                      | CCL4 volume mixing ratio                                                      | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_cfc113     | volume_mixing_ratio_cfc113                                    | CFC113 volume mixing ratio                                                    | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | faersw1           | aerosol_optical_depth_for_shortwave_bands_01-16               | aerosol optical depth for shortwave bands 01-16                               | none     |    3 | real             | kind_phys | out    | F        |
!! | faersw2           | aerosol_single_scattering_albedo_for_shortwave_bands_01-16    | aerosol single scattering albedo for shortwave bands 01-16                    | frac     |    3 | real             | kind_phys | out    | F        |
!! | faersw3           | aerosol_asymmetry_parameter_for_shortwave_bands_01-16         | aerosol asymmetry parameter for shortwave bands 01-16                         | none     |    3 | real             | kind_phys | out    | F        |
!! | faerlw1           | aerosol_optical_depth_for_longwave_bands_01-16                | aerosol optical depth for longwave bands 01-16                                | none     |    3 | real             | kind_phys | out    | F        |
!! | faerlw2           | aerosol_single_scattering_albedo_for_longwave_bands_01-16     | aerosol single scattering albedo for longwave bands 01-16                     | frac     |    3 | real             | kind_phys | out    | F        |
!! | faerlw3           | aerosol_asymmetry_parameter_for_longwave_bands_01-16          | aerosol asymmetry parameter for longwave bands 01-16                          | none     |    3 | real             | kind_phys | out    | F        |
!! | aerodp            | atmosphere_optical_thickness_due_to_ambient_aerosol_particles | vertical integrated optical depth for various aerosol species                 | none     |    2 | real             | kind_phys | out    | F        |
!! | clouds1           | total_cloud_fraction                                          | layer total cloud fraction                                                    | frac     |    2 | real             | kind_phys | out    | F        |
!! | clouds2           | cloud_liquid_water_path                                       | layer cloud liquid water path                                                 | g m-2    |    2 | real             | kind_phys | out    | F        |
!! | clouds3           | mean_effective_radius_for_liquid_cloud                        | mean effective radius for liquid cloud                                        | micron   |    2 | real             | kind_phys | out    | F        |
!! | clouds4           | cloud_ice_water_path                                          | layer cloud ice water path                                                    | g m-2    |    2 | real             | kind_phys | out    | F        |
!! | clouds5           | mean_effective_radius_for_ice_cloud                           | mean effective radius for ice cloud                                           | micron   |    2 | real             | kind_phys | out    | F        |
!! | clouds6           | cloud_rain_water_path                                         | cloud rain water path                                                         | g m-2    |    2 | real             | kind_phys | out    | F        |
!! | clouds7           | mean_effective_radius_for_rain_drop                           | mean effective radius for rain drop                                           | micron   |    2 | real             | kind_phys | out    | F        |
!! | clouds8           | cloud_snow_water_path                                         | cloud snow water path                                                         | g m-2    |    2 | real             | kind_phys | out    | F        |
!! | clouds9           | mean_effective_radius_for_snow_flake                          | mean effective radius for snow flake                                          | micron   |    2 | real             | kind_phys | out    | F        |
!! | cldsa             | cloud_area_fraction_for_radiation                             | fraction of clouds for low, middle,high, total and BL                         | frac     |    2 | real             | kind_phys | out    | F        |
!! | mtopa             | model_layer_number_at_cloud_top                               | vertical indices for low, middle and high cloud tops                          | index    |    2 | integer          |           | out    | F        |
!! | mbota             | model_layer_number_at_cloud_base                              | vertical indices for low, middle and high cloud bases                         | index    |    2 | integer          |           | out    | F        |
!! | de_lgth           | cloud_decorrelation_length                                    | cloud decorrelation length                                                    | km       |    1 | real             | kind_phys | out   | F        |
!! | alb1d             | surface_albedo_perturbation                                   | surface albedo perturbation                                                   | frac     |    1 | real             | kind_phys | out    | F        |
!! | errmsg            | ccpp_error_message                                            | error message for error handling in CCPP                                      | none     |    0 | character        | len=*     | out    | F        |
!! | errflg            | ccpp_error_flag                                               | error flag for error handling in CCPP                                         | flag     |    0 | integer          |           | out    | F        |
!! | kdist_lw          | K_distribution_file_for_RRTMGP_LW_scheme                      | DDT containing spectral information for RRTMGP LW radiation scheme            | DDT      |    0 | ty_gas_optics_rrtmgp_type |           | in    | F        |
!!
      ! Attention - the output arguments lm, im, lmk, lmp must not be set
      ! in the CCPP version - they are defined in the interstitial_create routine
      ! #########################################################################################
      subroutine GFS_rrtmgp_pre_run (Model, Grid, Sfcprop, Statein, Tbd, Cldprop, Coupling,     & ! IN
           Radtend,                                                                             & ! INOUT
           lm, im, lmk, lmp, kdist_lw,                                                          & ! IN
           kd, kt, kb, raddt, delp, dz, plvl, plyr, tlvl, tlyr, tsfg, tsfa, qlyr, olyr,         & ! OUT
           gasvmr_co2, gasvmr_n2o, gasvmr_ch4, gasvmr_o2, gasvmr_co, gasvmr_cfc11,              & ! OUT
           gasvmr_cfc12, gasvmr_cfc22, gasvmr_ccl4, gasvmr_cfc113, faersw1, faersw2,  faersw3,  & ! OUT
           faerlw1, faerlw2, faerlw3, aerodp, clouds1, clouds2, clouds3, clouds4, clouds5,      & ! OUT
           clouds6, clouds7, clouds8, clouds9, cldsa, mtopa, mbota, de_lgth, alb1d,             & ! OUT
           errmsg, errflg)

        use physparam
        use machine, only: &
             kind_phys                   ! Working type
        use GFS_typedefs, only:        &
             GFS_statein_type,         & ! Prognostic state data in from dycore
             GFS_stateout_type,        & ! Prognostic state or tendencies return to dycore
             GFS_sfcprop_type,         & ! Surface fields
             GFS_coupling_type,        & ! Fields to/from coupling with other components (e.g. land/ice/ocean/etc.)
             GFS_control_type,         & ! Model control parameters
             GFS_grid_type,            & ! Grid and interpolation related data
             GFS_tbd_type,             & ! To-Be-Determined data that doesn't fit in any one container
             GFS_cldprop_type,         & ! Cloud fields needed by radiation from physics
             GFS_radtend_type,         & ! Radiation tendencies needed in physics
             GFS_diag_type               ! Fields targetted for diagnostic output
        use physcons, only:            &
             eps   => con_eps,         & ! Rd/Rv
             epsm1 => con_epsm1,       & ! Rd/Rv-1
             fvirt => con_fvirt,       & ! Rv/Rd-1
             rog   => con_rog,         & ! Rd/g
             rocp  => con_rocp           ! Rd/cp
        use radcons, only: &
             itsfc,                    & ! Flag for LW sfc. temp.
             ltp,                      & ! 1-add extra-top layer; 0-no extra layer
             lextop,                   & ! ltp > 0
             qmin,qme5, qme6, epsq       ! Minimum vlaues for varius calculations
        use funcphys, only:            &
             fpvs                        ! Function ot compute sat. vapor pressure over liq.
        use module_radiation_astronomy,only: &
             coszmn                      ! Function to compute cos(SZA)
        use module_radiation_gases,    only: &
             NF_VGAS,                  & ! Number of active gas species
             getgases,                 & ! Routine to setup trace gases
             getozn                      ! Routine to setup ozone
        use module_radiation_aerosols, only: &
             NF_AESW,                  & ! Number of optical-fields in SW output (3=tau+g+omega)
             NF_AELW,                  & ! Number of optical-fields in LW output (3=tau+g+omega)
             setaer,                   & ! Routine to compute aerosol radiative properties (tau,g,omega)
             NSPC1                       ! Number of species for vertically integrated aerosol optical-depth
        use module_radiation_clouds, only: &
             NF_CLDS,                  & ! Number of fields in "clouds" array (e.g. (cloud(1)=lwp,clouds(2)=ReffLiq,...)
             progcld1,                 & ! Zhao/Moorthi's prognostic cloud scheme
             progcld3,                 & ! Zhao/Moorthi's prognostic cloud+pdfcld
             progcld4,                 & ! GFDL cloud scheme
             progcld5,                 & ! Thompson / WSM6 cloud micrphysics scheme
             progclduni                  ! Unified cloud-scheme
        use surface_perturbation, only: & 
             cdfnor                      ! Routine to compute CDF (used to compute percentiles)
        ! RRTMGP stuff
        use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp_type
        use GFS_rrtmgp_lw, only:       &
             nBandsLW!,                 & ! Number of LW bands in RRTMGP
             !kdist_lw                    ! DDT contining LW spectral information
        use GFS_rrtmgp_sw, only:       &
             nBandsSW,                 & ! Number of SW bands in RRTMGP
             kdist_sw                    ! DDT contining SW spectral information
        implicit none

        ! Inputs
        type(GFS_control_type),   intent(in)    :: Model
        type(GFS_grid_type),      intent(in)    :: Grid
        type(GFS_sfcprop_type),   intent(in)    :: Sfcprop
        type(GFS_statein_type),   intent(in)    :: Statein
        type(GFS_radtend_type),   intent(inout) :: Radtend
        type(GFS_tbd_type),       intent(in)    :: Tbd
        type(GFS_cldprop_type),   intent(in)    :: Cldprop
        type(GFS_coupling_type),  intent(in)    :: Coupling
        integer,                  intent(in)    :: im, lm, lmk, lmp
        type(ty_gas_optics_rrtmgp_type),intent(in)   :: kdist_lw

        ! Outputs
        integer,         intent(out) :: kd, kt, kb
        real(kind_phys), intent(out) :: raddt
        real(kind_phys), dimension(size(Grid%xlon,1)), intent(out) :: &
             tsfg, tsfa
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP), intent(out) ::  &
             delp, dz, plyr, tlyr, qlyr, olyr, gasvmr_co2, gasvmr_n2o, gasvmr_ch4,     &
             gasvmr_o2, gasvmr_co, gasvmr_cfc11, gasvmr_cfc12, gasvmr_cfc22,           &
             gasvmr_ccl4, gasvmr_cfc113, clouds1, clouds2, clouds3, clouds4, clouds5,  &
             clouds6, clouds7, clouds8, clouds9
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+1+LTP), intent(out) :: &
             plvl, tlvl
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,nBandsSW), intent(out) :: &
             faersw1, faersw2, faersw3
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,nBandsLW), intent(out) :: &
             faerlw1, faerlw2, faerlw3
        real(kind_phys), dimension(size(Grid%xlon,1),NSPC1), intent(out) :: &
             aerodp

        real(kind_phys), dimension(size(Grid%xlon,1),5), intent(out) :: cldsa
        integer,         dimension(size(Grid%xlon,1),3), intent(out) :: mbota,mtopa
        real(kind_phys), dimension(size(Grid%xlon,1)), intent(out) :: de_lgth,alb1d
        character(len=*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        ! Local variables
        integer :: me, nfxr, ntrac, ntcw, ntiw, ncld, ntrw, ntsw, ntgl,i, j, k, k1, k2, lsk, &
             lv, n, itop, ibtc, LP1, lla, llb, lya, lyb
        real(kind_phys) :: es, qs, delt, tem0d
        real(kind_phys), dimension(size(Grid%xlon,1)) :: cvt1, cvb1, tem1d, tskn
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP) :: htswc, htlwc,   &
             gcice, grain, grime, htsw0, htlw0, rhly, tvly,qstl, vvel, clw, ciw, prslk1,     &
             tem2da, cldcov, deltaq, cnvc, cnvw, effrl, effri, effrr, effrs
        real (kind_phys) :: clwmin, clwm, clwt, onemrh, value, tem1, tem2, tem3
        real (kind_phys), parameter :: xrc3 = 100.
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP+1) :: tem2db
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,Model%ncnd) :: ccnd
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,2:Model%ntrac) :: tracer1
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NF_CLDS) :: clouds
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NF_VGAS) :: gasvmr
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,nBandsSW,NF_AESW)::faersw
        real(kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,nBandsLW,NF_AELW)::faerlw

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (.not. (Model%lsswr .or. Model%lslwr)) return

        ! Define some commonly used integers
        me    = Model%me    ! MPI rank designator
        NFXR  = Model%nfxr  ! second dimension for fluxr diagnostic variable (radiation)
        NTRAC = Model%ntrac ! Number of tracers
        ntcw  = Model%ntcw  ! Tracer index for cloud condensate (or liquid water)
        ntiw  = Model%ntiw  ! Tracer index for ice 
        ncld  = Model%ncld  ! Cloud scheme
        ntrw  = Model%ntrw  ! Tracer index for rain 
        ntsw  = Model%ntsw  ! Tracer index for snow 
        ntgl  = Model%ntgl  ! Tracer index for groupel
        LP1   = LM + 1      ! num of in/out levels

        ! Set local /level/layer indexes corresponding to in/out variables
        if ( lextop ) then
           if ( ivflip == 1 ) then     ! vertical from sfc upward
              kd = 0                   ! index diff between in/out and local
              kt = 1                   ! index diff between lyr and upper bound
              kb = 0                   ! index diff between lyr and lower bound
              lla = LMK                ! local index at the 2nd level from top
              llb = LMP                ! local index at toa level
              lya = LM                 ! local index for the 2nd layer from top
              lyb = LP1                ! local index for the top layer
           else                        ! vertical from toa downward
              kd = 1                   ! index diff between in/out and local
              kt = 0                   ! index diff between lyr and upper bound
              kb = 1                   ! index diff between lyr and lower bound
              lla = 2                  ! local index at the 2nd level from top
              llb = 1                  ! local index at toa level
              lya = 2                  ! local index for the 2nd layer from top
              lyb = 1                  ! local index for the top layer
           endif                       ! end if_ivflip_block
        else
           kd = 0
           if ( ivflip == 1 ) then     ! vertical from sfc upward
              kt = 1                   ! index diff between lyr and upper bound
              kb = 0                   ! index diff between lyr and lower bound
           else                        ! vertical from toa downward
              kt = 0                   ! index diff between lyr and upper bound
              kb = 1                   ! index diff between lyr and lower bound
           endif                       ! end if_ivflip_block
        endif                          ! end if_lextop_block

        ! Radiation time step (output)
        raddt = min(Model%fhswr, Model%fhlwr)

        ! Setup surface ground temperature and ground/air skin temperature if required.
        if ( itsfc == 0 ) then            ! use same sfc skin-air/ground temp
           tskn(1:IM) = Sfcprop%tsfc(1:IM)
           tsfg(1:IM) = Sfcprop%tsfc(1:IM)
!           do i = 1, IM
!              tskn(i) = Sfcprop%tsfc(i)
!              tsfg(i) = Sfcprop%tsfc(i)
!           enddo
        else                              ! use diff sfc skin-air/ground temp
           tskn(1:IM) = Sfcprop%tsfc(1:IM)
           tsfg(1:IM) = Sfcprop%tsfc(1:IM)
!           do i = 1, IM
!              tskn(i) = Sfcprop%tsfc(i)
!              tsfg(i) = Sfcprop%tsfc(i)
!           enddo
        endif

        ! Prepare atmospheric profiles for radiation input.
        lsk = 0
        if (ivflip == 0 .and. lm < Model%levs) lsk = Model%levs - lm
      
        ! Copy over state fields into fields, compute some needed quantities.
        do k = 1, LM
           k1 = k + kd
           k2 = k + lsk
           do i = 1, IM
              plvl(i,k1+kb) = Statein%prsi(i,k2+kb)
              plyr(i,k1)    = Statein%prsl(i,k2)
              tlyr(i,k1)    = Statein%tgrs(i,k2)
              prslk1(i,k1)  = Statein%prslk(i,k2)
              
              ! Compute relative humidity.
              es  = min( Statein%prsl(i,k2),  fpvs( Statein%tgrs(i,k2) ) )  ! fpvs and prsl in pa
              qs  = max( QMIN, eps * es / (Statein%prsl(i,k2) + epsm1*es) )
              rhly(i,k1) = max( 0.0, min( 1.0, max(QMIN, Statein%qgrs(i,k2,1))/qs ) )
              qstl(i,k1) = qs
           enddo
           plvl(i,k1+kb) = kdist_lw%get_press_min()
        enddo
        
        ! Recast remaining all tracers (except sphum) forcing them all to be positive
        do j = 2, NTRAC
           do k = 1, LM
              k1 = k + kd
              k2 = k + lsk
              tracer1(:,k1,j) = max(0.0, Statein%qgrs(:,k2,j))
           enddo
        enddo
        
        ! Input data from toa to sfc
        if (ivflip == 0) then                                
           plvl(1:IM-1,1+kd) = Statein%prsi(1:IM-1,1)
           plvl(IM,1+kd)     = kdist_lw%get_press_min()
           if (lsk /= 0) then
              plvl(1:IM-1,1+kd) = 0.5 * (plvl(1:IM-1,2+kd) + plvl(1:IM-1,1+kd))
              plvl(IM,1+kd)     = kdist_lw%get_press_min()
           endif
           ! Input data from sfc to top
        else                                                
           plvl(1:IM-1,LP1+kd) = Statein%prsi(1:IM-1,LP1+lsk)
           plvl(IM,LP1+kd)     = kdist_lw%get_press_min()
           if (lsk /= 0) then
              plvl(1:IM-1,LM+kd) = 0.5 * (plvl(1:IM-1,LP1+kd) + plvl(1:IM-1,LM+kd))
              plvl(IM,LM+kd)     = kdist_lw%get_press_min()
           endif
        endif
        
        if ( lextop ) then                 ! values for extra top layer
           do i = 1, IM
              plvl(i,llb) =  kdist_lw%get_press_min()
              if ( plvl(i,lla) <=  kdist_lw%get_press_min() ) plvl(i,lla) = 2.0* kdist_lw%get_press_min()
              plyr(i,lyb)   = 0.5 * plvl(i,lla)
              tlyr(i,lyb)   = tlyr(i,lya)
              prslk1(i,lyb) = (plyr(i,lyb)*0.001) ** rocp ! plyr in Pa
              rhly(i,lyb)   = rhly(i,lya)
              qstl(i,lyb)   = qstl(i,lya)
           enddo
           
           ! note: may need to take care the top layer amount
           tracer1(:,lyb,:) = tracer1(:,lya,:)
        endif
        
        ! Get layer ozone mass mixing ratio 
        if (Model%ntoz > 0) then 
           do k=1,lmk
              do i=1,im
                 olyr(i,k) = max( QMIN, tracer1(i,k,Model%ntoz) )
              enddo
           enddo
           ! Use climatological ozone data
        else                               
           call getozn (prslk1, Grid%xlat, IM, LMK, olyr) 
        endif

        ! Compute cosine of zenith angle (only when SW is called)
        if (Model%lsswr) then
           call coszmn (Grid%xlon, Grid%sinlat, Grid%coslat, Model%solhr, IM, me, &
                Radtend%coszen, Radtend%coszdg)
        endif

        ! Call getgases(), to set up non-prognostic gas volume mixing ratios (gasvmr).
        call getgases (plvl/100., Grid%xlon, Grid%xlat, IM, LMK, gasvmr)
        
        ! Assign to gasvmr_XXXX
        gasvmr_co2    (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,1)
        gasvmr_n2o    (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,2)
        gasvmr_ch4    (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,3)
        gasvmr_o2     (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,4)
        gasvmr_co     (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,5)
        gasvmr_cfc11  (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,6)
        gasvmr_cfc12  (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,7)
        gasvmr_cfc22  (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,8)
        gasvmr_ccl4   (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,9)
        gasvmr_cfc113 (1:IM,1:LMK)  = gasvmr(1:IM,1:LMK,10)
 
        ! Get temperature at layer interface, and layer moisture.
        tem2da(1:IM,2:LMK) = log( plyr(1:IM,2:LMK) )
        tem2db(1:IM,2:LMK) = log( plvl(1:IM,2:LMK) )

        if (ivflip == 0) then              ! input data from toa to sfc
           do i = 1, IM
              tem1d (i)   = QME6
              tem2da(i,1) = log( plyr(i,1) )
              tem2db(i,1) = log( max( kdist_lw%get_press_min(), plvl(i,1)) )
              tem2db(i,LMP) = log( plvl(i,LMP) )
              tsfa  (i)   = tlyr(i,LMK)                  ! sfc layer air temp
              tlvl(i,1)   = tlyr(i,1)
              tlvl(i,LMP) = tskn(i)
           enddo
           do k = 1, LM
              k1 = k + kd
              do i = 1, IM
                 qlyr(i,k1) = max( tem1d(i), Statein%qgrs(i,k,1) )
                 tem1d(i)   = min( QME5, qlyr(i,k1) )
                 tvly(i,k1) = Statein%tgrs(i,k) * (1.0 + fvirt*qlyr(i,k1)) ! virtual T (K)
                 delp(i,k1) = plvl(i,k1+1) - plvl(i,k1)
              enddo
           enddo

           if ( lextop ) then
              do i = 1, IM
                 qlyr(i,lyb) = qlyr(i,lya)
                 tvly(i,lyb) = tvly(i,lya)
                 delp(i,lyb) = plvl(i,lla) - plvl(i,llb)
              enddo
           endif

           do k = 2, LMK
              do i = 1, IM
                 tlvl(i,k) = tlyr(i,k) + (tlyr(i,k-1) - tlyr(i,k))  * (tem2db(i,k)   - tem2da(i,k))  / &
                      (tem2da(i,k-1) - tem2da(i,k))
              enddo
           enddo

           ! Comput lvel height and layer thickness (km)
           tem0d = 0.001 * rog
           do i = 1, IM
              do k = 1, LMK
                 dz(i,k) = tem0d * (tem2db(i,k+1) - tem2db(i,k)) * tvly(i,k)
              enddo
           enddo

        else                               ! input data from sfc to toa
           do i = 1, IM
              tem1d (i)   = QME6
              tem2da(i,1) = log( plyr(i,1) )
              tem2db(i,1) = log( plvl(i,1) )
              tem2db(i,LMP) = log( max( kdist_lw%get_press_min(), plvl(i,LMP)) )
              tsfa  (i)   = tlyr(i,1)                    ! sfc layer air temp
              tlvl(i,1)   = tskn(i)
              tlvl(i,LMP) = tlyr(i,LMK)
           enddo
           
           do k = LM, 1, -1
              do i = 1, IM
                 qlyr(i,k) = max( tem1d(i), Statein%qgrs(i,k,1) )
                 tem1d(i)  = min( QME5, qlyr(i,k) )
                 tvly(i,k) = Statein%tgrs(i,k) * (1.0 + fvirt*qlyr(i,k)) ! virtual T (K)
                 delp(i,k) = plvl(i,k) - plvl(i,k+1)
              enddo
           enddo
           
           if ( lextop ) then
              do i = 1, IM
                 qlyr(i,lyb) = qlyr(i,lya)
                 tvly(i,lyb) = tvly(i,lya)
                 delp(i,lyb) = plvl(i,lla) - plvl(i,llb)
              enddo
           endif
           
           do k = 1, LMK-1
              do i = 1, IM
                 tlvl(i,k+1) = tlyr(i,k) + (tlyr(i,k+1) - tlyr(i,k)) * (tem2db(i,k+1) - tem2da(i,k))  / &
                      (tem2da(i,k+1) - tem2da(i,k))
              enddo
           enddo

           ! Compute level height and layer thickness (km)
           tem0d = 0.001 * rog
           do i = 1, IM
              do k = LMK, 1, -1
                 dz(i,k) = tem0d * (tem2db(i,k) - tem2db(i,k+1)) * tvly(i,k)
              enddo
           enddo
        endif                              ! end_if_ivflip

        ! Call module_radiation_aerosols::setaer(),to setup aerosols property profile for radiation.
        call setaer (plvl, plyr, prslk1, tvly, rhly, Sfcprop%slmsk,  tracer1, Grid%xlon, &
                     Grid%xlat, IM, LMK, LMP, Model%lsswr, Model%lslwr, faersw, faerlw, aerodp)

        ! Store aerosol optical properties
        ! SW. 
        ! For RRTMGP SW the bands are now ordered from [IR(band) -> nIR -> UV], in RRTMG the 
        ! band ordering was [nIR -> UV -> IR(band)]
        faersw1(1:IM,1:LMK,1)          = faersw(1:IM,1:LMK,nBandsSW,1)
        faersw2(1:IM,1:LMK,1)          = faersw(1:IM,1:LMK,nBandsSW,2)
        faersw3(1:IM,1:LMK,1)          = faersw(1:IM,1:LMK,nBandsSW,3)
        faersw1(1:IM,1:LMK,2:nBandsSW) = faersw(1:IM,1:LMK,1:nBandsSW-1,1)
        faersw2(1:IM,1:LMK,2:nBandsSW) = faersw(1:IM,1:LMK,1:nBandsSW-1,2)
        faersw3(1:IM,1:LMK,2:nBandsSW) = faersw(1:IM,1:LMK,1:nBandsSW-1,3)
        ! LW
        faerlw1(1:IM,1:LMK,1:nBandsLW) = faerlw(1:IM,1:LMK,1:nBandsLW,1)
        faerlw2(1:IM,1:LMK,1:nBandsLW) = faerlw(1:IM,1:LMK,1:nBandsLW,2)
        faerlw3(1:IM,1:LMK,1:nBandsLW) = faerlw(1:IM,1:LMK,1:nBandsLW,3)
   
        !  Obtain cloud information for radiation calculations
        !    (clouds,cldsa,mtopa,mbota)
        !   for  prognostic cloud:
        !    - For Zhao/Moorthi's prognostic cloud scheme,
        !      call module_radiation_clouds::progcld1()
        !    - For Zhao/Moorthi's prognostic cloud+pdfcld,
        !      call module_radiation_clouds::progcld3()
        !      call module_radiation_clouds::progclduni() for unified cloud and ncld=2
        ccnd = 0.0_kind_phys
        if (Model%ncnd == 1) then                                       ! Zhao_Carr_Sundqvist
           ccnd(1:IM,1:LMK,1) = tracer1(1:IM,1:LMK,ntcw)                ! -liquid water/ice
        elseif (Model%ncnd == 2) then                                   ! MG
           ccnd(1:IM,1:LMK,1) = tracer1(1:IM,1:LMK,ntcw)                ! -liquid water
           ccnd(1:IM,1:LMK,2) = tracer1(1:IM,1:LMK,ntiw)                ! -ice water
        elseif (Model%ncnd == 4) then                                   ! MG2
           ccnd(1:IM,1:LMK,1) = tracer1(1:IM,1:LMK,ntcw)                ! -liquid water
           ccnd(1:IM,1:LMK,2) = tracer1(1:IM,1:LMK,ntiw)                ! -ice water
           ccnd(1:IM,1:LMK,3) = tracer1(1:IM,1:LMK,ntrw)                ! -rain water
           ccnd(1:IM,1:LMK,4) = tracer1(1:IM,1:LMK,ntsw)                ! -snow water
        elseif (Model%ncnd == 5) then                                   ! GFDL MP, Thompson, MG3
           ccnd(1:IM,1:LMK,1) = tracer1(1:IM,1:LMK,ntcw)                ! -liquid water
           ccnd(1:IM,1:LMK,2) = tracer1(1:IM,1:LMK,ntiw)                ! -ice water
           ccnd(1:IM,1:LMK,3) = tracer1(1:IM,1:LMK,ntrw)                ! -rain water
           ccnd(1:IM,1:LMK,4) = tracer1(1:IM,1:LMK,ntsw) + &            ! -snow + grapuel
                                tracer1(1:IM,1:LMK,ntgl) 
        endif
        where(ccnd < epsq) ccnd = 0.0

        if (Model%imp_physics == 11 ) then
           if (.not. Model%lgfdlmprad) then
              ccnd(:,:,1) =               tracer1(:,1:LMK,ntcw)
              ccnd(:,:,1) = ccnd(:,:,1) + tracer1(:,1:LMK,ntrw)
              ccnd(:,:,1) = ccnd(:,:,1) + tracer1(:,1:LMK,ntiw)
              ccnd(:,:,1) = ccnd(:,:,1) + tracer1(:,1:LMK,ntsw)
              ccnd(:,:,1) = ccnd(:,:,1) + tracer1(:,1:LMK,ntgl)
           endif
           do k=1,LMK
              do i=1,IM
                 if (ccnd(i,k,1) < EPSQ ) ccnd(i,k,1) = 0.0
              enddo
           enddo
        endif

        ! Add suspended convective cloud water to grid-scale cloud water
        ! only for cloud fraction & radiation computation it is to enhance 
        ! cloudiness due to suspended convec cloud water for zhao/moorthi's 
        ! (imp_phys=99) & ferrier's (imp_phys=5) microphysics schemes
        if ((Model%num_p3d == 4) .and. (Model%npdf3d == 3)) then       ! same as Model%imp_physics = 99
           deltaq(1:im,1+kd:lm+kd) = Tbd%phy_f3d(1:im,1:lm,5)
           cnvw  (1:im,1+kd:lm+kd) = Tbd%phy_f3d(1:im,1:lm,6)
           cnvc  (1:im,1+kd:lm+kd) = Tbd%phy_f3d(1:im,1:lm,7)
        elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then ! same as MOdel%imp_physics=98
           deltaq(1:im,1+kd:lm+kd) = 0.0
           cnvw  (1:im,1+kd:lm+kd) = Tbd%phy_f3d(1:im,1:lm,Model%num_p3d+1)
           cnvc  (1:im,1+kd:lm+kd) = 0.0
        else                                                           ! all the rest
           deltaq(1:im,1:lmk) = 0.0
           cnvw  (1:im,1:lmk) = 0.0
           cnvc  (1:im,1:lmk) = 0.0
        endif
        
        if (lextop) then
           cldcov(1:im,lyb) = cldcov(1:im,lya)
           deltaq(1:im,lyb) = deltaq(1:im,lya)
           cnvw  (1:im,lyb) = cnvw  (1:im,lya)
           cnvc  (1:im,lyb) = cnvc  (1:im,lya)
           if (Model%effr_in) then
              effrl(1:im,lyb) = effrl(1:im,lya)
              effri(1:im,lyb) = effri(1:im,lya)
              effrr(1:im,lyb) = effrr(1:im,lya)
              effrs(1:im,lyb) = effrs(1:im,lya)
           endif
        endif
        
        if (Model%imp_physics == 99) then
           ccnd(1:IM,1:LMK,1) = ccnd(1:IM,1:LMK,1) + cnvw(1:IM,1:LMK)
        endif
        
        if (Model%imp_physics == 10) then
           ccnd(1:IM,1:LMK,1) = ccnd(1:IM,1:LMK,1) + cnvw(1:IM,1:LMK) + ccnd(1:IM,1:LMK,2)
        endif

        ! DJS2019: START        
        ! Compute layer cloud fraction.
        clwmin = 0.0
        cldcov(:,:) = 0.0
        if (.not. Model%lmfshal) then
           do k = 1, LMK
              do i = 1, IM
                 clwt = 1.0e-6 * (plyr(i,k)*0.1)
                 if (ccnd(i,k,1) > 0.) then
                    onemrh= max( 1.e-10, 1.0-rhly(i,k) )
                    clwm  = clwmin / max( 0.01, plyr(i,k)*0.1 )
                    tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
                    tem1  = 2000.0 / tem1
                    value = max( min( tem1*(ccnd(i,k,1)-clwm), 50.0 ), 0.0 )
                    tem2  = sqrt( sqrt(rhly(i,k)) )
                    cldcov(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
                 endif
              enddo
           enddo
        else
           do k = 1, LMK
              do i = 1, IM
                 clwt = 1.0e-6 * (plyr(i,k)*0.1)
                 if (ccnd(i,k,1) .gt. 0) then
                    onemrh= max( 1.e-10, 1.0-rhly(i,k) )
                    clwm  = clwmin / max( 0.01, plyr(i,k)*0.1 )
                    tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
                    if (Model%lmfdeep2) then
                       tem1  = xrc3 / tem1
                    else
                       tem1  = 100.0 / tem1
                    endif
                    value = max( min( tem1*(ccnd(i,k,1)-clwm), 50.0 ), 0.0 )
                    tem2  = sqrt( sqrt(rhly(i,k)) )
                    cldcov(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
                 endif
              enddo
           enddo
        endif
        ! DJS2019: END

        if (Model%uni_cld) then
           if (Model%effr_in) then
              cldcov(1:im,1+kd:lm+kd) = Tbd%phy_f3d(1:im,1:lm,Model%indcld)
              effrl(1:im,1+kd:lm+kd)  = Tbd%phy_f3d(1:im,1:lm,2)
              effri(1:im,1+kd:lm+kd)  = Tbd%phy_f3d(1:im,1:lm,3)
              effrr(1:im,1+kd:lm+kd)  = Tbd%phy_f3d(1:im,1:lm,4)
              effrs(1:im,1+kd:lm+kd)  = Tbd%phy_f3d(1:im,1:lm,5)
           else
              do k=1,lm
                 k1 = k + kd
                 do i=1,im
                    !cldcov(i,k1) = Tbd%phy_f3d(i,k,Model%indcld)
                    !if (tracer1(i,k,ntcw) .gt. 0 .or. tracer1(i,k,ntiw) .gt. 0) then
                    !   cldcov(i,k1) = 0.1
                    !else
                    !   cldcov(i,k1) = 0.0
                    !endif
                 enddo
              enddo
           endif
        elseif (Model%imp_physics == Model%imp_physics_gfdl) then                          ! GFDL MP
           cldcov(1:IM,1+kd:LM+kd) = tracer1(1:IM,1:LM,Model%ntclamt)
        else                                                           ! neither of the other two cases
           ! cldcov = 0.0
        endif

        ! MICROPHYSICS
        ! *) zhao/moorthi's prognostic cloud scheme or unified cloud and/or with MG microphysics
        if (Model%imp_physics == 99 .or. Model%imp_physics == 10) then           
           if (Model%uni_cld .and. Model%ncld >= 2) then
              call progclduni (plyr/100., plvl/100., tlyr, tvly, ccnd, Model%ncnd,          & ! IN
                               Grid%xlat, Grid%xlon, Sfcprop%slmsk, dz, delp/100.,IM,       & ! IN
                               LMK, LMP, cldcov, effrl, effri, effrr, effrs, Model%effr_in, & ! IN
                               clouds, cldsa, mtopa, mbota, de_lgth)                          ! OUT
           else
              call progcld1 (plyr/100. ,plvl/100., tlyr, tvly, qlyr, qstl, rhly,            & ! IN
                             ccnd(1:IM,1:LMK,1), Grid%xlat,Grid%xlon,Sfcprop%slmsk, dz,     & ! IN
                             delp/100., IM, LMK, LMP, Model%uni_cld, Model%lmfshal,         & ! IN
                             Model%lmfdeep2, cldcov, effrl, effri, effrr, effrs,            & ! IN
                             Model%effr_in,                                                 & ! IN
                             clouds, cldsa, mtopa, mbota, de_lgth)                            ! OUT 
           endif
        ! *) zhao/moorthi's prognostic cloud+pdfcld
        elseif(Model%imp_physics == 98) then
           call progcld3 (plyr/100., plvl/100., tlyr, tvly, qlyr, qstl, rhly,               & ! IN
                          ccnd(1:IM,1:LMK,1), cnvw, cnvc, Grid%xlat, Grid%xlon,             & ! IN
                          Sfcprop%slmsk, dz, delp/100., im, lmk, lmp, deltaq, Model%sup,    & ! IN
                          Model%kdt, me,                                                    & ! IN
                          clouds, cldsa, mtopa, mbota, de_lgth)                               ! OUT
        ! *) GFDL cloud scheme
        elseif (Model%imp_physics == 11) then 
           if (.not.Model%lgfdlmprad) then
              call progcld4 (plyr/100., plvl/100., tlyr, tvly, qlyr, qstl, rhly,            & ! IN
                           ccnd(1:IM,1:LMK,1), cnvw, cnvc,Grid%xlat, Grid%xlon,             & ! IN
                           Sfcprop%slmsk, cldcov, dz, delp/100., im, lmk, lmp,              & ! IN
                           clouds, cldsa, mtopa, mbota, de_lgth)                              ! OUT
           else
              call progclduni (plyr/100., plvl/100., tlyr, tvly, ccnd, Model%ncnd,          & ! IN
                              Grid%xlat, Grid%xlon, Sfcprop%slmsk, dz,delp/100., IM, LMK,   & ! IN
                              LMP, cldcov, effrl, effri, effrr, effrs, Model%effr_in,       & ! IN
                              clouds, cldsa, mtopa, mbota, de_lgth)                           ! OUT
           endif
        ! *) Thompson / WSM6 cloud micrphysics scheme
        elseif(Model%imp_physics == 8 .or. Model%imp_physics == 6) then
           if (Model%kdt == 1) then
              Tbd%phy_f3d(:,:,1) = 10.
              Tbd%phy_f3d(:,:,2) = 50.
              Tbd%phy_f3d(:,:,3) = 250.
           endif

           call progcld5 (plyr/100., plvl/100., tlyr, qlyr, qstl, rhly, tracer1, Grid%xlat, & ! IN
                         Grid%xlon,Sfcprop%slmsk,dz,delp/100., ntrac-1, ntcw-1, ntiw-1,     & ! IN
                         ntrw-1, ntsw-1, ntgl-1, im, lmk, lmp, Model%uni_cld, Model%lmfshal,& ! IN
                         Model%lmfdeep2, cldcov(:,1:LMK),Tbd%phy_f3d(:,:,1),                & ! IN
                         Tbd%phy_f3d(:,:,2), Tbd%phy_f3d(:,:,3),                            & ! IN
                         clouds,cldsa,mtopa,mbota, de_lgth)                                   ! OUT
        endif                            ! end if_imp_physics

        ! CCPP
        clouds1(1:IM,1:LMK)  = clouds(1:IM,1:LMK,1)
        clouds2(1:IM,1:LMK)  = clouds(1:IM,1:LMK,2)
        clouds3(1:IM,1:LMK)  = clouds(1:IM,1:LMK,3)
        clouds4(1:IM,1:LMK)  = clouds(1:IM,1:LMK,4)
        clouds5(1:IM,1:LMK)  = clouds(1:IM,1:LMK,5)
        clouds6(1:IM,1:LMK)  = clouds(1:IM,1:LMK,6)
        clouds7(1:IM,1:LMK)  = clouds(1:IM,1:LMK,7)
        clouds8(1:IM,1:LMK)  = clouds(1:IM,1:LMK,8)
        clouds9(1:IM,1:LMK)  = clouds(1:IM,1:LMK,9)

        ! mg, sfc-perts
        !  ---  scale random patterns for surface perturbations with
        !  perturbation size
        !  ---  turn vegetation fraction pattern into percentile pattern
        alb1d(:) = 0.
        if (Model%do_sfcperts) then
           if (Model%pertalb(1) > 0.) then
              do i=1,im
                 call cdfnor(Coupling%sfc_wts(i,5),alb1d(i))
              enddo
           endif
        endif
        ! mg, sfc-perts

      end subroutine GFS_rrtmgp_pre_run

!> \section arg_table_GFS_rrtmgp_pre_finalize Argument Table
!!
      subroutine GFS_rrtmgp_pre_finalize ()
        close(58)
      end subroutine GFS_rrtmgp_pre_finalize

!! @}
      end module GFS_rrtmgp_pre
