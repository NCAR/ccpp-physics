!> \file GFS_rrtmgp_pre.f90
!! This file contains
module GFS_rrtmgp_pre
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
  use module_radiation_surface,  only: &
       setemis,                  & ! Routine to compute surface-emissivity
       NF_ALBD,                  & ! Number of surface albedo categories (4; nir-direct, nir-diffuse, uvvis-direct, uvvis-diffuse)
       setalb                      ! Routine to compute surface albedo
  ! RRTMGP types
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use rrtmgp_aux,            only: check_error_msg

  real(kind_phys), parameter :: &
       amd   = 28.9644_kind_phys,  & ! Molecular weight of dry-air     (g/mol)
       amw   = 18.0154_kind_phys,  & ! Molecular weight of water vapor (g/mol)
       amo3  = 47.9982_kind_phys,  & ! Modelular weight of ozone       (g/mol)
       amdw  = amd/amw,            & ! Molecular weight of dry air / water vapor
       amdo3 = amd/amo3              ! Molecular weight of dry air / ozone
  
  public GFS_rrtmgp_pre_run,GFS_rrtmgp_pre_init,GFS_rrtmgp_pre_finalize
  
contains

  subroutine GFS_rrtmgp_pre_init ()
  end subroutine GFS_rrtmgp_pre_init

!> \section arg_table_GFS_rrtmgp_pre_run Argument Table
!! | local_name         | standard_name                                                 | long_name                                                                     | units    | rank |  type                 |   kind    | intent | optional |
!! |--------------------|---------------------------------------------------------------|-------------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model              | GFS_control_type_instance                                     | Fortran DDT containing FV3-GFS model control parameters                       | DDT      |    0 | GFS_control_type      |           | in     | F        |
!! | Grid               | GFS_grid_type_instance                                        | Fortran DDT containing FV3-GFS grid and interpolation related data            | DDT      |    0 | GFS_grid_type         |           | in     | F        |
!! | Sfcprop            | GFS_sfcprop_type_instance                                     | Fortran DDT containing FV3-GFS surface fields                                 | DDT      |    0 | GFS_sfcprop_type      |           | in     | F        |
!! | Statein            | GFS_statein_type_instance                                     | Fortran DDT containing FV3-GFS prognostic state data in from dycore           | DDT      |    0 | GFS_statein_type      |           | in     | F        |
!! | Tbd                | GFS_tbd_type_instance                                         | Fortran DDT containing FV3-GFS data not yet assigned to a defined container   | DDT      |    0 | GFS_tbd_type          |           | in     | F        |
!! | Coupling           | GFS_coupling_type_instance                                    | Fortran DDT containing FV3-GFS fields needed for coupling                     | DDT      |    0 | GFS_coupling_type     |           | in     | F        |
!! | Radtend            | GFS_radtend_type_instance                                     | Fortran DDT containing FV3-GFS radiation tendencies                           | DDT      |    0 | GFS_radtend_type      |           | inout  | F        |
!! | ncol               | horizontal_loop_extent                                        | horizontal loop extent                                                        | count    |    0 | integer               |           | in     | F        |
!! | lw_gas_props       | coefficients_for_lw_gas_optics                                | DDT containing spectral information for RRTMGP LW radiation scheme            | DDT      |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | sw_gas_props       | coefficients_for_sw_gas_optics                                | DDT containing spectral information for RRTMGP SW radiation scheme            | DDT      |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | raddt              | time_step_for_radiation                                       | radiation time step                                                           | s        |    0 | real                  | kind_phys | out    | F        |
!! | p_lay              | air_pressure_at_layer_for_RRTMGP_in_hPa                       | air pressure at vertical layer for radiation calculation                      | hPa      |    2 | real                  | kind_phys | out    | F        |
!! | p_lev              | air_pressure_at_interface_for_RRTMGP_in_hPa                   | air pressure at vertical interface for radiation calculation                  | hPa      |    2 | real                  | kind_phys | out    | F        |
!! | t_lay              | air_temperature_at_layer_for_RRTMGP                           | air temperature at vertical layer for radiation calculation                   | K        |    2 | real                  | kind_phys | out    | F        |
!! | t_lev              | air_temperature_at_interface_for_RRTMGP                       | air temperature  at vertical interface for radiation calculation              | K        |    2 | real                  | kind_phys | out    | F        |
!! | tsfg               | surface_ground_temperature_for_radiation                      | surface ground temperature for radiation                                      | K        |    1 | real                  | kind_phys | out    | F        |
!! | tsfa               | surface_air_temperature_for_radiation                         | lowest model layer air temperature for radiation                              | K        |    1 | real                  | kind_phys | out    | F        |
!! | cld_frac           | total_cloud_fraction                                          | layer total cloud fraction                                                    | frac     |    2 | real                  | kind_phys | out    | F        |
!! | cld_lwp            | cloud_liquid_water_path                                       | layer cloud liquid water path                                                 | g m-2    |    2 | real                  | kind_phys | out    | F        |
!! | cld_reliq          | mean_effective_radius_for_liquid_cloud                        | mean effective radius for liquid cloud                                        | micron   |    2 | real                  | kind_phys | out    | F        |
!! | cld_iwp            | cloud_ice_water_path                                          | layer cloud ice water path                                                    | g m-2    |    2 | real                  | kind_phys | out    | F        |
!! | cld_reice          | mean_effective_radius_for_ice_cloud                           | mean effective radius for ice cloud                                           | micron   |    2 | real                  | kind_phys | out    | F        |
!! | cld_swp            | cloud_snow_water_path                                         | layer cloud snow water path                                                   | g m-2    |    2 | real                  | kind_phys | out    | F        |
!! | cld_resnow         | mean_effective_radius_for_snow_flake                          | mean effective radius for snow cloud                                          | micron   |    2 | real                  | kind_phys | out    | F        |
!! | cld_rwp            | cloud_rain_water_path                                         | layer cloud rain water path                                                   | g m-2    |    2 | real                  | kind_phys | out    | F        |
!! | cld_rerain         | mean_effective_radius_for_rain_drop                           | mean effective radius for rain cloud                                          | micron   |    2 | real                  | kind_phys | out    | F        |
!! | faerlw             | aerosol_optical_properties_for_longwave_bands_01-16           | aerosol optical properties for longwave bands 01-16                           | various  |    4 | real                  | kind_phys | out    | F        |
!! | faersw             | aerosol_optical_properties_for_shortwave_bands_01-16          | aerosol optical properties for shortwave bands 01-16                          | various  |    4 | real                  | kind_phys | out    | F        |
!! | mtopa              | model_layer_number_at_cloud_top                               | vertical indices for low, middle and high cloud tops                          | index    |    2 | integer               |           | out    | F        |
!! | mbota              | model_layer_number_at_cloud_base                              | vertical indices for low, middle and high cloud bases                         | index    |    2 | integer               |           | out    | F        |
!! | cldsa              | cloud_area_fraction_for_radiation                             | fraction of clouds for low, middle, high, total and BL                        | frac     |    2 | real                  | kind_phys | out    | F        |
!! | aerodp             | atmosphere_optical_thickness_due_to_ambient_aerosol_particles | vertical integrated optical depth for various aerosol species                 | none     |    2 | real                  | kind_phys | out    | F        |
!! | alb1d              | surface_albedo_perturbation                                   | surface albedo perturbation                                                   | frac     |    1 | real                  | kind_phys | out    | F        |
!! | gas_concentrations | Gas_concentrations_for_RRTMGP_suite                           | DDT containing gas concentrations for RRTMGP radiation scheme                 | DDT      |    0 | ty_gas_concs          |           | out    | F        |
!! | nday               | daytime_points_dimension                                      | daytime points dimension                                                      | count    |    0 | integer               |           | out    | F        |
!! | idxday             | daytime_points                                                | daytime points                                                                | index    |    1 | integer               |           | out    | F        |
!! | errmsg             | ccpp_error_message                                            | error message for error handling in CCPP                                      | none     |    0 | character             | len=*     | out    | F        |
!! | errflg             | ccpp_error_flag                                               | error flag for error handling in CCPP                                         | flag     |    0 | integer               |           | out    | F        |
!!
  ! Attention - the output arguments lm, im, lmk, lmp must not be set
  ! in the CCPP version - they are defined in the interstitial_create routine
  subroutine GFS_rrtmgp_pre_run (Model, Grid, Statein, Coupling, Radtend, Sfcprop, Tbd, & ! IN
       ncol, lw_gas_props, sw_gas_props,                                                & ! IN
       raddt, p_lay, t_lay, p_lev, t_lev, tsfg, tsfa, alb1d, cld_frac, cld_lwp,         & ! OUT
       cld_reliq, cld_iwp, cld_reice, cld_swp, cld_resnow, cld_rwp, cld_rerain, faerlw, & ! OUT
       faersw, cldsa, mtopa, mbota, aerodp, nday, idxday, gas_concentrations, errmsg, errflg)
    
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                ! Fortran DDT containing FV3-GFS model control parameters
    type(GFS_grid_type), intent(in) :: &
         Grid                 ! Fortran DDT containing FV3-GFS grid and interpolation related data 
    type(GFS_statein_type), intent(in) :: &
         Statein              ! Fortran DDT containing FV3-GFS prognostic state data in from dycore    
    type(GFS_coupling_type), intent(in) :: &
         Coupling             ! Fortran DDT containing FV3-GFS fields to/from coupling with other components 
    type(GFS_radtend_type), intent(inout) :: &
         Radtend              ! Fortran DDT containing FV3-GFS radiation tendencies 
    type(GFS_sfcprop_type), intent(in) :: &
         Sfcprop              ! Fortran DDT containing FV3-GFS surface fields
    type(GFS_tbd_type), intent(in) :: &
         Tbd                  ! Fortran DDT containing FV3-GFS data not yet assigned to a defined container
    integer, intent(in)    :: &
         ncol                 ! Number of horizontal grid points
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props,      & ! RRTMGP DDT containing spectral information for LW calculation
         sw_gas_props         ! RRTMGP DDT containing spectral information for SW calculation

    ! Outputs
    real(kind_phys), dimension(ncol,Model%levs), intent(out) :: &
         p_lay,             & ! Pressure at model-layer
         t_lay                ! Temperature at model layer
    real(kind_phys), dimension(ncol,Model%levs+1), intent(out) :: &
         p_lev,             & ! Pressure at model-interface
         t_lev                ! Temperature at model-interface
    real(kind_phys), intent(out) :: &
         raddt                ! Radiation time-step
    real(kind_phys), dimension(ncol), intent(out) :: &
         tsfg,              & ! Ground temperature
         tsfa                 ! Skin temperature
    integer, intent(out)   :: &
         nday                 ! Number of daylit points
    integer, dimension(ncol), intent(out) :: &
         idxday               ! Indices for daylit points
    real(kind_phys), dimension(ncol), intent(out) :: &
         alb1d                ! Surface albedo pertubation
    type(ty_gas_concs),intent(out) :: &
         gas_concentrations   ! RRTMGP DDT containing gas volumne mixing ratios
    character(len=*), intent(out) :: &
         errmsg               ! Error message
    integer, intent(out) :: &  
         errflg               ! Error flag
    real(kind_phys), dimension(ncol,Model%levr+LTP),intent(out) :: &
         cld_frac,          & ! Total cloud fraction
         cld_lwp,           & ! Cloud liquid water path
         cld_reliq,         & ! Cloud liquid effective radius
         cld_iwp,           & ! Cloud ice water path
         cld_reice,         & ! Cloud ice effecive radius
         cld_swp,           & ! Cloud snow water path
         cld_resnow,        & ! Cloud snow effective radius
         cld_rwp,           & ! Cloud rain water path
         cld_rerain           ! Cloud rain effective radius
    real(kind_phys), dimension(ncol,Model%levs,sw_gas_props%get_nband(),NF_AESW), intent(out) ::&
         faersw               ! Aerosol radiative properties in each SW band.
    real(kind_phys), dimension(ncol,Model%levs,lw_gas_props%get_nband(),NF_AELW), intent(out) ::&
         faerlw               ! Aerosol radiative properties in each LW band.
    integer,dimension(ncol,3),intent(out) :: &
         mbota,             & ! Vertical indices for cloud tops
         mtopa                ! Vertical indices for cloud bases
    real(kind_phys), dimension(ncol,5), intent(out) :: &
         cldsa                ! Fraction of clouds for low, middle, high, total and BL 
    real(kind_phys), dimension(ncol,NSPC1), intent(out) :: &
         aerodp               ! Vertical integrated optical depth for various aerosol species  

    ! Local variables
    integer :: i, j, k, iCol, iBand, iSFC, iTOA, iLay
    logical :: top_at_1
    real(kind_phys),dimension(NCOL,Model%levs) :: vmr_o3, vmr_h2o
    real(kind_phys) :: es, qs
    real(kind_phys), dimension(ncol)  :: de_lgth
    real(kind_phys), dimension(ncol, NF_ALBD) :: sfcalb
    real(kind_phys), dimension(ncol, Model%levs) :: relhum, qs_lay, q_lay, deltaZ, tv_lay,&
         deltaP, o3_lay
    real(kind_phys), dimension(ncol, Model%levs, 2:Model%ntrac) :: tracer
    real(kind_phys), dimension(ncol, Model%levs, NF_VGAS) :: gas_vmr
    real(kind_phys), dimension(ncol, Model%levs, NF_CLDS) :: clouds
    real(kind_phys), dimension(ncol, Model%levs, sw_gas_props%get_nband(), NF_AESW)::faersw2

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. (Model%lsswr .or. Model%lslwr)) return
    
    ! #######################################################################################
    ! What is vertical ordering?
    ! #######################################################################################
    top_at_1 = (Statein%prsi(1,1) .lt.  Statein%prsi(1, Model%levs))
    if (top_at_1) then 
       iSFC = Model%levs
       iTOA = 1
    else
       iSFC = 1
       iTOA = Model%levs
    endif

    ! #######################################################################################
    ! Compute some fields needed by RRTMGP
    ! #######################################################################################
    ! Copy state fields over for use in RRTMGP
    p_lev(1:NCOL,iSFC:iTOA) = Statein%prsi(1:NCOL,1:Model%levs)
    p_lev(1:NCOL,iTOA+1)    = spread(lw_gas_props%get_press_min(),dim=1, ncopies=NCOL)
    p_lay(1:NCOL,iSFC:iTOA) = Statein%prsl(1:NCOL,1:Model%levs)
    t_lay(1:NCOL,iSFC:iTOA) = Statein%tgrs(1:NCOL,1:Model%levs)

    ! Compute layer pressure thicknes
    deltaP = p_lev(:,iSFC:iTOA)-p_lev(:,iSFC+1:iTOA+1)

    ! Compute temperature at layer-interfaces
    t_lev(1:NCOL,iSFC) = Sfcprop%tsfc(1:NCOL)
    do iCol=1,NCOL
       do iLay=iSFC+1,iTOA
          t_lev(iCol,iLay) = (t_lay(iCol,iLay)+t_lay(iCol,iLay-1))/2._kind_phys
       enddo
       t_lev(iCol,iTOA+1) = lw_gas_props%get_temp_min()
    enddo

    ! Compute a bunch of thermodynamic fields needed by the macrophysics schemes. Relative humidity, 
    ! saturation mixing-ratio, vapor mixing-ratio, virtual temperature, layer thickness,...
    do iCol=1,NCOL
       do iLay=iSFC,iTOA
          es                = min( Statein%prsl(iCol,iLay),  fpvs( Statein%tgrs(iCol,iLay) ) )  ! fpvs and prsl in pa
          qs                = max( QMIN, eps * es / (Statein%prsl(iCol,iLay) + epsm1*es) )
          relhum(iCol,iLay) = max( 0._kind_phys, min( 1._kind_phys, max(QMIN, Statein%qgrs(iCol,iLay,1))/qs ) )
          qs_lay(iCol,iLay) = qs
          q_lay(iCol,iLay)  = max( 1.e-6, Statein%qgrs(iCol,iLay,1) )
          tv_lay(iCol,iLay) = Statein%tgrs(iCol,iLay) * (1._kind_phys + fvirt*q_lay(iCol,iLay)) 
          deltaZ(iCol,iLay) = (rog*0.001) * (log(p_lev(iCol,iLay)) - log(p_lev(iCol,iLay+1))) * tv_lay(iCol,iLay)
       enddo
    enddo

    ! #######################################################################################
    ! Get layer ozone mass mixing ratio 
    ! #######################################################################################
    ! First recast remaining all tracers (except sphum) forcing them all to be positive
    do j = 2, model%NTRAC
       tracer(1:NCOL,1:Model%levs,j) = max(0.0, Statein%qgrs(1:NCOL,1:Model%levs,j))
    enddo

    if (Model%ntoz > 0) then 
       do iLay=iSFC,iTOA
          do iCol=1,NCOL
             o3_lay(iCol,iLay) = max( QMIN, tracer(iCol,iLay,Model%ntoz) )
          enddo
       enddo
    ! OR Use climatological ozone data
    else                               
       call getozn (Statein%prslk(1:NCOL,iSFC:iTOA), Grid%xlat, NCOL, Model%levs, o3_lay) 
    endif

    ! #######################################################################################
    ! Set gas concentrations for RRTMGP
    ! #######################################################################################
    ! Call getgases(), to set up non-prognostic gas volume mixing ratios (gas_vmr).
    call getgases (p_lev/100., Grid%xlon, Grid%xlat, NCOL, Model%levs, gas_vmr)

    ! Compute volume mixing-ratios for ozone (mmr) and specific-humidity.
    vmr_h2o = merge((q_lay/(1-q_lay))*amdw, 0., q_lay  .ne. 1.)
    vmr_o3  = merge(o3_lay*amdo3,           0., o3_lay .gt. 0.)
    !
    call gas_concentrations%reset()
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr('o2',  gas_vmr(:,:,4)))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr('co2', gas_vmr(:,:,1)))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr('ch4', gas_vmr(:,:,3)))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr('n2o', gas_vmr(:,:,2)))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr('h2o', vmr_h2o))
    call check_error_msg('GFS_rrtmgp_pre_run',gas_concentrations%set_vmr('o3',  vmr_o3))

    ! #######################################################################################
    ! Radiation time step (output) (Is this really needed?)
    ! #######################################################################################
    raddt = min(Model%fhswr, Model%fhlwr)

    ! #######################################################################################
    ! Compute cosine of zenith angle (only when SW is called)
    ! #######################################################################################
    if (Model%lsswr) then
       call coszmn (Grid%xlon, Grid%sinlat, Grid%coslat, Model%solhr, NCOL, Model%me, &
            Radtend%coszen, Radtend%coszdg)
    endif

    ! #######################################################################################
    ! Call module_radiation_aerosols::setaer(),to setup aerosols property profile for both 
    ! LW and SW radiation.
    ! #######################################################################################
    call setaer(p_lev, p_lay, Statein%prslk(1:NCOL,iSFC:iTOA), tv_lay, relhum,              &
         Sfcprop%slmsk,  tracer, Grid%xlon, Grid%xlat, NCOL, Model%levs, Model%levs+1,      &
         Model%lsswr, Model%lslwr, faersw2, faerlw, aerodp)
    
    ! Store aerosol optical properties
    ! SW. 
    ! For RRTMGP SW the bands are now ordered from [IR(band) -> nIR -> UV], in RRTMG the 
    ! band ordering was [nIR -> UV -> IR(band)]
    faersw(1:NCOL,1:Model%levs,1,1)                      = faersw2(1:NCOL,1:Model%levs,sw_gas_props%get_nband(),1)
    faersw(1:NCOL,1:Model%levs,1,2)                      = faersw2(1:NCOL,1:Model%levs,sw_gas_props%get_nband(),2)
    faersw(1:NCOL,1:Model%levs,1,3)                      = faersw2(1:NCOL,1:Model%levs,sw_gas_props%get_nband(),3)
    faersw(1:NCOL,1:Model%levs,2:sw_gas_props%get_nband(),1) = faersw2(1:NCOL,1:Model%levs,1:sw_gas_props%get_nband()-1,1)
    faersw(1:NCOL,1:Model%levs,2:sw_gas_props%get_nband(),2) = faersw2(1:NCOL,1:Model%levs,1:sw_gas_props%get_nband()-1,2)
    faersw(1:NCOL,1:Model%levs,2:sw_gas_props%get_nband(),3) = faersw2(1:NCOL,1:Model%levs,1:sw_gas_props%get_nband()-1,3)

    ! Setup surface ground temperature and ground/air skin temperature if required.
    tsfg(1:NCOL) = Sfcprop%tsfc(1:NCOL)
    tsfa(1:NCOL) = Sfcprop%tsfc(1:NCOL)

    ! #######################################################################################
    ! Cloud microphysics
    ! #######################################################################################
    call cloud_microphysics(Model, Tbd, Grid, Sfcprop, ncol, tracer, p_lay, t_lay,    &
         p_lev, tv_lay, relhum, qs_lay, q_lay, deltaZ, deltaP, &
         clouds, cldsa, mbota, mtopa, de_lgth)

    ! Copy output cloud fields
    cld_frac   = clouds(:,:,1)
    cld_lwp    = clouds(:,:,2)
    cld_reliq  = clouds(:,:,3)
    cld_iwp    = clouds(:,:,4)
    cld_reice  = clouds(:,:,5)   
    cld_rwp    = clouds(:,:,6)  
    cld_rerain = clouds(:,:,7)  
    cld_swp    = clouds(:,:,8)  
    cld_resnow = clouds(:,:,9)  

    ! #######################################################################################
    ! mg, sfc-perts
    !  ---  scale random patterns for surface perturbations with perturbation size
    !  ---  turn vegetation fraction pattern into percentile pattern
    ! #######################################################################################
    alb1d(:) = 0.
    if (Model%do_sfcperts) then
       if (Model%pertalb(1) > 0.) then
          do i=1,ncol
             call cdfnor(Coupling%sfc_wts(i,5),alb1d(i))
          enddo
       endif
    endif    

    ! #######################################################################################
    ! Call module_radiation_surface::setemis(),to setup surface emissivity for LW radiation.
    ! #######################################################################################
    if (Model%lslwr) then
       call setemis (Grid%xlon, Grid%xlat, Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr,     &
            Sfcprop%zorl, tsfg, tsfa, Sfcprop%hprim, NCOL,  Radtend%semis)
       do iBand=1,lw_gas_props%get_nband()
          Radtend%sfc_emiss_byband(iBand,1:NCOL) = Radtend%semis(1:NCOL)
       enddo
    endif

    ! #######################################################################################
    ! For SW, gather daylit points, compute surface albedo in each band,
    ! #######################################################################################
    if (Model%lsswr) then
       ! Check for daytime points for SW radiation.
       nday = 0
       idxday = 0
       do i = 1, NCOL
          if (Radtend%coszen(i) >= 0.0001) then
             nday = nday + 1
             idxday(nday) = i
          endif
       enddo
       
       ! Call module_radiation_surface::setalb() to setup surface albedo.
       call setalb (Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr,&    !  ---  inputs:
                    Sfcprop%snoalb, Sfcprop%zorl, Radtend%coszen,&
                    tsfg, tsfa, Sfcprop%hprim, Sfcprop%alvsf,    &
                    Sfcprop%alnsf, Sfcprop%alvwf, Sfcprop%alnwf, &
                    Sfcprop%facsf, Sfcprop%facwf, Sfcprop%fice,  &
                    Sfcprop%tisfc, NCOL,                         &
                    alb1d, Model%pertalb,                        &    !  mg, sfc-perts
                    sfcalb)                                           !  ---  outputs
       
       ! Approximate mean surface albedo from vis- and nir-  diffuse values.
       Radtend%sfalb(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
    else
       nday   = 0
       idxday = 0
       sfcalb = 0.0
    endif
      
    ! Spread across all SW bands
    do iBand=1,sw_gas_props%get_nband()
       Radtend%sfc_alb_nir_dir(iBand,1:NCOL)   = sfcalb(1:NCOL,1)
       Radtend%sfc_alb_nir_dif(iBand,1:NCOL)   = sfcalb(1:NCOL,2)
       Radtend%sfc_alb_uvvis_dir(iBand,1:NCOL) = sfcalb(1:NCOL,3)
       Radtend%sfc_alb_uvvis_dif(iBand,1:NCOL) = sfcalb(1:NCOL,4)
    enddo

  end subroutine GFS_rrtmgp_pre_run
  
!> \section arg_table_GFS_rrtmgp_pre_finalize Argument Table
!!
  subroutine GFS_rrtmgp_pre_finalize ()
  end subroutine GFS_rrtmgp_pre_finalize

  ! #######################################################################################
  ! Subroutine cloud_microphysics()
  ! #######################################################################################
  subroutine cloud_microphysics(Model, Tbd, Grid, Sfcprop, ncol, tracer, p_lay, t_lay,    &
       p_lev, tv_lay, relhum, qs_lay, q_lay, deltaZ, deltaP, &
       clouds, cldsa, mbota, mtopa, de_lgth)
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                ! Fortran DDT containing FV3-GFS model control parameters
    type(GFS_tbd_type), intent(in) :: &
         Tbd                  ! Fortran DDT containing FV3-GFS data not yet assigned to a defined container
    type(GFS_grid_type), intent(in) :: &
         Grid                 ! Fortran DDT containing FV3-GFS grid and interpolation related data 
    type(GFS_sfcprop_type), intent(in) :: &
         Sfcprop              ! Fortran DDT containing FV3-GFS surface fields
    integer, intent(in) :: &
         ncol ! Number of horizontal gridpoints
    real(kind_phys), dimension(ncol, Model%levs, 2:Model%ntrac),intent(in) :: &
         tracer               !
    real(kind_phys), dimension(ncol,Model%levs), intent(in) :: &
         p_lay,             & !
         t_lay,             & !
         tv_lay,            & !
         relhum,            & !
         qs_lay,            & !
         q_lay,             & !
         deltaZ,            & !
         deltaP
    real(kind_phys), dimension(ncol,Model%levs+1), intent(in) :: &
         p_lev                !

    ! Outputs
    real(kind_phys), dimension(ncol, Model%levs, NF_CLDS),intent(out) :: clouds
    integer,dimension(ncol,3), intent(out) :: mbota, mtopa
    real(kind_phys), dimension(ncol), intent(out)  :: de_lgth
    real(kind_phys), dimension(ncol, 5), intent(out) :: cldsa

    ! Local variables
    real(kind_phys), dimension(ncol, Model%levs, Model%ncnd) :: cld_condensate
    integer :: i,k
    real(kind_phys) :: clwmin, clwm, clwt, onemrh, value, tem1, tem2
    real(kind_phys), parameter :: xrc3 = 100.
    real(kind_phys), dimension(ncol, Model%levs) :: delta_q, cnv_w, cnv_c, effr_l, effr_i, effr_r, effr_s, cldcov

    ! #######################################################################################
    !  Obtain cloud information for radiation calculations
    !    (clouds,cldsa,mtopa,mbota)
    !   for  prognostic cloud:
    !    - For Zhao/Moorthi's prognostic cloud scheme,
    !      call module_radiation_clouds::progcld1()
    !    - For Zhao/Moorthi's prognostic cloud+pdfcld,
    !      call module_radiation_clouds::progcld3()
    !      call module_radiation_clouds::progclduni() for unified cloud and ncld=2
    ! #######################################################################################
    cld_condensate = 0.0_kind_phys
    if (Model%ncnd == 1) then                                                                    ! Zhao_Carr_Sundqvist
       cld_condensate(1:NCOL,1:Model%levs,1) = tracer(1:NCOL,1:Model%levs,Model%ntcw)            ! -liquid water/ice
    elseif (Model%ncnd == 2) then                                                                ! MG
       cld_condensate(1:NCOL,1:Model%levs,1) = tracer(1:NCOL,1:Model%levs,Model%ntcw)            ! -liquid water
       cld_condensate(1:NCOL,1:Model%levs,2) = tracer(1:NCOL,1:Model%levs,Model%ntiw)            ! -ice water
    elseif (Model%ncnd == 4) then                                                                ! MG2
       cld_condensate(1:NCOL,1:Model%levs,1) = tracer(1:NCOL,1:Model%levs,Model%ntcw)            ! -liquid water
       cld_condensate(1:NCOL,1:Model%levs,2) = tracer(1:NCOL,1:Model%levs,Model%ntiw)            ! -ice water
       cld_condensate(1:NCOL,1:Model%levs,3) = tracer(1:NCOL,1:Model%levs,Model%ntrw)            ! -rain water
       cld_condensate(1:NCOL,1:Model%levs,4) = tracer(1:NCOL,1:Model%levs,Model%ntsw)            ! -snow water
    elseif (Model%ncnd == 5) then                                                                ! GFDL MP, Thompson, MG3
       cld_condensate(1:NCOL,1:Model%levs,1) = tracer(1:NCOL,1:Model%levs,Model%ntcw)            ! -liquid water
       cld_condensate(1:NCOL,1:Model%levs,2) = tracer(1:NCOL,1:Model%levs,Model%ntiw)            ! -ice water
       cld_condensate(1:NCOL,1:Model%levs,3) = tracer(1:NCOL,1:Model%levs,Model%ntrw)            ! -rain water
       cld_condensate(1:NCOL,1:Model%levs,4) = tracer(1:NCOL,1:Model%levs,Model%ntsw) + &        ! -snow + grapuel
                                               tracer(1:NCOL,1:Model%levs,Model%ntgl) 
    endif
    where(cld_condensate < epsq) cld_condensate = 0.0
    
    ! For GFDL microphysics scheme...
    if (Model%imp_physics == 11 ) then
       if (.not. Model%lgfdlmprad) then
          cld_condensate(:,:,1) =                         tracer(:,1:Model%levs,Model%ntcw)
          cld_condensate(:,:,1) = cld_condensate(:,:,1) + tracer(:,1:Model%levs,Model%ntrw)
          cld_condensate(:,:,1) = cld_condensate(:,:,1) + tracer(:,1:Model%levs,Model%ntiw)
          cld_condensate(:,:,1) = cld_condensate(:,:,1) + tracer(:,1:Model%levs,Model%ntsw)
          cld_condensate(:,:,1) = cld_condensate(:,:,1) + tracer(:,1:Model%levs,Model%ntgl)
       endif
       do k=1,Model%levs
          do i=1,NCOL
             if (cld_condensate(i,k,1) < EPSQ ) cld_condensate(i,k,1) = 0.0
          enddo
       enddo
    endif

    ! Add suspended convective cloud water to grid-scale cloud water
    ! only for cloud fraction & radiation computation it is to enhance 
    ! cloudiness due to suspended convec cloud water for zhao/moorthi's 
    ! (imp_phys=99) & ferrier's (imp_phys=5) microphysics schemes
    if ((Model%num_p3d == 4) .and. (Model%npdf3d == 3)) then       ! same as Model%imp_physics = 99
       delta_q(1:ncol,1:Model%levs) = Tbd%phy_f3d(1:ncol,1:Model%levs,5)
       cnv_w  (1:ncol,1:Model%levs) = Tbd%phy_f3d(1:ncol,1:Model%levs,6)
       cnv_c  (1:ncol,1:Model%levs) = Tbd%phy_f3d(1:ncol,1:Model%levs,7)
    elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then ! same as MOdel%imp_physics=98
       delta_q(1:ncol,1:Model%levs) = 0.0
       cnv_w  (1:ncol,1:Model%levs) = Tbd%phy_f3d(1:ncol,1:Model%levs,Model%num_p3d+1)
       cnv_c  (1:ncol,1:Model%levs) = 0.0
    else                                                           ! all the rest
       delta_q(1:ncol,1:Model%levs) = 0.0
       cnv_w  (1:ncol,1:Model%levs) = 0.0
       cnv_c  (1:ncol,1:Model%levs) = 0.0
    endif

    ! For zhao/moorthi's prognostic cloud scheme, add in convective cloud water to liquid-cloud water
    if (Model%imp_physics == 99) then
       cld_condensate(1:NCOL,1:Model%levs,1) = cld_condensate(1:NCOL,1:Model%levs,1) + cnv_w(1:NCOL,1:Model%levs)
    endif
    
    ! For MG prognostic cloud scheme, add in convective cloud water to liquid-and-ice-cloud condensate
    if (Model%imp_physics == 10) then
       cld_condensate(1:NCOL,1:Model%levs,1) = cld_condensate(1:NCOL,1:Model%levs,1) + cnv_w(1:NCOL,1:Model%levs) + cld_condensate(1:NCOL,1:Model%levs,2)
    endif

    if (Model%uni_cld) then
       if (Model%effr_in) then
          cldcov(1:ncol,1:Model%levs)  = Tbd%phy_f3d(1:ncol,1:Model%levs,Model%indcld)
          effr_l(1:ncol,1:Model%levs)  = Tbd%phy_f3d(1:ncol,1:Model%levs,2)
          effr_i(1:ncol,1:Model%levs)  = Tbd%phy_f3d(1:ncol,1:Model%levs,3)
          effr_r(1:ncol,1:Model%levs)  = Tbd%phy_f3d(1:ncol,1:Model%levs,4)
          effr_s(1:ncol,1:Model%levs)  = Tbd%phy_f3d(1:ncol,1:Model%levs,5)
       else
          do k=1,model%levs
             do i=1,ncol
                cldcov(i,k) = Tbd%phy_f3d(i,k,Model%indcld)
                if (tracer(i,k,model%ntcw) .gt. 0 .or. tracer(i,k,model%ntiw) .gt. 0) then
                   cldcov(i,k) = 0.1
                else
                   cldcov(i,k) = 0.0
                endif
             enddo
          enddo
       endif
    elseif (Model%imp_physics == Model%imp_physics_gfdl) then                          ! GFDL MP
       cldcov(1:NCOL,1:Model%levs) = tracer(1:NCOL,1:Model%levs,Model%ntclamt)
    else                                                           ! neither of the other two cases
       cldcov = 0.0
    endif

    ! #######################################################################################
    ! This is a hack to get the first-column in a file to contain a cloud.
    ! #######################################################################################
    ! DJS2019: START        
    ! Compute layer cloud fraction.
    clwmin = 0.0
    cldcov(:,:) = 0.0
    if (.not. Model%lmfshal) then
       do k = 1, Model%levs
          do i = 1, NCOL
             clwt = 1.0e-6 * (p_lay(i,k)*0.1)
             if (cld_condensate(i,k,1) > 0.) then
                onemrh= max( 1.e-10, 1.0-relhum(i,k) )
                clwm  = clwmin / max( 0.01, p_lay(i,k)*0.1 )
                tem1  = min(max(sqrt(sqrt(onemrh*qs_lay(i,k))),0.0001),1.0)
                tem1  = 2000.0 / tem1
                value = max( min( tem1*(cld_condensate(i,k,1)-clwm), 50.0 ), 0.0 )
                tem2  = sqrt( sqrt(relhum(i,k)) )
                cldcov(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
             endif
          enddo
       enddo
    else
       do k = 1, Model%levs
          do i = 1, NCOL
             clwt = 1.0e-6 * (p_lay(i,k)*0.1)
             if (cld_condensate(i,k,1) .gt. 0) then
                onemrh= max( 1.e-10, 1.0-relhum(i,k) )
                clwm  = clwmin / max( 0.01, p_lay(i,k)*0.1 )
                tem1  = min(max((onemrh*qs_lay(i,k))**0.49,0.0001),1.0)  !jhan
                if (Model%lmfdeep2) then
                   tem1  = xrc3 / tem1
                else
                   tem1  = 100.0 / tem1
                endif
                value = max( min( tem1*(cld_condensate(i,k,1)-clwm), 50.0 ), 0.0 )
                tem2  = sqrt( sqrt(relhum(i,k)) )
                cldcov(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
             endif
          enddo
       enddo
    endif
    ! DJS2019: END

    ! #######################################################################################
    ! MICROPHYSICS
    ! #######################################################################################
    ! *) zhao/moorthi's prognostic cloud scheme or unified cloud and/or with MG microphysics
    if (Model%imp_physics == 99 .or. Model%imp_physics == 10) then           
       if (Model%uni_cld .and. Model%ncld >= 2) then
          call progclduni(           &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               cld_condensate,       & ! IN  - Cloud condensate amount (Model%ncnd types)     ()
               Model%ncnd,           & ! IN  - Number of cloud condensate types               ()
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               deltaZ,               & ! IN  - Layer thickness                                (m)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               cldcov,               & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
               effr_l,               & ! IN  - Liquid-water effective radius                  (microns)
               effr_i,               & ! IN  - Ice-water effective radius                     (microns)
               effr_r,               & ! IN  - Rain-water effective radius                    (microns)
               effr_s,               & ! IN  - Snow-water effective radius                    (microns)
               Model%effr_in,        & ! IN  - Logical, if .true. use input effective radii
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       else
          call progcld1 (            &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               q_lay,                & ! IN  - Specific humidity at layer center              (kg/kg)
               qs_lay,               & ! IN  - Saturation specific humidity at layer center   (kg/kg)
               relhum,               & ! IN  - Relative humidity at layer center              (1)
               cld_condensate(:,:,1),& ! IN  - Cloud condensate amount                        ()
                                       !       (Zhao: liq+convective; MG: liq+ice+convective) 
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               deltaZ,               & ! IN  - Layer thickness                                (m)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               Model%uni_cld,        & ! IN  - True for cloud fraction from shoc
               Model%lmfshal,        & ! IN  - True for mass flux shallow convection
               Model%lmfdeep2,       & ! IN  - True for mass flux deep convection
               cldcov,               & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
               effr_l,               & ! IN  - Liquid-water effective radius                  (microns)
               effr_i,               & ! IN  - Ice-water effective radius                     (microns)
               effr_r,               & ! IN  - Rain-water effective radius                    (microns)
               effr_s,               & ! IN  - Snow-water effective radius                    (microns)
               Model%effr_in,        & ! IN  - Logical, if .true. use input effective radii
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       endif
       ! *) zhao/moorthi's prognostic cloud+pdfcld
    elseif(Model%imp_physics == 98) then
       call progcld3 (               &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               q_lay,                & ! IN  - Specific humidity at layer center              (kg/kg)
               qs_lay,               & ! IN  - Saturation specific humidity at layer center   (kg/kg)
               relhum,               & ! IN  - Relative humidity at layer center              (1)
               cld_condensate(:,:,1),& ! IN  - Cloud condensate amount (only h20)             ()
               cnv_w,                & ! IN  - Layer convective cloud condensate
               cnv_c,                & ! IN  - Layer convective cloud cover
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               deltaZ,               & ! IN  - Layer thickness                                (m)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               delta_q,              & ! IN  - Total water distribution width
               Model%sup,            & ! IN  - ??? Supersaturation?
               Model%kdt,            & ! IN  - ??? 
               Model%me,             & ! IN  - ??? NOT USED IN PROGCLD3()
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       ! *) GFDL cloud scheme
    elseif (Model%imp_physics == 11) then 
       if (.not.Model%lgfdlmprad) then
          call progcld4 (            &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               q_lay,                & ! IN  - Specific humidity at layer center              (kg/kg)
               qs_lay,               & ! IN  - Saturation specific humidity at layer center   (kg/kg)
               relhum,               & ! IN  - Relative humidity at layer center              (1)
               cld_condensate(:,:,1),& ! IN  - Cloud condensate amount (only h20)             ()
               cnv_w,                & ! IN  - Layer convective cloud condensate
               cnv_c,                & ! IN  - Layer convective cloud cover
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               cldcov,               & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
               deltaZ,               & ! IN  - Layer thickness                                (m)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       else
          call progclduni(           &
               p_lay/100.,           & ! IN  - Pressure at model layer centers                (mb)
               p_lev/100.,           & ! IN  - Pressure at model interfaces                   (mb)
               t_lay,                & ! IN  - Temperature at layer centers                   (K)
               tv_lay,               & ! IN  - Virtual temperature at layer centers           (K)
               cld_condensate,       & ! IN  - Cloud condensate amount (Model%ncnd types)     ()
               Model%ncnd,           & ! IN  - Number of cloud condensate types               ()
               Grid%xlat,            & ! IN  - Latitude                                       (radians)
               Grid%xlon,            & ! IN  - Longitude                                      (radians)
               Sfcprop%slmsk,        & ! IN  - Land/Sea mask                                  ()
               deltaZ,               & ! IN  - Layer thickness                                (m)
               deltaP/100.,          & ! IN  - Layer thickness                                (hPa)
               NCOL,                 & ! IN  - Number of horizontal gridpoints
               MODEL%LEVS,           & ! IN  - Number of model layers
               MODEL%LEVS+1,         & ! IN  - Number of model levels
               cldcov,               & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
               effr_l,               & ! IN  - Liquid-water effective radius                  (microns)
               effr_i,               & ! IN  - Ice-water effective radius                     (microns)
               effr_r,               & ! IN  - Rain-water effective radius                    (microns)
               effr_s,               & ! IN  - Snow-water effective radius                    (microns)
               Model%effr_in,        & ! IN  - Logical, if .true. use input effective radii
               clouds,               & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
               cldsa,                & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
               mtopa,                & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
               mbota,                & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
               de_lgth)                ! OUT - clouds decorrelation length (km)
       endif
       ! *) Thompson / WSM6 cloud micrphysics scheme
    elseif(Model%imp_physics == 8 .or. Model%imp_physics == 6) then
 
       call progcld5 ( & ! IN
            p_lay/100.,               & ! IN  - Pressure at model layer centers                (mb)
            p_lev/100.,               & ! IN  - Pressure at model interfaces                   (mb)
            t_lay,                    & ! IN  - Temperature at layer centers                   (K)
            q_lay,                    & ! IN  - Specific humidity at layer center              (kg/kg)
            qs_lay,                   & ! IN  - Saturation specific humidity at layer center   (kg/kg)
            relhum,                   & ! IN  - Relative humidity at layer center              (1)
            tracer,                   & ! IN  - Cloud condensate amount in layer by type       ()
            Grid%xlat,                & ! IN  - Latitude                                       (radians)
            Grid%xlon,                & ! IN  - Longitude                                      (radians)
            Sfcprop%slmsk,            & ! IN  - Land/Sea mask                                  ()
            deltaZ,                   & ! IN  - Layer thickness                                (m)
            deltaP/100.,              & ! IN  - Layer thickness                                (hPa)
            Model%ntrac-1,            & ! IN  - Number of tracers
            Model%ntcw-1,             & ! IN  - Tracer index for cloud condensate (or liquid water)
            Model%ntiw-1,             & ! IN  - Tracer index for ice 
            Model%ntrw-1,             & ! IN  - Tracer index for rain 
            Model%ntsw-1,             & ! IN  - Tracer index for snow 
            Model%ntgl-1,             & ! IN  - Tracer index for groupel
            NCOL,                     & ! IN  - Number of horizontal gridpoints
            MODEL%LEVS,               & ! IN  - Number of model layers
            MODEL%LEVS+1,             & ! IN  - Number of model levels
            Model%uni_cld,            & ! IN  - True for cloud fraction from shoc
            Model%lmfshal,            & ! IN  - True for mass flux shallow convection
            Model%lmfdeep2,           & ! IN  - True for mass flux deep convection
            cldcov(:,1:Model%levs),   & ! IN  - Layer cloud fraction (used if uni_cld=.true.)
            Tbd%phy_f3d(:,:,1),       & ! IN  - Liquid-water effective radius                  (microns)
            Tbd%phy_f3d(:,:,2),       & ! IN  - Ice-water effective radius                     (microns)
            Tbd%phy_f3d(:,:,3),       & ! IN  - LSnow-water effective radius                   (microns)
            clouds,                   & ! OUT - Cloud properties                               (NCOL,Model%levs,NF_CLDS)
            cldsa,                    & ! OUT - fraction of clouds for low, mid, hi, tot, bl   (NCOL,5)
            mtopa,                    & ! OUT - vertical indices for low, mid, hi cloud tops   (NCOL,3)
            mbota,                    & ! OUT - vertical indices for low, mid, hi cloud bases  (NCOL,3)
            de_lgth)                    ! OUT - clouds decorrelation length (km)
    endif ! end if_imp_physics
     
  end subroutine cloud_microphysics
end module GFS_rrtmgp_pre
