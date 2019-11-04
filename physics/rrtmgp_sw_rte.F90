! ###########################################################################################
! ###########################################################################################
module rrtmgp_sw_rte
  use machine,                 only: kind_phys
  use GFS_typedefs,            only: GFS_control_type, GFS_radtend_type, GFS_statein_type, GFS_interstitial_type
  use mo_rte_kind,             only: wl
  use mo_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,         only: ty_cloud_optics
  use mo_optical_props,        only: ty_optical_props_2str
  use mo_rte_sw,               only: rte_sw
  use mo_gas_concentrations,   only: ty_gas_concs
  use mo_fluxes_byband,        only: ty_fluxes_byband
  use module_radsw_parameters, only: cmpfsw_type
  use rrtmgp_aux,              only: check_error_msg

  public rrtmgp_sw_rte_init, rrtmgp_sw_rte_run, rrtmgp_sw_rte_finalize

contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_rte_init
  ! #########################################################################################
  subroutine rrtmgp_sw_rte_init()
  end subroutine rrtmgp_sw_rte_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_rte_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_rte_run
!! \htmlinclude rrtmgp_sw_rte.html
!!
  subroutine rrtmgp_sw_rte_run(Model, Interstitial, Radtend, Statein, ncol, sw_gas_props, p_lay, t_lay,   &
       p_lev, gas_concentrations, sw_optical_props_clrsky, sw_optical_props_clouds,         &
       sw_optical_props_aerosol, lsswr, nday, idxday, hsw0, hswb, scmpsw,                   &
       fluxswUP_allsky, fluxswDOWN_allsky, fluxswUP_clrsky, fluxswDOWN_clrsky, errmsg, errflg)

    ! Inputs
    type(GFS_control_type),   intent(in)    :: &
         Model
    type(GFS_interstitial_type), intent(in) :: &
         Interstitial
    type(GFS_radtend_type),   intent(in)    :: &
         Radtend
    type(GFS_statein_type), intent(in) :: &
         Statein                 ! Fortran DDT containing FV3-GFS prognostic state data in from dycore 
    integer, intent(in) :: &
         ncol,                 & ! Number of horizontal gridpoints
         nday                    ! Number of daytime points
    integer, intent(in), dimension(nday) :: &
         idxday                  ! Index array for daytime points
    real(kind_phys), dimension(ncol,Model%levs), intent(in) :: &
         p_lay,                & ! Pressure @ model layer-centers         (hPa)
         t_lay                   ! Temperature                            (K)
    real(kind_phys), dimension(ncol,Model%levs+1), intent(in) :: &
         p_lev                   ! Pressure @ model layer-interfaces      (hPa)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props                ! DDT containing SW spectral information
    type(ty_optical_props_2str),intent(in) :: &
         sw_optical_props_clrsky, & ! RRTMGP DDT: longwave clear-sky radiative properties 
         sw_optical_props_clouds, & ! RRTMGP DDT: longwave cloud radiative properties 
         sw_optical_props_aerosol ! RRTMGP DDT: longwave aerosol radiative properties

    type(ty_gas_concs),intent(in) :: &
         gas_concentrations      ! RRTMGP DDT: trace gas concentrations   (vmr)
    logical, intent(in) :: &
         lsswr                   ! Flag to calculate SW irradiances

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg
    real(kind_phys), dimension(ncol,Model%levs+1), intent(out) :: &
         fluxswUP_allsky,   & ! All-sky flux                    (W/m2)
         fluxswDOWN_allsky, & ! All-sky flux                    (W/m2)
         fluxswUP_clrsky,   & ! Clear-sky flux                  (W/m2)
         fluxswDOWN_clrsky    ! All-sky flux                    (W/m2)

    ! Inputs (optional) (NOTE. We only need the optional arguments to know what fluxes to output, HR's are computed later)
    real(kind_phys), dimension(ncol,Model%levs), optional, intent(inout) :: &
         hsw0             ! Clear-sky heating rate            (K/sec)
    real(kind_phys), dimension(ncol,Model%levs,sw_gas_props%get_nband()), intent(inout), optional :: &
         hswb             ! All-sky heating rate, by band     (K/sec)
    ! Outputs (optional)
    type(cmpfsw_type), dimension(ncol), intent(inout),optional :: &
         scmpsw           ! 2D surface fluxes, components:
                          ! uvbfc - total sky downward uv-b flux at  (W/m2)
                          ! uvbf0 - clear sky downward uv-b flux at  (W/m2)
                          ! nirbm - downward nir direct beam flux    (W/m2)
                          ! nirdf - downward nir diffused flux       (W/m2)
                          ! visbm - downward uv+vis direct beam flux (W/m2)
                          ! visdf - downward uv+vis diffused flux    (W/m2)

    ! Local variables
    type(ty_fluxes_byband) :: &
         flux_allsky, & ! All-sky flux                      (W/m2)
         flux_clrsky    ! Clear-sky flux                    (W/m2)
    real(kind_phys), dimension(nday,Model%levs+1),target :: &
         fluxSW_up_allsky, fluxSW_up_clrsky, fluxSW_dn_allsky, fluxSW_dn_clrsky, fluxSW_dn_dir_allsky
    real(kind_phys), dimension(nday,Model%levs+1,sw_gas_props%get_nband()),target :: &
         fluxSWBB_up_allsky, fluxSWBB_dn_allsky
    real(kind_phys), dimension(ncol,Model%levs) :: vmrTemp
    logical :: l_ClrSky_HR=.false., l_AllSky_HR_byband=.false., l_scmpsw=.false., top_at_1
    integer :: iGas,iSFC,iTOA
    type(ty_optical_props_2str)  :: &
         sw_optical_props_clouds_daylit,  & ! RRTMGP DDT: longwave cloud radiative properties 
         sw_optical_props_clrsky_daylit, & ! RRTMGP DDT: longwave clear-sky radiative properties 
         sw_optical_props_aerosol_daylit   ! RRTMGP DDT: longwave aerosol radiative properties
    type(ty_gas_concs) :: &
         gas_concentrations_daylit    ! RRTMGP DDT: trace gas concentrations   (vmr)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg  = 0

    if (.not. lsswr) return

    ! Vertical ordering?
    top_at_1 = (Statein%prsi(1,1) .lt.  Statein%prsi(1, Model%levs))
    if (top_at_1) then 
       iSFC = Model%levs+1
       iTOA = 1
    else
       iSFC = 1
       iTOA = Model%levs+1
    endif

    ! Are any optional outputs requested? Need to know now to compute correct fluxes.
    l_ClrSky_HR        = present(hsw0)
    l_AllSky_HR_byband = present(hswb)
    l_scmpsw           = present(scmpsw)
    if ( l_scmpsw ) then
       scmpsw = cmpfsw_type (0., 0., 0., 0., 0., 0.)
    endif
    fluxswUP_allsky(:,:)   = 0._kind_phys
    fluxswDOWN_allsky(:,:) = 0._kind_phys
    fluxswUP_clrsky(:,:)   = 0._kind_phys
    fluxswDOWN_clrsky(:,:) = 0._kind_phys

    if (nDay .gt. 0) then

       ! Subset the cloud and aerosol radiative properties over daylit points.
       ! Cloud optics [nDay,Model%levs,nGpts]
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clouds_daylit%alloc_2str(nday, Model%levs, sw_gas_props))
       sw_optical_props_clouds_daylit%tau    = sw_optical_props_clouds%tau(idxday,:,:)
       sw_optical_props_clouds_daylit%ssa    = sw_optical_props_clouds%ssa(idxday,:,:)
       sw_optical_props_clouds_daylit%g      = sw_optical_props_clouds%g(idxday,:,:)
       ! Aerosol optics [nDay,Model%levs,nBands]
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_aerosol_daylit%alloc_2str(nday, Model%levs, sw_gas_props%get_band_lims_wavenumber()))
       sw_optical_props_aerosol_daylit%tau = sw_optical_props_aerosol%tau(idxday,:,:)
       sw_optical_props_aerosol_daylit%ssa = sw_optical_props_aerosol%ssa(idxday,:,:)
       sw_optical_props_aerosol_daylit%g   = sw_optical_props_aerosol%g(idxday,:,:)
       ! Clear-sky optics [nDay,Model%levs,nGpts]
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clrsky_daylit%alloc_2str(nday, Model%levs, sw_gas_props))
       sw_optical_props_clrsky_daylit%tau = sw_optical_props_clrsky%tau(idxday,:,:)
       sw_optical_props_clrsky_daylit%ssa = sw_optical_props_clrsky%ssa(idxday,:,:)
       sw_optical_props_clrsky_daylit%g   = sw_optical_props_clrsky%g(idxday,:,:)
      
       ! Similarly, subset the gas concentrations.
       do iGas=1,Model%nGases
          call check_error_msg('rrtmgp_sw_rte_run',gas_concentrations%get_vmr(trim(Model%active_gases_array(iGas)),vmrTemp))
          call check_error_msg('rrtmgp_sw_rte_run',gas_concentrations_daylit%set_vmr(trim(Model%active_gases_array(iGas)),vmrTemp(idxday,:)))
       enddo

       ! Initialize RRTMGP DDT containing 2D(3D) fluxes
       flux_allsky%flux_up     => fluxSW_up_allsky
       flux_allsky%flux_dn     => fluxSW_dn_allsky
       flux_allsky%flux_dn_dir => fluxSW_dn_dir_allsky
       flux_clrsky%flux_up     => fluxSW_up_clrsky
       flux_clrsky%flux_dn     => fluxSW_dn_clrsky
       ! Only calculate fluxes by-band, only when heating-rate profiles by band are requested.
       if (l_AllSky_HR_byband) then
          flux_allsky%bnd_flux_up => fluxSWBB_up_allsky
          flux_allsky%bnd_flux_dn => fluxSWBB_dn_allsky
       endif

       ! Compute clear-sky fluxes (if requested)
       ! Clear-sky fluxes are gas+aerosol
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_aerosol_daylit%increment(sw_optical_props_clrsky_daylit))
       if (l_ClrSky_HR) then
          call check_error_msg('rrtmgp_sw_rte_run',rte_sw(               &
               sw_optical_props_clrsky_daylit,     & ! IN  - optical-properties
               top_at_1,                           & ! IN  - veritcal ordering flag
               Radtend%coszen(idxday),             & ! IN  - Cosine of solar zenith angle
               Interstitial%toa_src_sw(idxday,:),       & ! IN  - incident solar flux at TOA
               Interstitial%sfc_alb_nir_dir(:,idxday),  & ! IN  - Shortwave surface albedo (direct)
               Interstitial%sfc_alb_nir_dif(:,idxday),  & ! IN  - Shortwave surface albedo (diffuse)
               flux_clrsky))                         ! OUT - Fluxes, clear-sky, 3D (nCol,Model%levs,nBand) 
          ! Store fluxes
          fluxswUP_clrsky(idxday,:)   = flux_clrsky%flux_up
          fluxswDOWN_clrsky(idxday,:) = flux_clrsky%flux_dn
       endif

       ! Compute all-sky fluxes
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clouds_daylit%increment(sw_optical_props_clrsky_daylit))
       call check_error_msg('rrtmgp_sw_rte_run',rte_sw(               &
            sw_optical_props_clrsky_daylit,     & ! IN  - optical-properties
            top_at_1,                           & ! IN  - veritcal ordering flag
            Radtend%coszen(idxday),             & ! IN  - Cosine of solar zenith angle
            Interstitial%toa_src_sw(idxday,:),       & ! IN  - incident solar flux at TOA
            Interstitial%sfc_alb_nir_dir(:,idxday),  & ! IN  - Shortwave surface albedo (direct)
            Interstitial%sfc_alb_nir_dif(:,idxday),  & ! IN  - Shortwave surface albedo (diffuse)
            flux_allsky))                         ! OUT - Fluxes, clear-sky, 3D (nCol,Model%levs,nBand) 
       ! Store fluxes
       fluxswUP_allsky(idxday,:)   = flux_allsky%flux_up
       fluxswDOWN_allsky(idxday,:) = flux_allsky%flux_dn
       if ( l_scmpsw ) then
          scmpsw(idxday)%nirbm = flux_allsky%flux_dn_dir(idxday,iSFC) !Interstitial%sfc_alb_nir_dir(iSFC,idxday)
          scmpsw(idxday)%nirdf = flux_allsky%flux_dn(idxday,iSFC)  - flux_allsky%flux_dn_dir(idxday,iSFC) !Interstitial%sfc_alb_nir_dif(iSFC,idxday)
       endif
    endif
  end subroutine rrtmgp_sw_rte_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_rte_finalize
  ! #########################################################################################
  subroutine rrtmgp_sw_rte_finalize()
  end subroutine rrtmgp_sw_rte_finalize

end module rrtmgp_sw_rte
