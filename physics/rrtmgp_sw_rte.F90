module rrtmgp_sw_rte
  use machine,                 only: kind_phys
  use GFS_typedefs,            only: GFS_control_type, GFS_radtend_type, GFS_statein_type
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
  subroutine rrtmgp_sw_rte_run(Model, Radtend, Statein, ncol, sw_gas_props, p_lay, t_lay,   &
       p_lev, gas_concentrations, sw_optical_props_clrsky, sfc_alb_nir_dir, sfc_alb_nir_dif,&
       sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, toa_src_sw, sw_optical_props_clouds,           &
       sw_optical_props_aerosol, lsswr, nday, idxday, hsw0, hswb, rrtmgp_nGases, active_gases_array, scmpsw, fluxswUP_allsky,  &
       fluxswDOWN_allsky, fluxswUP_clrsky, fluxswDOWN_clrsky, errmsg, errflg)

    ! Inputs

    integer, intent(in) :: &
         rrtmgp_nGases       ! Number of trace gases active in RRTMGP
    character(len=*),dimension(rrtmgp_nGases), intent(in) :: &
         active_gases_array  ! Character array containing trace gases to include in RRTMGP

    type(GFS_control_type), intent(in)    :: &
         Model
    type(GFS_radtend_type), intent(in)    :: &
         Radtend
    type(GFS_statein_type), intent(in) :: &
         Statein                    ! DDT: FV3-GFS prognostic state data in from dycore 
    integer, intent(in) :: &
         ncol,                    & ! Number of horizontal gridpoints
         nday                       ! Number of daytime points
    integer, intent(in), dimension(ncol) :: &
         idxday                     ! Index array for daytime points
    real(kind_phys), dimension(ncol,Model%levs), intent(in) :: &
         p_lay,                   & ! Pressure @ model layer-centers (Pa)
         t_lay                      ! Temperature (K)
    real(kind_phys), dimension(ncol,Model%levs+1), intent(in) :: &
         p_lev                      ! Pressure @ model layer-interfaces (Pa)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props               ! RRTMGP DDT: SW spectral information
    real(kind_phys), dimension(sw_gas_props%get_nband(),ncol), intent(in) :: &
         sfc_alb_nir_dir,         & ! Surface albedo (direct) 
         sfc_alb_nir_dif,         & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir,       & ! Surface albedo (direct)
         sfc_alb_uvvis_dif          ! Surface albedo (diffuse)
    real(kind_phys), dimension(ncol,sw_gas_props%get_ngpt()), intent(in) :: &
         toa_src_sw                 ! TOA incident spectral flux (W/m2)
    type(ty_optical_props_2str),intent(inout) :: &
         sw_optical_props_clrsky, & ! RRTMGP DDT: longwave clear-sky radiative properties 
         sw_optical_props_clouds, & ! RRTMGP DDT: longwave cloud radiative properties 
         sw_optical_props_aerosol   ! RRTMGP DDT: longwave aerosol radiative properties
    type(ty_gas_concs),intent(in) :: &
         gas_concentrations         ! RRTMGP DDT: trace gas concentrations   (vmr)
    logical, intent(in) :: &
         lsswr                      ! Flag to calculate SW irradiances

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg
    real(kind_phys), dimension(ncol,Model%levs+1), intent(inout) :: &
         fluxswUP_allsky,         & ! RRTMGP upward all-sky flux profiles (W/m2)
         fluxswDOWN_allsky,       & ! RRTMGP downward all-sky flux profiles (W/m2)
         fluxswUP_clrsky,         & ! RRTMGP upward clear-sky flux profiles (W/m2)
         fluxswDOWN_clrsky          ! RRTMGP downward clear-sky flux profiles (W/m2)

    ! Inputs (optional) (NOTE. We only need the optional arguments to know what fluxes to output, HR's are computed later)
    real(kind_phys), dimension(ncol,Model%levs), optional, intent(inout) :: &
         hsw0                       ! Clear-sky heating rate (K/sec)
    real(kind_phys), dimension(ncol,Model%levs,sw_gas_props%get_nband()), intent(inout), optional :: &
         hswb                       ! All-sky heating rate, by band (K/sec)
    ! Outputs (optional)
    type(cmpfsw_type), dimension(ncol), intent(inout),optional :: &
         scmpsw                     ! 2D surface fluxes, components:
                                    ! uvbfc - total sky downward uv-b flux (W/m2)
                                    ! uvbf0 - clear sky downward uv-b flux (W/m2)
                                    ! nirbm - downward nir direct beam flux (W/m2)
                                    ! nirdf - downward nir diffused flux (W/m2)
                                    ! visbm - downward uv+vis direct beam flux (W/m2)
                                    ! visdf - downward uv+vis diffused flux (W/m2)

    ! Local variables
    real(kind_phys), dimension(sw_gas_props%get_nband(),nday) :: &
         sfc_alb_dir,sfc_alb_dif
    type(ty_fluxes_byband) :: &
         flux_allsky, & ! All-sky flux (W/m2)
         flux_clrsky    ! Clear-sky flux (W/m2)
    real(kind_phys), dimension(nday,Model%levs+1,sw_gas_props%get_nband()),target :: &
         fluxSW_up_allsky, fluxSW_up_clrsky, fluxSW_dn_allsky, fluxSW_dn_clrsky, fluxSW_dn_dir_allsky
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
    if (nDay .gt. 0) then

       ! Vertical ordering?
       top_at_1 = (p_lev(1,1) .lt. p_lev(1, Model%levs))
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

       ! Initialize fluxes
       fluxswUP_allsky(:,:)   = 0._kind_phys
       fluxswDOWN_allsky(:,:) = 0._kind_phys
       fluxswUP_clrsky(:,:)   = 0._kind_phys
       fluxswDOWN_clrsky(:,:) = 0._kind_phys
       
       ! Subset the gas concentrations, only need daylit points.
       do iGas=1,rrtmgp_nGases
          call check_error_msg('rrtmgp_sw_rte_run',&
               gas_concentrations%get_vmr(trim(active_gases_array(iGas)),vmrTemp))
          call check_error_msg('rrtmgp_sw_rte_run',&
               gas_concentrations_daylit%set_vmr(trim(active_gases_array(iGas)),vmrTemp(idxday(1:nday),:)))
       enddo

       ! Initialize RRTMGP DDT containing 2D(3D) fluxes
       flux_allsky%bnd_flux_up     => fluxSW_up_allsky
       flux_allsky%bnd_flux_dn     => fluxSW_dn_allsky
       flux_allsky%bnd_flux_dn_dir => fluxSW_dn_dir_allsky
       flux_clrsky%bnd_flux_up     => fluxSW_up_clrsky
       flux_clrsky%bnd_flux_dn     => fluxSW_dn_clrsky

       ! In RRTMG, the near-IR and uv-visible surface albedos are averaged.
       sfc_alb_dir = 0.5_kind_phys*(sfc_alb_nir_dir(:,idxday(1:nday)) + sfc_alb_uvvis_dir(:,idxday(1:nday)))
       sfc_alb_dif = 0.5_kind_phys*(sfc_alb_nir_dif(:,idxday(1:nday)) + sfc_alb_uvvis_dif(:,idxday(1:nday)))

       ! Compute clear-sky fluxes (if requested)
       ! Clear-sky fluxes (gas+aerosol)
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_aerosol%increment(sw_optical_props_clrsky))
       ! Delta-scale optical properties
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clrsky%delta_scale())
       if (l_ClrSky_HR) then
          call check_error_msg('rrtmgp_sw_rte_run',rte_sw(     &
               sw_optical_props_clrsky,           & ! IN  - optical-properties
               top_at_1,                          & ! IN  - veritcal ordering flag
               Radtend%coszen(idxday(1:nday)),    & ! IN  - Cosine of solar zenith angle
               toa_src_sw(idxday(1:nday),:),      & ! IN  - incident solar flux at TOA
               sfc_alb_dir,                       & ! IN  - Shortwave surface albedo (direct)
               sfc_alb_dif,                       & ! IN  - Shortwave surface albedo (diffuse)
               flux_clrsky))                        ! OUT - Fluxes, clear-sky, 3D (nCol,Model%levs,nBand) 
          ! Store fluxes
          fluxswUP_clrsky(idxday(1:nday),:)   = sum(flux_clrsky%bnd_flux_up,dim=3)
          fluxswDOWN_clrsky(idxday(1:nday),:) = sum(flux_clrsky%bnd_flux_dn,dim=3)
       endif

       ! Compute all-sky fluxes
       ! All-sky fluxes (clear-sky + clouds)
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clouds%increment(sw_optical_props_clrsky))
       ! Delta-scale optical properties
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clouds%delta_scale())
       call check_error_msg('rrtmgp_sw_rte_run',rte_sw(     &
            sw_optical_props_clrsky,           & ! IN  - optical-properties
            top_at_1,                          & ! IN  - veritcal ordering flag
            Radtend%coszen(idxday(1:nday)),    & ! IN  - Cosine of solar zenith angle
            toa_src_sw(idxday(1:nday),:),      & ! IN  - incident solar flux at TOA
            sfc_alb_dir,                       & ! IN  - Shortwave surface albedo (direct)
            sfc_alb_dif,                       & ! IN  - Shortwave surface albedo (diffuse)
            flux_allsky))                        ! OUT - Fluxes, clear-sky, 3D (nCol,Model%levs,nBand) 
       ! Store fluxes
       fluxswUP_allsky(idxday(1:nday),:)   = sum(flux_allsky%bnd_flux_up,dim=3)
       fluxswDOWN_allsky(idxday(1:nday),:) = sum(flux_allsky%bnd_flux_dn,dim=3)
       if ( l_scmpsw ) then
          scmpsw(idxday(1:nday))%nirbm = sum(flux_allsky%bnd_flux_dn_dir(idxday(1:nday),iSFC,:),dim=2)
          scmpsw(idxday(1:nday))%nirdf = sum(flux_allsky%bnd_flux_dn(idxday(1:nday),iSFC,:),dim=2)  - &
               sum(flux_allsky%bnd_flux_dn_dir(idxday(1:nday),iSFC,:),dim=2)
       endif
    endif
  end subroutine rrtmgp_sw_rte_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_rte_finalize
  ! #########################################################################################
  subroutine rrtmgp_sw_rte_finalize()
  end subroutine rrtmgp_sw_rte_finalize

end module rrtmgp_sw_rte
