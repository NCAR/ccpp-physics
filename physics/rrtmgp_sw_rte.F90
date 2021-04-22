module rrtmgp_sw_rte
  use machine,                 only: kind_phys
  use mo_rte_kind,             only: wl
  use mo_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,         only: ty_cloud_optics
  use mo_optical_props,        only: ty_optical_props_2str
  use mo_rte_sw,               only: rte_sw
  use mo_gas_concentrations,   only: ty_gas_concs
  use mo_fluxes_byband,        only: ty_fluxes_byband
  use module_radsw_parameters, only: cmpfsw_type
  use rrtmgp_aux,              only: check_error_msg
  use rrtmgp_sw_gas_optics,    only: sw_gas_props
  implicit none

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
  subroutine rrtmgp_sw_rte_run(doSWrad, doSWclrsky, nCol, nLev, nDay, idxday, coszen, p_lay, &
       t_lay, p_lev, sw_optical_props_clrsky, sfc_alb_nir_dir, sfc_alb_nir_dif,              &
       sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, toa_src_sw, sw_optical_props_clouds,            &
       sw_optical_props_aerosol, scmpsw, fluxswUP_allsky, fluxswDOWN_allsky, fluxswUP_clrsky,&
       fluxswDOWN_clrsky, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doSWrad,                 & ! Flag to calculate SW irradiances
         doSWclrsky                 ! Compute clear-sky fluxes?
    integer, intent(in) :: &
         nCol,                    & ! Number of horizontal gridpoints
         nday,                    & ! Number of daytime points
         nLev                       ! Number of vertical levels
    integer, intent(in), dimension(ncol) :: &
         idxday                     ! Index array for daytime points
    real(kind_phys),intent(in), dimension(ncol) :: &
         coszen                     ! Cosize of SZA
    real(kind_phys), dimension(ncol,NLev), intent(in) :: &
         p_lay,                   & ! Pressure @ model layer-centers (Pa)
         t_lay                      ! Temperature (K)
    real(kind_phys), dimension(ncol,NLev+1), intent(in) :: &
         p_lev                      ! Pressure @ model layer-interfaces (Pa)
    type(ty_optical_props_2str),intent(inout) :: &
         sw_optical_props_clrsky    ! RRTMGP DDT: shortwave clear-sky radiative properties 
   type(ty_optical_props_2str),intent(in) :: &
         sw_optical_props_clouds, & ! RRTMGP DDT: shortwave cloud radiative properties 
         sw_optical_props_aerosol   ! RRTMGP DDT: shortwave aerosol radiative properties
    real(kind_phys), dimension(sw_gas_props%get_nband(),ncol), intent(in) :: &
         sfc_alb_nir_dir,         & ! Surface albedo (direct) 
         sfc_alb_nir_dif,         & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir,       & ! Surface albedo (direct)
         sfc_alb_uvvis_dif          ! Surface albedo (diffuse)
    real(kind_phys), dimension(ncol,sw_gas_props%get_ngpt()), intent(in) :: &
         toa_src_sw                 ! TOA incident spectral flux (W/m2)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                     ! CCPP error message
    integer, intent(out) :: &
         errflg                     ! CCPP error flag
    real(kind_phys), dimension(ncol,NLev+1), intent(inout) :: &
         fluxswUP_allsky,         & ! RRTMGP upward all-sky flux profiles (W/m2)
         fluxswDOWN_allsky,       & ! RRTMGP downward all-sky flux profiles (W/m2)
         fluxswUP_clrsky,         & ! RRTMGP upward clear-sky flux profiles (W/m2)
         fluxswDOWN_clrsky          ! RRTMGP downward clear-sky flux profiles (W/m2)

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
    real(kind_phys), dimension(nday,NLev+1,sw_gas_props%get_nband()),target :: &
         fluxSW_up_allsky, fluxSW_up_clrsky, fluxSW_dn_allsky, fluxSW_dn_clrsky, fluxSW_dn_dir_allsky
    real(kind_phys), dimension(ncol,NLev) :: vmrTemp
    logical :: l_scmpsw=.false., top_at_1
    integer :: iGas,iSFC,iTOA,iBand

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg  = 0

    if (.not. doSWrad) return

    ! Initialize output fluxes
    fluxswUP_allsky(:,:)   = 0._kind_phys
    fluxswDOWN_allsky(:,:) = 0._kind_phys
    fluxswUP_clrsky(:,:)   = 0._kind_phys
    fluxswDOWN_clrsky(:,:) = 0._kind_phys

    if (nDay .gt. 0) then

       ! Vertical ordering?
       top_at_1 = (p_lev(1,1) .lt. p_lev(1, NLev))
       if (top_at_1) then 
          iSFC = NLev+1
          iTOA = 1
       else
          iSFC = 1
          iTOA = NLev+1
       endif
       
       ! Are any optional outputs requested? Need to know now to compute correct fluxes.
       l_scmpsw           = present(scmpsw)
       if ( l_scmpsw ) then
          scmpsw = cmpfsw_type (0., 0., 0., 0., 0., 0.)
       endif
       
       ! Initialize RRTMGP DDT containing 2D(3D) fluxes
       fluxSW_up_allsky(:,:,:)     = 0._kind_phys
       fluxSW_dn_allsky(:,:,:)     = 0._kind_phys
       fluxSW_dn_dir_allsky(:,:,:) = 0._kind_phys
       fluxSW_up_clrsky(:,:,:)     = 0._kind_phys
       fluxSW_dn_clrsky(:,:,:)     = 0._kind_phys
       flux_allsky%bnd_flux_up     => fluxSW_up_allsky
       flux_allsky%bnd_flux_dn     => fluxSW_dn_allsky
       flux_allsky%bnd_flux_dn_dir => fluxSW_dn_dir_allsky
       flux_clrsky%bnd_flux_up     => fluxSW_up_clrsky
       flux_clrsky%bnd_flux_dn     => fluxSW_dn_clrsky

       !  *Note* Legacy RRTMG code. May need to revisit
       do iBand=1,sw_gas_props%get_nband()
          if (iBand .lt. 10) then
             sfc_alb_dir(iBand,:) = sfc_alb_nir_dir(iBand,idxday(1:nday))
             sfc_alb_dif(iBand,:) = sfc_alb_nir_dif(iBand,idxday(1:nday))
          endif
          if (iBand .eq. 10) then
             sfc_alb_dir(iBand,:) = 0.5_kind_phys*(sfc_alb_nir_dir(iBand,idxday(1:nday)) + sfc_alb_uvvis_dir(iBand,idxday(1:nday)))
             sfc_alb_dif(iBand,:) = 0.5_kind_phys*(sfc_alb_nir_dif(iBand,idxday(1:nday)) + sfc_alb_uvvis_dif(iBand,idxday(1:nday)))
          endif
          if (iBand .gt. 10) then
             sfc_alb_dir(iBand,:) = sfc_alb_uvvis_dir(iBand,idxday(1:nday))
             sfc_alb_dif(iBand,:) = sfc_alb_uvvis_dif(iBand,idxday(1:nday))
          endif
       enddo

       ! Compute clear-sky fluxes (if requested)
       ! Clear-sky fluxes (gas+aerosol)
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_aerosol%increment(sw_optical_props_clrsky))
       ! Delta-scale optical properties
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clrsky%delta_scale())
       if (doSWclrsky) then
          call check_error_msg('rrtmgp_sw_rte_run',rte_sw(     &
               sw_optical_props_clrsky,      & ! IN  - optical-properties
               top_at_1,                     & ! IN  - veritcal ordering flag
               coszen(idxday(1:nday)),       & ! IN  - Cosine of solar zenith angle
               toa_src_sw(idxday(1:nday),:), & ! IN  - incident solar flux at TOA
               sfc_alb_dir,                  & ! IN  - Shortwave surface albedo (direct)
               sfc_alb_dif,                  & ! IN  - Shortwave surface albedo (diffuse)
               flux_clrsky))                   ! OUT - Fluxes, clear-sky, 3D (nCol,NLev,nBand) 
          ! Store fluxes
          fluxswUP_clrsky(idxday(1:nday),:)   = sum(flux_clrsky%bnd_flux_up,dim=3)
          fluxswDOWN_clrsky(idxday(1:nday),:) = sum(flux_clrsky%bnd_flux_dn,dim=3)
       endif
       
       ! Compute all-sky fluxes
       ! All-sky fluxes (clear-sky + clouds)
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clouds%increment(sw_optical_props_clrsky))
       ! Delta-scale optical properties
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_clrsky%delta_scale())
       call check_error_msg('rrtmgp_sw_rte_run',rte_sw(     &
            sw_optical_props_clrsky,      & ! IN  - optical-properties
            top_at_1,                     & ! IN  - veritcal ordering flag
            coszen(idxday(1:nday)),       & ! IN  - Cosine of solar zenith angle
            toa_src_sw(idxday(1:nday),:), & ! IN  - incident solar flux at TOA
            sfc_alb_dir,                  & ! IN  - Shortwave surface albedo (direct)
            sfc_alb_dif,                  & ! IN  - Shortwave surface albedo (diffuse)
            flux_allsky))                   ! OUT - Fluxes, clear-sky, 3D (nCol,NLev,nBand) 
       ! Store fluxes
       fluxswUP_allsky(idxday(1:nday),:)   = sum(flux_allsky%bnd_flux_up,dim=3)
       fluxswDOWN_allsky(idxday(1:nday),:) = sum(flux_allsky%bnd_flux_dn,dim=3)
       if ( l_scmpsw ) then
          scmpsw(idxday(1:nday))%nirbm = sum(flux_allsky%bnd_flux_dn_dir(1:nday,iSFC,:),dim=2)
          scmpsw(idxday(1:nday))%nirdf = sum(flux_allsky%bnd_flux_dn(1:nday,iSFC,:),dim=2)  - &
               sum(flux_allsky%bnd_flux_dn_dir(1:nday,iSFC,:),dim=2)
       endif
    endif
  end subroutine rrtmgp_sw_rte_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_rte_finalize
  ! #########################################################################################
  subroutine rrtmgp_sw_rte_finalize()
  end subroutine rrtmgp_sw_rte_finalize

end module rrtmgp_sw_rte
