module rrtmgp_sw_rte
  use machine,                 only: kind_phys
  use mo_optical_props,        only: ty_optical_props_2str
  use mo_rte_sw,               only: rte_sw
  use mo_fluxes_byband,        only: ty_fluxes_byband
  use module_radsw_parameters, only: cmpfsw_type
  use radiation_tools,         only: check_error_msg
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
  subroutine rrtmgp_sw_rte_run(doSWrad, doSWclrsky, nCol, nLev, nDay, idxday, coszen, p_lay,&
       t_lay, top_at_1, doGP_sgs_cnv, doGP_sgs_mynn, iSFC, sfc_alb_nir_dir, sfc_alb_nir_dif,&
       sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, toa_src_sw, sw_optical_props_clrsky,           &
       sw_optical_props_clouds, sw_optical_props_precipByBand,                              &
       sw_optical_props_cnvcloudsByBand, sw_optical_props_MYNNcloudsByBand,                 &
       sw_optical_props_aerosol, scmpsw, fluxswUP_allsky, fluxswDOWN_allsky,                &
       fluxswUP_clrsky, fluxswDOWN_clrsky, errmsg, errflg)
    
    ! Inputs
    logical, intent(in) :: &
         top_at_1,                          & ! Vertical ordering flag
         doGP_sgs_mynn,                     & ! Flag for MYNN-EDMF PBL cloud scheme
         doGP_sgs_cnv,                      & ! Flag for sgs convective clouds scheme
         doSWrad,                           & ! Flag to calculate SW irradiances
         doSWclrsky                           ! Compute clear-sky fluxes?
    integer, intent(in) :: &
         nCol,                              & ! Number of horizontal gridpoints
         nday,                              & ! Number of daytime points
         nLev,                              & ! Number of vertical levels
         iSFC                                 ! Vertical index for surface-level
    integer, intent(in), dimension(:) :: &
         idxday                               ! Index array for daytime points
    real(kind_phys),intent(in), dimension(:) :: &
         sfc_alb_nir_dir,                   & ! Surface albedo (direct)
         sfc_alb_nir_dif,                   & ! Surface albedo (diffuse)
         sfc_alb_uvvis_dir,                 & ! Surface albedo (direct)
         sfc_alb_uvvis_dif,                 & ! Surface albedo (diffuse)
         coszen                               ! Cosize of SZA
    real(kind_phys), dimension(:,:), intent(in) :: &
         p_lay,                             & ! Pressure @ model layer-centers (Pa)
         t_lay,                             & ! Temperature (K)
         toa_src_sw                           ! TOA incident spectral flux (W/m2)
    type(ty_optical_props_2str),intent(inout) :: &
         sw_optical_props_clrsky              ! RRTMGP DDT: shortwave clear-sky radiative properties 
    type(ty_optical_props_2str),intent(in) :: &
         sw_optical_props_clouds,           & ! RRTMGP DDT: shortwave cloud optical properties 
         sw_optical_props_cnvcloudsByBand,  & ! RRTMGP DDT: shortwave convecive cloud optical properties
         sw_optical_props_MYNNcloudsByBand, & ! RRTMGP DDT: shortwave MYNN-EDMF PBL cloud optical properties
         sw_optical_props_precipByBand,     & ! RRTMGP DDT: shortwave precipitation optical properties
         sw_optical_props_aerosol             ! RRTMGP DDT: shortwave aerosol optical properties

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                     ! CCPP error message
    integer, intent(out) :: &
         errflg                     ! CCPP error flag
    real(kind_phys), dimension(:,:), intent(inout) :: &
         fluxswUP_allsky,         & ! RRTMGP upward all-sky flux profiles (W/m2)
         fluxswDOWN_allsky,       & ! RRTMGP downward all-sky flux profiles (W/m2)
         fluxswUP_clrsky,         & ! RRTMGP upward clear-sky flux profiles (W/m2)
         fluxswDOWN_clrsky          ! RRTMGP downward clear-sky flux profiles (W/m2)
    type(cmpfsw_type), dimension(:), intent(inout) :: &
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
    integer :: iBand, iDay,ibd
    real(kind_phys), dimension(2,sw_gas_props%get_nband()) :: bandlimits
    real(kind_phys), dimension(2), parameter :: nIR_uvvis_bnd = (/12850,16000/)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg  = 0

    if (.not. doSWrad) return

    if (nDay .gt. 0) then

       ! Initialize RRTMGP DDT containing 2D(3D) fluxes
       flux_allsky%bnd_flux_up     => fluxSW_up_allsky
       flux_allsky%bnd_flux_dn     => fluxSW_dn_allsky
       flux_allsky%bnd_flux_dn_dir => fluxSW_dn_dir_allsky
       flux_clrsky%bnd_flux_up     => fluxSW_up_clrsky
       flux_clrsky%bnd_flux_dn     => fluxSW_dn_clrsky

       ! Use near-IR albedo for bands with wavenumbers extending to 12850cm-1
       ! Use uv-vis albedo for bands with wavenumbers greater than 16000cm-1
       ! For overlapping band, average near-IR and us-vis albedos.
       bandlimits = sw_gas_props%get_band_lims_wavenumber()
       do iBand=1,sw_gas_props%get_nband()
          if (bandlimits(1,iBand) .lt. nIR_uvvis_bnd(1)) then
             sfc_alb_dir(iBand,:) = sfc_alb_nir_dir(idxday(1:nday))
             sfc_alb_dif(iBand,:) = sfc_alb_nir_dif(idxday(1:nday))
          endif
          if (bandlimits(1,iBand) .eq. nIR_uvvis_bnd(1)) then
             sfc_alb_dir(iBand,:) = 0.5_kind_phys*(sfc_alb_nir_dir(idxday(1:nday)) + sfc_alb_uvvis_dir(idxday(1:nday)))
             sfc_alb_dif(iBand,:) = 0.5_kind_phys*(sfc_alb_nir_dif(idxday(1:nday)) + sfc_alb_uvvis_dif(idxday(1:nday)))
             ibd = iBand
          endif
          if (bandlimits(1,iBand) .ge. nIR_uvvis_bnd(2)) then
             sfc_alb_dir(iBand,:) = sfc_alb_uvvis_dir(idxday(1:nday))
             sfc_alb_dif(iBand,:) = sfc_alb_uvvis_dif(idxday(1:nday))
          endif
       enddo

       !
       ! Compute clear-sky fluxes (if requested)
       !

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

       !
       ! Compute all-sky fluxes
       !

       ! Include convective cloud?
       if (doGP_sgs_cnv) then
          call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_cnvcloudsByBand%increment(sw_optical_props_clrsky))
       endif

       ! Include MYNN-EDMF PBL cloud?
       if (doGP_sgs_mynn) then
          call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_MYNNcloudsByBand%increment(sw_optical_props_clrsky))
       endif

       ! All-sky fluxes (clear-sky + clouds + precipitation)
       call check_error_msg('rrtmgp_sw_rte_run',sw_optical_props_precipByBand%increment(sw_optical_props_clrsky))
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
       do iDay=1,nDay
          ! Near IR
          scmpsw(idxday(iDay))%nirbm = sum(flux_allsky%bnd_flux_dn_dir(iDay,iSFC,1:ibd-1))  + &
                                           flux_allsky%bnd_flux_dn_dir(iDay,iSFC,ibd)/2.
          scmpsw(idxday(iDay))%nirdf = (sum(flux_allsky%bnd_flux_dn(iDay,iSFC,1:ibd-1))     + &
                                            flux_allsky%bnd_flux_dn(iDay,iSFC,ibd)/2.)      - &
                                       (sum(flux_allsky%bnd_flux_dn_dir(iDay,iSFC,1:ibd-1)) + &
                                            flux_allsky%bnd_flux_dn_dir(iDay,iSFC,ibd)/2.)
          ! UV-VIS
          scmpsw(idxday(iDay))%visbm = sum(flux_allsky%bnd_flux_dn_dir(iDay,iSFC,ibd+1:sw_gas_props%get_nband()))  + &
                                           flux_allsky%bnd_flux_dn_dir(iDay,iSFC,ibd)/2.
          scmpsw(idxday(iDay))%visdf = (sum(flux_allsky%bnd_flux_dn(iDay,iSFC,ibd+1:sw_gas_props%get_nband()))     + &
                                            flux_allsky%bnd_flux_dn(iDay,iSFC,ibd)/2. )                            - &
                                       (sum(flux_allsky%bnd_flux_dn_dir(iDay,iSFC,ibd+1:sw_gas_props%get_nband())) + &
                                            flux_allsky%bnd_flux_dn_dir(iDay,iSFC,ibd)/2.)
       enddo
    else
       fluxswUP_allsky(:,:)   = 0._kind_phys
       fluxswDOWN_allsky(:,:) = 0._kind_phys
       fluxswUP_clrsky(:,:)   = 0._kind_phys
       fluxswDOWN_clrsky(:,:) = 0._kind_phys
       scmpsw                 = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )       
    endif

  end subroutine rrtmgp_sw_rte_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_rte_finalize
  ! #########################################################################################
  subroutine rrtmgp_sw_rte_finalize()
  end subroutine rrtmgp_sw_rte_finalize

end module rrtmgp_sw_rte
