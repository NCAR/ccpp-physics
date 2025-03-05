!> \file GFS_suite_interstitial_rad_reset.f90
!! ##########################################################################################
!!
!!  Contains code to reset radiation-related interstitial variables in the GFS physics suite. 
!!
!! ##########################################################################################
module GFS_suite_interstitial_rad_reset
  use machine, only: kind_phys
  use module_radsw_parameters, only: profsw_type, cmpfsw_type
  use module_radlw_parameters, only: proflw_type
  implicit none
contains
!> \section arg_table_GFS_suite_interstitial_rad_reset_run Argument Table
!! \htmlinclude GFS_suite_interstitial_rad_reset_run.html
!!
  subroutine GFS_suite_interstitial_rad_reset_run (do_RRTMGP, clear_val, aerodp, alb1d,     &
       alpha, cldsa, cldtaulw, cldtausw, cloud_frac, cloud_lwp, cloud_reliq, cloud_iwp,     &
       cloud_reice, cloud_rwp, cloud_rerain, cloud_swp, cloud_resnow, de_lgth, delr, dzlyr, &
       faerlw, faersw, gasvmr_co2, gasvmr_n2o, gasvmr_ch4, gasvmr_o2, gasvmr_co,            &
       gasvmr_cfc11, gasvmr_cfc12, gasvmr_cfc22, gasvmr_ccl4, gasvmr_cfc113, htlwc, htlw0,  &
       htswc, htsw0, idxday, kb, kd, kt, mbota, mtopa, nday, olyr, plvl, plyr, qlyr, raddt, &
       sfcalb_nIR_dir, sfcalb_nIR_dif, sfcalb_uvvis_dir, sfcalb_uvvis_dif, tlvl, tlyr, tsfa,&
       tsfg, scmpsw, tracer, tv_lay, relhum, qs_lay, q_lay, deltaZ, deltaZc, deltaP, p_lev, &
       p_lay, t_lev, t_lay, cloud_overlap_param, precip_overlap_param, fluxlwUP_allsky,     &
       fluxlwDOWN_allsky, fluxlwUP_clrsky, fluxlwDOWN_clrsky, fluxswUP_allsky,              &
       fluxswDOWN_allsky, fluxswUP_clrsky, fluxswDOWN_clrsky, aerosolslw, aerosolssw,       &
       precip_frac, cld_cnv_frac, cnv_cloud_overlap_param, cld_cnv_lwp, cld_cnv_reliq,      &
       cld_cnv_iwp, cld_cnv_reice, cld_pbl_lwp, cld_pbl_reliq, cld_pbl_iwp, cld_pbl_reice,  &
       sfc_emiss_byband, sec_diff_byband, sfc_alb_nir_dir, sfc_alb_nir_dif,                 &
       sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, toa_src_sw, toa_src_lw, vmr_o2, vmr_h2o,       &
       vmr_o3, vmr_ch4, vmr_n2o, vmr_co2, flxprf_lw, flxprf_sw, errmsg, errflg)
    !
    logical,            intent(in   ) :: do_RRTMGP
    real(kind_phys),    intent(in   ) :: clear_val
    real (kind_phys),   intent(inout) :: aerodp(:,:), alb1d(:), cldsa(:,:),                 &
         cldtaulw(:,:), cldtausw(:,:), cloud_frac(:,:), cloud_lwp(:,:), cloud_reliq(:,:),   &
         cloud_iwp(:,:), cloud_reice(:,:), cloud_rwp(:,:), cloud_rerain(:,:),               &
         cloud_swp(:,:), cloud_resnow(:,:), de_lgth(:), delr(:,:), dzlyr(:,:),              &
         faerlw(:,:,:,:), faersw(:,:,:,:), gasvmr_co2(:,:), gasvmr_n2o(:,:),                &
         gasvmr_ch4(:,:), gasvmr_o2(:,:), gasvmr_co(:,:), gasvmr_cfc11(:,:),                &
         gasvmr_cfc12(:,:), gasvmr_cfc22(:,:), gasvmr_ccl4(:,:), gasvmr_cfc113(:,:),        &
         htlwc(:,:), htlw0(:,:), htswc(:,:), htsw0(:,:), olyr(:,:), plvl(:,:), plyr(:,:),   &
         qlyr(:,:), raddt, sfcalb_nIR_dir(:), sfcalb_nIR_dif(:), sfcalb_uvvis_dir(:),       &
         sfcalb_uvvis_dif(:), tlvl(:,:), tlyr(:,:), tsfa(:), tsfg(:)
    type (cmpfsw_type), intent(inout) :: scmpsw(:)
    real(kind_phys),    intent(inout), optional :: alpha(:,:), tracer(:,:,:),               &
         tv_lay(:,:), relhum(:,:), qs_lay(:,:), q_lay(:,:), deltaZ(:,:), deltaZc(:,:),      &
         deltaP(:,:), p_lev(:,:), p_lay(:,:), t_lev(:,:), t_lay(:,:),                       &
         cloud_overlap_param(:,:), precip_overlap_param(:,:), fluxlwUP_allsky(:,:),         &
         fluxlwDOWN_allsky(:,:), fluxlwUP_clrsky(:,:), fluxlwDOWN_clrsky(:,:),              &
         fluxswUP_allsky(:,:), fluxswDOWN_allsky(:,:), fluxswUP_clrsky(:,:),                &
         fluxswDOWN_clrsky(:,:), aerosolslw(:,:,:,:), aerosolssw(:,:,:,:), precip_frac(:,:),&
         cld_cnv_frac(:,:), cnv_cloud_overlap_param(:,:), cld_cnv_lwp(:,:),                 &
         cld_cnv_reliq(:,:), cld_cnv_iwp(:,:), cld_cnv_reice(:,:), cld_pbl_lwp(:,:),        &
         cld_pbl_reliq(:,:), cld_pbl_iwp(:,:), cld_pbl_reice(:,:), sfc_emiss_byband(:,:),   &
         sec_diff_byband(:,:), sfc_alb_nir_dir(:,:), sfc_alb_nir_dif(:,:),                  &
         sfc_alb_uvvis_dir(:,:), sfc_alb_uvvis_dif(:,:), toa_src_sw(:,:), toa_src_lw(:,:),  &
         vmr_o2(:,:), vmr_h2o(:,:), vmr_o3(:,:), vmr_ch4(:,:), vmr_n2o(:,:), vmr_co2(:,:)
    type(proflw_type),  intent(inout) :: flxprf_lw(:,:)
    type(profsw_type),  intent(inout) :: flxprf_sw(:,:)
    integer,            intent(inout) :: idxday(:), kb, kd, kt, mbota(:,:), mtopa(:,:), nday
    character(len=*),   intent(out  ) :: errmsg
    integer,            intent(out  ) :: errflg
    
    ! Initialize CCPP error logging
    errmsg = ''
    errflg = 0

    ! Reset fields
    aerodp           = clear_val
    alb1d            = clear_val
    if (.not. do_RRTMGP) then
       alpha         = clear_val
    end if
    cldsa            = clear_val
    cldtaulw         = clear_val
    cldtausw         = clear_val
    cloud_frac       = clear_val
    cloud_lwp        = clear_val
    cloud_reliq      = clear_val
    cloud_iwp        = clear_val
    cloud_reice      = clear_val
    cloud_rwp        = clear_val
    cloud_rerain     = clear_val
    cloud_swp        = clear_val
    cloud_resnow     = clear_val
    de_lgth          = clear_val
    delr             = clear_val
    dzlyr            = clear_val
    faerlw           = clear_val
    faersw           = clear_val
    gasvmr_co2       = clear_val
    gasvmr_n2o       = clear_val
    gasvmr_ch4       = clear_val
    gasvmr_o2        = clear_val
    gasvmr_co        = clear_val
    gasvmr_cfc11     = clear_val
    gasvmr_cfc12     = clear_val
    gasvmr_cfc22     = clear_val
    gasvmr_ccl4      = clear_val
    gasvmr_cfc113    = clear_val
    htlwc            = clear_val
    htlw0            = clear_val
    htswc            = clear_val
    htsw0            = clear_val
    idxday           = 0
    kb               = 0
    kd               = 0
    kt               = 0
    mbota            = 0
    mtopa            = 0
    nday             = 0
    olyr             = clear_val
    plvl             = clear_val
    plyr             = clear_val
    qlyr             = clear_val
    raddt            = clear_val
    sfcalb_nIR_dir   = clear_val
    sfcalb_nIR_dif   = clear_val
    sfcalb_uvvis_dir = clear_val
    sfcalb_uvvis_dif = clear_val
    tlvl             = clear_val
    tlyr             = clear_val
    tsfa             = clear_val
    tsfg             = clear_val
    
    ! Shortwave surface fluxes.
    scmpsw%uvbfc     = clear_val
    scmpsw%uvbf0     = clear_val
    scmpsw%nirbm     = clear_val
    scmpsw%nirdf     = clear_val
    scmpsw%visbm     = clear_val
    scmpsw%visdf     = clear_val

    ! RRMTGP
    if (do_RRTMGP) then
       tracer               = clear_val
       tv_lay               = clear_val
       relhum               = clear_val
       qs_lay               = clear_val
       q_lay                = clear_val
       deltaZ               = clear_val
       deltaZc              = clear_val
       deltaP               = clear_val
       p_lev                = clear_val
       p_lay                = clear_val
       t_lev                = clear_val
       t_lay                = clear_val
       cloud_overlap_param  = clear_val
       precip_overlap_param = clear_val
       fluxlwUP_allsky      = clear_val
       fluxlwDOWN_allsky    = clear_val
       fluxlwUP_clrsky      = clear_val
       fluxlwDOWN_clrsky    = clear_val
       fluxswUP_allsky      = clear_val
       fluxswDOWN_allsky    = clear_val
       fluxswUP_clrsky      = clear_val
       fluxswDOWN_clrsky    = clear_val
       aerosolslw           = clear_val
       aerosolssw           = clear_val
       precip_frac          = clear_val
       cld_cnv_frac         = clear_val
       cnv_cloud_overlap_param  = clear_val
       cld_cnv_lwp          = clear_val
       cld_cnv_reliq        = clear_val
       cld_cnv_iwp          = clear_val
       cld_cnv_reice        = clear_val
       cld_pbl_lwp          = clear_val
       cld_pbl_reliq        = clear_val
       cld_pbl_iwp          = clear_val
       cld_pbl_reice        = clear_val
       sfc_emiss_byband     = clear_val
       sec_diff_byband      = clear_val
       sfc_alb_nir_dir      = clear_val
       sfc_alb_nir_dif      = clear_val
       sfc_alb_uvvis_dir    = clear_val
       sfc_alb_uvvis_dif    = clear_val
       toa_src_sw           = clear_val
       toa_src_lw           = clear_val
       vmr_o2               = clear_val
       vmr_h2o              = clear_val
       vmr_o3               = clear_val
       vmr_ch4              = clear_val
       vmr_n2o              = clear_val
       vmr_co2              = clear_val
       flxprf_lw%upfxc      = clear_val
       flxprf_lw%dnfxc      = clear_val
       flxprf_lw%upfx0      = clear_val
       flxprf_lw%dnfx0      = clear_val
       flxprf_sw%upfxc      = clear_val
       flxprf_sw%dnfxc      = clear_val
       flxprf_sw%upfx0      = clear_val
       flxprf_sw%dnfx0      = clear_val
    end if
    !
  end subroutine GFS_suite_interstitial_rad_reset_run
  !
end module GFS_suite_interstitial_rad_reset
