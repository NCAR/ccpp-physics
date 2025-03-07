!> \file GFS_suite_interstitial_phys_reset.F90
!!
!! ##########################################################################################
!!
!!  Contains code to reset physics-related interstitial variables in the GFS physics suite.
!!
!! ##########################################################################################
module GFS_suite_interstitial_phys_reset
  use machine, only: kind_phys
  implicit none
contains
!> \section arg_table_GFS_suite_interstitial_phys_reset_run Argument Table
!! \htmlinclude GFS_suite_interstitial_phys_reset_run.html
!!
  subroutine GFS_suite_interstitial_phys_reset_run(clear_val, huge, levs, kdt, ldiag_ugwp, &
       do_ugwp_v0, do_ugwp_v0_nst_only, do_ugwp_v1, gwd_opt, imp_physics,                  &
       imp_physics_gfdl, imp_physics_thompson, imp_physics_nssl, imp_physics_mg, lsm,      &
       ilsm_noahmp, dtp, avg_max_length, nsfullradar_diag,                                 &
       adjsfculw_land, adjsfculw_ice, adjsfculw_water, adjnirbmd, adjnirbmu, adjnirdfd,    &
       adjnirdfu, adjvisbmd, adjvisbmu, adjvisdfu, adjvisdfd, bexp1d, cd, cd_ice, cd_land, &
       cd_water, cdq, cdq_ice, cdq_land, cdq_water, chh_ice, chh_land, chh_water, cld1d,   &
       cldf, clw, clx, cmm_ice, cmm_land, cmm_water, cnvc, cnvw, ctei_r, ctei_rml, cumabs, &
       dd_mf, del, del_gz, dlength, dqdt, dqsfc1, drain, dt_mf, dtdt, dtsfc1, dtzm, dudt,  &
       dusfcg, dusfc1, dvdftra, dvdt, dvsfcg, dvsfc1, elvmax, ep1d, ep1d_ice, ep1d_land,   &
       ep1d_water, evap_ice, evap_land, evap_water, pah, ecan, etran, edir,                &
       ffhh_ice, ffhh_land, ffhh_water, fh2, fh2_ice, fh2_land, fh2_water, flag_cice,      &
       flag_guess, flag_iter, flag_lakefreeze, ffmm_ice, ffmm_land, ffmm_water, fm10,      &
       fm10_ice, fm10_land, fm10_water, frland, fscav, fswtr, gabsbdlw, gabsbdlw_ice,      &
       gabsbdlw_land, gabsbdlw_water, gamma, gamq, gamt, gflx, gflx_ice, gflx_land,        &
       gflx_water, gwdcu, gwdcv, zvfun, hffac, hflxq, hflx_ice,  hflx_land, hflx_water,    &
       dry, icy, lake, ocean, islmsk, islmsk_cice, wet, kbot, kcnv, kinver, kpbl, ktop,    &
       oa4, oc, prcpmp, prnum, qss_ice, qss_land, qss_water, raincd, raincs, rainmcadj,    &
       rainp, rb, rb_ice, rb_land, rb_water, rhc, runoff, save_q, save_t, save_tcp,        &
       save_u, save_v, sigma, sigmaf, sigmafrac, sigmatot, snowc, snohf, snowmt,           &
       stress, stress_ice, stress_land, stress_water, theta, tprcp_ice, tprcp_land,        &
       tprcp_water, tseal, tsfc_water, tsurf_ice, tsurf_land, tsurf_water,                 &
       uustar_ice, uustar_land, uustar_water, vdftra, vegf1d, lndp_vgf, wcbmax, wind,      &
       work1, work2, work3, xcosz, xlai1d, xmu, z01d, zt1d, ztmax_ice, ztmax_land,         &
       ztmax_water, tau_mtb, tau_ogw, tau_tofd, tau_ngw, tau_oss, dudt_mtb, dudt_tms,      &
       zmtb, zlwb, zogw, zngw, dudt_ngw, dvdt_ngw, dtdt_ngw, kdis_ngw, varss, ocss, oa4ss, &
       clxss, graupelmp, icemp, rainmp, snowmp, ncgl, ncpr, ncps, qgl, qrn, qsnw, qlcn,    &
       qicn, w_upi, cf_upi, cnv_mfd, cnv_dqldt, clcn, cnv_fice, cnv_ndrop, cnv_nice,       &
       t2mmp, q2mp, max_hourly_reset, ext_diag_thompson_reset, fullradar_diag,             &
       errmsg, errflg)

    ! Inputs
    real (kind_phys), intent(in   ) :: clear_val, huge, avg_max_length, nsfullradar_diag,  &
         dtp
    integer,          intent(in   ) :: levs, kdt, gwd_opt, imp_physics, imp_physics_gfdl,  &
         imp_physics_thompson, imp_physics_nssl, imp_physics_mg, lsm, ilsm_noahmp
    logical,          intent(in   ) :: ldiag_ugwp, do_ugwp_v0, do_ugwp_v0_nst_only,        &
         do_ugwp_v1
    ! Outputs
    real(kind_phys),  intent(inout) :: adjsfculw_land(:), adjsfculw_ice(:),                &
         adjsfculw_water(:), adjnirbmd(:), adjnirbmu(:), adjnirdfd(:), adjnirdfu(:),       &
         adjvisbmd(:), adjvisbmu(:), adjvisdfu(:), adjvisdfd(:), bexp1d(:), cd(:),         &
         cd_ice(:), cd_land(:), cd_water(:), cdq(:), cdq_ice(:), cdq_land(:), cdq_water(:),&
         chh_ice(:), chh_land(:), chh_water(:), cldf(:), cld1d(:), clw(:,:,:),             &
         clx(:,:), cmm_ice(:), cmm_land(:), cmm_water(:), cnvc(:,:), cnvw(:,:), ctei_r(:), &
         ctei_rml(:), cumabs(:), dd_mf(:,:), del(:,:),  del_gz(:,:), dlength(:),           &
         dqdt(:,:,:), dqsfc1(:), drain(:), dtdt(:,:), dtsfc1(:), dtzm(:), dt_mf(:,:),      &
         dudt(:,:), dusfcg(:), dusfc1(:), dvdftra(:,:,:), dvdt(:,:), dvsfcg(:), dvsfc1(:), &
         elvmax(:), ep1d(:), ep1d_ice(:), ep1d_land(:), ep1d_water(:), evap_ice(:),        &
         evap_land(:), evap_water(:), pah(:), ecan(:), etran(:), edir(:),                  &
         ffhh_ice(:), ffhh_land(:), ffhh_water(:), fh2(:), fh2_ice(:), fh2_land(:),        &
         fh2_water(:), ffmm_ice(:), ffmm_land(:), ffmm_water(:), fm10(:), fm10_ice(:),     &
         fm10_land(:), fm10_water(:), frland(:), fscav(:), fswtr(:), gabsbdlw(:),          &
         gabsbdlw_ice(:), gabsbdlw_land(:), gabsbdlw_water(:), gamma(:), gamq(:), gamt(:), &
         gflx(:), gflx_ice(:), gflx_land(:), gflx_water(:), gwdcu(:,:), gwdcv(:,:),        &
         zvfun(:), hffac(:), hflxq(:), hflx_ice(:), hflx_land(:), hflx_water(:), oa4(:,:), &
         oc(:), prcpmp(:), prnum(:,:), qss_ice(:), qss_land(:), qss_water(:), raincd(:),   &
         raincs(:), rainmcadj(:), rainp(:,:), rb(:), rb_ice(:), rb_land(:), rb_water(:),   &
         rhc(:,:), runoff(:), save_q(:,:,:), save_t(:,:), save_tcp(:,:), save_u(:,:),      &
         save_v(:,:), sigma(:), sigmaf(:), sigmafrac(:,:), sigmatot(:,:),                  &
         snowc(:), snohf(:), snowmt(:), stress(:), stress_ice(:), stress_land(:),          &
         stress_water(:), theta(:), tprcp_ice(:), tprcp_land(:), tprcp_water(:),           &
         tseal(:), tsfc_water(:), tsurf_ice(:), tsurf_land(:), tsurf_water(:),             &
         uustar_ice(:), uustar_land(:), uustar_water(:), vdftra(:,:,:),                    &
         vegf1d(:), lndp_vgf, wcbmax(:), wind(:), work1(:), work2(:), work3(:), xcosz(:),  &
         xlai1d(:), xmu(:), z01d(:), zt1d(:), ztmax_ice(:), ztmax_land(:), ztmax_water(:), &
         tau_oss(:), tau_tofd(:), tau_mtb(:), tau_ogw(:), tau_ngw(:), zngw(:), zmtb(:),    &
         zlwb(:), zogw(:), dudt_mtb(:,:), dudt_tms(:,:)
    real(kind_phys),  intent(inout), optional :: clcn(:,:), cnv_dqldt(:,:), cnv_fice(:,:), &
         cnv_mfd(:,:), cnv_ndrop(:,:), cnv_nice(:,:), graupelmp(:), icemp(:), ncgl(:,:),   &
         ncpr(:,:), ncps(:,:), q2mp(:), qgl(:,:), qicn(:,:), qlcn(:,:), qrn(:,:),          &
         qsnw(:,:), rainmp(:), snowmp(:), t2mmp(:), w_upi(:,:), dudt_ngw(:,:),             &
         dvdt_ngw(:,:), dtdt_ngw(:,:) ,kdis_ngw(:,:), varss(:), ocss(:), oa4ss(:,:),       &
         clxss(:,:), cf_upi(:,:)
    logical, intent(inout) :: flag_cice(:), flag_guess(:), flag_iter(:),flag_lakefreeze(:),&
         dry(:), icy(:), lake(:), ocean(:), wet(:), fullradar_diag, max_hourly_reset,      &
         ext_diag_thompson_reset
    integer, intent(inout) :: islmsk(:), islmsk_cice(:), kbot(:), kcnv(:), kinver(:),      &
         kpbl(:), ktop(:)
    !
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    errmsg = ''
    errflg = 0

    adjsfculw_land  = clear_val
    adjsfculw_ice   = clear_val
    adjsfculw_water = clear_val
    adjnirbmd       = clear_val
    adjnirbmu       = clear_val
    adjnirdfd       = clear_val
    adjnirdfu       = clear_val
    adjvisbmd       = clear_val
    adjvisbmu       = clear_val
    adjvisdfu       = clear_val
    adjvisdfd       = clear_val
    bexp1d          = clear_val
    cd              = clear_val
    cd_ice          = huge
    cd_land         = huge
    cd_water        = huge
    cdq             = clear_val
    cdq_ice         = huge
    cdq_land        = huge
    cdq_water       = huge
    chh_ice         = huge
    chh_land        = huge
    chh_water       = huge
    cld1d           = clear_val
    cldf            = clear_val
    clw             = clear_val
    clw(:,:,2)      = -999.9
    clx             = clear_val
    cmm_ice         = huge
    cmm_land        = huge
    cmm_water       = huge
    cnvc            = clear_val
    cnvw            = clear_val
    ctei_r          = clear_val
    ctei_rml        = clear_val
    cumabs          = clear_val
    dd_mf           = clear_val
    del             = clear_val
    del_gz          = clear_val
    dlength         = clear_val
    dqdt            = clear_val
    dqsfc1          = clear_val
    drain           = clear_val
    dt_mf           = clear_val
    dtdt            = clear_val
    dtsfc1          = clear_val
    dtzm            = clear_val
    dudt            = clear_val
    dusfcg          = clear_val
    dusfc1          = clear_val
    dvdftra         = clear_val
    dvdt            = clear_val
    dvsfcg          = clear_val
    dvsfc1          = clear_val
    elvmax          = clear_val
    ep1d            = clear_val
    ep1d_ice        = huge
    ep1d_land       = huge
    ep1d_water      = huge
    evap_ice        = huge
    evap_land       = huge
    evap_water      = huge
    pah             = clear_val
    ecan            = clear_val
    etran           = clear_val
    edir            = clear_val
    ffhh_ice        = huge
    ffhh_land       = huge
    ffhh_water      = huge
    fh2             = clear_val
    fh2_ice         = huge
    fh2_land        = huge
    fh2_water       = huge
    flag_cice       = .false.
    flag_guess      = .false.
    flag_iter       = .true.
    flag_lakefreeze = .false.
    ffmm_ice        = huge
    ffmm_land       = huge
    ffmm_water      = huge
    fm10            = clear_val
    fm10_ice        = huge
    fm10_land       = huge
    fm10_water      = huge
    frland          = clear_val
    fscav           = clear_val
    fswtr           = clear_val
    gabsbdlw        = clear_val
    gabsbdlw_ice    = clear_val
    gabsbdlw_land   = clear_val
    gabsbdlw_water  = clear_val
    gamma           = clear_val
    gamq            = clear_val
    gamt            = clear_val
    gflx            = clear_val
    gflx_ice        = clear_val
    gflx_land       = clear_val
    gflx_water      = clear_val
    gwdcu           = clear_val
    gwdcv           = clear_val
    zvfun           = clear_val
    hffac           = clear_val
    hflxq           = clear_val
    hflx_ice        = huge
    hflx_land       = huge
    hflx_water      = huge
    dry             = .false.
    icy             = .false.
    lake            = .false.
    ocean           = .false.
    islmsk          = 0
    islmsk_cice     = 0
    wet             = .false.
    kbot            = levs
    kcnv            = 0
    kinver          = levs
    kpbl            = 0
    ktop            = 1
    oa4             = clear_val
    oc              = clear_val
    prcpmp          = clear_val
    prnum           = clear_val
    qss_ice         = huge
    qss_land        = huge
    qss_water       = huge
    raincd          = clear_val
    raincs          = clear_val
    rainmcadj       = clear_val
    rainp           = clear_val
    rb              = clear_val
    rb_ice          = huge
    rb_land         = huge
    rb_water        = huge
    rhc             = clear_val
    runoff          = clear_val
    save_q          = clear_val
    save_t          = clear_val
    save_tcp        = clear_val
    save_u          = clear_val
    save_v          = clear_val
    sigma           = clear_val
    sigmaf          = clear_val
    sigmafrac       = clear_val
    sigmatot        = clear_val
    snowc           = clear_val
    snohf           = clear_val
    snowmt          = clear_val
    stress          = clear_val
    stress_ice      = huge
    stress_land     = huge
    stress_water    = huge
    theta           = clear_val
    tprcp_ice       = huge
    tprcp_land      = huge
    tprcp_water     = huge
    tseal           = clear_val
    tsfc_water      = huge
    tsurf_ice       = huge
    tsurf_land      = huge
    tsurf_water     = huge
    uustar_ice      = huge
    uustar_land     = huge
    uustar_water    = huge
    vdftra          = clear_val
    vegf1d          = clear_val
    lndp_vgf        = clear_val
    wcbmax          = clear_val
    wind            = huge
    work1           = clear_val
    work2           = clear_val
    work3           = clear_val
    xcosz           = clear_val
    xlai1d          = clear_val
    xmu             = clear_val
    z01d            = clear_val
    zt1d            = clear_val
    ztmax_ice       = clear_val
    ztmax_land      = clear_val
    ztmax_water     = clear_val
      
    ! UGWP common
    tau_mtb         = clear_val
    tau_ogw         = clear_val
    tau_tofd        = clear_val
    tau_ngw         = clear_val
    tau_oss         = clear_val
    dudt_mtb        = clear_val
    dudt_tms        = clear_val
    zmtb            = clear_val
    zlwb            = clear_val
    zogw            = clear_val
    zngw            = clear_val

    ! CIRES UGWP v1
    if (ldiag_ugwp .or. do_ugwp_v0 .or. do_ugwp_v0_nst_only .or. do_ugwp_v1) then
       dudt_ngw        = clear_val
       dvdt_ngw        = clear_val
       dtdt_ngw        = clear_val
       kdis_ngw        = clear_val
    end if

    !-- GSL drag suite
    if (gwd_opt==3 .or. gwd_opt==33 .or. gwd_opt==2 .or. gwd_opt==22) then
       varss           = clear_val
       ocss            = clear_val
       oa4ss           = clear_val
       clxss           = clear_val
    end if

    ! Reset fields that are conditional on physics choices
    if (imp_physics == imp_physics_gfdl      .or. &
        imp_physics == imp_physics_thompson  .or. & 
        imp_physics == imp_physics_nssl ) then
       graupelmp = clear_val
       icemp     = clear_val
       rainmp    = clear_val
       snowmp    = clear_val
    else if (imp_physics == imp_physics_mg) then
       ncgl      = clear_val
       ncpr      = clear_val
       ncps      = clear_val
       qgl       = clear_val
       qrn       = clear_val
       qsnw      = clear_val
       qlcn      = clear_val
       qicn      = clear_val
       w_upi     = clear_val
       cf_upi    = clear_val
       cnv_mfd   = clear_val
       cnv_dqldt = clear_val
       clcn      = clear_val
       cnv_fice  = clear_val
       cnv_ndrop = clear_val
       cnv_nice  = clear_val
    end if
    if (lsm == ilsm_noahmp) then
       t2mmp     = clear_val
       q2mp      = clear_val
    end if

    ! Set flag for resetting maximum hourly output fields
    max_hourly_reset = mod(kdt-1, nint(avg_max_length/dtp)) == 0
    ! Use same logic in UFS to reset Thompson extended diagnostics
    ext_diag_thompson_reset = max_hourly_reset

    ! Frequency flag for computing the full radar reflectivity (water coated ice) 
    if (nsfullradar_diag<0) then
       fullradar_diag = .true.
    else
       fullradar_diag = (kdt == 1 .or. mod(kdt, nint(nsfullradar_diag/dtp)) == 0) 
    end if
      
  end subroutine GFS_suite_interstitial_phys_reset_run
    
end module GFS_suite_interstitial_phys_reset
