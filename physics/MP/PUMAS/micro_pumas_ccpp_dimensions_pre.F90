module micro_pumas_ccpp_dimensions_pre

  implicit none

contains

  !> \section arg_table_micro_pumas_ccpp_dimensions_pre_run Argument Table
  !! \htmlinclude micro_pumas_ccpp_dimensions_pre_run.html
  subroutine micro_pumas_ccpp_dimensions_pre_run(ncol, nlev, nlevp1,                &
                             micro_ncol, micro_nlev, micro_nlevp1,                  &
                             airT_in, micro_airT, airq_in, micro_airq,              &
                             cldliq_in, micro_cldliq,                               &
                             cldice_in, micro_cldice,                               &
                             numliq_in, micro_numliq,                               &
                             numice_in, micro_numice,                               &
                             rainliq_in, micro_rainliq,                             &
                             snowice_in, micro_snowice,                             &
                             numrain_in, micro_numrain,                             &
                             numsnow_in, micro_numsnow,                             &
                             graupice_in, micro_graupice,                           &
                             numgraup_in, micro_numgraup,                           &
                             relvar_in, micro_relvar,                               &
                             accre_enhan_in, micro_accre_enhan,                     &
                             pmid_in, micro_pmid,                                   &
                             pdel_in, micro_pdel,                                   &
                             pint_in, micro_pint,                                   &
                             strat_cldfrc_in, micro_strat_cldfrc,                   &
                             strat_liq_cldfrc_in, micro_strat_liq_cldfrc,           &
                             strat_ice_cldfrc_in, micro_strat_ice_cldfrc,           & 
                             qsatfac_in, micro_qsatfac,                             &
                             naai_in, micro_naai,                                   &
                             npccn_in, micro_npccn,                                 &
                             rndst_in, micro_rndst,                                 &
                             nacon_in, micro_nacon,                                 &
                             snowice_tend_external_in, micro_snowice_tend_external, &
                             numsnow_tend_external_in, micro_numsnow_tend_external, &
                             effi_external_in, micro_effi_external,                 &
                             frzimm_in, micro_frzimm,                               &
                             frzcnt_in, micro_frzcnt,                               &
                             frzdep_in, micro_frzdep,                               &
                             errmsg, errcode)

    !External dependencies:
    use ccpp_kinds,        only: kind_phys

    !Host model dimensions/parameters:
    integer,         intent(in) :: ncol
    integer,         intent(in) :: nlev
    integer,         intent(in) :: nlevp1
    integer,         intent(in) :: micro_ncol         !Number of horizontal microphysics columns (count)
    integer,         intent(in) :: micro_nlev         !Number of microphysics vertical layers (count)
    integer,         intent(in) :: micro_nlevp1       !Number of microphysics vertical interfaces (count)

    ! Air temperature (K)
    real(kind_phys), intent(in)  :: airT_in(:, :)
    real(kind_phys), intent(out) :: micro_airT(:, :)
    ! Water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: airq_in(:, :)
    real(kind_phys), intent(out) :: micro_airq(:, :)
    ! Cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: cldliq_in(:, :)
    real(kind_phys), intent(out) :: micro_cldliq(:, :)
    ! Cloud ice mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: cldice_in(:, :)
    real(kind_phys), intent(out) :: micro_cldice(:, :)
    ! Mass number concentration of cloud liquid water wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numliq_in(:, :)
    real(kind_phys), intent(out) :: micro_numliq(:, :)
    ! Mass number concentration of cloud ice wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numice_in(:, :)
    real(kind_phys), intent(out) :: micro_numice(:, :)
    ! Rain mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: rainliq_in(:, :)
    real(kind_phys), intent(out) :: micro_rainliq(:, :)
    ! Snow mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: snowice_in(:, :)
    real(kind_phys), intent(out) :: micro_snowice(:, :)
    ! Mass number concentration of rain wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numrain_in(:, :)
    real(kind_phys), intent(out) :: micro_numrain(:, :)
    ! Mass number concentration of snow wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numsnow_in(:, :)
    real(kind_phys), intent(out) :: micro_numsnow(:, :)
    ! Graupel mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: graupice_in(:, :)
    real(kind_phys), intent(out) :: micro_graupice(:, :)
    ! Mass number concentration of graupel wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numgraup_in(:, :)
    real(kind_phys), intent(out) :: micro_numgraup(:, :)
    ! Relative variance of cloud water (1)
    real(kind_phys), intent(in)  :: relvar_in(:, :)
    real(kind_phys), intent(out) :: micro_relvar(:, :)
    ! Accretion enhancement factor (1)
    real(kind_phys), intent(in)  :: accre_enhan_in(:, :)
    real(kind_phys), intent(out) :: micro_accre_enhan(:, :)
    ! Air pressure (Pa)
    real(kind_phys), intent(in)  :: pmid_in(:, :)
    real(kind_phys), intent(out) :: micro_pmid(:, :)
    ! Air pressure thickness (Pa)
    real(kind_phys), intent(in)  :: pdel_in(:, :)
    real(kind_phys), intent(out) :: micro_pdel(:, :)
    ! Air pressure at interfaces (Pa)
    real(kind_phys), intent(in)  :: pint_in(:, :)
    real(kind_phys), intent(out) :: micro_pint(:, :)
    ! Stratiform cloud area fraction (fraction)
    real(kind_phys), intent(in)  :: strat_cldfrc_in(:, :)
    real(kind_phys), intent(out) :: micro_strat_cldfrc(:, :)
    ! Stratiform cloud liquid area fraction (fraction)
    real(kind_phys), intent(in)  :: strat_liq_cldfrc_in(:, :)
    real(kind_phys), intent(out) :: micro_strat_liq_cldfrc(:, :)
    ! Stratiform cloud ice area fraction (fraction)
    real(kind_phys), intent(in)  :: strat_ice_cldfrc_in(:, :)
    real(kind_phys), intent(out) :: micro_strat_ice_cldfrc(:, :)
    ! Subgrid cloud water saturation scaling factor (1)
    real(kind_phys), intent(in)  :: qsatfac_in(:, :)
    real(kind_phys), intent(out) :: micro_qsatfac(:, :)
    ! Tendency of activated ice nuclei mass number concentration (kg-1 s-1)
    real(kind_phys), intent(in)  :: naai_in(:, :)
    real(kind_phys), intent(out) :: micro_naai(:, :)
    ! Tendency of activated cloud condensation nuclei mass number concentration (kg-1 s-1)
    real(kind_phys), intent(in)  :: npccn_in(:, :)
    real(kind_phys), intent(out) :: micro_npccn(:, :)
    ! Dust radii by size bin  (m)
    real(kind_phys), intent(in)  :: rndst_in(:, :, :)
    real(kind_phys), intent(out) :: micro_rndst(:, micro_nlev, :)
    ! Dust number concentration by size bin (m-3)
    real(kind_phys), intent(in)  :: nacon_in(:, :, :)
    real(kind_phys), intent(out) :: micro_nacon(:, micro_nlev, :)
    ! Tendency of snow mixing ratio wrt moist air and condensed water from external microphysics (kg kg-1 s-1)
    real(kind_phys), intent(in)  :: snowice_tend_external_in(:, :)
    real(kind_phys), intent(out) :: micro_snowice_tend_external(:, :)
    ! Tendency of mass number concentration of snow wrt moist air and condensed water from external microphysics (kg-1 s-1)
    real(kind_phys), intent(in)  :: numsnow_tend_external_in(:, :)
    real(kind_phys), intent(out) :: micro_numsnow_tend_external(:, :)
    ! Effective radius of stratiform cloud ice particle from external microphysics (m)
    real(kind_phys), intent(in)  :: effi_external_in(:, :)
    real(kind_phys), intent(out) :: micro_effi_external(:, :)
    ! Tendency of cloud liquid droplet number concentration due to immersion freezing (cm-3)
    real(kind_phys), intent(in)  :: frzimm_in(:, :)
    real(kind_phys), intent(out) :: micro_frzimm(:, :)
    ! Tendency of cloud liquid droplet number concentration due to contact freezing (cm-3)
    real(kind_phys), intent(in)  :: frzcnt_in(:, :)
    real(kind_phys), intent(out) :: micro_frzcnt(:, :)
    ! Tendency of cloud ice number concentration due to deposition nucleation (cm-3)
    real(kind_phys), intent(in)  :: frzdep_in(:, :)
    real(kind_phys), intent(out) :: micro_frzdep(:, :)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    !Initialize error message and error code:
    errmsg  = ''
    errcode = 0

!+ IH 
! For now we just use nocls = micro_ncol, but we need to constrain the vertical extent for the microphysical fields.
! Therefore micro_xxx(:ncol,:) = xxx(:,::)
!- IH
    micro_airT(:ncol,:) = airT_in(:,::)
    micro_airq(:ncol,:) = airq_in(:,::)
    micro_cldliq(:ncol,:) = cldliq_in(:,::)
    micro_cldice(:ncol,:) = cldice_in(:,::)
    micro_numliq(:ncol,:) = numliq_in(:,::)
    micro_numice(:ncol,:) = numice_in(:,::)
    micro_rainliq(:ncol,:) = rainliq_in(:,::)
    micro_snowice(:ncol,:) = snowice_in(:,::)
    micro_numrain(:ncol,:) = numrain_in(:,::)
    micro_numsnow(:ncol,:) = numsnow_in(:,::)
    micro_graupice(:ncol,:) = graupice_in(:,::)
    micro_numgraup(:ncol,:) = numgraup_in(:,::)
    micro_relvar(:ncol,:) = relvar_in(:,::)
    micro_accre_enhan(:ncol,:) = accre_enhan_in(:,::)
    micro_pmid(:ncol,:) = pmid_in(:,::)
    micro_pdel(:ncol,:) = pdel_in(:,::)
    micro_pint(:ncol,:) = pint_in(:,:micro_nlevp1)
    micro_strat_cldfrc(:ncol,:) = strat_cldfrc_in(:,::)
    micro_strat_liq_cldfrc(:ncol,:) = strat_liq_cldfrc_in(:,::)
    micro_strat_ice_cldfrc(:ncol,:) = strat_ice_cldfrc_in(:,::)
    micro_qsatfac(:ncol,:) = qsatfac_in(:,::)
    micro_naai(:ncol,:) = naai_in(:,::)
    micro_npccn(:ncol,:) = npccn_in(:,::)
    micro_rndst(:ncol,:,:) = rndst_in(:,:micro_nlev,:)
    micro_nacon(:ncol,:,:) = nacon_in(:,:micro_nlev,:)
    micro_snowice_tend_external(:ncol,:) = snowice_tend_external_in(:,::)
    micro_numsnow_tend_external(:ncol,:) = numsnow_tend_external_in(:,::)
    micro_effi_external(:ncol,:) = effi_external_in(:,::)
    micro_frzcnt(:ncol,:) = frzcnt_in(:,::)
    micro_frzdep(:ncol,:) = frzdep_in(:,::)


  end subroutine micro_pumas_ccpp_dimensions_pre_run

end module micro_pumas_ccpp_dimensions_pre
