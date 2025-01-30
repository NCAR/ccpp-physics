! ########################################################################################
!>\file mp_pumas_pre.F90
!!

!> This module contains the pre-processing step prior to calling the PUMAS microphysics.
!
! ########################################################################################
module mp_pumas_pre
  use machine, only: kind_phys, kind_dbl_prec

  implicit none
  public mp_pumas_pre_init, mp_pumas_pre_run, mp_pumas_pre_finalize
contains
  ! ######################################################################################
  !> \section arg_table_mp_pumas_pre_init Argument Table
  !! \htmlinclude mp_pumas_pre_init.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_pre_init(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_pre_init

  ! ######################################################################################
  !> \section arg_table_mp_pumas_pre_run Argument Table
  !! \htmlinclude mp_pumas_pre_run.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_pre_run(micro_ncol, micro_nlev, micro_nlevp1, micro_dust_nbins, &
       micro_airT, micro_airq, micro_cldliq, micro_cldice, micro_numliq, &
       micro_numice, micro_rainliq, micro_snowice, micro_numrain, micro_numsnow, micro_graupice,&
       micro_numgraup, micro_relvar, micro_accre_enhan, micro_pmid, micro_pdel, micro_pint,     &
       micro_strat_cldfrc, micro_strat_liq_cldfrc, micro_strat_ice_cldfrc, micro_qsatfac,       &
       micro_naai, micro_npccn, micro_rndst, micro_nacon, micro_snowice_tend_external,          &
       micro_numsnow_tend_external, micro_effi_external, micro_frzimm, micro_frzcnt,            &
       micro_frzdep, errmsg, errflg)
    ! Inputs
    integer, intent(in) :: micro_ncol, micro_nlev, micro_nlevp1, micro_dust_nbins

    ! Outputs
    real(kind_phys), dimension(:,:), intent(out) :: micro_airT, micro_airq, micro_cldliq,&
         micro_cldice, micro_numliq,  micro_numice, micro_rainliq, micro_snowice,        &
         micro_numrain, micro_numsnow, micro_graupice, micro_numgraup, micro_relvar,     &
         micro_accre_enhan, micro_pmid, micro_pdel, micro_pint, micro_strat_cldfrc,      &
         micro_strat_liq_cldfrc, micro_strat_ice_cldfrc, micro_qsatfac, micro_naai,      &
         micro_npccn, micro_snowice_tend_external, micro_numsnow_tend_external,          &
         micro_effi_external, micro_frzimm, micro_frzcnt, micro_frzdep
    real(kind_phys), dimension(:,:,:), intent(out) :: micro_rndst, micro_nacon

    ! CCPP error reporting
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_pre_run

  ! ######################################################################################
  !> \section arg_table_mp_pumas_pre_finalize Argument Table
  !! \htmlinclude mp_pumas_pre_finalize.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_pre_finalize(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_pre_finalize

end module mp_pumas_pre
