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
  subroutine mp_pumas_pre_run(tgrs, prsl, prsi, qgrs, &
       ntqv, ntcw, ntiw, ntlnc, ntinc, ntrnc, ntsnc, ntrw, ntsw, ntgl, nthl, ntwa, ntgnc, &
       micro_ncol, micro_nlev, micro_nlevp1, micro_dust_nbins, &
       micro_airT, micro_airq, micro_cldliq, micro_cldice, micro_numliq, &
       micro_numice, micro_rainliq, micro_snowice, micro_numrain, micro_numsnow, micro_graupice,&
       micro_numgraup, micro_relvar, micro_accre_enhan, micro_pmid, micro_pdel, micro_pint,     &
       micro_strat_cldfrc, micro_strat_liq_cldfrc, micro_strat_ice_cldfrc, micro_qsatfac,       &
       micro_rndst, micro_nacon, micro_snowice_tend_external,          &
       micro_numsnow_tend_external, micro_effi_external, micro_frzimm, micro_frzcnt,            &
       micro_frzdep, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: micro_ncol, micro_nlev, micro_nlevp1, micro_dust_nbins
    integer, intent(in) :: ntqv, ntcw, ntiw, ntlnc, ntinc, ntrnc, ntsnc, ntrw, ntsw, ntgl, nthl, ntwa, ntgnc
    real(kind_phys), dimension(:,:),   intent(in) :: tgrs, prsl, prsi
    real(kind_phys), dimension(:,:,:), intent(in) :: qgrs

    ! Outputs
    real(kind_phys), dimension(:,:), intent(out) :: micro_airT, micro_airq, micro_cldliq,&
         micro_cldice, micro_numliq,  micro_numice, micro_rainliq, micro_snowice,        &
         micro_numrain, micro_numsnow, micro_graupice, micro_numgraup, micro_relvar,     &
         micro_accre_enhan, micro_pmid, micro_pdel, micro_pint, micro_strat_cldfrc,      &
         micro_strat_liq_cldfrc, micro_strat_ice_cldfrc, micro_qsatfac, &
          micro_snowice_tend_external, micro_numsnow_tend_external,          &
         micro_effi_external, micro_frzimm, micro_frzcnt, micro_frzdep
    real(kind_phys), dimension(:,:,:), intent(out) :: micro_rndst, micro_nacon

    ! CCPP error reporting
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    ! Sub sample input state by n-subcolumns
    micro_airT(:, 1:micro_nlev)     = tgrs(:, 1:micro_nlev)
    micro_airq(:, 1:micro_nlev)     = qgrs(:, 1:micro_nlev, ntqv)
    micro_cldliq(:, 1:micro_nlev)   = qgrs(:, 1:micro_nlev, ntcw)
    micro_cldice(:, 1:micro_nlev)   = qgrs(:, 1:micro_nlev, ntiw)
    micro_rainliq(:, 1:micro_nlev)  = qgrs(:, 1:micro_nlev, ntrw)
    micro_snowice(:, 1:micro_nlev)  = qgrs(:, 1:micro_nlev, ntsw)
    micro_graupice(:, 1:micro_nlev) = qgrs(:, 1:micro_nlev, ntgl)
    micro_numliq(:, 1:micro_nlev)   = qgrs(:, 1:micro_nlev, ntlnc)
    micro_numice(:, 1:micro_nlev)   = qgrs(:, 1:micro_nlev, ntinc)
    micro_numrain(:, 1:micro_nlev)  = qgrs(:, 1:micro_nlev, ntrnc)
    micro_numsnow(:, 1:micro_nlev)  = qgrs(:, 1:micro_nlev, ntsnc)
    micro_numgraup(:, 1:micro_nlev) = qgrs(:, 1:micro_nlev, ntgnc)
    micro_pmid(:, 1:micro_nlev)     = prsl(:, 1:micro_nlev)
    micro_pint(:, 1:micro_nlev)     = prsi(:, 1:micro_nlevp1)
    
    ! microphysics relative variance of cloud water
    micro_relvar(:, 1:micro_nlev) = 0._kind_phys

    ! microphysics accretion enhancement factor
    micro_accre_enhan(:, 1:micro_nlev) = 0._kind_phys

    ! Pressure thickness (compute)
    !micro_pdel(:, 1:micro_nlev) = 

    ! microphysics stratiform cloud area fraction
    micro_strat_cldfrc(:, 1:micro_nlev) = 0._kind_phys

    ! microphysics stratiform cloud liquid area fraction
    micro_strat_liq_cldfrc(:, 1:micro_nlev) = 0._kind_phys

    ! microphysics stratiform cloud ice area fraction
    micro_strat_ice_cldfrc(:, 1:micro_nlev) = 0._kind_phys

    ! microphysics subgrid cloud water saturation scaling factor
    micro_qsatfac(:, 1:micro_nlev) = 1._kind_phys

    ! microphysics tendency of snow mixing ratio wrt moist air and condensed water from external microphysics
    micro_snowice_tend_external(:, 1:micro_nlev) = 0._kind_phys

    ! microphysics tendency of mass number concentration of snow wrt moist air and condensed water from external microphysics
    micro_numsnow_tend_external(:, 1:micro_nlev) = 0._kind_phys

    ! microphysics effective radius of stratiform cloud ice particle from external microphysics
    micro_effi_external(:, 1:micro_nlev) = 0._kind_phys

    ! microphysics tendency of cloud liquid droplet number concentration due to immersion freezing
    micro_frzimm(:, 1:micro_nlev) = 0._kind_phys

    ! microphysics tendency of cloud liquid droplet number concentration due to contact freezing
    micro_frzcnt(:, 1:micro_nlev) = 0._kind_phys

    ! microphysics tendency of cloud ice number concentration due to deposition nucleation
    micro_frzdep(:, 1:micro_nlev) = 0._kind_phys
    
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
