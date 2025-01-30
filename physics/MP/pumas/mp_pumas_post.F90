! ########################################################################################
!>\file mp_pumas_post.F90
!!

!> This module contains the post-processing step after calling the PUMAS microphysics.
!
! ########################################################################################
module mp_pumas_post
  use machine, only: kind_phys, kind_dbl_postc
  use proc_rates_type
  implicit none
  public mp_pumas_post_init, mp_pumas_post_run, mp_pumas_post_finalize
contains
  ! ######################################################################################
  !> \section arg_table_mp_pumas_post_init Argument Table
  !! \htmlinclude mp_pumas_post_init.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_post_init(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_post_init

  ! ######################################################################################
  !> \section arg_table_mp_pumas_post_run Argument Table
  !! \htmlinclude mp_pumas_post_run.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_post_run(micro_prect, micro_preci, micro_proc_rates, micro_qcsinksum_rate1ord, &
             micro_airT_tend, micro_airq_tend, micro_cldliq_tend, micro_cldice_tend, &
             micro_numliq_tend, micro_numice_tend, micro_rainliq_tend, micro_snowice_tend,&
             micro_numrain_tend, micro_numsnow_tend, micro_graupice_tend, micro_numgraup_tend,&
             micro_effc, micro_effc_fn, micro_effi, micro_sadice, micro_sadsnow,&
             micro_prec_evap, micro_am_evap_st, micro_prec_prod, micro_cmeice,&
             micro_deffi, micro_pgamrad, micro_lamcrad, micro_snowice_in_prec,&
             micro_scaled_diam_snow, micro_graupice_in_prec, micro_numgraup_vol_in_prec,&
             micro_scaled_diam_graup, micro_lflx, micro_iflx, micro_gflx,&
             micro_rflx, micro_sflx, micro_rainliq_in_prec, micro_reff_rain,&
             micro_reff_snow, micro_reff_grau, micro_numrain_vol_in_prec,&
             micro_numsnow_vol_in_prec, micro_refl, micro_arefl, micro_areflz, &
             micro_frefl, micro_csrfl, micro_acsrfl, micro_fcsrfl, micro_refl10cm,&
             micro_reflz10cm, micro_rercld, micro_ncai, micro_ncal, micro_rainliq,&
             micro_snowice, micro_numrain_vol, micro_numsnow_vol, micro_diam_rain,&
             micro_diam_snow, micro_graupice, micro_numgraup_vol, micro_diam_graup,&
             micro_lflx, micro_iflx, micro_gflx, micro_rflx, micro_sflx,&
             micro_rainliq_in_prec, micro_reff_rain, micro_reff_snow, micro_reff_grau, &
             micro_numrain_vol_in_prec, micro_numsnow_vol_in_prec, micro_refl, micro_arefl, &
             micro_areflz, micro_frefl, micro_csrfl, micro_acsrfl, micro_fcsrfl, &
             micro_refl10cm, micro_reflz10cm, micro_rercld, micro_ncai, micro_ncal, &
             micro_rainliq, micro_snowice, micro_numrain_vol, micro_numsnow_vol, &
             micro_diam_rain, micro_diam_snow, micro_graupice, micro_numgraup_vol, &
             micro_diam_graup, micro_freq_graup, micro_freq_snow, micro_freq_rain, &
             micro_frac_ice, micro_frac_cldliq_tend, micro_rain_evap, micro_proc_rates, &
             errmsg, errflg)

    ! Inputs
    real(kind_phys), dimension(:),   intent(in) :: micro_prect, micro_preci
    real(kind_phys), dimension(:,:), intent(in) :: &
         micro_qcsinksum_rate1ord, &
         micro_airT_tend, &
         micro_airq_tend, &
         micro_cldliq_tend, &
         micro_cldice_tend, &
         micro_numliq_tend, &
         micro_numice_tend, &
         micro_rainliq_tend, &
         micro_snowice_tend, &
         micro_numrain_tend, &
         micro_numsnow_tend, &
         micro_graupice_tend, &
         micro_numgraup_tend, &
         micro_effc, &
         micro_effc_fn, &
         micro_effi, &
         micro_sadice, &
         micro_sadsnow, &
         micro_prec_evap, &
         micro_am_evap_st, &
         micro_prec_prod, &
         micro_cmeice, &
         micro_deffi, &
         micro_pgamrad, &
         micro_lamcrad, &
         micro_snowice_in_prec, &
         micro_scaled_diam_snow, &
         micro_graupice_in_prec, &
         micro_numgraup_vol_in_prec, &
         micro_scaled_diam_graup, &
         micro_lflx, &
         micro_iflx, &
         micro_gflx, &
         micro_rflx, &
         micro_sflx, &
         micro_rainliq_in_prec, &
         micro_reff_rain, &
         micro_reff_snow, &
         micro_reff_grau, &
         micro_numrain_vol_in_prec, &
         micro_numsnow_vol_in_prec, &
         micro_refl, &
         micro_arefl, &
         micro_areflz, &
         micro_frefl, &
         micro_csrfl, &
         micro_acsrfl, &
         micro_fcsrfl, &
         micro_refl10cm, &
         micro_reflz10cm, &
         micro_rercld, &
         micro_ncai, &
         micro_ncal, &
         micro_rainliq, &
         micro_snowice, &
         micro_numrain_vol, &
         micro_numsnow_vol, &
         micro_diam_rain, &
         micro_diam_snow, &
         micro_graupice, &
         micro_numgraup_vol, &
         micro_diam_graup, &
         micro_freq_graup, &
         micro_freq_snow, &
         micro_freq_rain, &
         micro_frac_ice, &
         micro_frac_cldliq_tend, &
         micro_rain_evap, &
         micro_proc_rates

    type(proc_rates_type), intent(inout) ::  micro_proc_rates
   
    ! CCPP error handling
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_post_run

  ! ######################################################################################
  !> \section arg_table_mp_pumas_post_finalize Argument Table
  !! \htmlinclude mp_pumas_post_finalize.html
  !!
  ! ######################################################################################
  subroutine mp_pumas_post_finalize(errmsg, errflg)
    character(len=*), intent(  out) :: errmsg
    integer,          intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine mp_pumas_post_finalize

end module mp_pumas_post
