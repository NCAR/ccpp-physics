module rrtmgp_lw_pre
  use machine, only: &
       kind_phys                   ! Working type
  use module_radiation_surface,  only: &
       setemis                     ! Routine to compute surface-emissivity
  use mo_gas_optics_rrtmgp,  only: &
       ty_gas_optics_rrtmgp
  use rrtmgp_lw_gas_optics, only: lw_gas_props

  implicit none

  public rrtmgp_lw_pre_run,rrtmgp_lw_pre_init,rrtmgp_lw_pre_finalize

contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_pre_init
  ! #########################################################################################
  subroutine rrtmgp_lw_pre_init ()
  end subroutine rrtmgp_lw_pre_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_pre_run
  ! #########################################################################################
!> \section arg_table_rrtmgp_lw_pre_run
!! \htmlinclude rrtmgp_lw_pre_run.html
!!
  subroutine rrtmgp_lw_pre_run (doLWrad, nCol, xlon, xlat, slmsk, zorl, snowd, sncovr, &
       tsfg, tsfa, hprime, sfc_emiss_byband, emiss, semis, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doLWrad          ! Logical flag for longwave radiation call
    integer, intent(in) :: &
         nCol             ! Number of horizontal grid points
    real(kind_phys), dimension(nCol), intent(in) :: &
         xlon,          & ! Longitude
         xlat,          & ! Latitude
         slmsk,         & ! Land/sea/sea-ice mask
         zorl,          & ! Surface roughness length (cm)
         snowd,         & ! water equivalent snow depth (mm)
         sncovr,        & ! Surface snow are fraction (1)
         tsfg,          & ! Surface ground temperature for radiation (K)
         tsfa,          & ! Lowest model layer air temperature for radiation (K)
         hprime           ! Standard deviation of subgrid orography
    real(kind_phys), dimension(:), intent(in) :: &
         emiss            ! Surface emissivity from Noah MP

    ! Outputs 
    real(kind_phys), dimension(lw_gas_props%get_nband(),ncol), intent(out) :: &
         sfc_emiss_byband ! Surface emissivity in each band
    character(len=*), intent(out) :: &
         errmsg           ! Error message
    integer, intent(out) :: &  
         errflg           ! Error flag
    real(kind_phys), dimension(nCol), intent(out) :: &
         semis

    ! Local variables
    integer :: iBand

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (.not. doLWrad) return

    ! #######################################################################################
    ! Call module_radiation_surface::setemis(),to setup surface emissivity for LW radiation.
    ! #######################################################################################
    call setemis (xlon, xlat, slmsk, snowd, sncovr, zorl, tsfg, tsfa, hprime, emiss, nCol, semis)

    ! Assign same emissivity to all bands
    do iBand=1,lw_gas_props%get_nband()
       sfc_emiss_byband(iBand,:) = semis
    enddo

  end subroutine rrtmgp_lw_pre_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_pre_finalize
  ! #########################################################################################
  subroutine rrtmgp_lw_pre_finalize ()
  end subroutine rrtmgp_lw_pre_finalize

end module rrtmgp_lw_pre
