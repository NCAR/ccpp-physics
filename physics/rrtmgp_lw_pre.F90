!> \file rrtmgp_lw_pre.F90
!!
!> \defgroup rrtmgp_lw_pre rrtmgp_lw_pre.F90
!!
!! \brief RRTMGP Longwave pre-processing routine.
!!
module rrtmgp_lw_pre
  use machine, only: &
       kind_phys                   ! Working type
  use mo_gas_optics_rrtmgp,  only: &
       ty_gas_optics_rrtmgp
  use rrtmgp_lw_gas_optics, only: lw_gas_props

  implicit none

  public rrtmgp_lw_pre_run

contains

!>\defgroup rrtmgp_lw_pre_mode GFS RRTMGP-LW Pre Module
!> \section arg_table_rrtmgp_lw_pre_run
!! \htmlinclude rrtmgp_lw_pre_run.html
!!
!> \ingroup rrtmgp_lw_pre
!!
!! \brief
!!
!! \section rrtmgp_lw_pre_run
  subroutine rrtmgp_lw_pre_run (doLWrad, semis, sfc_emiss_byband, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doLWrad
    real(kind_phys), dimension(:), intent(in) :: &
         semis

    ! Outputs
    real(kind_phys), dimension(:,:), intent(inout) :: &
         sfc_emiss_byband ! Surface emissivity in each band
    character(len=*), intent(out) :: &
         errmsg           ! Error message
    integer, intent(out) :: &  
         errflg           ! Error flag

    ! Local variables
    integer :: iBand

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return

    ! Assign same emissivity to all bands
    do iBand=1,lw_gas_props%get_nband()
       sfc_emiss_byband(iBand,:) = semis
    enddo

  end subroutine rrtmgp_lw_pre_run
  
end module rrtmgp_lw_pre
