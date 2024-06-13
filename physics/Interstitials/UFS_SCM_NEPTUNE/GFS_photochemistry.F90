! #########################################################################################
!> \file GFS_photochemistry.f90
!!
! #########################################################################################
module GFS_photochemistry
  use machine,        only: kind_phys
  use module_ozphys,  only: ty_ozphys
  use module_h2ophys, only: ty_h2ophys
  implicit none
contains
! #########################################################################################
!> \section arg_table_GFS_photochemistry_run Argument Table
!! \htmlinclude GFS_photochemistry_run.html
!!
! #########################################################################################
  subroutine GFS_photochemistry_run (dtp, ozphys, oz_phys_2015, oz_phys_2006, con_1ovg,   &
       prsl, dp, ozpl, h2o_phys, h2ophys, h2opl, h2o0, oz0, gt0, do3_dt_prd, do3_dt_ozmx, &
       do3_dt_temp, do3_dt_ohoz, errmsg, errflg)
    
    ! Inputs
    real(kind=kind_phys), intent(in ) :: &
         dtp,          & ! Model timestep
         con_1ovg        ! Physical constant (1./gravity)
    real(kind=kind_phys), intent(in ), dimension(:,:) :: &
         prsl,         & ! Air pressure (Pa)
         dp              ! Pressure thickness (Pa)
    real(kind=kind_phys), intent(in ), dimension(:,:,:) :: &
         ozpl,         & ! Ozone data for current model timestep
         h2opl           ! h2o data for curent model timestep
    logical,              intent(in) :: &
         oz_phys_2015, & ! Do ozone photochemistry? (2015)
         oz_phys_2006, & ! Do ozone photochemistry? (2006)
         h2o_phys        ! Do h2o photochemistry?
    type(ty_ozphys), intent(in) :: &
         ozphys          ! DDT with ozone photochemistry scheme.
    type(ty_h2ophys), intent(in) :: &
         h2ophys         ! DDT with h2o photochemistry scheme.

    ! Outputs (optional)
    real(kind=kind_phys), intent(inout), dimension(:,:), optional :: &
         do3_dt_prd,   & ! Physics tendency: production and loss effect
         do3_dt_ozmx,  & ! Physics tendency: ozone mixing ratio effect
         do3_dt_temp,  & ! Physics tendency: temperature effect
         do3_dt_ohoz     ! Physics tendency: overhead ozone effect

    ! Outputs
    real(kind=kind_phys), intent(out), dimension(:,:) :: &
         oz0,          & ! Update ozone concentration.
         h2o0            ! Updated h2o concentration.
    real(kind=kind_phys), intent(inout), dimension(:,:) :: &
         gt0             ! Updated temperature
    
    character(len=*), intent(out) :: &
         errmsg          ! CCPP Error message.
    integer,  intent(out) :: &
         errflg          ! CCPP Error flag.

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (oz_phys_2015) then
       call ozphys%run_o3prog_2015(con_1ovg, dtp, prsl, gt0, dp, ozpl, oz0, do3_dt_prd,    &
            do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz)
    endif
    if (oz_phys_2006) then
       call ozphys%run_o3prog_2006(con_1ovg, dtp, prsl, gt0, dp, ozpl, oz0, do3_dt_prd,    &
            do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz)
    endif
    if (h2o_phys) then
       call h2ophys%run(dtp, prsl, h2opl, h2o0)
    endif

  end subroutine GFS_photochemistry_run

end module GFS_photochemistry
