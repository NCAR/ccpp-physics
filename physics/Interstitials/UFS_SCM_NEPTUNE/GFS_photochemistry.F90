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
!> \section arg_table_GFS_photochemistry_init Argument Table
!! \htmlinclude GFS_photochemistry_init.html
!!
! #########################################################################################
  subroutine GFS_photochemistry_init(oz_phys_2006, oz_phys_2015, h2o_phys, errmsg, errflg)
    logical, intent(in) :: &
         oz_phys_2015, & ! Do ozone photochemistry? (2015)
         oz_phys_2006, & ! Do ozone photochemistry? (2006)
         h2o_phys        ! Do stratospheric h2o photochemistry?
    character(len=*), intent(out) :: &
         errmsg          ! CCPP Error message.
    integer,  intent(out) :: &
         errflg          ! CCPP Error flag.

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! If no photchemical scheme is on, but SDF has this module, report an error?
    if ((.not. oz_phys_2006) .and. (.not. oz_phys_2015) .and. (.not. h2o_phys)) then
       write (errmsg,'(*(a))') 'Logic error: One of [oz_phys_2006, oz_phys_2015, or h2o_phys] must == .true. '
       errflg = 1
       return
    endif
    
    ! Only one ozone scheme can be on. Otherwise, return and report error.
    if (oz_phys_2006 .and. oz_phys_2015) then
       write (errmsg,'(*(a))') 'Logic error: Only one ozone scheme can be enabled at a time'
       errflg = 1
       return
    endif
    
  end subroutine GFS_photochemistry_init

! #########################################################################################
!> \section arg_table_GFS_photochemistry_run Argument Table
!! \htmlinclude GFS_photochemistry_run.html
!!
! #########################################################################################
  subroutine GFS_photochemistry_run (dtp, ozphys, oz_phys_2015, oz_phys_2006, con_1ovg,   &
       prsl, dp, ozpl, h2o_phys, h2ophys, h2opl, h2o0, oz0, gt0, do3_dt_prd, do3_dt_ozmx, &
       do3_dt_temp, do3_dt_ohoz, dqv_dt_prd, dqv_dt_qvmx, errmsg, errflg)
    
    ! Inputs
    real(kind=kind_phys), intent(in) :: &
         dtp,          & ! Model timestep
         con_1ovg        ! Physical constant (1./gravity)
    real(kind=kind_phys), intent(in), dimension(:,:) :: &
         prsl,         & ! Air pressure (Pa)
         dp,           & ! Pressure thickness (Pa)
         gt0             ! Air temperature (K)
    real(kind=kind_phys), intent(in), dimension(:,:,:) :: &
         ozpl,         & ! Ozone data for current model timestep.
         h2opl           ! h2o data for curent model timestep.
    logical,              intent(in) :: &
         oz_phys_2015, & ! Do ozone photochemistry? (2015)
         oz_phys_2006, & ! Do ozone photochemistry? (2006)
         h2o_phys        ! Do stratospheric h2o photochemistry?
    type(ty_ozphys), intent(in) :: &
         ozphys          ! DDT with ozone photochemistry scheme/data.
    type(ty_h2ophys), intent(in) :: &
         h2ophys         ! DDT with h2o photochemistry scheme/data.

    ! Outputs (optional)
    real(kind=kind_phys), intent(inout), dimension(:,:), optional :: &
         do3_dt_prd,   & ! Physics tendency: production and loss effect
         do3_dt_ozmx,  & ! Physics tendency: ozone mixing ratio effect
         do3_dt_temp,  & ! Physics tendency: temperature effect
         do3_dt_ohoz,  & ! Physics tendency: overhead ozone effect
         dqv_dt_prd,   & ! Physics tendency: Climatological net production effect
         dqv_dt_qvmx     ! Physics tendency: specific humidity effect

    ! Outputs
    real(kind=kind_phys), intent(inout), dimension(:,:) :: &
         oz0,          & ! Update ozone concentration.
         h2o0            ! Updated h2o concentration.
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
       call h2ophys%run(dtp, prsl, h2opl, h2o0, dqv_dt_prd, dqv_dt_qvmx)
    endif

  end subroutine GFS_photochemistry_run

end module GFS_photochemistry
