! ###########################################################################################
!> \file GFS_physics_diagnostics.F90
!!
! ###########################################################################################
module GFS_physics_diagnostics
  use machine, only : kind_phys, kind_dbl_prec, kind_sngl_prec
  implicit none
  public GFS_physics_diagnostics_init, GFS_physics_diagnostics_run
contains

! ########################################################################################### 
! SUBROUTINE GFS_physics_diagnostics_init
! ###########################################################################################
!! \section arg_table_GFS_physics_diagnostics_init Argument Table
!! \htmlinclude GFS_physics_diagnostics_init.html
!!
  subroutine GFS_physics_diagnostics_init(errmsg, errflg)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg         ! CCPP error message
    integer, intent(out) :: &
         errflg         ! CCPP error flag

  end subroutine GFS_physics_diagnostics_init

! ###########################################################################################
! SUBROUTINE GFS_physics_diagnostics_run
! ###########################################################################################
!! \section arg_table_GFS_physics_diagnostics_run Argument Table
!! \htmlinclude GFS_physics_diagnostics_run.html
!!
  subroutine GFS_physics_diagnostics_run(nCol, nLev, ntoz, dtidx, ip_prod_loss, ip_ozmix,   &
       ip_temp, ip_overhead_ozone, do3_dt_prd, do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz, dtend,&
       errmsg, errflg)
    ! Inputs
    integer, intent(in) :: &
         nCol,           & ! Horizontal dimension
         nLev,           & ! Number of vertical layers
         ntoz,           & ! Index for ozone mixing ratio
         ip_prod_loss,   & ! Index for process in diagnostic tendency output
         ip_ozmix,       & ! Index for process in diagnostic tendency output
         ip_temp,        & ! Index for process in diagnostic tendency output
         ip_overhead_ozone ! Index for process in diagnostic tendency output    
    integer, intent(in), dimension(:,:) :: &
         dtidx             ! Bookkeeping indices for GFS diagnostic tendencies

    ! Inputs (optional)
    real(kind=kind_phys), intent(in), dimension(:,:), pointer, optional :: &
         do3_dt_prd,  & ! Physics tendency: production and loss effect
         do3_dt_ozmx, & ! Physics tendency: ozone mixing ratio effect
         do3_dt_temp, & ! Physics tendency: temperature effect
         do3_dt_ohoz    ! Physics tendency: overhead ozone effect

    ! Outputs
    real(kind=kind_phys), intent(inout), dimension(:,:,:) :: &
         dtend          ! Diagnostic tendencies for state variables
    character(len=*), intent(out) :: &
         errmsg         ! CCPP error message
    integer, intent(out) :: &
         errflg         ! CCPP error flag

    ! Locals
    integer :: idtend
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! #######################################################################################
    !
    ! Ozone physics diagnostics
    !
    ! #######################################################################################
    idtend = dtidx(100+ntoz,ip_prod_loss)
    if (idtend >= 1 .and. associated(do3_dt_prd)) then  
       dtend(:,:,idtend) = dtend(:,:,idtend) + do3_dt_prd
    endif
    !
    idtend = dtidx(100+ntoz,ip_ozmix)
    if (idtend >= 1 .and. associated(do3_dt_ozmx)) then
       dtend(:,:,idtend) = dtend(:,:,idtend) + do3_dt_ozmx
    endif
    !
    idtend = dtidx(100+ntoz,ip_temp)
    if (idtend >= 1 .and. associated(do3_dt_temp)) then
       dtend(:,:,idtend) = dtend(:,:,idtend) + do3_dt_temp
    endif
    !
    idtend = dtidx(100+ntoz,ip_overhead_ozone)
    if (idtend >= 1 .and. associated(do3_dt_ohoz)) then
       dtend(:,:,idtend) = dtend(:,:,idtend) + do3_dt_ohoz
    endif

  end subroutine GFS_physics_diagnostics_run

end module GFS_physics_diagnostics
