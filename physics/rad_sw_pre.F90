! ######################################################################################
!>\file rad_sw_pre.f90
!!
!! This file gathers the sunlit points for the shortwave radiation schemes.
!!
!> \defgroup rad_sw_pre GFS radiation pre routine.
!! @{ 
!!
! ###################################################################################### 
module rad_sw_pre
contains

  ! ###################################################################################
!> \section arg_table_rad_sw_pre_run Argument Table
!! \htmlinclude rad_sw_pre_run.html
!!
!! \section rad_sw_pre_run
!! @{
  ! ###################################################################################
  subroutine rad_sw_pre_run (im, lsswr, coszen, nday, idxday, errmsg, errflg)
    use machine,  only: kind_phys
    implicit none
    
    ! Inputs
    integer,                      intent(in)    :: im
    logical,                      intent(in)    :: lsswr
    real(kind_phys), dimension(:), intent(in)   :: coszen
    
    ! Outputs
    integer,                      intent(out)   :: nday
    integer, dimension(:),        intent(out)   :: idxday
    character(len=*),             intent(out)   :: errmsg
    integer,                      intent(out)   :: errflg
    
    ! Local variables
    integer :: i
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (lsswr) then
       ! Check for daytime points for SW radiation.
       nday = 0
       idxday = 0
       do i = 1, IM
          if (coszen(i) >= 0.0001) then
             nday = nday + 1
             idxday(nday) = i
          endif
       enddo
    else
       nday   = 0
       idxday = 0
    endif
    
  end subroutine rad_sw_pre_run
!! @}
end module rad_sw_pre
