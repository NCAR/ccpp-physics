module rrtmgp_aux
  use machine, only: &
       kind_phys                   ! Working type
  implicit none

  real(kind_phys) :: &
       rrtmgp_minP, & ! Minimum pressure allowed in RRTMGP
       rrtmgp_minT    ! Minimum temperature allowed in RRTMGP
contains
  ! #########################################################################################
  ! SUBROUTINE check_error_msg
  ! #########################################################################################
  subroutine check_error_msg(routine_name, error_msg)
    character(len=*), intent(in) :: &
         error_msg, routine_name
    
    if(error_msg /= "") then
       print*,"ERROR("//trim(routine_name)//"): "
       print*,trim(error_msg)
       return
    end if
  end subroutine check_error_msg  
end module rrtmgp_aux
