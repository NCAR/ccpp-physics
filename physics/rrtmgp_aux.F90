module rrtmgp_aux
  implicit none
contains
  !
  subroutine rrtmgp_aux_init()
  end subroutine rrtmgp_aux_init
  !
  subroutine rrtmgp_aux_run()
  end subroutine rrtmgp_aux_run
  !
  subroutine rrtmgp_aux_finalize()
  end subroutine rrtmgp_aux_finalize

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
