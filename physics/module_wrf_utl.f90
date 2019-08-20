module module_wrf_utl
 implicit none
contains

subroutine wrf_error_fatal(string)
  implicit none
  character(len=*), intent(in) :: string
  print*, string
  stop
end subroutine wrf_error_fatal

subroutine wrf_message(msg)
  implicit none
  character(len=*), intent(in) :: msg
  write(*,'(A)') msg
end subroutine wrf_message

logical function wrf_dm_on_monitor() result (return_value)
  implicit none
  return_value = .TRUE.
end function wrf_dm_on_monitor

subroutine wrf_dm_bcast_real(rval, ival)
  implicit none
  real, intent(in) :: rval
  integer, intent(in) :: ival
end subroutine wrf_dm_bcast_real

subroutine wrf_dm_bcast_integer(ival1, ival2)
  implicit none
  real, intent(in) :: ival1
  integer, intent(in) :: ival2
end subroutine wrf_dm_bcast_integer

subroutine wrf_dm_bcast_string(sval, ival)
  implicit none
  character(len=*), intent(in) :: sval
  integer, intent(in) :: ival
end subroutine wrf_dm_bcast_string

subroutine wrf_debug( level , str ) 
  implicit none 
  character*(*) str 
  integer , intent (in) :: level 
  call wrf_message( str ) 
  return 
end subroutine wrf_debug 

end module module_wrf_utl

