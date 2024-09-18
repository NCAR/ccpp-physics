! This module is the NRL replacement for the parse_tracers module
! in GFDL_parse_tracers.F90. It translates the tracer names used
! in GFS_typedefs (stemming from the FV3 dycore at GFDL) into
! the correct names used at NRL so that all tracer indices defined
! in GFS_typedefs either retain the original value from NRL (if
! in use) or are correctly set to NO_TRACER (if not in use).
module parse_tracers

  integer, parameter :: NO_TRACER = -99

  public get_tracer_index, NO_TRACER

CONTAINS

  function get_tracer_index (tracer_names, gfdl_name, me, master, debug)

    character(len=32), intent(in) :: tracer_names(:)
    character(len=*),  intent(in) :: gfdl_name
    integer,           intent(in) :: me
    integer,           intent(in) :: master
    logical,           intent(in) :: debug
    !--- local variables
    character(len=32) :: name
    integer :: get_tracer_index
    integer :: i

    get_tracer_index = NO_TRACER

    do i=1, size(tracer_names)
       select case(trim(gfdl_name))
         ! Should not make it to here, since ntqv is hardcoded to 1
         !case '???'
         !  name = 'water_vapor_mixing_ratio'
         case ('o3mr')
           name = 'ozone_mixing_ratio'
         case ('liq_wat')
           name = 'cloud_water_mixing_ratio'
         case ('ice_wat')
           name = 'cloud_ice_mixing_ratio'
         case ('rainwat')
           name = 'rain_mixing_ratio'
         case ('snowwat')
           name = 'snow_mixing_ratio'
         case ('graupel')
           name = 'graupel_mixing_ratio'
         case ('ice_nc')
           name = 'number_concentration_of_cloud_ic'
         case ('rain_nc')
           name = 'number_concentration_of_rain'
         case ('sgs_tke')
           name = 'turbulent_kinetic_energy'
         case default
           name = trim(gfdl_name)
           if (debug .and. (me == master)) then
             print *,' PE ',me,' no translation found for tracer with GFDL name '//trim(gfdl_name)//' - use as is'
           endif
       end select
       if (trim(name) == trim(tracer_names(i))) then
           get_tracer_index = i
           exit
       endif
    enddo

    if (debug .and. (me == master)) then
      if (get_tracer_index == NO_TRACER) then
        print *,' PE ',me,' tracer with name '//trim(gfdl_name)//'='//trim(name)//' not found'
      else
        print *,' PE ',me,' tracer FOUND: '//trim(gfdl_name)//'='//trim(name)
      endif
    endif

    return

  end function get_tracer_index

end module parse_tracers
