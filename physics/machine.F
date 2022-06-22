      module machine

!! \section arg_table_machine
!! \htmlinclude machine.html
!!

      implicit none
      
      integer, parameter :: kind_io4  = 4, kind_io8  = 8 , kind_ior = 8 &
     &,                     kind_evod = 8, kind_dbl_prec = 8            &
     &,                     kind_sngl_prec = 4, kind_INTEGER = 4        &
     &,                     kind_LOGICAL = 4

#ifdef SINGLE_PREC
      integer, parameter :: kind_rad  = kind_sngl_prec                  &
     &,                     kind_phys = kind_sngl_prec                  &
     &,                     kind_grid = kind_dbl_prec                   &! atmos_cubed_sphere requres kind_grid=8
     &,                     kind_REAL = kind_sngl_prec                   ! used in cmp_comm
#else
      integer, parameter :: kind_rad  = kind_dbl_prec                   &
     &,                     kind_phys = kind_dbl_prec                   &
     &,                     kind_grid = kind_dbl_prec                   &! atmos_cubed_sphere requres kind_grid=8
     &,                     kind_REAL = kind_dbl_prec                    ! used in cmp_comm
#endif

#ifdef OVERLOAD_R4
      integer, parameter :: kind_dyn  = 4
#else
      integer, parameter :: kind_dyn  = 8
#endif

!
      real(kind=kind_evod), parameter :: mprec = 1.e-12           ! machine precision to restrict dep
      real(kind=kind_evod), parameter :: grib_undef = 9.99e20     ! grib undefine value
!
      end module machine
