!>\file h2o_def.f 
!! This file contains array definition in H2O scheme.

!>\ingroup mod_GFS_phys_time_vary
!! This module defines arrays in H2O scheme.
      module h2o_def

!> \section arg_table_h2o_def
!! \htmlinclude h2o_def.html
!!

      use machine , only : kind_phys
      implicit none

      integer, parameter :: kh2opltc=29

      integer latsh2o, levh2o, timeh2o,  h2o_coeff
      real (kind=kind_phys), allocatable :: h2o_lat(:), h2o_pres(:)
     &,                                     h2o_time(:)
      real (kind=kind_phys), allocatable :: h2oplin(:,:,:,:)

      end module h2o_def
