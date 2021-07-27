!>\file ozne_def.f
!! This file contains the ozone array definition used in ozone physics.

!>\ingroup mod_GFS_phys_time_vary
!! This module defines arrays in Ozone scheme.
      module ozne_def

!> \section arg_table_ozne_def
!! \htmlinclude ozne_def.html
!!

      use machine , only : kind_phys
      implicit none
      
      integer, parameter :: kozpl=28, kozc=48

      integer latsozp, levozp, timeoz, latsozc, levozc, timeozc
     &,       oz_coeff
      real (kind=kind_phys) blatc, dphiozc
      real (kind=kind_phys), allocatable :: oz_lat(:), oz_pres(:)
     &,                                     oz_time(:)
      real (kind=kind_phys), allocatable :: ozplin(:,:,:,:)

      end module ozne_def
