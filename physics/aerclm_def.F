!>\file aerclm_def.F
!! This file contains aerosol climatology definition in MG microphysics

!>\ingroup mod_GFS_phys_time_vary
!! This module defines aerosol arrays in MG microphysics.
      module aerclm_def
      use machine , only : kind_phys
      implicit none
     
! only read monthly merra2 data for m-1, m, m+1 
      integer, parameter :: levsaer=45, latsaer=91, lonsaer=144
      integer, parameter :: lmerra=72,  ntrcaerm=15, timeaer=12

      integer              :: ntrcaer
      character*10         :: specname(ntrcaerm)
      real (kind=kind_phys):: aer_lat(latsaer), aer_lon(lonsaer)
     &                       ,aer_time(13)
      real (kind=4),         allocatable, dimension(:,:,:,:,:) :: aerin
      real (kind=kind_phys), allocatable, dimension(:,:,:,:) :: aer_pres

      data aer_time/15.5, 45.,  74.5,  105., 135.5, 166., 196.5,     
     &             227.5, 258., 288.5, 319., 349.5, 380.5/

      data specname /'DU001','DU002','DU003','DU004','DU005', 
     &               'SS001','SS002','SS003','SS004','SS005','SO4',
     &               'BCPHOBIC','BCPHILIC','OCPHILIC','OCPHOBIC'/

      end module aerclm_def
