      module aerclm_def
      use machine , only : kind_phys, kind_io4
      implicit none

      integer, parameter   :: levsaer=72, ntrcaerm=15, timeaer=2
      integer              :: latsaer, lonsaer, ntrcaer, levsw
      integer              :: n1sv, n2sv
      integer              :: iamin, iamax, jamin, jamax

      character*10         :: specname(ntrcaerm)
      real (kind=kind_phys):: aer_time(13)

      real (kind=kind_phys), allocatable, dimension(:) :: aer_lat
      real (kind=kind_phys), allocatable, dimension(:) :: aer_lon
      real (kind=kind_io4),  allocatable, dimension(:,:,:,:) :: aer_pres
      real (kind=kind_io4),  allocatable, dimension(:,:,:,:,:) :: aerin

      data aer_time/15.5, 45.,  74.5,  105., 135.5, 166., 196.5,
     &             227.5, 258., 288.5, 319., 349.5, 380.5/

      data specname /'DU001','DU002','DU003','DU004','DU005',
     &               'SS001','SS002','SS003','SS004','SS005','SO4',
     &               'BCPHOBIC','BCPHILIC','OCPHOBIC','OCPHILIC'/

      end module aerclm_def
