      module date_def
      use machine,   ONLY: kind_phys
      implicit none
      
!jw      integer idate(4)
!jw      real(kind=kind_phys) fhour,shour,thour,z00
       real(kind=kind_phys) shour,thour,z00
       real(kind=kind_phys),target :: fhour, zhour
       integer,target :: idate(4),idate7(7)
!
      REAL(KIND=KIND_PHYS) ,ALLOCATABLE :: spdmax(:)

      end module date_def
