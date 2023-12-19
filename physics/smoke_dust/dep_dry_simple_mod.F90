!>\file dep_dry_simple_mod.F90
!! This file is for the dry depostion driver.

module dep_dry_simple_mod

  use machine ,        only : kind_phys

  implicit none

  private

  public :: dry_dep_driver_simple

contains

    subroutine dry_dep_driver_simple(rmol,ust,ndvel,ddvel,                &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
!----------------------------------------------------------------------
  IMPLICIT NONE

      INTEGER,      INTENT(IN   ) :: ndvel,                               &
                                  ids,ide, jds,jde, kds,kde,              &
                                  ims,ime, jms,jme, kms,kme,              &
                                  its,ite, jts,jte, kts,kte
      REAL(kind_phys),  DIMENSION( ims:ime , jms:jme )        ,           &
          INTENT(INOUT) :: ust, rmol

      REAL(kind_phys), PARAMETER :: kpart=500.
      REAL(kind_phys) :: dvpart

!
! Output array
      REAL(kind_phys), DIMENSION( its:ite, jts:jte, ndvel ), INTENT(INOUT) ::   ddvel


      integer :: i,j,k,nv
!
! necessary for aerosols (module dependent)
!
! .. Intrinsic Functions ..
      INTRINSIC max, min

! compute dry deposition velocities = ddvel

        ddvel(:,:,:) = 0.0
        do nv = 1, ndvel
         do j = jts, jte
          do i = its, ite
           dvpart = ust(i,j)/kpart

            IF (rmol(i,j)<0.) THEN       ! UNSTABLE LAYERING CORRECTION
               dvpart = dvpart*(1.+(-300.*rmol(i,j))**0.66667)
            ENDIF

            ddvel(i,j,nv) = MIN(0.50,dvpart)        ! m/s
          enddo
         enddo
        enddo

end subroutine dry_dep_driver_simple

end module dep_dry_simple_mod
