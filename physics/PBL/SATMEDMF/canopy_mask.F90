   module canopy_mask_mod

   use machine , only : kind_phys

   implicit none

! Vertical arrays
   integer :: nkt
   integer, parameter :: nkc = 3 ! # of canopy layers for shading effects

   public :: nkt  ! # of resolved model layers plus canopy layers

   public :: canopy_mask_init, canopy_mask_run

   contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   subroutine canopy_mask_init(im, ix, km, &! nkt,    &
           claie, cfch, cfrt, cclu, cpopu, &  !in:
           FRT_mask,                       &  ! out
           errmsg,errflg)

   implicit none

! Horizontal arrays
   integer  :: im, ix, km  ! horizontal & vertical domain specifications

   real(kind=kind_phys) :: claie(im), cfch(im), cfrt(im), &
                            cclu(im),cpopu(im)
   real(kind=kind_phys) :: FRT_mask(im)

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!...local variables

! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

!...Allocate and initialize new canopy arrays

! Initializations

   FRT_mask(:)=0.0

   nkt= km + nkc   ! # of resolved model layers plus canopy layers

   return
   end subroutine canopy_mask_init

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   subroutine canopy_mask_run (im, ix, km, &  !in:
           claie, cfch, cfrt, cclu, cpopu, &  !in:
           FRT_mask,                       &  !out:
           errmsg,errflg)

   implicit none

!...Arguments:

! Horizontal arrays
   integer  :: im, ix, km  ! horizontal & vertical domain specifications

   real(kind=kind_phys) :: claie(im), cfch(im),  cfrt(im), &
                                      cclu(im), cpopu(im)
   real(kind=kind_phys) :: FRT_mask(im)

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!...local variables

   integer i,is,k,n

! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   do i=1,im

      !NOT a Continuos forest canopy
      if (    claie(i) .LT. 0.1                      &
         .OR. cfch (i) .LT. 0.5                      &
!IVAI: modified contiguous canopy condition
!        .OR. MAX(0.0, 1.0 - cfrt(i)) .GT. 0.5
         .OR. MAX(0.0, 1.0 - cfrt(i)) .GT. 0.75      &
         .OR. cpopu(i) .GT. 10000.0                  &
         .OR. (EXP(-0.5*claie(i)*cclu(i)) .GT. 0.45  &
         .AND. cfch(i)  .LT. 18.) ) THEN

         FRT_mask(i) = -1.0

      ! Continuous forest canopy
      ELSE

         FRT_mask(i) = 1.0

      END IF ! Forest Canopy Mask

   end do ! i=1,im

   return
   end subroutine canopy_mask_run

   end module canopy_mask_mod
