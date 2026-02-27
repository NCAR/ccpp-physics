   module canopy_mask_mod

   use machine , only : kind_phys

   implicit none

   public :: canopy_mask_init, canopy_mask_run

   contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   subroutine canopy_mask_init(im, km, nkc, nkt, &
           claie, cfch, cfrt, cclu, cpopu, &  !in:
           FRT_mask,                       &  ! out
           errmsg,errflg)

   implicit none

! Horizontal arrays
   integer  :: im, km  ! horizontal & vertical domain specifications
   integer, intent(in)  :: nkc, nkt

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

   return
   end subroutine canopy_mask_init

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   subroutine canopy_mask_run (im, km, nkc, nkt, &  !in:
           claie, cfch, cfrt, cclu, cpopu, &  !in:
           FRT_mask,                       &  !out:
           errmsg,errflg)

   implicit none

!...Arguments:

! Horizontal arrays
   integer  :: im, km  ! horizontal & vertical domain specifications
   integer, intent(in) :: nkc, nkt

   real(kind=kind_phys) :: claie(im), cfch(im),  cfrt(im), &
                                      cclu(im), cpopu(im)
   real(kind=kind_phys) :: FRT_mask(im)

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!...local variables

   integer i

! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   do i=1,im

      !NOT a Continuous forest canopy
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

   end do

   return
   end subroutine canopy_mask_run

   end module canopy_mask_mod
