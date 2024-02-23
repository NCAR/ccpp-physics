!>\file sfc_data.F
!! This file contains an data surface scheme.

!> This module contains the CCPP-compliant CDEPS data scheme 
!! scheme when the model is using data provided by CDEPS.
      module sfc_data
      implicit none
      private
      public :: sfc_data_run

      contains

!>\defgroup gfs_data_main Simple Wrapper for CDEPS inline 
!! This subroutine pass CDEPS inline provided data to other schemes
!! by overwriting the specific part
!>@{
!! \section arg_table_sfc_data_run Argument Table
!! \htmlinclude sfc_data_run.html
!!
!!>\section gen_sfc_data CDEPS Inline data scheme
      subroutine sfc_data_run                                           &
!...................................
!  ---  inputs:
     &     ( im, use_data, tsfco_dat, mask_dat, tsfc_wat,               &
     &       errmsg, errflg,                                            &
     &       naux2d, aux2d)

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_data                                                      !
!       inputs:                                                         !
!        ( im, use_data, tsfco_dat                                      !
!       outputs:                                                        !
!          tsfc_wat )                                                   !
!                                                                       !
!                                                                       !
!  subprograms/functions called: fpvs                                   !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!    feb  2024  -- u. turuncoglu added initially to interact with       !
!                  fv3 cdeps inline                                     !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horizontal dimension                     1    !
!     use_data - logical, =T if data provided by cdeps inline      1    !
!     tsfco_dat- real, sea surface temperature ( k )               im   !
!                                                                       !
!  outputs:                                                             !
!     tsfc_wat - real, sea surface temperature ( k )               im   !
!                                                                       !
! ===================================================================== !
!
      use machine , only : kind_phys
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: im
      logical, intent(in) :: use_data 
      real (kind=kind_phys), dimension(:), intent(in) :: tsfco_dat
      real (kind=kind_phys), dimension(:), intent(in) :: mask_dat

!  ---  outputs:
      real (kind=kind_phys), dimension(:), intent(inout) :: tsfc_wat
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      integer, intent(in) :: naux2d
      real(kind=kind_phys), intent(inout) :: aux2d(:,:)

!  ---  locals:
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check coupling from component land to atmosphere
      if (.not. use_data) return

      aux2d(:,1) = tsfco_dat(:)
      aux2d(:,2) = mask_dat(:)
      aux2d(:,3) = tsfc_wat(:)

      ! Fill variables
      !do i = 1, im
      !   if (tsfco_dat(i)
      !end do
!
      return
!...................................
      end subroutine sfc_data_run
!-----------------------------------
!>@}
      end module sfc_data
