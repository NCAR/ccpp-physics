!>\file rrfs_smoke_postpbl.F90
!! This file is CCPP RRFS smoke postpbl driver
!! Haiqin.Li@noaa.gov 03/2022

 module rrfs_smoke_postpbl

   use machine ,        only : kind_phys

   implicit none

   private

   public :: rrfs_smoke_postpbl_run

contains

!> \section arg_table_rrfs_smoke_postpbl_run Argument Table
!! \htmlinclude rrfs_smoke_postpbl_run.html
!!
    subroutine rrfs_smoke_postpbl_run(ite, kte, ntsmoke, ntdust, ntcoarsepm, ntrac,      &
                   qgrs, chem3d, rrfs_sd, errmsg, errflg)

    implicit none


    integer,        intent(in) :: ite,kte,ntsmoke,ntdust,ntcoarsepm,ntrac

    integer, parameter :: its=1,kts=1

    real(kind_phys), dimension(:,:,:), intent(inout) :: qgrs
    real(kind_phys), dimension(:,:,:), intent(inout), optional :: chem3d
    logical, intent(in) :: rrfs_sd
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

!>-- local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    if (.not. rrfs_sd) return

    !--- put smoke stuff back into tracer array

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntsmoke)= chem3d(i,k,1)
       qgrs(i,k,ntdust )= chem3d(i,k,2)
       qgrs(i,k,ntcoarsepm)= chem3d(i,k,3)
     enddo
    enddo

    return

   end subroutine rrfs_smoke_postpbl_run

  end module rrfs_smoke_postpbl
