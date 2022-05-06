!>\file rrfs_smoke_postpbl.F90
!! This file is CCPP RRFS smoke postpbl driver
!! Haiqin.Li@noaa.gov 03/2022

 module rrfs_smoke_postpbl

   use machine ,        only : kind_phys
   use rrfs_smoke_config

   implicit none

   private

   public :: rrfs_smoke_postpbl_run

contains

!>\defgroup rrfs_smoke_postpbl GSD Chem emission driver Module  
!> \ingroup gsd_chem_group
!! This is the GSD Chem emission driver Module
!! \section arg_table_rrfs_smoke_postpbl_run Argument Table
!! \htmlinclude rrfs_smoke_postpbl_run.html
!!
!>\section rrfs_smoke_postpbl GSD Chemistry Scheme General Algorithm
!> @{
    subroutine rrfs_smoke_postpbl_run(ite, kte, ntsmoke, ntdust, ntrac,      &
                   qgrs, chem3d, errmsg, errflg)

    implicit none


    integer,        intent(in) :: ite,kte,ntsmoke,ntdust,ntrac

    integer, parameter :: its=1,kts=1

    real(kind_phys), dimension(:,:,:), intent(inout) :: qgrs
    real(kind_phys), dimension(:,:,:), intent(inout) :: chem3d
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

!>-- local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    !--- put smoke stuff back into tracer array

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntsmoke)= chem3d(i,k,1)
       qgrs(i,k,ntdust )= chem3d(i,k,2)
     enddo
    enddo

 end subroutine rrfs_smoke_postpbl_run

!> @}
  end module rrfs_smoke_postpbl
