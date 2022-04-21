!>\file rrfs_smoke_postpbl.F90
!! This file is CCPP RRFS smoke postpbl driver
!! Haiqin.Li@noaa.gov 03/2022

 module rrfs_smoke_postpbl

   use machine ,        only : kind_phys
   use rrfs_smoke_config

   implicit none

   private

   public :: rrfs_smoke_postpbl_init, rrfs_smoke_postpbl_run, rrfs_smoke_postpbl_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine rrfs_smoke_postpbl_init()
      end subroutine rrfs_smoke_postpbl_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_rrfs_smoke_postpbl_finalize Argument Table
!!
      subroutine rrfs_smoke_postpbl_finalize()
      end subroutine rrfs_smoke_postpbl_finalize

!> \defgroup gsd_chem_group GSD Chem emission driver Module
!! This is the gsd chemistry
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

    real(kind_phys), dimension(ite,kte,ntrac), intent(inout) :: qgrs
    real(kind_phys), dimension(ite,kte,    2), intent(inout) :: chem3d
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
