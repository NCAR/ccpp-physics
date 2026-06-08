!> \file m_micro_pre.F90
!! This file contains subroutines that prepare data for the Morrison-Gettelman microphysics scheme
!! as part of the GFS physics suite.
      module m_micro_pre

      implicit none

      contains

!! \section arg_table_m_micro_pre_run Argument Table
!! \htmlinclude m_micro_pre_run.html
!!
      subroutine m_micro_pre_run (im, levs, do_shoc, skip_macro, cld_shoc,     &
        cld_frc_MG, clcn, errmsg, errflg )

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs
      logical, intent(in) :: do_shoc
      logical, intent(inout) :: skip_macro

      real(kind=kind_phys), intent(in), optional :: cld_shoc(:,:)
      real(kind=kind_phys), intent(inout) :: cld_frc_MG(:,:)

      real(kind=kind_phys), intent(in) :: clcn(:,:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      skip_macro = do_shoc
      if (do_shoc) then
        do k=1,levs
          do i=1,im
            cld_frc_MG(i,k) = cld_shoc(i,k)
          enddo
        enddo 
      end if

      ! add convective cloud fraction
      do k = 1,levs
        do i = 1,im
          cld_frc_MG(i,k) = min(1.0, cld_frc_MG(i,k) + clcn(i,k))
        enddo
      enddo

      end subroutine m_micro_pre_run

      end module m_micro_pre
