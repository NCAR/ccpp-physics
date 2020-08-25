!>\file rrtmg_lw_post
!!This file contains
      module rrtmg_lw_post 
      contains

!>\defgroup rrtmg_lw_post GFS RRTMG scheme post
!! @{
!> \section arg_table_rrtmg_lw_post_init Argument Table
!!
      subroutine rrtmg_lw_post_init()
      end subroutine rrtmg_lw_post_init

! PGI compiler does not accept lines longer than 264 characters, remove during pre-processing
#ifndef __PGI
!> \section arg_table_rrtmg_lw_post_run Argument Table
!! \htmlinclude rrtmg_lw_post_run.html
!!
#endif
      subroutine rrtmg_lw_post_run (Model, Grid, Radtend, Coupling,   &
                 im, ltp, lm, kd, tsfa, htlwc, htlw0, errmsg, errflg)
    
      use machine,                   only: kind_phys
      use GFS_typedefs,              only: GFS_coupling_type,          &
                                           GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type
      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_grid_type),            intent(in)    :: Grid
      type(GFS_radtend_type),         intent(inout) :: Radtend
      integer,                        intent(in)    :: im, ltp, LM, kd
      real(kind=kind_phys), dimension(size(Grid%xlon,1), lm+LTP), intent(in) :: htlwc, htlw0
      real(kind=kind_phys), dimension(size(Grid%xlon,1)),         intent(in) :: tsfa
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! local variables
      integer :: k1, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lslwr) then
!> -# Save calculation results
!>  - Save surface air temp for diurnal adjustment at model t-steps

        Radtend%tsflw (:) = tsfa(:)

        do k = 1, LM
          k1 = k + kd
            Radtend%htrlw(1:im,k) = htlwc(1:im,k1)
        enddo
        ! --- repopulate the points above levr
        if (lm < Model%levs) then
          do k = lm+1,Model%levs
            Radtend%htrlw (1:im,k) = Radtend%htrlw (1:im,LM)
          enddo
        endif

        if (Model%lwhtr) then
          do k = 1, lm
            k1 = k + kd
            Radtend%lwhc(1:im,k) = htlw0(1:im,k1)
          enddo
          ! --- repopulate the points above levr
          if (lm < Model%levs) then
            do k = lm+1,Model%levs
              Radtend%lwhc(1:im,k) = Radtend%lwhc(1:im,LM)
            enddo
          endif
        endif

! --- radiation fluxes for other physics processes
        Coupling%sfcdlw(:) = Radtend%sfcflw(:)%dnfxc

      endif                                ! end_if_lslwr

      end subroutine rrtmg_lw_post_run

!> \section arg_table_rrtmg_lw_post_finalize Argument Table
!!
      subroutine rrtmg_lw_post_finalize ()
      end subroutine rrtmg_lw_post_finalize

!! @}
      end module rrtmg_lw_post
