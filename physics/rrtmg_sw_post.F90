!>\file rrtmg_sw_post
!! This file contains
      module rrtmg_sw_post
      contains

!>\defgroup rrtmg_sw_post GFS RRTMG scheme post
!! @{
!> \section arg_table_rrtmg_sw_post_init Argument Table
!!
      subroutine rrtmg_sw_post_init ()
      end subroutine rrtmg_sw_post_init
! PGI compiler does not accept lines longer than 264 characters, remove during pre-processing
#ifndef __PGI
!> \section arg_table_rrtmg_sw_post_run Argument Table
!! \htmlinclude rrtmg_sw_post_run.html
!!
#endif
      subroutine rrtmg_sw_post_run (Model, Grid, Diag, Radtend, Coupling,  &
                 im, ltp, nday, lm, kd, htswc, htsw0,                      &
                 sfcalb1, sfcalb2, sfcalb3, sfcalb4, scmpsw, errmsg, errflg)

      use machine,                   only: kind_phys
      use module_radsw_parameters,   only: topfsw_type, sfcfsw_type,   &
                                           cmpfsw_type
      use GFS_typedefs,              only: GFS_coupling_type,          &
                                           GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type,           &
                                           GFS_diag_type

      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_grid_type),            intent(in)    :: Grid
      type(GFS_diag_type),            intent(inout) :: Diag
      integer,                        intent(in)    :: im, lm, kd, nday, ltp
      type(cmpfsw_type),    dimension(size(Grid%xlon,1)), intent(inout) :: scmpsw
      real(kind=kind_phys), dimension(Size(Grid%xlon,1), lm+LTP), intent(in) ::  htswc, htsw0
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: sfcalb1, sfcalb2, sfcalb3, sfcalb4
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer :: i, k1, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lsswr) then
        if (nday > 0) then
          do k = 1, LM
            k1 = k + kd
            Radtend%htrsw(1:im,k) = htswc(1:im,k1)
          enddo
          ! We are assuming that radiative tendencies are from bottom to top 
          ! --- repopulate the points above levr i.e. LM
          if (lm < Model%levs) then
            do k = lm+1,Model%levs
              Radtend%htrsw (1:im,k) = Radtend%htrsw (1:im,LM)
            enddo
          endif

          if (Model%swhtr) then
            do k = 1, lm
               k1 = k + kd
               Radtend%swhc(1:im,k) = htsw0(1:im,k1)
             enddo
             ! --- repopulate the points above levr i.e. LM
             if (lm < Model%levs) then
               do k = lm+1,Model%levs
                 Radtend%swhc(1:im,k) = Radtend%swhc(1:im,LM)
               enddo
             endif
          endif

!  --- surface down and up spectral component fluxes
!>  - Save two spectral bands' surface downward and upward fluxes for
!!    output.

          do i=1,im
            Coupling%nirbmdi(i) = scmpsw(i)%nirbm
            Coupling%nirdfdi(i) = scmpsw(i)%nirdf
            Coupling%visbmdi(i) = scmpsw(i)%visbm
            Coupling%visdfdi(i) = scmpsw(i)%visdf

            Coupling%nirbmui(i) = scmpsw(i)%nirbm * sfcalb1(i)
            Coupling%nirdfui(i) = scmpsw(i)%nirdf * sfcalb2(i)
            Coupling%visbmui(i) = scmpsw(i)%visbm * sfcalb3(i)
            Coupling%visdfui(i) = scmpsw(i)%visdf * sfcalb4(i)
          enddo

        else                   ! if_nday_block

          Radtend%htrsw(:,:) = 0.0

          Radtend%sfcfsw = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          Diag%topfsw    = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw         = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )

          do i=1,im
            Coupling%nirbmdi(i) = 0.0
            Coupling%nirdfdi(i) = 0.0
            Coupling%visbmdi(i) = 0.0
            Coupling%visdfdi(i) = 0.0

            Coupling%nirbmui(i) = 0.0
            Coupling%nirdfui(i) = 0.0
            Coupling%visbmui(i) = 0.0
            Coupling%visdfui(i) = 0.0
          enddo

          if (Model%swhtr) then
            Radtend%swhc(:,:) = 0
          endif

        endif                  ! end_if_nday

! --- radiation fluxes for other physics processes
        do i=1,im
          Coupling%sfcnsw(i) = Radtend%sfcfsw(i)%dnfxc - Radtend%sfcfsw(i)%upfxc
          Coupling%sfcdsw(i) = Radtend%sfcfsw(i)%dnfxc
        enddo

      endif                                ! end_if_lsswr

      end subroutine rrtmg_sw_post_run
 
!> \section arg_table_rrtmg_sw_post_finalize Argument Table
!!
      subroutine rrtmg_sw_post_finalize ()
      end subroutine rrtmg_sw_post_finalize
!! @}
      end module rrtmg_sw_post
