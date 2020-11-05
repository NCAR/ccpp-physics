!>\file rrtmg_sw_post
!! This file contains
      module rrtmg_sw_post
      contains

!>\defgroup rrtmg_sw_post GFS RRTMG scheme post
!! @{
      subroutine rrtmg_sw_post_init ()
      end subroutine rrtmg_sw_post_init

!> \section arg_table_rrtmg_sw_post_run Argument Table
!! \htmlinclude rrtmg_sw_post_run.html
!!
      subroutine rrtmg_sw_post_run (im, levr, levs, ltp, nday, lm, kd, lsswr,  &
                 swhtr, sfcalb1, sfcalb2, sfcalb3, sfcalb4, htswc, htsw0,      &
                 nirbmdi, nirdfdi, visbmdi, visdfdi, nirbmui, nirdfui, visbmui,&
                 visdfui, sfcdsw, sfcnsw, htrsw, swhc, scmpsw, sfcfsw, topfsw, &
                 errmsg, errflg)

      use machine,                   only: kind_phys
      use module_radsw_parameters,   only: topfsw_type, sfcfsw_type,   &
                                           cmpfsw_type

      implicit none

      integer,                              intent(in)    :: im, levr, levs,   &
                                                             ltp, nday, lm, kd   
      logical,                              intent(in)    :: lsswr, swhtr
      real(kind=kind_phys), dimension(im),  intent(in)    :: sfcalb1, sfcalb2, &
                                                             sfcalb3, sfcalb4
      real(kind=kind_phys), dimension(im, levr+LTP), intent(in) ::  htswc, htsw0
      
      real(kind=kind_phys), dimension(im),  intent(inout) :: nirbmdi, nirdfdi, &
                                                             visbmdi, visdfdi, &
                                                             nirbmui, nirdfui, &
                                                             visbmui, visdfui, &
                                                             sfcdsw,  sfcnsw
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: htrsw, swhc

      type(cmpfsw_type), dimension(im),     intent(inout) :: scmpsw
      type(sfcfsw_type), dimension(im),     intent(inout) :: sfcfsw
      type(topfsw_type), dimension(im),     intent(inout) :: topfsw

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer :: i, k1, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lsswr) then
        if (nday > 0) then
          do k = 1, LM
            k1 = k + kd
            htrsw(1:im,k) = htswc(1:im,k1)
          enddo
          ! We are assuming that radiative tendencies are from bottom to top 
          ! --- repopulate the points above levr i.e. LM
          if (lm < levs) then
            do k = lm+1, levs
              htrsw (1:im,k) = htrsw (1:im,LM)
            enddo
          endif

          if (swhtr) then
            do k = 1, lm
               k1 = k + kd
               swhc(1:im,k) = htsw0(1:im,k1)
             enddo
             ! --- repopulate the points above levr i.e. LM
             if (lm < levs) then
               do k = lm+1, levs
                 swhc(1:im,k) = swhc(1:im,LM)
               enddo
             endif
          endif

!  --- surface down and up spectral component fluxes
!>  - Save two spectral bands' surface downward and upward fluxes for
!!    output.

          do i=1,im
            nirbmdi(i) = scmpsw(i)%nirbm
            nirdfdi(i) = scmpsw(i)%nirdf
            visbmdi(i) = scmpsw(i)%visbm
            visdfdi(i) = scmpsw(i)%visdf

            nirbmui(i) = scmpsw(i)%nirbm * sfcalb1(i)
            nirdfui(i) = scmpsw(i)%nirdf * sfcalb2(i)
            visbmui(i) = scmpsw(i)%visbm * sfcalb3(i)
            visdfui(i) = scmpsw(i)%visdf * sfcalb4(i)
          enddo

        else                   ! if_nday_block

          htrsw(:,:) = 0.0

          sfcfsw = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          topfsw = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )

          do i=1,im
            nirbmdi(i) = 0.0
            nirdfdi(i) = 0.0
            visbmdi(i) = 0.0
            visdfdi(i) = 0.0

            nirbmui(i) = 0.0
            nirdfui(i) = 0.0
            visbmui(i) = 0.0
            visdfui(i) = 0.0
          enddo

          if (swhtr) then
            swhc(:,:) = 0
          endif

        endif                  ! end_if_nday

! --- radiation fluxes for other physics processes
        do i=1,im
          sfcnsw(i) = sfcfsw(i)%dnfxc - sfcfsw(i)%upfxc
          sfcdsw(i) = sfcfsw(i)%dnfxc
        enddo

      endif                                ! end_if_lsswr

      end subroutine rrtmg_sw_post_run
 
      subroutine rrtmg_sw_post_finalize ()
      end subroutine rrtmg_sw_post_finalize
!! @}
      end module rrtmg_sw_post
