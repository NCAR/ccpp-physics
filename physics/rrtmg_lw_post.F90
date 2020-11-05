!>\file rrtmg_lw_post
!!This file contains
      module rrtmg_lw_post 
      contains

!>\defgroup rrtmg_lw_post GFS RRTMG scheme post
!! @{
      subroutine rrtmg_lw_post_init()
      end subroutine rrtmg_lw_post_init

!> \section arg_table_rrtmg_lw_post_run Argument Table
!! \htmlinclude rrtmg_lw_post_run.html
!!
      subroutine rrtmg_lw_post_run (im, levs, ltp, lm, kd, lslwr, lwhtr,       &
                 tsfa, htlwc, htlw0, sfcflw, tsflw, sfcdlw, htrlw, lwhc,       &
                 errmsg, errflg)
    
      use machine,                   only: kind_phys
      use module_radlw_parameters,   only: sfcflw_type
      
      implicit none
      
      integer,                                     intent(in) :: im, levs, ltp, lm, kd
      logical,                                     intent(in) :: lslwr, lwhtr
      real(kind=kind_phys), dimension(im),         intent(in) ::  tsfa
      real(kind=kind_phys), dimension(im, LM+LTP), intent(in) ::  htlwc
      real(kind=kind_phys), dimension(im, LM+LTP), intent(in) ::  htlw0
      
      type(sfcflw_type), dimension(im),            intent(in) :: sfcflw
      
      real(kind=kind_phys), dimension(im),         intent(inout) ::  tsflw, sfcdlw
      real(kind=kind_phys), dimension(im, levs),   intent(inout) ::  htrlw, lwhc
      character(len=*),                            intent(out) :: errmsg
      integer,                                     intent(out) :: errflg
      
      ! local variables
      integer :: k1, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lslwr) then
!> -# Save calculation results
!>  - Save surface air temp for diurnal adjustment at model t-steps

        tsflw (:) = tsfa(:)

        do k = 1, LM
          k1 = k + kd
            htrlw(1:im,k) = htlwc(1:im,k1)
        enddo
        ! --- repopulate the points above levr
        if (lm < levs) then
          do k = lm+1, levs
            htrlw (1:im,k) = htrlw (1:im,LM)
          enddo
        endif

        if (lwhtr) then
          do k = 1, lm
            k1 = k + kd
            lwhc(1:im,k) = htlw0(1:im,k1)
          enddo
          ! --- repopulate the points above levr
          if (lm < levs) then
            do k = lm+1, levs
              lwhc(1:im,k) = lwhc(1:im,LM)
            enddo
          endif
        endif

! --- radiation fluxes for other physics processes
        sfcdlw(:) = sfcflw(:)%dnfxc

      endif                                ! end_if_lslwr

      end subroutine rrtmg_lw_post_run

      subroutine rrtmg_lw_post_finalize ()
      end subroutine rrtmg_lw_post_finalize

!! @}
      end module rrtmg_lw_post
