!> \file module_MYNNrad_post.F90
!!  Contains the post (interstitial) work after the call to the radiation schemes:
!!    1) Restores the original qc & qi

      MODULE mynnrad_post

      contains

      subroutine mynnrad_post_init ()
      end subroutine mynnrad_post_init

      subroutine mynnrad_post_finalize ()
      end subroutine mynnrad_post_finalize

!>\defgroup gsd_mynnrad_post GSD mynnrad_post_run Module
!>\ingroup gsd_mynn_edmf
!!  This interstitial code restores the original resolved-scale clouds (qc and qi).
#if 0
!! \section arg_table_mynnrad_post_run Argument Table
!! \htmlinclude mynnrad_post_run.html
!!
#endif
SUBROUTINE mynnrad_post_run(               &
     &     ix,im,levs,                     &
     &     flag_init,flag_restart,         &
     &     qc,qi,                          &
     &     qc_save, qi_save,               &
     &     errmsg, errflg                  )

! should be moved to inside the mynn:
      use machine , only : kind_phys

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 

      integer, intent(in)  :: ix, im, levs
      logical,          intent(in)  :: flag_init, flag_restart
      real(kind=kind_phys), dimension(im,levs), intent(out) :: qc, qi
      real(kind=kind_phys), dimension(im,levs), intent(in)  :: qc_save, qi_save
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variable
      integer              :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !write(0,*)"=============================================="
      !write(0,*)"in mynn rad post"

      if (flag_init .and. (.not. flag_restart)) then
        !write (0,*) 'Skip MYNNrad_post flag_init = ', flag_init
        return
      endif

     ! Add subgrid cloud information:
        do k = 1, levs
           do i = 1, im

               qc(i,k) = qc_save(i,k)
               qi(i,k) = qi_save(i,k)

           enddo
        enddo

       ! print*,"===Finished restoring the resolved-scale clouds"
       ! print*,"qc_save:",qc_save(1,1)," qc:",qc(1,1)

  END SUBROUTINE mynnrad_post_run

!###=================================================================

END MODULE mynnrad_post
