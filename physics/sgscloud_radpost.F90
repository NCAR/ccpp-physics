!> \file sgscloud_radpost.F90
!!  Contains the post (interstitial) work after the call to the radiation schemes:
!!    1) Restores the original qc & qi
      module sgscloud_radpost

      contains

!>\defgroup sgscloud_radpost_mod sgscloud_radpost_run Module
!>  This interstitial code restores the original resolved-scale clouds (qc and qi).
!> \section arg_table_sgscloud_radpost_run Argument Table
!! \htmlinclude sgscloud_radpost_run.html
!!
      subroutine sgscloud_radpost_run( &
           im,levs,                    &
           flag_init,flag_restart,     &
           qc,qi,qs,                   &
           qc_save,qi_save,qs_save,    &
           errmsg,errflg               )

! should be moved to inside the mynn:
      use machine , only : kind_phys

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 

      integer, intent(in)  :: im, levs
      logical,          intent(in)  :: flag_init, flag_restart
      real(kind=kind_phys), dimension(:,:), intent(inout) :: qc, qi, qs
      real(kind=kind_phys), dimension(:,:), intent(in)    :: qc_save, qi_save, qs_save
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variable
      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !write(0,*)"=============================================="
      !write(0,*)"in sgscloud rad post"

      if (flag_init .and. (.not. flag_restart)) then
        !write (0,*) 'Skip MYNNrad_post flag_init = ', flag_init
        return
      endif

      ! Add subgrid cloud information:
      do k = 1, levs
        do i = 1, im
          qc(i,k) = qc_save(i,k)
          qi(i,k) = qi_save(i,k)
          qs(i,k) = qs_save(i,k)
        enddo
      enddo

      ! print*,"===Finished restoring the resolved-scale clouds"
      ! print*,"qc_save:",qc_save(1,1)," qc:",qc(1,1)

      end subroutine sgscloud_radpost_run
      end module sgscloud_radpost
