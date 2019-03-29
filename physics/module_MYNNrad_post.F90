!> \file module_MYNNrad_post.F90
!!  Contains the post (interstitial) work after the call to the radiation schemes:
!!    1) Restores the original qc & qi

      MODULE mynnrad_post

      contains

      subroutine mynnrad_post_init ()
      end subroutine mynnrad_post_init

      subroutine mynnrad_post_finalize ()
      end subroutine mynnrad_post_finalize

!!
!> \brief This interstitial code restores the original resolved-scale clouds (qc and qi).
#if 0
!! \section arg_table_mynnrad_post_run Argument Table
!! | local_name          | standard_name                                                               | long_name                                                                  | units   | rank | type      |    kind   | intent | optional |
!! |---------------------|-----------------------------------------------------------------------------|----------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | ix                  | horizontal_dimension                                                        | horizontal dimension                                                       | count   |    0 | integer   |           | in     | F        |
!! | im                  | horizontal_loop_extent                                                      | horizontal loop extent                                                     | count   |    0 | integer   |           | in     | F        |
!! | levs                | vertical_dimension                                                          | vertical layer dimension                                                   | count   |    0 | integer   |           | in     | F        |
!! | qc                  | cloud_condensed_water_mixing_ratio                                          | moist (dry+vapor, no condensates) mixing ratio of cloud water (condensate) | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | qi                  | ice_water_mixing_ratio                                                      | moist (dry+vapor, no condensates) mixing ratio of ice water                | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | qc_save             | cloud_condensed_water_mixing_ratio_save    | moist (dry+vapor, no condensates) mixing ratio of cloud water (condensate) before entering a physics scheme | kg kg-1 |    2 | real      | kind_phys | in     | F        |
!! | qi_save             | ice_water_mixing_ratio_save                | moist (dry+vapor, no condensates) mixing ratio of ice water before entering a physics scheme                | kg kg-1 |    2 | real      | kind_phys | in     | F        |
!! | errmsg              | ccpp_error_message                                                          | error message for error handling in CCPP                                   | none    |    0 | character | len=*     | out    | F        |
!! | errflg              | ccpp_error_flag                                                             | error flag for error handling in CCPP                                      | flag    |    0 | integer   |           | out    | F        |
!!
#endif
!###===================================================================
SUBROUTINE mynnrad_post_run(               &
     &     ix,im,levs,                     &
     &     qc,qi,                          &
     &     qc_save, qi_save,               &
     &     errmsg, errflg                  )

! should be moved to inside the mynn:
      use machine , only : kind_phys

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 

      integer, intent(in)  :: ix, im, levs
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
