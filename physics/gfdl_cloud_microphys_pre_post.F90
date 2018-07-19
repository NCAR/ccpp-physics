module gfdl_cloud_microphys_pre

   use machine, only : kind_phys

   implicit none

   private

   public gfdl_cloud_microphys_pre_run, gfdl_cloud_microphys_pre_init, gfdl_cloud_microphys_pre_finalize

contains

   subroutine gfdl_cloud_microphys_pre_init()
   end subroutine gfdl_cloud_microphys_pre_init

   subroutine gfdl_cloud_microphys_pre_finalize()
   end subroutine gfdl_cloud_microphys_pre_finalize

!! \section arg_table_gfdl_cloud_microphys_pre_run Argument Table
!! | local_name       | standard_name                                              | long_name                                              | units      | rank | type      | kind      | intent| optional |
!! |------------------|------------------------------------------------------------|--------------------------------------------------------|------------|------|-----------|-----------|-------|----------|
!! | im               | horizontal_loop_extent                                     | horizontal loop extent                                 | count      |    0 | integer   |           | in    | F        |
!! | rain0            | lwe_thickness_of_stratiform_precipitation_amount_per_day   | stratiform rain over 24h period                        | mm         |    1 | real      | kind_phys | out   | F        |
!! | ice0             | lwe_thickness_of_ice_amount_per_day                        | ice fall over 24h period                               | mm         |    1 | real      | kind_phys | out   | F        |
!! | snow0            | lwe_thickness_of_snow_amount_per_day                       | snow fall over 24h period                              | mm         |    1 | real      | kind_phys | out   | F        |
!! | graupel0         | lwe_thickness_of_graupel_amount_per_day                    | graupel fall over 24h period                           | mm         |    1 | real      | kind_phys | out   | F        |
!! | errmsg           | error_message                                              | error message for error handling in CCPP               | none       |    0 | character | len=*     | out   | F        | 
!! | errflg           | error_flag                                                 | error flag for error handling in CCPP                  | flag       |    0 | integer   |           | out   | F        |
!!
   subroutine gfdl_cloud_microphys_pre_run(im, rain0, ice0, snow0, graupel0, errmsg, errflg)

      implicit none

      integer, intent(in) :: im
      real(kind_phys), dimension(1:im), intent(out) :: rain0
      real(kind_phys), dimension(1:im), intent(out) :: ice0
      real(kind_phys), dimension(1:im), intent(out) :: snow0
      real(kind_phys), dimension(1:im), intent(out) :: graupel0
      character(len=*),                 intent(out) :: errmsg
      integer,                          intent(out) :: errflg

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      rain0    = 0
      ice0     = 0
      snow0    = 0
      graupel0 = 0

   end subroutine gfdl_cloud_microphys_pre_run

end module gfdl_cloud_microphys_pre


module gfdl_cloud_microphys_post

   use machine, only : kind_phys

   implicit none

   private

   public gfdl_cloud_microphys_post_run, gfdl_cloud_microphys_post_init, gfdl_cloud_microphys_post_finalize

   ! DH* CLEANUP !!!
   real(kind=kind_phys), parameter :: con_p001= 0.001d0
   real(kind=kind_phys), parameter :: con_day = 86400.d0
   real(kind=kind_phys), parameter :: rainmin = 1.0e-13

contains

   subroutine gfdl_cloud_microphys_post_init()
   end subroutine gfdl_cloud_microphys_post_init

   subroutine gfdl_cloud_microphys_post_finalize()
   end subroutine gfdl_cloud_microphys_post_finalize

!! \section arg_table_gfdl_cloud_microphys_post_run Argument Table
!! | local_name       | standard_name                                                         | long_name                                              | units      | rank | type      | kind      | intent| optional |
!! |------------------|-----------------------------------------------------------------------|--------------------------------------------------------|------------|------|-----------|-----------|-------|----------|
!! | im               | horizontal_loop_extent                                                | horizontal loop extent                                 | count      |    0 | integer   |           | in    | F        |
!! | rain0            | lwe_thickness_of_stratiform_precipitation_amount_per_day              | stratiform rain over 24h period                        | mm         |    1 | real      | kind_phys | in    | F        |
!! | ice0             | lwe_thickness_of_ice_amount_per_day                                   | ice fall over 24h period                               | mm         |    1 | real      | kind_phys | in    | F        |
!! | snow0            | lwe_thickness_of_snow_amount_per_day                                  | snow fall over 24h period                              | mm         |    1 | real      | kind_phys | in    | F        |
!! | graupel0         | lwe_thickness_of_graupel_amount_per_day                               | graupel fall over 24h period                           | mm         |    1 | real      | kind_phys | in    | F        |
!! | rain1            | lwe_thickness_of_precipitation_amount_on_dynamics_timestep            | rainfall at this timestep                              | m          |    1 | real      | kind_phys | out   | F        |
!! | ice1             | lwe_thickness_of_ice_amount_on_dynamics_timestep                      | ice fall at this time step                             | m          |    1 | real      | kind_phys | out   | F        |
!! | snow1            | lwe_thickness_of_snow_amount_on_dynamics_timestep                     | snow fall at this time step                            | m          |    1 | real      | kind_phys | out   | F        |
!! | graupel1         | lwe_thickness_of_graupel_amount_on_dynamics_timestep                  | graupel fall at this time step                         | m          |    1 | real      | kind_phys | out   | F        |
!! | sr               | ratio_of_snowfall_to_rainfall                                         | snow ratio: ratio of snow to total precipitation       | frac       |    1 | real      | kind_phys | out   | F        |
!! | dtp              | time_step_for_physics                                                 | physics timestep                                       | s          |    0 | real      | kind_phys | in    | F        |
!! | errmsg           | error_message                                                         | error message for error handling in CCPP               | none       |    0 | character | len=*     | out   | F        | 
!! | errflg           | error_flag                                                            | error flag for error handling in CCPP                  | flag       |    0 | integer   |           | out   | F        |
!!
   subroutine gfdl_cloud_microphys_post_run(im, rain0, ice0, snow0, graupel0, &
                                                rain1, ice1, snow1, graupel1, &
                                            sr, dtp, errmsg, errflg)

      implicit none

      integer, intent(in) :: im
      real(kind_phys), dimension(1:im), intent(in)  :: rain0
      real(kind_phys), dimension(1:im), intent(in)  :: ice0
      real(kind_phys), dimension(1:im), intent(in)  :: snow0
      real(kind_phys), dimension(1:im), intent(in)  :: graupel0
      real(kind_phys), dimension(1:im), intent(out) :: rain1
      real(kind_phys), dimension(1:im), intent(out) :: ice1
      real(kind_phys), dimension(1:im), intent(out) :: snow1
      real(kind_phys), dimension(1:im), intent(out) :: graupel1
      real(kind_phys), dimension(1:im), intent(out) :: sr
      real(kind_phys),                  intent(in)  :: dtp
      character(len=*),                 intent(out) :: errmsg
      integer,                          intent(out) :: errflg

      real(kind=kind_phys) :: tem
      integer :: i

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      tem   = dtp*con_p001/con_day

      do i=1,im
        rain1(i)    = (rain0(i)+snow0(i)+ice0(i)+graupel0(i)) * tem
        ice1(i)     = ice0    (i) * tem
        snow1(i)    = snow0   (i) * tem
        graupel1(i) = graupel0(i) * tem
        if ( rain1(i) > rainmin ) then
          sr(i) = (snow0(i) + ice0(i)  + graupel0(i)) &
                      / (rain0(i) + snow0(i) + ice0(i) + graupel0(i))
        else
          sr(i) = 0.0
        endif
      enddo

   end subroutine gfdl_cloud_microphys_post_run

end module gfdl_cloud_microphys_post
