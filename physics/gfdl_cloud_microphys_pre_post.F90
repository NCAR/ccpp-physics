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


module gfdl_cloud_microphys_post_1

   use machine, only : kind_phys

   implicit none

   private

   public gfdl_cloud_microphys_post_1_run, gfdl_cloud_microphys_post_1_init, gfdl_cloud_microphys_post_1_finalize

   ! DH* CLEANUP !!!
   real(kind=kind_phys), parameter :: con_p001= 0.001d0
   real(kind=kind_phys), parameter :: con_day = 86400.d0
   real(kind=kind_phys), parameter :: rainmin = 1.0e-13

contains

   subroutine gfdl_cloud_microphys_post_1_init()
   end subroutine gfdl_cloud_microphys_post_1_init

   subroutine gfdl_cloud_microphys_post_1_finalize()
   end subroutine gfdl_cloud_microphys_post_1_finalize

!! \section arg_table_gfdl_cloud_microphys_post_1_run Argument Table
!! | local_name       | standard_name                                                         | long_name                                              | units      | rank | type      | kind      | intent| optional |
!! |------------------|-----------------------------------------------------------------------|--------------------------------------------------------|------------|------|-----------|-----------|-------|----------|
!! | im               | horizontal_loop_extent                                                | horizontal loop extent                                 | count      |    0 | integer   |           | in    | F        |
!! | rain0            | lwe_thickness_of_stratiform_precipitation_amount_per_day              | stratiform rain over 24h period                        | mm         |    1 | real      | kind_phys | in    | F        |
!! | ice0             | lwe_thickness_of_ice_amount_per_day                                   | ice fall over 24h period                               | mm         |    1 | real      | kind_phys | in    | F        |
!! | snow0            | lwe_thickness_of_snow_amount_per_day                                  | snow fall over 24h period                              | mm         |    1 | real      | kind_phys | in    | F        |
!! | graupel0         | lwe_thickness_of_graupel_amount_per_day                               | graupel fall over 24h period                           | mm         |    1 | real      | kind_phys | in    | F        |
!! | rain1            | lwe_thickness_of_stratiform_precipitation_amount_on_dynamics_timestep | stratiform rainfall amount on physics timestep         | m          |    1 | real      | kind_phys | out   | F        |
!! | ice1             | lwe_thickness_of_ice_amount_on_dynamics_timestep                      | ice fall at this time step                             | m          |    1 | real      | kind_phys | out   | F        |
!! | snow1            | lwe_thickness_of_snow_amount_on_dynamics_timestep                     | snow fall at this time step                            | m          |    1 | real      | kind_phys | out   | F        |
!! | graupel1         | lwe_thickness_of_graupel_amount_on_dynamics_timestep                  | graupel fall at this time step                         | m          |    1 | real      | kind_phys | out   | F        |
!! | sr               | ratio_of_snowfall_to_rainfall                                         | snow ratio: ratio of snow to total precipitation       | frac       |    1 | real      | kind_phys | out   | F        |
!! | dtp              | time_step_for_physics                                                 | physics timestep                                       | s          |    0 | real      | kind_phys | in    | F        |
!! | errmsg           | error_message                                                         | error message for error handling in CCPP               | none       |    0 | character | len=*     | out   | F        | 
!! | errflg           | error_flag                                                            | error flag for error handling in CCPP                  | flag       |    0 | integer   |           | out   | F        |
!!
   subroutine gfdl_cloud_microphys_post_1_run(im, rain0, ice0, snow0, graupel0, &
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

   end subroutine gfdl_cloud_microphys_post_1_run

end module gfdl_cloud_microphys_post_1


module gfdl_cloud_microphys_post_2

   use machine, only : kind_phys

   implicit none

   private

   public gfdl_cloud_microphys_post_2_run, gfdl_cloud_microphys_post_2_init, gfdl_cloud_microphys_post_2_finalize

contains

   subroutine gfdl_cloud_microphys_post_2_init()
   end subroutine gfdl_cloud_microphys_post_2_init

   subroutine gfdl_cloud_microphys_post_2_finalize()
   end subroutine gfdl_cloud_microphys_post_2_finalize

!! \section arg_table_gfdl_cloud_microphys_post_2_run Argument Table
!! | local_name       | standard_name                                                          | long_name                                              | units      | rank | type      | kind      | intent| optional |
!! |------------------|------------------------------------------------------------------------|--------------------------------------------------------|------------|------|-----------|-----------|-------|----------|
!! | im               | horizontal_loop_extent                                                 | horizontal loop extent                                 | count      |    0 | integer   |           | in    | F        |
!! | tprcp            | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep | total precipitation amount in each time step           | m          |    1 | real      | kind_phys | out   | F        |
!! | rain             | lwe_thickness_of_precipitation_amount_on_dynamics_timestep             | total rain at this time step                           | m          |    1 | real      | kind_phys | in    | F        |
!! | rainc            | lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep  | convective rain at this time step                      | m          |    1 | real      | kind_phys | in    | F        |
!! | srflag           | flag_for_precipitation_type                                            | snow/rain flag for precipitation                       | flag       |    1 | real      | kind_phys | out   | F        |
!! | tsfc             | surface_skin_temperature                                               | ocean surface skin temperature                         | K          |    1 | real      | kind_phys | in    | F        |
!! | rain0            | lwe_thickness_of_stratiform_precipitation_amount_per_day               | stratiform rain over 24h period                        | mm         |    1 | real      | kind_phys | in    | F        |
!! | ice0             | lwe_thickness_of_ice_amount_per_day                                    | ice fall over 24h period                               | mm         |    1 | real      | kind_phys | in    | F        |
!! | snow0            | lwe_thickness_of_snow_amount_per_day                                   | snow fall over 24h period                              | mm         |    1 | real      | kind_phys | in    | F        |
!! | graupel0         | lwe_thickness_of_graupel_amount_per_day                                | graupel fall over 24h period                           | mm         |    1 | real      | kind_phys | in    | F        |
!! | errmsg           | error_message                                                          | error message for error handling in CCPP               | none       |    0 | character | len=*     | out   | F        | 
!! | errflg           | error_flag                                                             | error flag for error handling in CCPP                  | flag       |    0 | integer   |           | out   | F        |
!!
   subroutine gfdl_cloud_microphys_post_2_run(im, tprcp, rain, rainc, srflag, tsfc, &
                                              rain0, ice0, snow0, graupel0,         &
                                              errmsg, errflg)

      implicit none

      ! Interface variables
      integer,                          intent(in)  :: im
      real(kind_phys), dimension(1:im), intent(out) :: tprcp
      real(kind_phys), dimension(1:im), intent(in)  :: rain
      real(kind_phys), dimension(1:im), intent(in)  :: rainc
      real(kind_phys), dimension(1:im), intent(out) :: srflag
      real(kind_phys), dimension(1:im), intent(in)  :: tsfc
      real(kind_phys), dimension(1:im), intent(in)  :: rain0
      real(kind_phys), dimension(1:im), intent(in)  :: ice0
      real(kind_phys), dimension(1:im), intent(in)  :: snow0
      real(kind_phys), dimension(1:im), intent(in)  :: graupel0
      character(len=*),                 intent(out) :: errmsg
      integer,                          intent(out) :: errflg

      ! Local variables
      integer :: i
      real(kind=kind_phys) :: crain, csnow

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! determine convective rain/snow by surface temperature
      ! determine large-scale rain/snow by rain/snow coming out directly from MP
      do i=1,im
         tprcp(i)  = max(0.0, rain(i) )     ! clu: rain -> tprcp
         srflag(i) = 0.                     ! clu: default srflag as 'rain' (i.e. 0)
         if (tsfc(i) .ge. 273.15) then
            crain = rainc(i)
            csnow = 0.0
         else
           crain = 0.0
           csnow = rainc(i)
         endif
         if ((snow0(i)+ice0(i)+graupel0(i)+csnow) > (rain0(i)+crain)) then
           srflag(i) = 1.                   ! clu: set srflag to 'snow' (i.e. 1)
         endif
      enddo

   end subroutine gfdl_cloud_microphys_post_2_run

end module gfdl_cloud_microphys_post_2

