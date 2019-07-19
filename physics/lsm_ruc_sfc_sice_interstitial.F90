module lsm_ruc_sfc_sice_pre

   use machine, only: kind_phys

   implicit none

   private

   public :: lsm_ruc_sfc_sice_pre_init, lsm_ruc_sfc_sice_pre_run, lsm_ruc_sfc_sice_pre_finalize

contains

   subroutine lsm_ruc_sfc_sice_pre_init ()
   end subroutine lsm_ruc_sfc_sice_pre_init

   subroutine lsm_ruc_sfc_sice_pre_finalize ()
   end subroutine lsm_ruc_sfc_sice_pre_finalize

#if 0
!> \section arg_table_lsm_ruc_sfc_sice_pre_run Argument Table
!! | local_name           | standard_name                                                                | long_name                                                       | units         | rank | type      |    kind   | intent | optional |
!! |----------------------|------------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | im                   | horizontal_loop_extent                                                       | horizontal loop extent                                          | count         |    0 | integer   |           | in     | F        |
!! | lsoil_ruc            | soil_vertical_dimension_for_land_surface_model                               | number of soil layers internal to land surface model            | count         |    0 | integer   |           | in     | F        |
!! | lsoil                | soil_vertical_dimension                                                      | soil vertical layer dimension                                   | count         |    0 | integer   |           | in     | F        |
!! | land                 | flag_nonzero_land_surface_fraction                                           | flag indicating presence of some land surface area fraction     | flag          |    1 | logical   |           | in     | F        |
!! | stc                  | soil_temperature                                                             | soil temperature                                                | K             |    2 | real      | kind_phys | inout  | F        |
!! | tslb                 | soil_temperature_for_land_surface_model                                      | soil temperature for land surface model                         | K             |    2 | real      | kind_phys | in     | F        |
!! | errmsg               | ccpp_error_message                                                           | error message for error handling in CCPP                        | none          |    0 | character | len=*     | out    | F        |
!! | errflg               | ccpp_error_flag                                                              | error flag for error handling in CCPP                           | flag          |    0 | integer   |           | out    | F        |
!!
#endif
   subroutine lsm_ruc_sfc_sice_pre_run(im, lsoil_ruc, lsoil, land, stc, tslb, errmsg, errflg)

      implicit none

      ! Interface variables
      integer, intent(in) :: im, lsoil_ruc, lsoil
      logical, dimension(im), intent(in) :: land
!  --- on Noah levels
      real (kind=kind_phys), dimension(im,lsoil), intent(inout) :: stc
!  --- on RUC levels
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(in) :: tslb

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        if (.not.land(i)) then
          do k=1,min(lsoil,lsoil_ruc)
            stc(i,k) = tslb(i,k)
          end do
        end if
      end do

      end subroutine lsm_ruc_sfc_sice_pre_run

end module lsm_ruc_sfc_sice_pre

module lsm_ruc_sfc_sice_post

   use machine, only: kind_phys

   implicit none

   private

   public :: lsm_ruc_sfc_sice_post_init, lsm_ruc_sfc_sice_post_run, lsm_ruc_sfc_sice_post_finalize

contains

   subroutine lsm_ruc_sfc_sice_post_init ()
   end subroutine lsm_ruc_sfc_sice_post_init

   subroutine lsm_ruc_sfc_sice_post_finalize ()
   end subroutine lsm_ruc_sfc_sice_post_finalize

#if 0
!> \section arg_table_lsm_ruc_sfc_sice_post_run Argument Table
!! | local_name           | standard_name                                                                | long_name                                                       | units         | rank | type      |    kind   | intent | optional |
!! |----------------------|------------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | im                   | horizontal_loop_extent                                                       | horizontal loop extent                                          | count         |    0 | integer   |           | in     | F        |
!! | lsoil_ruc            | soil_vertical_dimension_for_land_surface_model                               | number of soil layers internal to land surface model            | count         |    0 | integer   |           | in     | F        |
!! | lsoil                | soil_vertical_dimension                                                      | soil vertical layer dimension                                   | count         |    0 | integer   |           | in     | F        |
!! | land                 | flag_nonzero_land_surface_fraction                                           | flag indicating presence of some land surface area fraction     | flag          |    1 | logical   |           | in     | F        |
!! | stc                  | soil_temperature                                                             | soil temperature                                                | K             |    2 | real      | kind_phys | in     | F        |
!! | tslb                 | soil_temperature_for_land_surface_model                                      | soil temperature for land surface model                         | K             |    2 | real      | kind_phys | inout  | F        |
!! | errmsg               | ccpp_error_message                                                           | error message for error handling in CCPP                        | none          |    0 | character | len=*     | out    | F        |
!! | errflg               | ccpp_error_flag                                                              | error flag for error handling in CCPP                           | flag          |    0 | integer   |           | out    | F        |
!!
#endif
   subroutine lsm_ruc_sfc_sice_post_run(im, lsoil_ruc, lsoil, land, stc, tslb, errmsg, errflg)

      implicit none

      ! Interface variables
      integer, intent(in) :: im, lsoil_ruc, lsoil
      logical, dimension(im), intent(in) :: land
!  --- on Noah levels
      real (kind=kind_phys), dimension(im,lsoil), intent(in) :: stc
!  --- on RUC levels
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: tslb

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        if (.not.land(i)) then
          do k=1,min(lsoil,lsoil_ruc)
            tslb(i,k) = stc(i,k)
          end do
        end if
      end do

      end subroutine lsm_ruc_sfc_sice_post_run

end module lsm_ruc_sfc_sice_post