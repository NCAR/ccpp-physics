!> \file GFS_stochastics.f90
!! This file contains code previously in GFS_stochastics_driver.

!>\defgroup gfs_stoch GFS Stochastics Main Module
!! This module
!! @ {
    module GFS_stochastics

      contains

      subroutine GFS_stochastics_init ()
      end subroutine GFS_stochastics_init

      subroutine GFS_stochastics_finalize()
      end subroutine GFS_stochastics_finalize

!> \section arg_table_GFS_stochastics_run Argument Table
!! | local_name     | standard_name                                                             | long_name                                                    | units   | rank | type      |    kind   | intent | optional |
!! |----------------|---------------------------------------------------------------------------|--------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                    | horizontal loop extent                                       | count   |    0 | integer   |           | in     | F        |
!! | km             | vertical_dimension                                                        | number of vertical levels                                    | count   |    0 | integer   |           | in     | F        |
!! | do_sppt        | flag_for_stochastic_surface_physics_perturbations                         | flag for stochastic surface physics perturbations            | flag    |    0 | logical   |           | in     | F        |
!! | use_zmtnblck   | flag_for_mountain_blocking                                                | flag for mountain blocking                                   | flag    |    0 | logical   |           | in     | F        |
!! | do_shum        | flag_for_stochastic_shum_option                                           | flag for stochastic shum option                              | flag    |    0 | logical   |           | in     | F        |
!! | do_skeb        | flag_for_stochastic_skeb_option                                           | flag for stochastic skeb option                              | flag    |    0 | logical   |           | in     | F        |
!! | zmtnblck       | level_of_dividing_streamline                                              | level of the dividing streamline                             | none    |    1 | real      | kind_phys | in     | F        |
!! | sppt_wts       | weights_for_stochastic_sppt_perturbation                                  | weights for stochastic sppt perturbation                     | none    |    2 | real      | kind_phys | inout  | F        |
!! | skebu_wts      | weights_for_stochastic_skeb_perturbation_of_x_wind                        | weights for stochastic skeb perturbation of x wind           | none    |    2 | real      | kind_phys | in     | F        |
!! | skebv_wts      | weights_for_stochastic_skeb_perturbation_of_y_wind                        | weights for stochastic skeb perturbation of y wind           | none    |    2 | real      | kind_phys | in     | F        |
!! | shum_wts       | weights_for_stochastic_shum_perturbation                                  | weights for stochastic shum perturbation                     | none    |    2 | real      | kind_phys | in     | F        |
!! | sppt_wts_inv   | weights_for_stochastic_sppt_perturbation_flipped                          | weights for stochastic sppt perturbation, flipped            | none    |    2 | real      | kind_phys | inout  | F        |
!! | skebu_wts_inv  | weights_for_stochastic_skeb_perturbation_of_x_wind_flipped                | weights for stochastic skeb perturbation of x wind, flipped  | none    |    2 | real      | kind_phys | inout  | F        |
!! | skebv_wts_inv  | weights_for_stochastic_skeb_perturbation_of_y_wind_flipped                | weights for stochastic skeb perturbation of y wind, flipped  | none    |    2 | real      | kind_phys | inout  | F        |
!! | shum_wts_inv   | weights_for_stochastic_shum_perturbation_flipped                          | weights for stochastic shum perturbation, flipped            | none    |    2 | real      | kind_phys | inout  | F        |
!! | diss_est       | dissipation_estimate_of_air_temperature_at_model_layers                   | dissipation estimate model layer mean temperature            | K       |    2 | real      | kind_phys | in     | F        |
!! | ugrs           | x_wind                                                                    | zonal wind                                                   | m s-1   |    2 | real      | kind_phys | in     | F        |
!! | vgrs           | y_wind                                                                    | meridional wind                                              | m s-1   |    2 | real      | kind_phys | in     | F        |
!! | tgrs           | air_temperature                                                           | model layer mean temperature                                 | K       |    2 | real      | kind_phys | in     | F        |
!! | qgrs           | water_vapor_specific_humidity                                             | water vapor specific humidity                                | kg kg-1 |    2 | real      | kind_phys | in     | F        |
!! | gu0            | x_wind_updated_by_physics                                                 | zonal wind updated by physics                                | m s-1   |    2 | real      | kind_phys | inout  | F        |
!! | gv0            | y_wind_updated_by_physics                                                 | meridional wind updated by physics                           | m s-1   |    2 | real      | kind_phys | inout  | F        |
!! | gt0            | air_temperature_updated_by_physics                                        | temperature updated by physics                               | K       |    2 | real      | kind_phys | inout  | F        |
!! | gq0            | water_vapor_specific_humidity_updated_by_physics                          | water vapor specific humidity updated by physics             | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | dtdtr          | tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step | temp. change due to radiative heating per time step          | K       |    2 | real      | kind_phys | in     | F        |
!! | rain           | lwe_thickness_of_precipitation_amount_on_dynamics_timestep                | total rain at this time step                                 | m       |    1 | real      | kind_phys | in     | F        |
!! | rainc          | lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep     | convective rain at this time step                            | m       |    1 | real      | kind_phys | in     | F        |
!! | tprcp          | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep    | total precipitation amount in each time step                 | m       |    1 | real      | kind_phys | inout  | F        |
!! | totprcp        | accumulated_lwe_thickness_of_precipitation_amount                         | accumulated total precipitation                              | m       |    1 | real      | kind_phys | inout  | F        |
!! | cnvprcp        | cumulative_lwe_thickness_of_convective_precipitation_amount               | cumulative convective precipitation                          | m       |    1 | real      | kind_phys | inout  | F        |
!! | totprcpb       | accumulated_lwe_thickness_of_precipitation_amount_in_bucket               | accumulated total precipitation in bucket                    | m       |    1 | real      | kind_phys | inout  | F        |
!! | cnvprcpb       | cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket     | cumulative convective precipitation in bucket                | m       |    1 | real      | kind_phys | inout  | F        |
!! | cplflx         | flag_for_flux_coupling                                                    | flag controlling cplflx collection (default off)             | flag    |    0 | logical   |           | in     | F        |
!! | rain_cpl       | lwe_thickness_of_precipitation_amount_for_coupling                        | total rain precipitation                                     | m       |    1 | real      | kind_phys | inout  | F        |
!! | snow_cpl       | lwe_thickness_of_snow_amount_for_coupling                                 | total snow precipitation                                     | m       |    1 | real      | kind_phys | inout  | F        |
!! | drain_cpl      | tendency_of_lwe_thickness_of_precipitation_amount_for_coupling            | change in rain_cpl (coupling_type)                           | m       |    1 | real      | kind_phys | in     | F        |
!! | dsnow_cpl      | tendency_of_lwe_thickness_of_snow_amount_for_coupling                     | change in show_cpl (coupling_type)                           | m       |    1 | real      | kind_phys | in     | F        |
!! | errmsg         | ccpp_error_message                                                        | error message for error handling in CCPP                     | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                           | error flag for error handling in CCPP                        | flag    |    0 | integer   |           | out    | F        |
!!
!-------------------------------------------------------------------------
! GFS stochastic_driver
!-------------------------------------------------------------------------
!    routine called prior to radiation and physics steps to handle:
!      1) sets up various time/date variables
!      2) sets up various triggers
!      3) defines random seed indices for radiation (in a reproducible way)
!      5) interpolates coefficients for prognostic ozone calculation
!      6) performs surface data cycling via the GFS gcycle routine
!-------------------------------------------------------------------------
      subroutine GFS_stochastics_run (im, km, do_sppt, use_zmtnblck, do_shum, do_skeb,   &
                                      zmtnblck, sppt_wts, skebu_wts, skebv_wts, shum_wts,&
                                      sppt_wts_inv, skebu_wts_inv, skebv_wts_inv,        &
                                      shum_wts_inv, diss_est,                            &
                                      ugrs, vgrs, tgrs, qgrs, gu0, gv0, gt0, gq0, dtdtr, &
                                      rain, rainc, tprcp, totprcp, cnvprcp,              &
                                      totprcpb, cnvprcpb, cplflx,                        &
                                      rain_cpl, snow_cpl, drain_cpl, dsnow_cpl,          &
                                      errmsg, errflg)

         use machine,               only: kind_phys

         implicit none

         integer,                               intent(in)    :: im
         integer,                               intent(in)    :: km
         logical,                               intent(in)    :: do_sppt
         logical,                               intent(in)    :: use_zmtnblck
         logical,                               intent(in)    :: do_shum
         logical,                               intent(in)    :: do_skeb
         real(kind_phys), dimension(1:im),      intent(in)    :: zmtnblck
         ! sppt_wts only allocated if do_sppt == .true.
         real(kind_phys), dimension(:,:),       intent(inout) :: sppt_wts
         ! skebu_wts, skebv_wts only allocated if do_skeb == .true.
         real(kind_phys), dimension(:,:),       intent(in)    :: skebu_wts
         real(kind_phys), dimension(:,:),       intent(in)    :: skebv_wts
         ! shum_wts only allocated if do_shum == .true.
         real(kind_phys), dimension(:,:),       intent(in)    :: shum_wts
         ! inverse/flipped weights are always allocated
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: sppt_wts_inv
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: skebu_wts_inv
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: skebv_wts_inv
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: shum_wts_inv
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: diss_est
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: ugrs
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: vgrs
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: tgrs
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: qgrs
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: gu0
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: gv0
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: gt0
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: gq0
         ! dtdtr only allocated if do_sppt == .true.
         real(kind_phys), dimension(:,:),       intent(in)    :: dtdtr
         real(kind_phys), dimension(1:im),      intent(in)    :: rain
         real(kind_phys), dimension(1:im),      intent(in)    :: rainc
         real(kind_phys), dimension(1:im),      intent(inout) :: tprcp
         real(kind_phys), dimension(1:im),      intent(inout) :: totprcp
         real(kind_phys), dimension(1:im),      intent(inout) :: cnvprcp
         real(kind_phys), dimension(1:im),      intent(inout) :: totprcpb
         real(kind_phys), dimension(1:im),      intent(inout) :: cnvprcpb
         logical,                               intent(in)    :: cplflx
         ! rain_cpl, snow_cpl only allocated if cplflx == .true. or do_sppt == .true.
         real(kind_phys), dimension(:),         intent(inout) :: rain_cpl
         real(kind_phys), dimension(:),         intent(inout) :: snow_cpl
         ! drain_cpl, dsnow_cpl only allocated if do_sppt == .true.
         real(kind_phys), dimension(:),         intent(in)    :: drain_cpl
         real(kind_phys), dimension(:),         intent(in)    :: dsnow_cpl
         character(len=*),                      intent(out)   :: errmsg
         integer,                               intent(out)   :: errflg

         !--- local variables
         integer :: k, i
         real(kind=kind_phys) :: upert, vpert, tpert, qpert, qnew, sppt_vwt

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (do_sppt) then
           do k=1,km
             do i=1,im
               sppt_vwt=1.0
               if (zmtnblck(i).EQ.0.0) then
                  sppt_vwt=1.0
               else 
                  if (k.GT.zmtnblck(i)+2) then
                     sppt_vwt=1.0
                  endif
                  if (k.LE.zmtnblck(i)) then
                     sppt_vwt=0.0
                  endif
                  if (k.EQ.zmtnblck(i)+1) then
                     sppt_vwt=0.333333
                  endif
                  if (k.EQ.zmtnblck(i)+2) then
                     sppt_vwt=0.666667
                  endif
               endif
               if (use_zmtnblck)then
                  sppt_wts(i,k)=(sppt_wts(i,k)-1)*sppt_vwt+1.0
               endif
               sppt_wts_inv(i,km-k+1)=sppt_wts(i,k)

               upert = (gu0(i,k) - ugrs(i,k))   * sppt_wts(i,k)
               vpert = (gv0(i,k) - vgrs(i,k))   * sppt_wts(i,k)
               tpert = (gt0(i,k) - tgrs(i,k) - dtdtr(i,k)) * sppt_wts(i,k)
               qpert = (gq0(i,k) - qgrs(i,k)) * sppt_wts(i,k)

               gu0(i,k)  = ugrs(i,k)+upert
               gv0(i,k)  = vgrs(i,k)+vpert

               !negative humidity check
               qnew = qgrs(i,k)+qpert
               if (qnew >= 1.0e-10) then
                  gq0(i,k) = qnew
                  gt0(i,k) = tgrs(i,k) + tpert + dtdtr(i,k)
               endif
             enddo
           enddo
           ! instantaneous precip rate going into land model at the next time step
           tprcp(:) = sppt_wts(:,15)*tprcp(:)
           totprcp(:) = totprcp(:) + (sppt_wts(:,15) - 1 )*rain(:)
           ! acccumulated total and convective preciptiation
           cnvprcp(:) = cnvprcp(:) + (sppt_wts(:,15) - 1 )*rainc(:)
           ! bucket precipitation adjustment due to sppt
           totprcpb(:) = totprcpb(:) + (sppt_wts(:,15) - 1 )*rain(:)
           cnvprcpb(:) = cnvprcpb(:) + (sppt_wts(:,15) - 1 )*rainc(:)

            if (cplflx) then
               rain_cpl(:) = rain_cpl(:) + (sppt_wts(:,15) - 1.0)*drain_cpl(:)
               snow_cpl(:) = snow_cpl(:) + (sppt_wts(:,15) - 1.0)*dsnow_cpl(:)
            endif
         
         endif

         if (do_shum) then
           do k=1,km
             gq0(:,k) = gq0(:,k)*(1.0 + shum_wts(:,k))
             shum_wts_inv(:,km-k+1) = shum_wts(:,k)
           end do
         endif
         
         if (do_skeb) then
           do k=1,km
             gu0(:,k) = gu0(:,k)+skebu_wts(:,k)*(diss_est(:,k))
             gv0(:,k) = gv0(:,k)+skebv_wts(:,k)*(diss_est(:,k))
             skebu_wts_inv(:,km-k+1) = skebu_wts(:,k)
             skebv_wts_inv(:,km-k+1) = skebv_wts(:,k)
           enddo
         endif

      end subroutine GFS_stochastics_run

    end module GFS_stochastics
!! @}
