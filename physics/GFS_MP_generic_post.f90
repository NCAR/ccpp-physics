!> \file GFS_MP_generic_post.f90
!! This file contains the subroutines that calculates physics/diagnotics variables
!! after calling microphysics scheme:
!! - totprcp: precipitation rate at surface
!! - dt3dt(:,:,6): large scale condensate heating rate at model layers
!! - dq3dt(:,:,4): large scale condensate moistening rate at model layers
!! - pwat: column integrated precipitable water

      module GFS_MP_generic_post
      contains

!> \defgroup GFS_MP_generic_post GFS MP generic post
!! @{
!! \section arg_table_GFS_MP_generic_post_init Argument Table
!!
      subroutine GFS_MP_generic_post_init
      end subroutine GFS_MP_generic_post_init


!! \section arg_table_GFS_MP_generic_post_run Argument Table
!! | local var name | longname                                                   | description                                                    | units       | rank |  type   |   kind    | intent | optional |
!! |----------------|------------------------------------------------------------|----------------------------------------------------------------|-------------|------|---------|-----------|--------|----------|
!! |   im           | horizontal_loop_extent                                     | horizontal loop extent                                         | count       | 0    | integer |           | in     |  F       |
!! |   ix           | horizontal_dimension                                       | horizontal dimension                                           | count       | 0    | integer |           | in     |  F       |
!! |   levs         | vertical_dimension                                         | vertical layer dimension                                       | count       | 0    | integer |           | in     |  F       |
!! |   dtf          | time_step_for_dynamics                                     | dynamics time step                                             | s           | 0    | real    | kind_phys | in     |  F       |
!! |   del          | air_pressure_difference_between_midlayers                  | air pressure difference between midlayers                      | Pa          | 2    | real    | kind_phys | in     |  F       |
!! |   lssav        | flag_diagnostics                                           | logical flag for model physics diagnostics                     | flag        | 0    | logical |           | in     |  F       |
!! |   ldiag3d      | flag_diagnostics_3D                                        | logical flag for 3D diagnostics                                | flag        | 0    | logical |           | in     |  F       |
!! |   rain         | lwe_thickness_of_precipitation_amount_on_dynamics_timestep | total rainfall amount on dynamics timestep                     | m           | 1    | real    | kind_phys | in     |  F       |
!! |   frain        | dynamics_to_physics_timestep_ratio                         | dtf/dtp, dynamics to physics timestep ratio                    | none        | 0    | real    | kind_phys | in     |  F       |
!! |   ntcw         | index_for_liquid_cloud_condensate                          | cloud condensate index in tracer array(3)                      | index       | 0    | integer |           | in     |  F       |
!! |   ncld         | number_of_hydrometeors                                     | number_of_hydrometeors(1 for Z-C)                              | count       | 0    | integer |           | in     |  F       |
!! |   cwm          | cloud_condensed_water_specific_humidity_updated_by_physics | cloud condensed water specific humidity                        | kg kg-1     | 2    | real    | kind_phys | in     |  F       |
!! |   t            | air_temperature_updated_by_physics                         | layer mean air temperature                                     | K           | 2    | real    | kind_phys | in     |  F       |
!! |   q            | water_vapor_specific_humidity_updated_by_physics           | water vapor specific humidity                                  | kg kg-1     | 2    | real    | kind_phys | in     |  F       |
!! |   save_t       | air_temperature_save                                       | air temperature before entering a physics scheme               | K           | 2    | real    | kind_phys | in     |  F       |
!! |   save_qv      | water_vapor_specific_humidity_save                         | water vapor specific humidity before entering a physics scheme | kg kg-1     | 2    | real    | kind_phys | in     |  F       |
!! |   totprcp      | accumulated_lwe_thickness_of_precipitation_amount          | accumulated total precipitation amount                         | m           | 1    | real    | kind_phys | inout  |  F       |
!! |   dt3dt6       | large_scale_condensate_heating_rate_at_model_layers        | large scale condensate heating rate at model layers            | K s-1       | 2    | real    | kind_phys | inout  |  F       |
!! |   dq3dt4       | large_scale_condensate_moistening_rate_at_model_layers     | large scale condensate moistening rate at model layers         | kg kg-1 s-1 | 2    | real    | kind_phys | inout  |  F       |
!! |   pwat         | column_precipitable_water                                  | column integrated precipitable water                           | kg m-2      | 1    | real    | kind_phys | out    |  F       |
!!
      subroutine GFS_MP_generic_post_run(im,ix,levs,dtf,del,          &
                 lssav,ldiag3d,rain,frain,ntcw,ncld,cwm,              & !input
                 t,q,save_t,save_qv,                                  &
                 totprcp,dt3dt6,dq3dt4,pwat)     ! output

!
      use machine,               only: kind_phys
      use physcons,              only:  con_g

      implicit none
!
!     declare variables.
!
      integer,intent(in)   :: im, ix, levs, ntcw, ncld
      integer              :: ic,i,k
      logical              :: lssav, ldiag3d
      real(kind=kind_phys) :: frain, dtf, tem
      real(kind=kind_phys),dimension(im)           :: work1
      real(kind=kind_phys),dimension(im), intent(in)      :: rain
      real(kind=kind_phys),dimension(ix,levs), intent(in) :: t,q,       &
                                             cwm, del, save_t,          &
                                             save_qv
      real(kind=kind_phys),dimension(im), intent(inout)   :: totprcp
      real(kind=kind_phys),dimension(im), intent(out)     :: pwat
      real(kind=kind_phys),dimension(ix,levs), intent(inout)  ::        &
                                                        dt3dt6,dq3dt4
!     CONSTANT PARAMETERS
      real(kind=kind_phys), parameter :: onebg   = 1.0/con_g

      if (lssav) then
        do i = 1, im
           totprcp(i) = totprcp(i) + rain(i)
        enddo

        if (ldiag3d) then
          do i = 1, im
            do k = 1,levs
               !CCPP dt3dt6(i,k) = dt3dt6(i,k) + (t(i,k)-dtdt(i,k)) * frain
               !CCPP dq3dt4(i,k) = dq3dt4(i,k) + (q(i,k)-dqdt1(i,k)) * frain
              dt3dt6(i,k) = dt3dt6(i,k) + (t(i,k)-save_t(i,k)) * frain
              dq3dt4(i,k) = dq3dt4(i,k) + (q(i,k)-save_qv(i,k)) * frain

            enddo
          enddo
        endif
      endif

!  --- ...  calculate column precipitable water "pwat"
      tem = dtf * 0.03456 / 86400.0
      do i = 1, im
        pwat(i) = 0.0
        !tem = dtf * 0.03456 / 86400.0

        do k = 1, levs
           work1(i) = 0.0
             !if (ncld > 0) then
             !do ic = ntcw, ntcw+ncld-1
             !  work1(i) = work1(i) +  Stateout%gq0(i,k,ic)
           work1(i) = work1(i) + cwm(i,k)
             !enddo
             !endif
           pwat(i) = pwat(i) + del(i,k)*(q(i,k)+work1(i))
         enddo
         pwat(i) = pwat(i) * onebg

      enddo

      !deallocate (clw)

      end subroutine GFS_MP_generic_post_run

!! \section arg_table_GFS_MP_generic_post_finalize Argument Table
!!
      subroutine GFS_MP_generic_post_finalize
      end subroutine GFS_MP_generic_post_finalize
!! @}
      end module GFS_MP_generic_post
