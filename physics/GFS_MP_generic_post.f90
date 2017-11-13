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
!!\section arg_table_GFS_MP_generic_post_init Argument Table
!!
      subroutine GFS_MP_generic_post_init     
      end subroutine GFS_MP_generic_post_init


!!\section arg_table_GFS_MP_generic_post_run Argument Table
!!| local var name | longname                                               |description                                               | units       | rank |  type   |   kind    | intent | optional |
!!|----------------|--------------------------------------------------------|----------------------------------------------------------|-------------|------|---------|-----------|--------|----------|
!!|   im           | horizontal_loop_extent                                 | horizontal loop extent, start at 1                       | index       | 0    | integer |           | in     |  F       |
!!|   ix           | horizontal_dimension                                   | horizontal dimension                                     | index       | 0    | integer |           | in     |  F       |
!!|   levs         | vertical_dimension                                     | vertical layer dimension                                 | index       | 0    | integer |           | in     |  F       |           
!!|   dtf          | time_step_for_dynamics                                 | dynamics time step                                       | s           | 0    | real    | kind_phys | in     |  F       |
!!|   del          | air_pressure_difference_between_midlayers              | air pressure difference between midlayers                | Pa          | 2    | real    | kind_phys | in     |  F       |
!!|   lssav        | logical_flag_for_physics_diagnostics                   | logical flag for model physics diagnostics               | flag        | 0    | logical |           | in     |  F       |
!!|   ldiag3d      | logical_flag_for_3D_diagnostics                        | logical flag for 3D diagnostics                          | flag        | 0    | logical |           | in     |  F       |
!!|   rain         | total_rainfall_at_surface                              | instantaneous total precipitation at surface (APCP)      | m           | 1    | real    | kind_phys | in     |  F       |
!!|   frain        | factor_for_centered_difference_scheme                  | dtf/dtp; factor for centered difference scheme correction| none        | 0    | real    | kind_phys | in     |  F       |
!!|   ntcw         | index_for_liquid_cloud_condensate                      | cloud condensate index in tracer array(3)                | index       | 0    | integer |           | in     |  F       |
!!|   ncld         | choice_of_cloud_scheme                                 | choice of cloud scheme(1 for Z-C)                        | none        | 0    | integer |           | in     |  F       |
!!|   cwm          | cloud_condensed_water_specific_humidity                | cloud condensed water specific humidity                  | kg kg-1     | 2    | real    | kind_phys | in     |  F       |
!!|   t            | air_temperature_updated_by_physics                     | layer mean air temperature                               | K           | 2    | real    | kind_phys | in     |  F       |
!!|   q            | water_vapor_specific_humidity                          | water vapor specific humidity                            | kg kg-1     | 2    | real    | kind_phys | in     |  F       |
!!|   dtdt         | air_temperature_before_microphysics_scheme             | air temperature saved before micophysics scheme          | K           | 2    | real    | kind_phys | in     |  F       |
!!|   dqdt1        | specific_humidity_before_microphysics_scheme           | specific humidity saved before microphysics schme        | kg kg-1     | 2    | real    | kind_phys | in     |  F       |
!!|   totprcp      | precipitation_rate_at_surface                          | precipitation rate at surface                            | kg m-2 s-1  | 1    | real    | kind_phys | inout  |  F       |
!!|   dt3dt6       | large_scale_condensate_heating_rate_at_model_layers    | large scale condensate heating rate at model layers      | K s-1       | 2    | real    | kind_phys | inout  |  F       |
!!|   dq3dt4       | large_scale_condensate_moistening_rate_at_model_layers | large scale condensate moistening rate at model layers   | kg kg-1 s-1 | 2    | real    | kind_phys | inout  |  F       |
!!|   pwat         | column_precipitable_water                              | column integrated precipitable water                     | kg m-2      | 1    | real    | kind_phys | out    |  F       |
!!
      subroutine GFS_MP_generic_post_run(im, ix, levs,dtf,del,          &
                 lssav,ldiag3d,rain,frain,ntcw,ncld,cwm,                & !input
                 t,q,dtdt,dqdt1,                                        &
                 totprcp, dt3dt6,dq3dt4,pwat  )     ! output
     
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
                                             cwm, del, dtdt, dqdt1
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
               dt3dt6(i,k) = dt3dt6(i,k) + (t(i,k)-dtdt(i,k)) * frain
               dq3dt4(i,k) = dq3dt4(i,k) + (q(i,k)-dqdt1(i,k)) * frain
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

!!\setction arg_table_GFS_MP_generic_post_finalize Argument Table
!!
      subroutine GFS_MP_generic_post_finalize
      end subroutine GFS_MP_generic_post_finalize
!! @}
      end module GFS_MP_generic_post

