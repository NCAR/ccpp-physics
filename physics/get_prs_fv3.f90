module get_prs_fv3

   use machine,  only: kind_phys
   use physcons, only: con_fvirt

!--- public declarations
   public get_prs_fv3_init, get_prs_fv3_run, get_prs_fv3_finalize

!--- local variables
   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
   real(kind=kind_phys), parameter :: half = 0.5_kind_phys

contains


!! \section arg_table_get_prs_fv3_init Argument Table
!!
   subroutine get_prs_fv3_init()
   end subroutine get_prs_fv3_init


!! \section arg_table_get_prs_fv3_run Argument Table
!! | local var name | longname                                             | description                                          | units      | rank | type    | kind      | intent | optional |
!! |----------------|------------------------------------------------------|------------------------------------------------------|------------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                               | horizontal loop extent                               | index      | 0    | integer | default   | in     | F        |
!! | ix             | horizontal_dimension                                 | horizontal dimension                                 | index      | 0    | integer | default   | in     | F        |
!! | iy             | horizontal_loop_extent                               | horizontal dimension                                 | index      | 0    | integer | default   | in     | F        |
!! | km             | vertical_dimension                                   | number of vertical layers                            | index      | 0    | integer | default   | in     | F        |
!! | A              | tendency_of_y_wind_due_to_physics                    | meridional wind tendency due to physics              | m s-2      | 2    | real    | kind_phys | inout  | F        |
!! | B              | tendency_of_x_wind_due_to_physics                    | zonal wind tendency due to physics                   | m s-2      | 2    | real    | kind_phys | inout  | F        |
!! | C              | tendency_of_air_temperature_due_to_physics           | air temperature tendency due to physics              | K s-1      | 2    | real    | kind_phys | inout  | F        |
!! | u1             | x_wind                                               | zonal wind                                           | m s-1      | 2    | real    | kind_phys | in     | F        |
!! | v1             | y_wind                                               | meridional wind                                      | m s-1      | 2    | real    | kind_phys | in     | F        |
!! | dt             | time_step_for_physics                                | physics time step                                    | s          | 0    | real    | kind_phys | in     | F        |
!! | cp             | specific_heat_of_dry_air_at_constant_pressure        | specific heat of dry air at constant pressure        | J kg-1 K-1 | 0    | real    | kind_phys | in     | F        |
!! | levr           | number_of_vertical_layers_for_radiation_calculations | number of vertical layers for radiation calculations | index      | 0    | integer | default   | in     | F        |
!! | pgr            | surface_air_pressure                                 | surface pressure                                     | Pa         | 1    | real    | kind_phys | in     | F        |
!! | prsl           | air_pressure                                         | mid-layer pressure                                   | Pa         | 2    | real    | kind_phys | in     | F        |
!! | prslrd0        | pressure_cutoff_for_rayleigh_damping                 | pressure level above which to apply Rayleigh damping | Pa         | 0    | real    | kind_phys | in     | F        |
!! | ral_ts         | time_scale_for_rayleigh_damping                      | time scale for Rayleigh damping                      | d          | 0    | real    | kind_phys | in     | F        |
!!
   subroutine get_prs_fv3_run(ix, levs, ntrac, phii, prsi, tgrs, qgrs, del, del_gz)
     integer, intent(in) :: ix, levs, ntrac
     real(kind=kind_phys), dimension(ix,levs+1),     intent(in)    :: phii
     real(kind=kind_phys), dimension(ix,levs+1),     intent(in)    :: prsi
     real(kind=kind_phys), dimension(ix,levs),       intent(in)    :: tgrs
     real(kind=kind_phys), dimension(ix,levs,ntrac), intent(in)    :: qgrs
     real(kind=kind_phys), dimension(ix,levs),       intent(inout) :: del
     real(kind=kind_phys), dimension(ix,levs+1),     intent(inout) :: del_gz

! SJL: Adjust the geopotential height hydrostatically in a way consistent with FV3 discretization
! del_gz is a temp array recording the old info before (t,q) are adjusted
     do k=1,levs
       do i=1,ix
            del(i,k) = prsi(i,k) - prsi(i,k+1)
         del_gz(i,k) = (phii(i,k+1) - phii(i,k)) /                    &
                        (tgrs(i,k)*(1.+con_fvirt*max(zero,qgrs(i,k,1))))
       enddo
     enddo

   end subroutine get_prs_fv3_run


!! \section arg_table_get_prs_fv3_finalize Argument Table
!!
   subroutine get_prs_fv3_finalize()
   end subroutine get_prs_fv3_finalize


end module get_prs_fv3



module get_phi_fv3

   use machine,  only: kind_phys
   use physcons, only: con_fvirt

!--- public declarations
   public get_phi_fv3_init, get_phi_fv3_run, get_phi_fv3_finalize

!--- local variables
   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
   real(kind=kind_phys), parameter :: half = 0.5_kind_phys

contains

!! \section arg_table_get_phi_fv3_init Argument Table
!!
   subroutine get_phi_fv3_init()
   end subroutine get_phi_fv3_init


!! \section arg_table_get_phi_fv3_run Argument Table
!! | local var name | longname                                             | description                                          | units      | rank | type    | kind      | intent | optional |
!! |----------------|------------------------------------------------------|------------------------------------------------------|------------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                               | horizontal loop extent                               | index      | 0    | integer | default   | in     | F        |
!! | ix             | horizontal_dimension                                 | horizontal dimension                                 | index      | 0    | integer | default   | in     | F        |
!! | iy             | horizontal_loop_extent                               | horizontal dimension                                 | index      | 0    | integer | default   | in     | F        |
!! | km             | vertical_dimension                                   | number of vertical layers                            | index      | 0    | integer | default   | in     | F        |
!! | A              | tendency_of_y_wind_due_to_physics                    | meridional wind tendency due to physics              | m s-2      | 2    | real    | kind_phys | inout  | F        |
!! | B              | tendency_of_x_wind_due_to_physics                    | zonal wind tendency due to physics                   | m s-2      | 2    | real    | kind_phys | inout  | F        |
!! | C              | tendency_of_air_temperature_due_to_physics           | air temperature tendency due to physics              | K s-1      | 2    | real    | kind_phys | inout  | F        |
!! | u1             | x_wind                                               | zonal wind                                           | m s-1      | 2    | real    | kind_phys | in     | F        |
!! | v1             | y_wind                                               | meridional wind                                      | m s-1      | 2    | real    | kind_phys | in     | F        |
!! | dt             | time_step_for_physics                                | physics time step                                    | s          | 0    | real    | kind_phys | in     | F        |
!! | cp             | specific_heat_of_dry_air_at_constant_pressure        | specific heat of dry air at constant pressure        | J kg-1 K-1 | 0    | real    | kind_phys | in     | F        |
!! | levr           | number_of_vertical_layers_for_radiation_calculations | number of vertical layers for radiation calculations | index      | 0    | integer | default   | in     | F        |
!! | pgr            | surface_air_pressure                                 | surface pressure                                     | Pa         | 1    | real    | kind_phys | in     | F        |
!! | prsl           | air_pressure                                         | mid-layer pressure                                   | Pa         | 2    | real    | kind_phys | in     | F        |
!! | prslrd0        | pressure_cutoff_for_rayleigh_damping                 | pressure level above which to apply Rayleigh damping | Pa         | 0    | real    | kind_phys | in     | F        |
!! | ral_ts         | time_scale_for_rayleigh_damping                      | time scale for Rayleigh damping                      | d          | 0    | real    | kind_phys | in     | F        |
!!
   subroutine get_phi_fv3_run(ix, levs, ntrac, gt0, gq0, del_gz, phii, phil)
     integer, intent(in) :: ix, levs, ntrac
     real(kind=kind_phys), dimension(ix,levs),       intent(in)    :: gt0
     real(kind=kind_phys), dimension(ix,levs,ntrac), intent(in)    :: gq0
     real(kind=kind_phys), dimension(ix,levs+1),     intent(inout) :: del_gz
     real(kind=kind_phys), dimension(ix,levs+1),     intent(inout) :: phii
     real(kind=kind_phys), dimension(ix,levs),       intent(inout) :: phil

! SJL: Adjust the heighz hydrostatically in a way consistent with FV3 discretization
     do i=1,ix
        phii(i,1) = zero
     enddo
     do k=1,levs
       do i=1,ix
         del_gz(i,k) = del_gz(i,k)*gt0(i,k) *                          &
     &                 (1.+con_fvirt*max(zero,gq0(i,k,1)))
         phii(i,k+1) = phii(i,k) + del_gz(i,k)
         phil(i,k)   = half*(phii(i,k) + phii(i,k+1))
       enddo
     enddo

   end subroutine get_phi_fv3_run


!! \section arg_table_get_phi_fv3_finalize Argument Table
!!
   subroutine get_phi_fv3_finalize()
   end subroutine get_phi_fv3_finalize


end module get_phi_fv3


