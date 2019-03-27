!>  \file cires_ugwp_post.F90
!! This file contains
module cire_ugwp_post

contains

!>\defgroup cire_ugwp_post CIRES UGWP Scheme Post
!! @{
!> \section arg_table_cire_ugwp_post_init Argument Table
!!
    subroutine cire_ugwp_post_init ()
    end subroutine cire_ugwp_post_init

!>@brief The subroutine initializes the CIRES UGWP
#if 0
!> \section arg_table_cire_ugwp_post_run Argument Table
!! | local_name       | standard_name                                                                  | long_name                                                    | units     | rank |  type     |   kind    | intent | optional |
!! |------------------|--------------------------------------------------------------------------------|--------------------------------------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | Model            | GFS_control_type_instance                                                      | Fortran DDT containing FV3-GFS model control parameters      | DDT       | 0    | GFS_control_type  |   | in     | F        |
!! | Diag             | GFS_diag_type_instance                                                         | Fortran DDT containing FV3-GFS diagnotics data               | DDT       | 0    | GFS_diag_type     |   | inout  | F        |
!! | Stateout         | GFS_stateout_type_instance                                                     | instance of derived type GFS_stateout_type                   | DDT       | 0    | GFS_stateout_type |   | in     | F        |
!! | dtp              | time_step_for_physics                                                          | physics timestep                                             | s         | 0    | real      | kind_phys | in     | F        |
!! | dtf              | time_step_for_dynamics                                                         | dynamics timestep                                            | s         | 0    | real      | kind_phys | none   | F        |
!! | im               | horizontal_loop_extent                                                         | horizontal loop extent                                       | count     | 0    | integer   |           | in     | F        |
!! | levs             | vertical_dimension                                                             | number of vertical levels                                    | count     | 0    | integer   |           | in     | F        |
!! | dusfcg           | instantaneous_x_stress_due_to_gravity_wave_drag                                | zonal surface stress due to orographic gravity wave drag     | Pa        | 1    | real      | kind_phys | none   | F        |
!! | dvsfcg           | instantaneous_y_stress_due_to_gravity_wave_drag                                | meridional surface stress due to orographic gravity wave drag| Pa        | 1    | real      | kind_phys | none   | F        |
!! | gw_dudt          | tendency_of_x_wind_due_to_ugwp                                                 | zonal wind tendency due to UGWP                              | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | gw_dvdt          | tendency_of_y_wind_due_to_ugwp                                                 | meridional wind tendency due to UGWP                         | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | gw_dtdt          | tendency_of_air_temperature_due_to_ugwp                                        | air temperature tendency due to UGWP                         | K s-1     | 2    | real      | kind_phys | none   | F        |
!! | gw_kdis          | eddy_mixing_due_to_ugwp                                                        | eddy mixing due to UGWP                                      | m2 s-1    | 2    | real      | kind_phys | none   | F        |
!! | tau_tofd         | instantaneous_momentum_flux_due_to_turbulent_orographic_form_drag              | momentum flux or stress due to TOFD                          | Pa        | 2    | real      | kind_phys | none   | F        | 
!! | tau_mtb          | instantaneous_momentum_flux_due_to_mountain_blocking_drag                      | momentum flux or stress due to mountain blocking drag        | Pa        | 2    | real      | kind_phys | none   | F        | 
!! | tau_ogw          | instantaneous_momentum_flux_due_to_orographic_gravity_wave_drag                | momentum flux or stress due to orographic gravity wave drag  | Pa        | 2    | real      | kind_phys | none   | F        | 
!! | tau_ngw          | instantaneous_momentum_flux_due_to_nonstationary_gravity_wave                  | momentum flux or stress due to nonstationary gravity waves   | Pa        | 2    | real      | kind_phys | none   | F        | 
!! | zmtb             | height_of_mountain_blocking                                                    | height of mountain blocking drag                             | m         | 1    | real      | kind_phys | none   | F        |
!! | zlwb             | height_of_low_level_wave_breaking                                              | height of low level wave breaking                            | m         | 1    | real      | kind_phys | none   | F        |   
!! | zogw             | height_of_launch_level_of_orographic_gravity_wave                              | height of launch level of orographic gravity wave            | m         | 1    | real      | kind_phys | none   | F        | 
!! | du3dt_mtb        | instantaneous_change_in_x_wind_due_to_mountain_blocking_drag                   | instantaneous change in x wind due to mountain blocking drag | m s-2     | 2    | real      | kind_phys | none   | F        | 
!! | du3dt_ogw        | instantaneous_change_in_x_wind_due_to_orographic_gravity_wave_drag             | instantaneous change in x wind due to orographic gw drag     | m s-2     | 2    | real      | kind_phys | none   | F        | 
!! | du3dt_tms        | instantaneous_change_in_x_wind_due_to_turbulent_orographic_form_drag           | instantaneous change in x wind due to TOFD                   | m s-2     | 2    | real      | kind_phys | none   | F        | 
!! | errmsg           | ccpp_error_message                                                             | error message for error handling in CCPP                     | none      | 0    | character | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                                | error flag for error handling in CCPP                        | flag      | 0    | integer   |           | out    | F        |
!!
#endif


     subroutine cire_ugwp_post_run (Model, Diag, Stateout, dtp, dtf, im, levs,
     &    dusfcg, dvsfcg,
! diag (tendencies due to ugwp) 
     &    gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,
! COORDE diag
     &    tau_tofd, tau_mtb, tau_ogw, tau_ngw,
     &    zmtb, zlwb, zogw,
     &    du3dt_mtb,du3dt_ogw, du3dt_tms,
     &    errmsg, errflg)

        use machine,                only: kind_phys
        use GFS_typedefs,           only: GFS_stateout_type,   &
                                          GFS_control_type,   &
                                          GFS_diag_type

        implicit none

        ! Interface variables
        type(GFS_control_type),   intent(in)    :: Model
        type(GFS_diag_type),      intent(inout) :: Diag
        type(GFS_stateout_type),  intent(in)    :: Stateout

        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtp, dtf

        real(kind=kind_phys), intent(out), dimension(im)      :: dusfcg, dvsfcg
        real(kind=kind_phys), intent(out), dimension(im)      :: zmtb, zlwb, zogw
        real(kind=kind_phys), intent(out), dimension(im)      :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
        real(kind=kind_phys), intent(out), dimension(im, levs):: gw_dudt, gw_dvdt, gw_dtdt, gw_kdis
        real(kind=kind_phys), intent(out), dimension(im, levs):: du3dt_mtb, du3dt_ogw, du3dt_tms

        character(len=*),        intent(out) :: errmsg
        integer,                 intent(out) :: errflg

        ! local variables
        integer :: i, k
        real(kind=kind_phys), parameter :: ftausec = 86400.
        real(kind=kind_phys) :: fdaily, frain
        ! Valery Yudin 2018: process-oriented diagnostics for 3D-fields in UGWP for COORDE
        ! PdXdt: local tendency for X after each physics chain
        ! TdXdt: total tendency for X due to ALL GFS_physics except radiance
        real(kind=kind_phys), dimension(im, levs) :: Pdtdt, Pdudt, Pdvdt, Tdtdt, Tdudt, Tdvdt


        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (.not. (Model%ldiag_ugwp)) return


        !WL* not sure if we need to put them in a _pre?
        Pdudt = 0.
        Pdvdt = 0.
        Pdtdt = 0.

        Tdudt = 0.
        Tdvdt = 0.
        Tdtdt = 0.
        !*WL

        fdaily  = dtp/ftausec
        if (Model%ldiag_ugwp) then
    
            Diag%zmtb =  Diag%zmtb + fdaily *zmtb
            Diag%zlwb =  Diag%zlwb + fdaily *zlwb
            Diag%zogw =  Diag%zogw + fdaily *zogw
    
            Diag%dugwd = Diag%dugwd + fdaily*dusfcg*ftausec
            Diag%dvgwd = Diag%dvgwd + fdaily*dvsfcg*ftausec
        
            Diag%tau_tofd  = Diag%tau_tofd + fdaily *tau_tofd
            Diag%tau_mtb   = Diag%tau_mtb + fdaily *tau_mtb
            Diag%tau_ogw   = Diag%tau_ogw + fdaily *tau_ogw
            Diag%tau_ngw   = Diag%tau_ngw + fdaily *tau_ngw
            Diag%du3dt_mtb = Diag%du3dt_mtb + fdaily *ax_mtb
            Diag%du3dt_tms = Diag%du3dt_tms + fdaily *ax_tms
            Diag%du3dt_ogw = Diag%du3dt_ogw + fdaily *ax_ogw

            Tdudt =     Tdudt + gw_dudt* fdaily
            Tdvdt =     Tdvdt + gw_dvdt* fdaily 
            Tdtdt =     Tdtdt + gw_dvdt* fdaily                             
                
         endif          


        ! Not sure how to add this part: tendencies from multiple parameterizations
        ! Standard accum-Update before "moist physics" by "PBL + GWP + RF" as imposed in GSM
        !                by accumulating   dxdt_pbl +dxdt_gwp + dxdt_rf
        ! do k=1,levs
        !   do i=1,im
        !       Stateout%gt0(i,k)  = Statein%tgrs(i,k) + (dtdt(i,k)+gw_dtdt(i,k)) * dtp
        !       Stateout%gu0(i,k)  = Statein%ugrs(i,k) + (dudt(i,k)+gw_dudt(i,k)) * dtp
        !       Stateout%gv0(i,k)  = Statein%vgrs(i,k) + (dvdt(i,k)+gw_dvdt(i,k)) * dtp
        !   enddo
        ! enddo
        ! Stateout%gq0(1:im,:,:) = Statein%qgrs(1:im,:,:) + dqdt(1:im,:,:) * dtp


        !=======================================================================         
        !     above: updates of the state by UGWP oro-GWS and RF-damp
        !  Diag%tav_ugwp & Diag%uav_ugwp(i,k)-Updated U-T state before moist/micro physics
        !================================================================================ 

        !frain=factor for centered difference scheme correction of rain amount.
        frain = dtf / dtp
        if (Model%ldiag_ugwp) then
            do k=1,levs
                do i=1,im
                    Diag%tav_ugwp(i,k) = Diag%tav_ugwp(i,k) + Stateout%gt0(i,k) * fdaily
                    Diag%uav_ugwp(i,k) = Diag%uav_ugwp(i,k) + Stateout%gu0(i,k) * fdaily
                    !Diag%vav_ogw(i,k)  = Diag%vav_ogw(i,k) +  Stateout%gv0(i,k) * fdaily    
                    PdUdt(i,k) =(Stateout%gu0(i,k)-dudt(i,k)) * frain/dtp
                    PdVdt(i,k) =(Stateout%gv0(i,k)-dVdt(i,k)) * frain/dtp 
                    PdTdt(i,k) =(Stateout%gt0(i,k)-dTdt(i,k)) * frain/dtp 
            
                    !WL* don't quite understand what they are (haven't defined them in typedefs)
                    !Diag%du3dt_moist(i,k) = Diag%du3dt_moist(i,k) + PdUdt(i,k)
                    !Diag%dv3dt_moist(i,k) = Diag%dv3dt_moist(i,k) + PdVdt(i,k)
                    !Diag%dt3dt_moist(i,k) = Diag%dt3dt_moist(i,k) + PdTdt(i,k)
                    !*WL

                    ! Attention : frain and increments
                    Tdudt(i,k) =  Tdudt(i,k) + PdUdt(i,k)*fdaily
                    Tdvdt(i,k) =  Tdvdt(i,k) + PdVdt(i,k)*fdaily  
                    Tdtdt(i,k) =  Tdtdt(i,k) + PdTdt(i,k)*fdaily                              
                enddo
            enddo
      endif



      end subroutine cire_ugwp_post_run

!> \section arg_table_cire_ugwp_post_finalize Argument Table
!!
      subroutine cire_ugwp_post_finalize ()
      end subroutine cire_ugwp_post_finalize

!! @}
end module cire_ugwp_post
