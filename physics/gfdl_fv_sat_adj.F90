!***********************************************************************
!*                   GNU Lesser General Public License                 
!*
!* This file is part of the GFDL Cloud Microphysics.
!*
!* The GFDL Cloud Microphysics is free software: you can 
!8 redistribute it and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or 
!* (at your option) any later version.
!*
!* The GFDL Cloud Microphysics is distributed in the hope it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty 
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the GFDL Cloud Microphysics.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'fv_cmp' implements the fast procesesses in the GFDL
!! microphysics
!!>@author Shian-Jiann Lin, Linjiong Zhou
! Fast saturation adjustment is part of the gfdl cloud microphysics
! =======================================================================

module fv_sat_adj 
! Modules Included:
! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>rvgas, rdgas, grav, hlv, hlf, cp_air</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td> r_grid</td>
!   </tr>
!   <tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>is_master</td>
!   </tr>
!   <tr>
!     <td>gfdl_cloud_microphys_mod</td>
!     <td>ql_gen, qi_gen, qi0_max, ql_mlt, ql0_max, qi_lim, qs_mlt,
!         tau_r2g, tau_smlt, tau_i2s, tau_v2l, tau_l2v, tau_imlt, tau_l2r,
!         rad_rain, rad_snow, rad_graupel, dw_ocean, dw_land</td>
!   </tr>
! </table>
    
    use constants_mod, only: rvgas, rdgas, grav, hlv, hlf, cp_air
    !use fv_mp_mod, only: is_master
    !use fv_arrays_mod, only: r_grid
    use machine,   only: r8 => kind_phys
    use gfdl_cloud_microphys_mod, only: ql_gen, qi_gen, qi0_max, ql_mlt, ql0_max, qi_lim, qs_mlt
    use gfdl_cloud_microphys_mod, only: icloud_f, sat_adj0, t_sub, cld_min
    use gfdl_cloud_microphys_mod, only: tau_r2g, tau_smlt, tau_i2s, tau_v2l, tau_l2v, tau_imlt, tau_l2r
    use gfdl_cloud_microphys_mod, only: rad_rain, rad_snow, rad_graupel, dw_ocean, dw_land
    use qs_init, only:wqs2_vect, iqs2, iqs1, wqs1
    
    implicit none
    
    private
    
    public fv_sat_adj_run
    
    ! real, parameter :: cp_air = cp_air ! 1004.6, heat capacity of dry air at constant pressure, come from constants_mod
    real, parameter :: cp_vap = 4.0 * rvgas !< 1846.0, heat capacity of water vapor at constant pressure
    real, parameter :: cv_air = cp_air - rdgas !< 717.55, heat capacity of dry air at constant volume
    real, parameter :: cv_vap = 3.0 * rvgas !< 1384.5, heat capacity of water vapor at constant volume
    
    ! http: / / www.engineeringtoolbox.com / ice - thermal - properties - d_576.html
    ! c_ice = 2050.0 at 0 deg c
    ! c_ice = 1972.0 at - 15 deg c
    ! c_ice = 1818.0 at - 40 deg c
    ! http: / / www.engineeringtoolbox.com / water - thermal - properties - d_162.html
    ! c_liq = 4205.0 at 4 deg c
    ! c_liq = 4185.5 at 15 deg c
    ! c_liq = 4178.0 at 30 deg c
    
    ! real, parameter :: c_ice = 2106.0 ! ifs: heat capacity of ice at 0 deg c
    ! real, parameter :: c_liq = 4218.0 ! ifs: heat capacity of liquid at 0 deg c
    real, parameter :: c_ice = 1972.0 !< gfdl: heat capacity of ice at - 15 deg c
    real, parameter :: c_liq = 4185.5 !< gfdl: heat capacity of liquid at 15 deg c
    
    real, parameter :: dc_vap = cp_vap - c_liq !< - 2339.5, isobaric heating / cooling
    real, parameter :: dc_ice = c_liq - c_ice !< 2213.5, isobaric heating / colling
    
    real, parameter :: tice = 273.16 !< freezing temperature
    real, parameter :: t_wfr = tice - 40. !< homogeneous freezing temperature
    
    real, parameter :: lv0 = hlv - dc_vap * tice !< 3.13905782e6, evaporation latent heat coefficient at 0 deg k
    real, parameter :: li00 = hlf - dc_ice * tice !< - 2.7105966e5, fusion latent heat coefficient at 0 deg k
    
    ! real (kind = r_grid), parameter :: e00 = 610.71 ! gfdl: saturation vapor pressure at 0 deg c
    real (r8), parameter :: e00 = 611.21 !< ifs: saturation vapor pressure at 0 deg c
    
    real (r8), parameter :: d2ice = dc_vap + dc_ice !< - 126, isobaric heating / cooling
    real (r8), parameter :: li2 = lv0 + li00 !< 2.86799816e6, sublimation latent heat coefficient at 0 deg k
    
    real, parameter :: lat2 = (hlv + hlf) ** 2 !< used in bigg mechanism
    
    real :: d0_vap !< the same as dc_vap, except that cp_vap can be cp_vap or cv_vap
    real :: lv00 !< the same as lv0, except that cp_vap can be cp_vap or cv_vap
    
    real, allocatable :: table (:), table2 (:), tablew (:), des2 (:), desw (:)
    
    logical :: mp_initialized = .false.

    
contains

    subroutine fv_sat_adj_init ()
    end subroutine fv_sat_adj_init

    subroutine fv_sat_adj_finalize ()
    end subroutine fv_sat_adj_finalize


!>@brief The subroutine 'fv_sat_adj' performs the fast processes in the GFDL microphysics.
!>@details This is designed for single-moment 6-class cloud microphysics schemes.
!! It handles the heat release due to in situ phase changes.
!!
!! \section arg_table_gfdl_fv_sat_adj_run Argument Table
!!| local_name     | standard_name                                                 | long_name                                                                              | units   | rank | type      |   kind    | intent | optional |
!!|----------------|---------------------------------------------------------------|----------------------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!!| mdt            | time_step_for_remapping                                       | remapping time step                                                                    |   s     |    0 | real      |           | in     | F        |
!!| zvir           | moist_gas_constant                                            | zvir=rv/rd-1.0                                                                         |         |    0 | real      |           | in     | F        |
!!| is             | starting_x_direction_index                                    | starting X direction index                                                             | count   |    0 | integer   |           | in     | F        |
!!| ie             | ending_x_direction_index                                      | ending X direction index                                                               | count   |    0 | integer   |           | in     | F        |
!!| js             | starting_y_direction_index                                    | starting Y direction index                                                             | count   |    0 | integer   |           | in     | F        |
!!| je             | ending_y_direction_index                                      | ending Y direction index                                                               | count   |    0 | integer   |           | in     | F        |
!!| ng             | number_of_ghost_zones                                         | number of ghost zones defined in fv_mp                                                 | count   |    0 | integer   |           | in     | F        |
!!| hydrostatic    | flag_for_hydrostatic_solver                                   | flag for use the hydrostatic or nonhydrostatic solver                                  |         |    0 | logical   |           | in     | F        |
!!| consv_te       | flag_for_fast_microphysics_conservation                       | defined in L2E, initerstitial                                                          |         |    0 | logical   |           | in     | F        |
!!| te0            | atmosphere_energy_content_at_Lagrangian_surface               | atmosphere total energy at Lagrangian surface                                          |  J m-2  |    2 | real      |           | out    | F        |                             
!!| qv             | water_vapor_specific_humidity_at_Lagrangian_surface           | water vapor specific humidity updated by fast physics at Lagrangian surface            | Kg Kg-1 |    2 | real      |           | inout  | F        |
!!| ql             | cloud_liquid_water_specific_humidity_at_Lagrangian_surface    | cloud liquid water specific humidity updated by fast physics at Lagrangian surface     | Kg Kg-1 |    2 | real      |           | inout  | F        |
!!| qi             | cloud_ice_specific_humidity_at_Lagrangian_surface             | cloud ice specific humidity updated by fast physics at Lagrangian surface              | Kg Kg-1 |    2 | real      |           | inout  | F        |
!!| qr             | cloud_rain_specific_humidity_at_Lagrangian_surface            | cloud rain specific humidity updated by fast physics at Lagrangian surface             | Kg Kg-1 |    2 | real      |           | inout  | F        |
!!| qs             | cloud_snow_specific_humidity_at_Lagrangian_surface            | cloud snow specific humidity updated by fast physics at Lagrangian surface             | Kg Kg-1 |    2 | real      |           | inout  | F        |   
!!| qg             | cloud_graupel_specific_humidity_at_Lagrangian_surface         | cloud graupel specific humidity updated by fast physics at Lagrangian surface          | Kg Kg-1 |    2 | real      |           | inout  | F        |
!!| hs             | surface_geopotential                                          | surface geopotential                                                                   | m2s-2   |    2 | real      |           | in     | F        |                         
!!| dpln           | log_pressure_depth_at_Lagrangian_surface                      | depth of logarithm pressure                                                            |   Pa    |    2 | real      |           | in     | F        |
!!| delz           | depth_at_Lagrangian_surface                                   | depth at Lagrangian_surface                                                            |  m      |    2 | real      |           | in     | F        | 
!!| pt             | virtual_temperature_at_Lagrangian_surface                     | virtual temperature at Lagrangian surface                                              |  K      |    2 | real      |           | inout  | F        |
!!| dp             | pressure_thickness_at_Lagrangian_surface                      | pressure thickness at Lagrangian surface                                               |  Pa     |    2 | real      |           | in     | F        |                  
!!| q_con          | cloud_condensed_water_specific_humidity_at_Lagrangian_surface | cloud condensed water specific humidity updated by fast physics at Lagrangian surface  | Kg Kg-1 |    2 | real      |           | inout  | F        |
!!| cappa          | cappa_moist_gas_constant_at_Lagrangian_surface                | cappa(i,j,k) = rdgas / ( rdgas +  cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )                   |         |    2 | real      |           | inout  | F        |
!!| area           | grid_area                                                     | Atm%gridstruct%area_64                                                                 |         |    2 | real      | r8        | in     | F        |
!!| dtdt           | tendency_of_air_temperature_at_Lagrangian_surface             | air temperature tendency due to fast physics at Lagrangian surface                     | K s-1   |    2 | real      |           | inout  | F        |
!!| out_dt         | flag_for_tendency_of_air_temperature_at_Lagrangian_surface    | flag for calculating tendency of air temperature due to fast physics                   |         |    0 | logical   |           | in     | F        |    
!!| last_step      | flag_for_the_last_step_of_k_split_remapping                   | flag for the last step of k-split remapping                                            |         |    0 | logical   |           | in     | F        |                                                                         |         |    0 | logical   |           | in     | F        |
!!| do_qa          | flag_for_inline_cloud_fraction_calculation                    | flag for the inline cloud fraction calculation                                         |         |    0 | logical   |           | in     | F        |
!!| qa             | cloud_fraction_at_Lagrangian_surface                          | cloud fraction at Lagrangian surface                                                   |         |    2 | real      |           | out    | F        |
!!
subroutine fv_sat_adj_run (mdt, zvir, is, ie, js, je, ng, hydrostatic, consv_te, &
        te0, qv, ql, qi, qr, qs, qg, hs, dpln, delz, pt, dp, q_con, cappa, &
        area, dtdt, out_dt, last_step, do_qa, qa)
    
    implicit none
    
    integer, intent (in) :: is, ie, js, je, ng
    
    logical, intent (in) :: hydrostatic, consv_te, out_dt, last_step, do_qa
    
    real, intent (in) :: zvir, mdt ! remapping time step
    
    real, intent (in), dimension (is - ng:ie + ng, js - ng:je + ng) :: dp, delz, hs
    real, intent (in), dimension (is:ie, js:je) :: dpln
    
    real, intent (inout), dimension (is - ng:ie + ng, js - ng:je + ng) :: pt, qv, ql, qi, qr, qs, qg
    real, intent (inout), dimension (is - ng:, js - ng:) :: q_con, cappa
    real, intent (inout), dimension (is:ie, js:je) :: dtdt
    
    real, intent (out), dimension (is - ng:ie + ng, js - ng:je + ng) :: qa, te0
    
    real (r8), intent (in), dimension (is - ng:ie + ng, js - ng:je + ng) :: area
    
    real, dimension (is:ie) :: wqsat, dq2dt, qpz, cvm, t0, pt1, qstar
    real, dimension (is:ie) :: icp2, lcp2, tcp2, tcp3
    real, dimension (is:ie) :: den, q_liq, q_sol, q_cond, src, sink, hvar
    real, dimension (is:ie) :: mc_air, lhl, lhi
    
    real :: qsw, rh
    real :: tc, qsi, dqsdt, dq, dq0, pidep, qi_crt, tmp, dtmp
    real :: tin, rqi, q_plus, q_minus
    real :: sdt, dt_bigg, adj_fac
    real :: fac_smlt, fac_r2g, fac_i2s, fac_imlt, fac_l2r, fac_v2l, fac_l2v
    real :: factor, qim, tice0, c_air, c_vap, dw
    
    integer :: i, j
    
    sdt = 0.5 * mdt ! half remapping time step
    dt_bigg = mdt ! bigg mechinism time step
    
    tice0 = tice - 0.01 ! 273.15, standard freezing temperature
    
    ! -----------------------------------------------------------------------
    ! define conversion scalar / factor
    ! -----------------------------------------------------------------------
    fac_i2s = 1. - exp (- mdt / tau_i2s)
    fac_v2l = 1. - exp (- sdt / tau_v2l)
    fac_r2g = 1. - exp (- mdt / tau_r2g)
    fac_l2r = 1. - exp (- mdt / tau_l2r)
    
    fac_l2v = 1. - exp (- sdt / tau_l2v)
    fac_l2v = min (sat_adj0, fac_l2v)
    
    fac_imlt = 1. - exp (- sdt / tau_imlt)
    fac_smlt = 1. - exp (- mdt / tau_smlt)
    
    ! -----------------------------------------------------------------------
    ! define heat capacity of dry air and water vapor based on hydrostatical property
    ! -----------------------------------------------------------------------
    
    if (hydrostatic) then
        c_air = cp_air
        c_vap = cp_vap
    else
        c_air = cv_air
        c_vap = cv_vap
    endif
    d0_vap = c_vap - c_liq
    lv00 = hlv - d0_vap * tice
    ! dc_vap = cp_vap - c_liq ! - 2339.5
    ! d0_vap = cv_vap - c_liq ! - 2801.0
    
    do j = js, je ! start j loop
        
        do i = is, ie
            q_liq (i) = ql (i, j) + qr (i, j)
            q_sol (i) = qi (i, j) + qs (i, j) + qg (i, j)
            qpz (i) = q_liq (i) + q_sol (i)
#ifdef USE_COND
            pt1 (i) = pt (i, j) / ((1 + zvir * qv (i, j)) * (1 - qpz (i)))
#else
            pt1 (i) = pt (i, j) / (1 + zvir * qv (i, j))
#endif
            t0 (i) = pt1 (i) ! true temperature
            qpz (i) = qpz (i) + qv (i, j) ! total_wat conserved in this routine
        enddo
        
        ! -----------------------------------------------------------------------
        ! define air density based on hydrostatical property
        ! -----------------------------------------------------------------------
        
        if (hydrostatic) then
            do i = is, ie
                den (i) = dp (i, j) / (dpln (i, j) * rdgas * pt (i, j))
            enddo
        else
            do i = is, ie
                den (i) = - dp (i, j) / (grav * delz (i, j)) ! moist_air density
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! define heat capacity and latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            mc_air (i) = (1. - qpz (i)) * c_air ! constant
            cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! fix energy conservation
        ! -----------------------------------------------------------------------
        if (consv_te) then
            if (hydrostatic) then
                do i = is, ie
                    te0 (i, j) = - c_air * t0 (i)
                enddo
            else
                do i = is, ie
#ifdef USE_COND
                    te0 (i, j) = - cvm (i) * t0 (i)
#else
                    te0 (i, j) = - c_air * t0 (i)
#endif
                enddo
            endif
        endif
        
        ! -----------------------------------------------------------------------
        ! fix negative cloud ice with snow
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qi (i, j) < 0.) then
                qs (i, j) = qs (i, j) + qi (i, j)
                qi (i, j) = 0.
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! melting of cloud ice to cloud water and rain
        ! -----------------------------------------------------------------------
        do i = is, ie
            if (qi (i, j) > 1.e-8 .and. pt1 (i) > tice) then
                sink (i) = min (qi (i, j), fac_imlt * (pt1 (i) - tice) / icp2 (i))
                qi (i, j) = qi (i, j) - sink (i)
                ! sjl, may 17, 2017
                ! tmp = min (sink (i), dim (ql_mlt, ql (i, j))) ! max ql amount
                ! ql (i, j) = ql (i, j) + tmp
                ! qr (i, j) = qr (i, j) + sink (i) - tmp
                ! sjl, may 17, 2017
                ql (i, j) = ql (i, j) + sink (i)
                q_liq (i) = q_liq (i) + sink (i)
                q_sol (i) = q_sol (i) - sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) - sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! fix negative snow with graupel or graupel with available snow
        ! -----------------------------------------------------------------------
        do i = is, ie
            if (qs (i, j) < 0.) then
                qg (i, j) = qg (i, j) + qs (i, j)
                qs (i, j) = 0.
            elseif (qg (i, j) < 0.) then
                tmp = min (- qg (i, j), max (0., qs (i, j)))
                qg (i, j) = qg (i, j) + tmp
                qs (i, j) = qs (i, j) - tmp
            endif
        enddo
        
        ! after this point cloud ice & snow are positive definite
        
        ! -----------------------------------------------------------------------
        ! fix negative cloud water with rain or rain with available cloud water
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (ql (i, j) < 0.) then
                tmp = min (- ql (i, j), max (0., qr (i, j)))
                ql (i, j) = ql (i, j) + tmp
                qr (i, j) = qr (i, j) - tmp
            elseif (qr (i, j) < 0.) then
                tmp = min (- qr (i, j), max (0., ql (i, j)))
                ql (i, j) = ql (i, j) - tmp
                qr (i, j) = qr (i, j) + tmp
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! enforce complete freezing of cloud water to cloud ice below - 48 c
        ! -----------------------------------------------------------------------

        do i = is, ie
            dtmp = tice - 48. - pt1 (i)
            if (ql (i, j) > 0. .and. dtmp > 0.) then
                sink (i) = min (ql (i, j), dtmp / icp2 (i))
                ql (i, j) = ql (i, j) - sink (i)
                qi (i, j) = qi (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
       
        do i = is, ie
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lhi (i) = li00 + dc_ice * pt1 (i)
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
            tcp3 (i) = lcp2 (i) + icp2 (i) * min (1., dim (tice, pt1 (i)) / 48.)
        enddo
        
        ! -----------------------------------------------------------------------
        ! condensation / evaporation between water vapor and cloud water
        ! -----------------------------------------------------------------------
        
        call wqs2_vect (is, ie, pt1, den, wqsat, dq2dt)
        
        adj_fac = sat_adj0
        do i = is, ie
            dq0 = (qv (i, j) - wqsat (i)) / (1. + tcp3 (i) * dq2dt (i))
            if (dq0 > 0.) then ! whole grid - box saturated
                src (i) = min (adj_fac * dq0, max (ql_gen - ql (i, j), fac_v2l * dq0))
            else ! evaporation of ql
                ! sjl 20170703 added ql factor to prevent the situation of high ql and rh<1
                ! factor = - min (1., fac_l2v * sqrt (max (0., ql (i, j)) / 1.e-5) * 10. * (1. - qv (i, j) / wqsat (i)))
                ! factor = - fac_l2v
                ! factor = - 1
                factor = - min (1., fac_l2v * 10. * (1. - qv (i, j) / wqsat (i))) ! the rh dependent factor = 1 at 90%
                src (i) = - min (ql (i, j), factor * dq0)
            endif
            qv (i, j) = qv (i, j) - src (i)
            ql (i, j) = ql (i, j) + src (i)
            q_liq (i) = q_liq (i) + src (i)
            cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
            pt1 (i) = pt1 (i) + src (i) * lhl (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lhi (i) = li00 + dc_ice * pt1 (i)
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
            tcp3 (i) = lcp2 (i) + icp2 (i) * min (1., dim (tice, pt1 (i)) / 48.)
        enddo
        
        if (last_step) then
            
            ! -----------------------------------------------------------------------
            ! condensation / evaporation between water vapor and cloud water, last time step
            ! enforce upper (no super_sat) & lower (critical rh) bounds
            ! final iteration:
            ! -----------------------------------------------------------------------
            
            call wqs2_vect (is, ie, pt1, den, wqsat, dq2dt)
            
            do i = is, ie
                dq0 = (qv (i, j) - wqsat (i)) / (1. + tcp3 (i) * dq2dt (i))
                if (dq0 > 0.) then ! remove super - saturation, prevent super saturation over water
                    src (i) = dq0
                else ! evaporation of ql
                    ! factor = - min (1., fac_l2v * sqrt (max (0., ql (i, j)) / 1.e-5) * 10. * (1. - qv (i, j) / wqsat (i))) ! the rh dependent factor = 1 at 90%
                    ! factor = - fac_l2v
                    ! factor = - 1
                    factor = - min (1., fac_l2v * 10. * (1. - qv (i, j) / wqsat (i))) ! the rh dependent factor = 1 at 90%
                    src (i) = - min (ql (i, j), factor * dq0)
                endif
                adj_fac = 1.
                qv (i, j) = qv (i, j) - src (i)
                ql (i, j) = ql (i, j) + src (i)
                q_liq (i) = q_liq (i) + src (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + src (i) * lhl (i) / cvm (i)
            enddo
            
            ! -----------------------------------------------------------------------
            ! update latend heat coefficient
            ! -----------------------------------------------------------------------
            
            do i = is, ie
                lhl (i) = lv00 + d0_vap * pt1 (i)
                lhi (i) = li00 + dc_ice * pt1 (i)
                lcp2 (i) = lhl (i) / cvm (i)
                icp2 (i) = lhi (i) / cvm (i)
            enddo
            
        endif
        
        ! -----------------------------------------------------------------------
        ! homogeneous freezing of cloud water to cloud ice
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = t_wfr - pt1 (i) ! [ - 40, - 48]
            if (ql (i, j) > 0. .and. dtmp > 0.) then
                sink (i) = min (ql (i, j), ql (i, j) * dtmp * 0.125, dtmp / icp2 (i))
                ql (i, j) = ql (i, j) - sink (i)
                qi (i, j) = qi (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! bigg mechanism (heterogeneous freezing of cloud water to cloud ice)
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            tc = tice0 - pt1 (i)
            if (ql (i, j) > 0.0 .and. tc > 0.) then
                sink (i) = 3.3333e-10 * dt_bigg * (exp (0.66 * tc) - 1.) * den (i) * ql (i, j) ** 2
                sink (i) = min (ql (i, j), tc / icp2 (i), sink (i))
                ql (i, j) = ql (i, j) - sink (i)
                qi (i, j) = qi (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! freezing of rain to graupel
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = (tice - 0.1) - pt1 (i)
            if (qr (i, j) > 1.e-7 .and. dtmp > 0.) then
                tmp = min (1., (dtmp * 0.025) ** 2) * qr (i, j) ! no limit on freezing below - 40 deg c
                sink (i) = min (tmp, fac_r2g * dtmp / icp2 (i))
                qr (i, j) = qr (i, j) - sink (i)
                qg (i, j) = qg (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! melting of snow to rain or cloud water
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = pt1 (i) - (tice + 0.1)
            if (qs (i, j) > 1.e-7 .and. dtmp > 0.) then
                tmp = min (1., (dtmp * 0.1) ** 2) * qs (i, j) ! no limter on melting above 10 deg c
                sink (i) = min (tmp, fac_smlt * dtmp / icp2 (i))
                tmp = min (sink (i), dim (qs_mlt, ql (i, j))) ! max ql due to snow melt
                qs (i, j) = qs (i, j) - sink (i)
                ql (i, j) = ql (i, j) + tmp
                qr (i, j) = qr (i, j) + sink (i) - tmp
                ! qr (i, j) = qr (i, j) + sink (i)
                q_liq (i) = q_liq (i) + sink (i)
                q_sol (i) = q_sol (i) - sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) - sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! autoconversion from cloud water to rain
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (ql (i, j) > ql0_max) then
                sink (i) = fac_l2r * (ql (i, j) - ql0_max)
                qr (i, j) = qr (i, j) + sink (i)
                ql (i, j) = ql (i, j) - sink (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
            tcp2 (i) = lcp2 (i) + icp2 (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! sublimation / deposition between water vapor and cloud ice
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            src (i) = 0.
            if (pt1 (i) < t_sub) then ! too cold to be accurate; freeze qv as a fix
                src (i) = dim (qv (i, j), 1.e-6)
            elseif (pt1 (i) < tice0) then
                qsi = iqs2 (pt1 (i), den (i), dqsdt)
                dq = qv (i, j) - qsi
                sink (i) = adj_fac * dq / (1. + tcp2 (i) * dqsdt)
                if (qi (i, j) > 1.e-8) then
                    pidep = sdt * dq * 349138.78 * exp (0.875 * log (qi (i, j) * den (i))) &
                         / (qsi * den (i) * lat2 / (0.0243 * rvgas * pt1 (i) ** 2) + 4.42478e4)
                else
                    pidep = 0.
                endif
                if (dq > 0.) then ! vapor - > ice
                    tmp = tice - pt1 (i)
                    qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den (i)
                    src (i) = min (sink (i), max (qi_crt - qi (i, j), pidep), tmp / tcp2 (i))
                else
                    pidep = pidep * min (1., dim (pt1 (i), t_sub) * 0.2)
                    src (i) = max (pidep, sink (i), - qi (i, j))
                endif
            endif
            qv (i, j) = qv (i, j) - src (i)
            qi (i, j) = qi (i, j) + src (i)
            q_sol (i) = q_sol (i) + src (i)
            cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
            pt1 (i) = pt1 (i) + src (i) * (lhl (i) + lhi (i)) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! virtual temp updated
        ! -----------------------------------------------------------------------
        
        do i = is, ie
#ifdef USE_COND
            q_con (i, j) = q_liq (i) + q_sol (i)
            tmp = 1. + zvir * qv (i, j)
            pt (i, j) = pt1 (i) * tmp * (1. - q_con (i, j))
            tmp = rdgas * tmp
            cappa (i, j) = tmp / (tmp + cvm (i))
#else
            pt (i, j) = pt1 (i) * (1. + zvir * qv (i, j))
#endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! fix negative graupel with available cloud ice
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qg (i, j) < 0.) then
                tmp = min (- qg (i, j), max (0., qi (i, j)))
                qg (i, j) = qg (i, j) + tmp
                qi (i, j) = qi (i, j) - tmp
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! autoconversion from cloud ice to snow
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            qim = qi0_max / den (i)
            if (qi (i, j) > qim) then
                sink (i) = fac_i2s * (qi (i, j) - qim)
                qi (i, j) = qi (i, j) - sink (i)
                qs (i, j) = qs (i, j) + sink (i)
            endif
        enddo
        
        if (out_dt) then
            do i = is, ie
                dtdt (i, j) = dtdt (i, j) + pt1 (i) - t0 (i)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! fix energy conservation
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            do i = is, ie
                if (hydrostatic) then
                    te0 (i, j) = dp (i, j) * (te0 (i, j) + c_air * pt1 (i))
                else
#ifdef USE_COND
                    te0 (i, j) = dp (i, j) * (te0 (i, j) + cvm (i) * pt1 (i))
#else
                    te0 (i, j) = dp (i, j) * (te0 (i, j) + c_air * pt1 (i))
#endif
                endif
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            lhl (i) = lv00 + d0_vap * pt1 (i)
            cvm (i) = mc_air (i) + (qv (i, j) + q_liq (i) + q_sol (i)) * c_vap
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! compute cloud fraction
        ! -----------------------------------------------------------------------
        
        if (do_qa .and. last_step) then
            
            ! -----------------------------------------------------------------------
            ! combine water species
            ! -----------------------------------------------------------------------
            
            if (rad_snow) then
                if (rad_graupel) then
                    do i = is, ie
                        q_sol (i) = qi (i, j) + qs (i, j) + qg (i, j)
                    enddo
                else
                    do i = is, ie
                        q_sol (i) = qi (i, j) + qs (i, j)
                    enddo
                endif
            else
                do i = is, ie
                    q_sol (i) = qi (i, j)
                enddo
            endif
            if (rad_rain) then
                do i = is, ie
                    q_liq (i) = ql (i, j) + qr (i, j)
                enddo
            else
                do i = is, ie
                    q_liq (i) = ql (i, j)
                enddo
            endif
            do i = is, ie
                q_cond (i) = q_sol (i) + q_liq (i)
            enddo
            
            ! -----------------------------------------------------------------------
            ! use the "liquid - frozen water temperature" (tin) to compute saturated specific humidity
            ! -----------------------------------------------------------------------
            
            do i = is, ie
                
                tin = pt1 (i) - (lcp2 (i) * q_cond (i) + icp2 (i) * q_sol (i)) ! minimum temperature
                ! tin = pt1 (i) - ((lv00 + d0_vap * pt1 (i)) * q_cond (i) + &
                ! (li00 + dc_ice * pt1 (i)) * q_sol (i)) / (mc_air (i) + qpz (i) * c_vap)
                
                ! -----------------------------------------------------------------------
                ! determine saturated specific humidity
                ! -----------------------------------------------------------------------
                
                if (tin <= t_wfr) then
                    ! ice phase:
                    qstar (i) = iqs1 (tin, den (i))
                elseif (tin >= tice) then
                    ! liquid phase:
                    qstar (i) = wqs1 (tin, den (i))
                else
                    ! mixed phase:
                    qsi = iqs1 (tin, den (i))
                    qsw = wqs1 (tin, den (i))
                    if (q_cond (i) > 1.e-6) then
                        rqi = q_sol (i) / q_cond (i)
                    else
                        ! mostly liquid water clouds at initial cloud development stage
                        rqi = ((tice - tin) / (tice - t_wfr))
                    endif
                    qstar (i) = rqi * qsi + (1. - rqi) * qsw
                endif
                
                ! higher than 10 m is considered "land" and will have higher subgrid variability
                dw = dw_ocean + (dw_land - dw_ocean) * min (1., abs (hs (i, j)) / (10. * grav))
                ! "scale - aware" subgrid variability: 100 - km as the base
                hvar (i) = min (0.2, max (0.01, dw * sqrt (sqrt (area (i, j)) / 100.e3)))
                
                ! -----------------------------------------------------------------------
                ! partial cloudiness by pdf:
                ! assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
                ! binary cloud scheme; qa = 0.5 if qstar (i) == qpz
                ! -----------------------------------------------------------------------
                
                
                rh = qpz (i) / qstar (i)
                
                ! -----------------------------------------------------------------------
                ! icloud_f = 0: bug - fixed
                ! icloud_f = 1: old fvgfs gfdl) mp implementation
                ! icloud_f = 2: binary cloud scheme (0 / 1)
                ! -----------------------------------------------------------------------
                
                if (rh > 0.75 .and. qpz (i) > 1.e-6) then
                    dq = hvar (i) * qpz (i)
                    q_plus = qpz (i) + dq
                    q_minus = qpz (i) - dq
                    if (icloud_f == 2) then
                        if (qpz (i) > qstar (i)) then
                            qa (i, j) = 1.
                        elseif (qstar (i) < q_plus .and. q_cond (i) > 1.e-6) then
                            qa (i, j) = ((q_plus - qstar (i)) / dq) ** 2
                            qa (i, j) = min (1., qa (i, j))
                        else
                            qa (i, j) = 0.
                        endif
                    else
                        if (qstar (i) < q_minus) then
                            qa (i, j) = 1.
                        else
                            if (qstar (i) < q_plus) then
                                if (icloud_f == 0) then
                                    qa (i, j) = (q_plus - qstar (i)) / (dq + dq)
                                else
                                    qa (i, j) = (q_plus - qstar (i)) / (2. * dq * (1. - q_cond (i)))
                                endif
                            else
                                qa (i, j) = 0.
                            endif
                            ! impose minimum cloudiness if substantial q_cond (i) exist
                            if (q_cond (i) > 1.e-6) then
                                qa (i, j) = max (cld_min, qa (i, j))
                            endif
                            qa (i, j) = min (1., qa (i, j))
                        endif
                    endif
                else
                    qa (i, j) = 0.
                endif
                
            enddo
            
        endif
        
    enddo ! end j loop
    
end subroutine fv_sat_adj_run

end module fv_sat_adj
