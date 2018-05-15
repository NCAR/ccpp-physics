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

module qs_init
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
    
    implicit none
    
    private
    
    public  qs_init_run,  wqs2_vect, iqs2, iqs1, wqs1, wqs2
    
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
          
subroutine qs_init_init ()
end subroutine qs_init_init

subroutine qs_init_finalize ()
end subroutine qs_init_finalize


! =======================================================================
! initialization
! prepare saturation water vapor pressure tables
! =======================================================================
!>@brief The subroutine 'qs_init' initializes lookup tables for the saturation mixing ratio.
!! \section arg_table_gfdl_qs_init_run Argument Table
!!| local_name     | standard_name                                                 | long_name                                                                              | units   | rank | type      |   kind    | intent | optional |
!!|----------------|---------------------------------------------------------------|----------------------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!!| kmp            | top_layer_index_for_gfdl_mp                                   | top_layer_inder_for_gfdl_mp                                                            |         |    0 | real      |           |    in  |   F      |
!!
subroutine qs_init_run (kmp)
    
    implicit none
    
    integer, intent (in) :: kmp
    
    integer, parameter :: length = 2621
    
    integer :: i
    
    if (mp_initialized) return
    
    !if (is_master ()) write (*, *) 'top layer for gfdl_mp = ', kmp
    
    ! generate es table (dt = 0.1 deg c)
    
    allocate (table (length))
    allocate (table2 (length))
    allocate (tablew (length))
    allocate (des2 (length))
    allocate (desw (length))
    
    call qs_table (length)
    call qs_table2 (length)
    call qs_tablew (length)
    
    do i = 1, length - 1
        des2 (i) = max (0., table2 (i + 1) - table2 (i))
        desw (i) = max (0., tablew (i + 1) - tablew (i))
    enddo
    des2 (length) = des2 (length - 1)
    desw (length) = desw (length - 1)
    
    mp_initialized = .true.
    
end subroutine qs_init_run

! =======================================================================
! saturation water vapor pressure table i
! 3 - phase table
! =======================================================================

subroutine qs_table (n)
    
    implicit none
    
    integer, intent (in) :: n
    
    real (r8) :: delt = 0.1
    real (r8) :: tmin, tem, esh20
    real (r8) :: wice, wh2o, fac0, fac1, fac2
    real (r8) :: esupc (200)
    
    integer :: i
    
    tmin = tice - 160.
    
    ! -----------------------------------------------------------------------
    ! compute es over ice between - 160 deg c and 0 deg c.
    ! -----------------------------------------------------------------------
    
    do i = 1, 1600
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * li2
        fac2 = (d2ice * log (tem / tice) + fac1) / rvgas
        table (i) = e00 * exp (fac2)
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute es over water between - 20 deg c and 102 deg c.
    ! -----------------------------------------------------------------------
    
    do i = 1, 1221
        tem = 253.16 + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
        esh20 = e00 * exp (fac2)
        if (i <= 200) then
            esupc (i) = esh20
        else
            table (i + 1400) = esh20
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! derive blended es over ice and supercooled water between - 20 deg c and 0 deg c
    ! -----------------------------------------------------------------------
    
    do i = 1, 200
        tem = 253.16 + delt * real (i - 1)
        wice = 0.05 * (tice - tem)
        wh2o = 0.05 * (tem - 253.16)
        table (i + 1400) = wice * table (i + 1400) + wh2o * esupc (i)
    enddo
    
end subroutine qs_table

! =======================================================================
! saturation water vapor pressure table ii
! 1 - phase table
! =======================================================================

subroutine qs_tablew (n)
    
    implicit none
    
    integer, intent (in) :: n
    
    real (r8) :: delt = 0.1
    real (r8) :: tmin, tem, fac0, fac1, fac2
    
    integer :: i
    
    tmin = tice - 160.
    
    ! -----------------------------------------------------------------------
    ! compute es over water
    ! -----------------------------------------------------------------------
    
    do i = 1, n
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
        tablew (i) = e00 * exp (fac2)
    enddo
    
end subroutine qs_tablew

! =======================================================================
! saturation water vapor pressure table iii
! 2 - phase table
! =======================================================================

subroutine qs_table2 (n)
    
    implicit none
    
    integer, intent (in) :: n
    
    real (r8) :: delt = 0.1
    real (r8) :: tmin, tem0, tem1, fac0, fac1, fac2
    
    integer :: i, i0, i1
    
    tmin = tice - 160.
    
    do i = 1, n
        tem0 = tmin + delt * real (i - 1)
        fac0 = (tem0 - tice) / (tem0 * tice)
        if (i <= 1600) then
            ! -----------------------------------------------------------------------
            ! compute es over ice between - 160 deg c and 0 deg c.
            ! -----------------------------------------------------------------------
            fac1 = fac0 * li2
            fac2 = (d2ice * log (tem0 / tice) + fac1) / rvgas
        else
            ! -----------------------------------------------------------------------
            ! compute es over water between 0 deg c and 102 deg c.
            ! -----------------------------------------------------------------------
            fac1 = fac0 * lv0
            fac2 = (dc_vap * log (tem0 / tice) + fac1) / rvgas
        endif
        table2 (i) = e00 * exp (fac2)
    enddo
    
    ! -----------------------------------------------------------------------
    ! smoother around 0 deg c
    ! -----------------------------------------------------------------------
    
    i0 = 1600
    i1 = 1601
    tem0 = 0.25 * (table2 (i0 - 1) + 2. * table (i0) + table2 (i0 + 1))
    tem1 = 0.25 * (table2 (i1 - 1) + 2. * table (i1) + table2 (i1 + 1))
    table2 (i0) = tem0
    table2 (i1) = tem1
    
end subroutine qs_table2

! =======================================================================
!>@brief the function 'wqs1' computes the 
!! saturated specific humidity for table ii
! =======================================================================
real function wqs1 (ta, den)

    implicit none

    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es, ap1, tmin

    integer :: it

    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs1 = es / (rvgas * ta * den)

end function wqs1

! =======================================================================
!>@brief the function 'wqs1' computes the  saturated specific humidity 
!! for table iii
! =======================================================================
real function iqs1 (ta, den)

    implicit none

    ! water - ice phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es, ap1, tmin

    integer :: it

    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs1 = es / (rvgas * ta * den)

end function iqs1

! =======================================================================
!>@brief The function 'wqs2'computes the gradient of saturated specific 
!! humidity for table ii
! =======================================================================
real function wqs2 (ta, den, dqdt)

    implicit none

    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real, intent (out) :: dqdt

    real :: es, ap1, tmin

    integer :: it

    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    ! finite diff, del_t = 0.1:
    dqdt = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta * den)

end function wqs2


! =======================================================================
!>@brief The function wqs2_vect computes the gradient of saturated 
!! specific humidity for table ii.
!! It is the same as "wqs2", but written as vector function.
! =======================================================================
subroutine wqs2_vect (is, ie, ta, den, wqsat, dqdt)

    implicit none

    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    integer, intent (in) :: is, ie

    real, intent (in), dimension (is:ie) :: ta, den

    real, intent (out), dimension (is:ie) :: wqsat, dqdt

    real :: es, ap1, tmin

    integer :: i, it

    tmin = tice - 160.

    do i = is, ie
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es = tablew (it) + (ap1 - it) * desw (it)
        wqsat (i) = es / (rvgas * ta (i) * den (i))
        it = ap1 - 0.5
        ! finite diff, del_t = 0.1:
        dqdt (i) = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta (i) * den (i))
    enddo

end subroutine wqs2_vect

! =======================================================================
!>@brief The function 'iqs2' computes the gradient of saturated specific 
!! humidity for table iii.
! =======================================================================
real function iqs2 (ta, den, dqdt)

    implicit none

    ! water - ice phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real, intent (out) :: dqdt

    real :: es, ap1, tmin

    integer :: it

    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    ! finite diff, del_t = 0.1:
    dqdt = 10. * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) / (rvgas * ta * den)

end function iqs2


end module qs_init
