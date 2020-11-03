!>\file GFS_rrtmg_post.f90
!! This file contains
       module GFS_rrtmg_post
       contains

!>\defgroup GFS_rrtmg_post GFS RRTMG Scheme Post
!! @{
       subroutine GFS_rrtmg_post_init ()
       end subroutine GFS_rrtmg_post_init

!> \section arg_table_GFS_rrtmg_post_run Argument Table
!! \htmlinclude GFS_rrtmg_post_run.html
!!
      subroutine GFS_rrtmg_post_run (im, km, kmp1, lm, ltp, kt, kb, kd, nspc1, &
              nfxr, nday, lsswr, lslwr, lssav, fhlwr, fhswr, raddt, coszen,    &
              coszdg, prsi, tgrs, aerodp, cldsa, mtopa, mbota, clouds1,        &
              cldtaulw, cldtausw, sfcflw, sfcfsw, topflw, topfsw, scmpsw,      &
              fluxr, errmsg, errflg)

      use machine,                             only: kind_phys
      use module_radsw_parameters,             only: topfsw_type, sfcfsw_type, &
                                                     cmpfsw_type
      use module_radlw_parameters,             only: topflw_type, sfcflw_type

      implicit none

      ! Interface variables
      integer,              intent(in) :: im, km, kmp1, lm, ltp, kt, kb, kd,   &
                                          nspc1, nfxr, nday
      logical,              intent(in) :: lsswr, lslwr, lssav
      real(kind=kind_phys), intent(in) :: raddt, fhlwr, fhswr
            
      real(kind=kind_phys), dimension(im),        intent(in) :: coszen, coszdg
      
      real(kind=kind_phys), dimension(im,kmp1),   intent(in) :: prsi
      real(kind=kind_phys), dimension(im,km),     intent(in) :: tgrs
      
      real(kind=kind_phys), dimension(im,NSPC1),  intent(in) :: aerodp
      real(kind=kind_phys), dimension(im,5),      intent(in) :: cldsa
      integer,              dimension(im,3),      intent(in) :: mbota, mtopa
      real(kind=kind_phys), dimension(im,lm+LTP), intent(in) :: clouds1
      real(kind=kind_phys), dimension(im,lm+LTP), intent(in) :: cldtausw
      real(kind=kind_phys), dimension(im,lm+LTP), intent(in) :: cldtaulw
      
      type(sfcflw_type), dimension(im), intent(in) :: sfcflw
      type(sfcfsw_type), dimension(im), intent(in) :: sfcfsw
      type(cmpfsw_type), dimension(im), intent(in) :: scmpsw
      type(topflw_type), dimension(im), intent(in) :: topflw
      type(topfsw_type), dimension(im), intent(in) :: topfsw
      
      real(kind=kind_phys), dimension(im,nfxr), intent(inout) :: fluxr

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! Local variables
      integer :: i, j, k, k1, itop, ibtc
      real(kind=kind_phys) :: tem0d, tem1, tem2

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (.not. (lsswr .or. lslwr)) return

!>  - For time averaged output quantities (including total-sky and
!!    clear-sky SW and LW fluxes at TOA and surface; conventional
!!    3-domain cloud amount, cloud top and base pressure, and cloud top
!!    temperature; aerosols AOD, etc.), store computed results in
!!    corresponding slots of array fluxr with appropriate time weights.

!  --- ...  collect the fluxr data for wrtsfc

      if (lssav) then
        if (lsswr) then
          do i=1,im
!            fluxr(i,34) = fluxr(i,34) + fhswr*aerodp(i,1)  ! total aod at 550nm
!            fluxr(i,35) = fluxr(i,35) + fhswr*aerodp(i,2)  ! DU aod at 550nm
!            fluxr(i,36) = fluxr(i,36) + fhswr*aerodp(i,3)  ! BC aod at 550nm
!            fluxr(i,37) = fluxr(i,37) + fhswr*aerodp(i,4)  ! OC aod at 550nm
!            fluxr(i,38) = fluxr(i,38) + fhswr*aerodp(i,5)  ! SU aod at 550nm
!            fluxr(i,39) = fluxr(i,39) + fhswr*aerodp(i,6)  ! SS aod at 550nm
             fluxr(i,34) = aerodp(i,1)  ! total aod at 550nm
             fluxr(i,35) = aerodp(i,2)  ! DU aod at 550nm
             fluxr(i,36) = aerodp(i,3)  ! BC aod at 550nm
             fluxr(i,37) = aerodp(i,4)  ! OC aod at 550nm
             fluxr(i,38) = aerodp(i,5)  ! SU aod at 550nm
             fluxr(i,39) = aerodp(i,6)  ! SS aod at 550nm
          enddo
        endif

!  ---  save lw toa and sfc fluxes
        if (lslwr) then
          do i=1,im
!  ---  lw total-sky fluxes
            fluxr(i,1 ) = fluxr(i,1 ) + fhlwr * topflw(i)%upfxc   ! total sky top lw up
            fluxr(i,19) = fluxr(i,19) + fhlwr * sfcflw(i)%dnfxc   ! total sky sfc lw dn
            fluxr(i,20) = fluxr(i,20) + fhlwr * sfcflw(i)%upfxc   ! total sky sfc lw up
!  ---  lw clear-sky fluxes
            fluxr(i,28) = fluxr(i,28) + fhlwr * topflw(i)%upfx0   ! clear sky top lw up
            fluxr(i,30) = fluxr(i,30) + fhlwr * sfcflw(i)%dnfx0   ! clear sky sfc lw dn
            fluxr(i,33) = fluxr(i,33) + fhlwr * sfcflw(i)%upfx0   ! clear sky sfc lw up
          enddo
        endif

!  ---  save sw toa and sfc fluxes with proper diurnal sw wgt. coszen=mean cosz over daylight
!       part of sw calling interval, while coszdg= mean cosz over entire interval
        if (lsswr) then
          do i = 1, IM
            if (coszen(i) > 0.) then
!  ---                                  sw total-sky fluxes
!                                       -------------------
              tem0d = fhswr * coszdg(i) / coszen(i)
              fluxr(i,2 ) = fluxr(i,2)  +  topfsw(i)%upfxc * tem0d  ! total sky top sw up
              fluxr(i,3 ) = fluxr(i,3)  +  sfcfsw(i)%upfxc * tem0d  ! total sky sfc sw up
              fluxr(i,4 ) = fluxr(i,4)  +  sfcfsw(i)%dnfxc * tem0d  ! total sky sfc sw dn
!  ---                                  sw uv-b fluxes
!                                       --------------
              fluxr(i,21) = fluxr(i,21) + scmpsw(i)%uvbfc * tem0d          ! total sky uv-b sw dn
              fluxr(i,22) = fluxr(i,22) + scmpsw(i)%uvbf0 * tem0d          ! clear sky uv-b sw dn
!  ---                                  sw toa incoming fluxes
!                                       ----------------------
              fluxr(i,23) = fluxr(i,23) + topfsw(i)%dnfxc * tem0d     ! top sw dn
!  ---                                  sw sfc flux components
!                                       ----------------------
              fluxr(i,24) = fluxr(i,24) + scmpsw(i)%visbm * tem0d          ! uv/vis beam sw dn
              fluxr(i,25) = fluxr(i,25) + scmpsw(i)%visdf * tem0d          ! uv/vis diff sw dn
              fluxr(i,26) = fluxr(i,26) + scmpsw(i)%nirbm * tem0d          ! nir beam sw dn
              fluxr(i,27) = fluxr(i,27) + scmpsw(i)%nirdf * tem0d          ! nir diff sw dn
!  ---                                  sw clear-sky fluxes
!                                       -------------------
              fluxr(i,29) = fluxr(i,29) + topfsw(i)%upfx0 * tem0d  ! clear sky top sw up
              fluxr(i,31) = fluxr(i,31) + sfcfsw(i)%upfx0 * tem0d  ! clear sky sfc sw up
              fluxr(i,32) = fluxr(i,32) + sfcfsw(i)%dnfx0 * tem0d  ! clear sky sfc sw dn
            endif
          enddo
        endif

!  ---  save total and boundary layer clouds

        if (lsswr .or. lslwr) then
          do i=1,im
            fluxr(i,17) = fluxr(i,17) + raddt * cldsa(i,4)
            fluxr(i,18) = fluxr(i,18) + raddt * cldsa(i,5)
          enddo

!  ---  save cld frac,toplyr,botlyr and top temp, note that the order
!       of h,m,l cloud is reversed for the fluxr output.
!  ---  save interface pressure (pa) of top/bot

          do j = 1, 3
            do i = 1, IM
              tem0d = raddt * cldsa(i,j)
              itop  = mtopa(i,j) - kd
              ibtc  = mbota(i,j) - kd
              fluxr(i, 8-j) = fluxr(i, 8-j) + tem0d
              fluxr(i,11-j) = fluxr(i,11-j) + tem0d * prsi(i,itop+kt)
              fluxr(i,14-j) = fluxr(i,14-j) + tem0d * prsi(i,ibtc+kb)
              fluxr(i,17-j) = fluxr(i,17-j) + tem0d * tgrs(i,itop)
            enddo
          enddo

!       Anning adds optical depth and emissivity output
          if (lsswr .and. (nday > 0)) then
            do j = 1, 3
              do i = 1, IM
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j) - kd
                ibtc  = mbota(i,j) - kd
                tem1 = 0.
                do k=ibtc,itop
                  tem1 = tem1 + cldtausw(i,k)      ! approx .55 um channel
                enddo
                fluxr(i,43-j) = fluxr(i,43-j) + tem0d * tem1
              enddo
            enddo
          endif

          if (lslwr) then
            do j = 1, 3
              do i = 1, IM
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j) - kd
                ibtc  = mbota(i,j) - kd
                tem2 = 0.
                do k=ibtc,itop
                  tem2 = tem2 + cldtaulw(i,k)      ! approx 10. um channel
                enddo
                fluxr(i,46-j) = fluxr(i,46-j) + tem0d * (1.0-exp(-tem2))
              enddo
            enddo
          endif

        endif

      endif                                ! end_if_lssav
!
      end subroutine GFS_rrtmg_post_run

      subroutine GFS_rrtmg_post_finalize ()
      end subroutine GFS_rrtmg_post_finalize

!! @}
      end module GFS_rrtmg_post
