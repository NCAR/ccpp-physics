!>\file sfc_ocean.F
!! This file contains an alternative GFS near-surface sea temperature
!! scheme when the model is initialized from GRIB2 data.

!> This module contains the CCPP-compliant GFS near-surface sea temperature
!! scheme when the model is initialized from GRIB2 data.
      module sfc_ocean
      implicit none
      private
      public :: sfc_ocean_run

      contains


!>\defgroup gfs_ocean_main GFS Simple Ocean Module
!! This subroutine calculates thermodynamical properties over
!! open water.
!>@{
!! \section arg_table_sfc_ocean_run Argument Table
!! \htmlinclude sfc_ocean_run.html
!!
!!>\section gen_sfc_ocean GFS Simple Ocean scheme General Algorithm
      subroutine sfc_ocean_run                                          &
!...................................
!  ---  inputs:
     &     ( im, hvap, cp, rd, eps, epsm1, rvrdm1, ps, u1, v1, t1, q1,  &
     &       tskin, cm, ch, lseaspray, fm, fm10,                        &
     &       prsl1, prslki, wet, use_lake_model, wind,                  &,  ! --- inputs
     &       flag_iter, use_med_flux, dqsfc_med, dtsfc_med,             &
     &       qsurf, cmm, chh, gflux, evap, hflx, ep,                    &   ! --- outputs
     &       errmsg, errflg                                             &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_ocean                                                     !
!       inputs:                                                         !
!        ( im, ps, u1, v1, t1, q1, tskin, cm, ch, lseaspray, fm, fm10,  !
!            prsl1, prslki, wet, use_lake_model, wind, flag_iter,       !
!            use_med_flux,                                              !
!       outputs:                                                        !
!            qsurf, cmm, chh, gflux, evap, hflx, ep )                   !
!                                                                       !
!                                                                       !
!  subprograms/functions called: fpvs                                   !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!         2005  -- created from the original progtm to account for      !
!                  ocean only                                           !
!    oct  2006  -- h. wei      added cmm and chh to the output          !
!    apr  2009  -- y.-t. hou   modified to match the modified gbphys.f  !
!                  reformatted the code and added program documentation !
!    sep  2009  -- s. moorthi removed rcl and made pa as pressure unit  !
!                  and furthur reformatted the code                     !
!    dec  2021  -- u. turuncoglu added support for receiving fluxes     !
!                  from mediator                                        !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horizontal dimension                     1    !
!     ps       - real, surface pressure                            im   !
!     u1, v1   - real, u/v component of surface layer wind (m/s)   im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     lseaspray- logical, .t. for parameterization for sea spray   1    !
!     fm       - real, a stability profile function for momentum   im   !
!     fm10     - real, a stability profile function for momentum   im   !
!                       at 10m                                          !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     wet      - logical, =T if any ocean/lak, =F otherwise        im   !
!     wind     - real, wind speed (m/s)                            im   !
!     flag_iter- logical,                                          im   !
!     use_med_flux - logical, =T to use fluxes coming from med     1    !
!     dqsfc_med- real, latent heat flux                            im   !
!     dtsfc_med- real, sensible heat flux                          im   !
!                                                                       !
!  outputs:                                                             !
!     qsurf    - real, specific humidity at sfc                    im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     gflux    - real, ground heat flux (zero for ocean)           im   !
!     evap     - real, evaporation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!     ep       - real, potential evaporation                       im   !
!                                                                       !
! ===================================================================== !
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
!
      implicit none

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: one  = 1.0_kind_phys,         &
     &           zero = 0.0_kind_phys, qmin = 1.0e-8_kind_phys
!  ---  inputs:
      integer, intent(in) :: im
      real (kind=kind_phys), intent(in) :: hvap, cp, rd,                &
     &                 eps, epsm1, rvrdm1

      real (kind=kind_phys), dimension(:), intent(in) :: ps, u1, v1,    &
     &      t1, q1, tskin, cm, ch, fm, fm10, prsl1, prslki, wind

! For sea spray effect
      logical, intent(in) :: lseaspray
!
      logical, dimension(:), intent(in) :: flag_iter, wet
      integer, dimension(:), intent(in) :: use_lake_model
!
      logical, intent(in) :: use_med_flux

! To receive fluxes from mediator
      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       dqsfc_med, dtsfc_med

!  ---  outputs:
      real (kind=kind_phys), dimension(:), intent(inout) :: qsurf,      &
     &       cmm, chh, gflux, evap, hflx, ep

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:

      real (kind=kind_phys) :: qss, rch, tem,
     &                         elocp, cpinv, hvapi
      real (kind=kind_phys), dimension(im) :: rho, q0

      integer :: i

      logical :: flag(im)
!
!  parameters for sea spray effect
!
      real (kind=kind_phys) :: f10m, u10m, v10m, ws10, ru10, qss1,
     &                         bb1, hflxs, evaps, ptem
!
!     real (kind=kind_phys), parameter :: alps=0.5, bets=0.5, gams=0.1,
!     real (kind=kind_phys), parameter :: alps=0.5, bets=0.5, gams=0.0,
!     real (kind=kind_phys), parameter :: alps=1.0, bets=1.0, gams=0.2,
      real (kind=kind_phys), parameter :: alps=0.75,bets=0.75,gams=0.15,
     &                       ws10cr=30., conlf=7.2e-9, consf=6.4e-8
!
!======================================================================================================
!===> ...  begin here
!
!  -- ...  initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      cpinv = one/cp
      hvapi = one/hvap
      elocp = hvap/cp
!
!> - Flag for open water
      do i = 1, im
        flag(i) = (wet(i) .and. flag_iter(i) .and. use_lake_model(i)/=1)
!> - Initialize variables. all units are supposedly m.k.s. unless specified
!           ps is in pascals, wind is wind speed,
!           rho is density, qss is sat. hum. at surface

        if ( flag(i) ) then
          if (use_med_flux) then
            q0(i)    = max( q1(i), qmin )
            rho(i)   = prsl1(i) / (rd*t1(i)*(one + rvrdm1*q0(i)))

            tem      = ch(i) * wind(i)
            cmm(i)   = cm(i) * wind(i)
            chh(i)   = rho(i) * tem

            hflx(i)  = dtsfc_med(i)
            evap(i)  = dqsfc_med(i)

            qsurf(i) = q1(i) + dqsfc_med(i) / (hvap*chh(i))
            gflux(i) = zero
          else
            q0(i)    = max( q1(i), qmin )
            rho(i)   = prsl1(i) / (rd*t1(i)*(one + rvrdm1*q0(i)))

            qss      = fpvs( tskin(i) )
            qss      = eps*qss / (ps(i) + epsm1*qss)

!  --- ...    rcp  = rho cp ch v

            rch      = rho(i) * cp * ch(i) * wind(i)
            tem      = ch(i) * wind(i)
            cmm(i)   = cm(i) * wind(i)
            chh(i)   = rho(i) * tem

!> - Calcualte sensible and latent heat flux over open water

            hflx(i)  = rch * (tskin(i) - t1(i) * prslki(i))

            evap(i)  = elocp * rch * (qss - q0(i))

            qsurf(i) = qss
            gflux(i) = zero
          endif
        endif
      enddo
!
!> - Include sea spray effects
!
      do i=1,im
        if(lseaspray .and. flag(i)) then
          f10m = fm10(i) / fm(i)
          u10m = f10m * u1(i)
          v10m = f10m * v1(i)
          ws10 = sqrt(u10m*u10m + v10m*v10m)
          ws10 = max(ws10,1.)
          ws10 = min(ws10,ws10cr)
          tem = .015 * ws10 * ws10
          ru10 = 1. - .087 * log(10./tem)
          qss1 = fpvs(t1(i))
          qss1 = eps * qss1 / (prsl1(i) + epsm1 * qss1)
          tem = rd * cp * t1(i) * t1(i)
          tem = 1. + eps * hvap * hvap * qss1 / tem
          bb1 = 1. / tem
          evaps = conlf * (ws10**5.4) * ru10 * bb1
          evaps = evaps * rho(i) * hvap * (qss1 - q0(i))
          evap(i) = evap(i) + alps * evaps
          hflxs = consf * (ws10**3.4) * ru10
          hflxs = hflxs * rho(i) * cp * (tskin(i) - t1(i))
          ptem = alps - gams
          hflx(i) = hflx(i) + bets * hflxs - ptem * evaps
        endif
      enddo
!
      do i = 1, im
        if ( flag(i) ) then
          tem      = 1.0 / rho(i)
          hflx(i)  = hflx(i) * tem * cpinv
          evap(i)  = evap(i) * tem * hvapi
          ep(i)    = evap(i)
        endif
      enddo
!
      return
!...................................
      end subroutine sfc_ocean_run
!-----------------------------------
!>@}
      end module sfc_ocean
