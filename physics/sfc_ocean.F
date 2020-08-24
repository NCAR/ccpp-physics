!>\file sfc_ocean.F
!! This file contains an alternative GFS near-surface sea temperature
!! scheme when the model is initialized from GRIB2 data.

!> This module contains the CCPP-compliant GFS near-surface sea temperature
!! scheme when the model is initialized from GRIB2 data.
      module sfc_ocean
      implicit none
      private
      public :: sfc_ocean_init, sfc_ocean_run, sfc_ocean_finalize

      contains

!! \section arg_table_sfc_ocean_init  Argument Table
!!
      subroutine sfc_ocean_init()
      end subroutine sfc_ocean_init

!! \section arg_table_sfc_ocean_finalize  Argument Table
!!
      subroutine sfc_ocean_finalize()
      end subroutine sfc_ocean_finalize

!>\defgroup gfs_ocean_main GFS Simple Ocean Scheme Module
!! This subroutine calculates thermodynamical properties over
!! open water.
!! \section arg_table_sfc_ocean_run Argument Table
!! \htmlinclude sfc_ocean_run.html
!!
      subroutine sfc_ocean_run                                          &
!...................................
!  ---  inputs:
     &     ( im, rd, eps, epsm1, rvrdm1, ps, t1, q1,                    &
     &       tskin, cm, ch, prsl1, prslki, wet, lake, wind,             &,  ! --- inputs
     &       flag_iter,                                                 &
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
!          ( im, ps, t1, q1, tskin, cm, ch,                             !
!!         ( im, ps, u1, v1, t1, q1, tskin, cm, ch,                     !
!            prsl1, prslki, wet, lake, wind, flag_iter,                 !
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
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horizontal dimension                     1    !
!     ps       - real, surface pressure                            im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     wet      - logical, =T if any ocean/lak, =F otherwise        im   !
!     wind     - real, wind speed (m/s)                            im   !
!     flag_iter- logical,                                          im   !
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
      real (kind=kind_phys), parameter :: one  = 1.0_kind_phys, zero = 0.0_kind_phys    &
     &,                                   qmin = 1.0e-8_kind_phys
!  ---  inputs:
      integer, intent(in) :: im
      real (kind=kind_phys), intent(in) :: rd, eps, epsm1, rvrdm1

      real (kind=kind_phys), dimension(im), intent(in) :: ps,           &
     &      t1, q1, tskin, cm, ch, prsl1, prslki, wind

      logical, dimension(im), intent(in) :: flag_iter, wet, lake

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: qsurf,     &
     &       cmm, chh, gflux, evap, hflx, ep

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:

      real (kind=kind_phys) :: q0, qss, rch, rho, tem

      integer :: i

!===> ...  begin here
!
!  -- ...  initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!  --- ...  flag for open water
      do i = 1, im

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specified
!           ps is in pascals, wind is wind speed,
!           rho is density, qss is sat. hum. at surface

        if (wet(i) .and. flag_iter(i) .and. .not. lake(i)) then
          q0       = max( q1(i), qmin )
          rho      = prsl1(i) / (rd*t1(i)*(one + rvrdm1*q0))

          qss      = fpvs( tskin(i) )
          qss      = eps*qss / (ps(i) + epsm1*qss)

!  --- ...    rcp  = rho cp ch v

          tem      = ch(i) * wind(i)
          cmm(i)   = cm(i) * wind(i)
          chh(i)   = rho * tem

!  --- ...  sensible and latent heat flux over open water

          hflx(i)  = tem * (tskin(i) - t1(i) * prslki(i))

          evap(i)  = tem * (qss - q0)

          ep(i)    = evap(i)
          qsurf(i) = qss
          gflux(i) = zero
        endif
      enddo
!
      return
!...................................
      end subroutine sfc_ocean_run
!-----------------------------------
      end module sfc_ocean
