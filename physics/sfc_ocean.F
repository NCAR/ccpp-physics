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

#if 0
!! \section arg_table_sfc_ocean_run Argument Table
!! \htmlinclude sfc_ocean_run.html
!!
#endif
      subroutine sfc_ocean_run                                          &
!...................................
!  ---  inputs:
     &     ( im, cp, rd, eps, epsm1, hvap, rvrdm1, ps, t1, q1,          &
     &       tskin, cm, ch, prsl1, prslki, wet, wind,                   &
     &       flag_iter,                                                 &
!  ---  outputs:
     &       qsurf, cmm, chh, gflux, evap, hflx, ep,                    &
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
!            prsl1, prslki, wet, wind, flag_iter,                       !
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

!  ---  inputs:
      integer, intent(in) :: im
      real (kind=kind_phys), intent(in) :: cp, rd, eps, epsm1, hvap,    &
     &      rvrdm1

      real (kind=kind_phys), dimension(im), intent(in) :: ps,           &
     &      t1, q1, tskin, cm, ch, prsl1, prslki, wind

      logical, dimension(im), intent(in) :: flag_iter, wet

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: qsurf,     &
     &       cmm, chh, gflux, evap, hflx, ep

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:

      real (kind=kind_phys) :: q0, qss, rch, rho, tem, cpinv,           &
     &      hvapi, elocp

      integer :: i

      logical :: flag(im)
!
!===> ...  begin here
      cpinv  = 1.0/cp
      hvapi  = 1.0/hvap
      elocp  = hvap/cp
!
!  -- ...  initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!  --- ...  flag for open water
      do i = 1, im
        flag(i) = (wet(i) .and. flag_iter(i))

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specified
!           ps is in pascals, wind is wind speed,
!           rho is density, qss is sat. hum. at surface

        if ( flag(i) ) then
          q0       = max( q1(i), 1.0e-8 )
          rho      = prsl1(i) / (rd*t1(i)*(1.0 + rvrdm1*q0))

          qss      = fpvs( tskin(i) )
          qss      = eps*qss / (ps(i) + epsm1*qss)

          evap(i)  = 0.0
          hflx(i)  = 0.0
          ep(i)    = 0.0
          gflux(i) = 0.0

!  --- ...    rcp  = rho cp ch v

          rch      = rho * cp * ch(i) * wind(i)
          cmm(i)   = cm(i) * wind(i)
          chh(i)   = rho * ch(i) * wind(i)

!  --- ...  sensible and latent heat flux over open water

          hflx(i)  = rch * (tskin(i) - t1(i) * prslki(i))

          evap(i)  = elocp*rch * (qss - q0)
          qsurf(i) = qss

          tem      = 1.0 / rho
          hflx(i)  = hflx(i) * tem * cpinv
          evap(i)  = evap(i) * tem * hvapi
        endif
      enddo
!
      return
!...................................
      end subroutine sfc_ocean_run
!-----------------------------------
      end module sfc_ocean
