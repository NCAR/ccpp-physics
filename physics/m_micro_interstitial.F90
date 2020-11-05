!> \file m_micro_interstitial.F90
!! This file contains subroutines that prepare data for and from the Morrison-Gettelman microphysics scheme
!! as part of the GFS physics suite.
      module m_micro_pre

      implicit none

      contains

      subroutine m_micro_pre_init()
      end subroutine m_micro_pre_init

! \brief Brief description of the subroutine
!!
!! \section arg_table_m_micro_pre_run Argument Table
!! \htmlinclude m_micro_pre_run.html
!!
      subroutine m_micro_pre_run (im, levs, do_shoc, skip_macro, fprcp, mg3_as_mg2, gq0_ice, gq0_water, gq0_rain,  &
        gq0_snow, gq0_graupel, gq0_rain_nc, gq0_snow_nc, gq0_graupel_nc, cld_shoc, cnvc, cnvw, tcr, tcrf, gt0,     &
        qrn, qsnw, qgl, ncpr, ncps, ncgl, cld_frc_MG, clw_water, clw_ice, clcn, errmsg, errflg )

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs, fprcp
      logical, intent(in) :: do_shoc, mg3_as_mg2
      logical, intent(inout) :: skip_macro
      real(kind=kind_phys), intent(in) :: tcr, tcrf

      real(kind=kind_phys), intent(in) ::                               &
          gq0_ice(:,:), gq0_water(:,:), gq0_rain(:,:), gq0_snow(:,:),   &
          gq0_graupel(:,:), gq0_rain_nc(:,:), gq0_snow_nc(:,:),         &
          gq0_graupel_nc(:,:), cld_shoc(:,:), cnvc(:,:), cnvw(:,:),     &
          gt0(:,:)

      real(kind=kind_phys), intent(inout) ::                              &
          qrn(:,:), qsnw(:,:), qgl(:,:), ncpr(:,:), ncps(:,:), ncgl(:,:), &
          cld_frc_MG(:,:)

      real(kind=kind_phys), intent(out) :: clw_ice(:,:), clw_water(:,:)

      real(kind=kind_phys), intent(in) :: clcn(:,:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k
      real(kind=kind_phys) :: tem

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !       Acheng used clw here for other code to run smoothly and minimum change
      !       to make the code work. However, the nc and clw should be treated
      !       in other procceses too.  August 28/2015; Hope that can be done next
      !       year. I believe this will make the physical interaction more reasonable
      !       Anning 12/5/2015 changed ntcw hold liquid only
      skip_macro = do_shoc
      if (do_shoc) then
        if (fprcp == 0) then
          do k=1,levs
            do i=1,im
              clw_ice(i,k)    = gq0_ice(i,k)
              clw_water(i,k)  = gq0_water(i,k)
              cld_frc_MG(i,k) = cld_shoc(i,k)
            enddo
          enddo
        else if ((abs(fprcp) == 1) .or. mg3_as_mg2) then
          do k=1,levs
            do i=1,im
              clw_ice(i,k)    = gq0_ice(i,k)
              clw_water(i,k)  = gq0_water(i,k)
              qrn(i,k)        = gq0_rain(i,k)
              qsnw(i,k)       = gq0_snow(i,k)
              ncpr(i,k)       = gq0_rain_nc(i,k)
              ncps(i,k)       = gq0_snow_nc(i,k)
              cld_frc_MG(i,k) = cld_shoc(i,k)
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              clw_ice(i,k)    = gq0_ice(i,k)
              clw_water(i,k)  = gq0_water(i,k)
              qrn(i,k)        = gq0_rain(i,k)
              qsnw(i,k)       = gq0_snow(i,k)
              qgl(i,k)        = gq0_graupel(i,k)
              ncpr(i,k)       = gq0_rain_nc(i,k)
              ncps(i,k)       = gq0_snow_nc(i,k)
              ncgl(i,k)       = gq0_graupel_nc(i,k)
              cld_frc_MG(i,k) = cld_shoc(i,k)
            enddo
          enddo
        end if
      else
        if (fprcp == 0 ) then
          do k=1,levs
            do i=1,im
              clw_ice(i,k)   = gq0_ice(i,k)
              clw_water(i,k) = gq0_water(i,k)
            enddo
          enddo
        elseif (abs(fprcp) == 1 .or. mg3_as_mg2) then
          do k=1,levs
            do i=1,im
              clw_ice(i,k)   = gq0_ice(i,k)
              clw_water(i,k) = gq0_water(i,k)
              qrn(i,k)       = gq0_rain(i,k)
              qsnw(i,k)      = gq0_snow(i,k)
              ncpr(i,k)      = gq0_rain_nc(i,k)
              ncps(i,k)      = gq0_snow_nc(i,k)
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              clw_ice(i,k)   = gq0_ice(i,k)
              clw_water(i,k) = gq0_water(i,k)
              qrn(i,k)       = gq0_rain(i,k)
              qsnw(i,k)      = gq0_snow(i,k)
              qgl(i,k)       = gq0_graupel(i,k)
              ncpr(i,k)      = gq0_rain_nc(i,k)
              ncps(i,k)      = gq0_snow_nc(i,k)
              ncgl(i,k)      = gq0_graupel_nc(i,k)
            enddo
          enddo
        endif
      end if

      ! add convective cloud fraction
      do k = 1,levs
        do i = 1,im
          cld_frc_MG(i,k) = min(1.0, cld_frc_MG(i,k) + clcn(i,k))
        enddo
      enddo

      end subroutine m_micro_pre_run

      subroutine m_micro_pre_finalize ()
      end subroutine m_micro_pre_finalize

      end module m_micro_pre

!> This module contains the CCPP-compliant MG microphysics
!! post intersititial codes.
      module m_micro_post

      implicit none

      contains

      subroutine m_micro_post_init()
      end subroutine m_micro_post_init

! \brief Brief description of the subroutine
!!
!! \section arg_table_m_micro_post_run Argument Table
!! \htmlinclude m_micro_post_run.html
!!
      subroutine m_micro_post_run(                                          &
        im, levs, fprcp, mg3_as_mg2, ncpr, ncps, ncgl, qrn, qsnw, qgl,      &
        gq0_ice, gq0_rain, gq0_snow, gq0_graupel, gq0_rain_nc, gq0_snow_nc, &
        gq0_graupel_nc, ice, snow, graupel, dtp, errmsg, errflg)

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs, fprcp
      logical, intent(in) :: mg3_as_mg2

      real(kind=kind_phys), intent(in   ) :: ncpr(1:im,1:levs)
      real(kind=kind_phys), intent(in   ) :: ncps(1:im,1:levs)
      real(kind=kind_phys), intent(in   ) :: ncgl(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: qrn(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: qsnw(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: qgl(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: gq0_ice(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: gq0_rain(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: gq0_snow(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: gq0_graupel(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: gq0_rain_nc(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: gq0_snow_nc(1:im,1:levs)
      real(kind=kind_phys), intent(inout) :: gq0_graupel_nc(1:im,1:levs)
      real(kind=kind_phys), intent(  out) :: ice(1:im)
      real(kind=kind_phys), intent(  out) :: snow(1:im)
      real(kind=kind_phys), intent(  out) :: graupel(1:im)
      real(kind=kind_phys), intent(in   ) :: dtp

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      real(kind=kind_phys), parameter :: qsmall   = 1.0d-20
      real(kind=kind_phys), parameter :: con_p001 = 0.001d0
      real(kind=kind_phys), parameter :: con_day  = 86400.0d0
      integer :: i, k
      real(kind=kind_phys) :: tem

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!     do k=1,levs
!     write(1000+me,*)' maxwatnca=',maxval(Stateout%gq0(1:im,k,ntlnc)),' k=',k,' kdt=',kdt
!     enddo
!     write(1000+me,*)' at latitude = ',lat
!     tx1 = 1000.0
!     call moist_bud(im,ix,ix,levs,me,kdt,con_g,tx1,del,rain1
!    &,                    txa, clw(1,1,2), clw(1,1,1)
!    &,           gq0(1,1,1),gq0(1,1,ntcw),gq0(1,1,ntcw+1),' m_micro  ')

!       if (lprnt) write(0,*) ' rain1=',rain1(ipr)*86400.0, &
!    &' rainc=',diag%rainc(ipr)*86400.0
!    &,' cn_prc=',cn_prc(ipr),' cn_snr=',cn_snr(ipr)
!       if(lprnt) write(0,*) ' aftgt0=',Stateout%gt0(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*) ' aftlsgq0=',stateout%gq0(ipr,:,1),' kdt=',kdt
!       if (lprnt) write(0,*)' clw1aft=',stateout%gq0(ipr,:,ntiw),' kdt=',kdt
!       if (ntgl > 0 .and. lprnt)  &
!                  write(0,*)' cgw1aft=',stateout%gq0(ipr,:,ntgl),' kdt=',kdt
!       if (lprnt) write(0,*)' cloudsm=',tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clw2aft=',stateout%gq0(ipr,:,ntcw),' kdt=',kdt
!       if (lprnt) write(0,*)' qrna=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwa=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglba',qgl(ipr,:),' kdt=',kdt

      tem = dtp * con_p001 / con_day
      if (abs(fprcp) == 1 .or. mg3_as_mg2) then
        do k=1,levs
          do i=1,im
            if (abs(qrn(i,k))  < qsmall) qrn(i,k)  = 0.0
            if (abs(qsnw(i,k)) < qsmall) qsnw(i,k) = 0.0
            gq0_rain(i,k)    = qrn(i,k)
            gq0_snow(i,k)    = qsnw(i,k)
            gq0_rain_nc(i,k) = ncpr(i,k)
            gq0_snow_nc(i,k) = ncps(i,k)
          enddo
        enddo
        do i=1,im
          ice(i)  = tem * gq0_ice(i,1)
          snow(i) = tem * qsnw(i,1)
        enddo
      elseif (fprcp > 1) then
        do k=1,levs
          do i=1,im
            if (abs(qrn(i,k))  < qsmall) qrn(i,k)  = 0.0
            if (abs(qsnw(i,k)) < qsmall) qsnw(i,k) = 0.0
            if (abs(qgl(i,k))  < qsmall) qgl(i,k)  = 0.0
            gq0_rain(i,k)       = qrn(i,k)
            gq0_snow(i,k)       = qsnw(i,k)
            gq0_graupel(i,k)    = qgl(i,k)
            gq0_rain_nc(i,k)    = ncpr(i,k)
            gq0_snow_nc(i,k)    = ncps(i,k)
            gq0_graupel_nc(i,k) = ncgl(i,k)
          enddo
        enddo
        do i=1,im
          ice(i)     = tem * gq0_ice(i,1)
          snow(i)    = tem * qsnw(i,1)
          graupel(i) = tem * qgl(i,1)
        enddo

      endif

!       if (lprnt) write(0,*)' cloudsm=',tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clw2aft=',stateout%gq0(ipr,:,ntcw),' kdt=',kdt
!       if (lprnt) write(0,*)' qrna=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwa=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglba',qgl(ipr,:),' kdt=',kdt
!


      end subroutine m_micro_post_run

      subroutine m_micro_post_finalize()
      end subroutine m_micro_post_finalize

      end module m_micro_post
