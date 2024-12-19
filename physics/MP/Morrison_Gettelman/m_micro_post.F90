!> \file m_micro_post.F90
!! This file contains subroutines that prepare data from the Morrison-Gettelman microphysics scheme
!! as part of the GFS physics suite.
      
      module m_micro_post

      implicit none

      contains

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

      real(kind=kind_phys), intent(in   ),optional :: ncpr(:,:)
      real(kind=kind_phys), intent(in   ),optional :: ncps(:,:)
      real(kind=kind_phys), intent(in   ),optional :: ncgl(:,:)
      real(kind=kind_phys), intent(inout),optional :: qrn(:,:)
      real(kind=kind_phys), intent(inout),optional :: qsnw(:,:)
      real(kind=kind_phys), intent(inout),optional :: qgl(:,:)
      real(kind=kind_phys), intent(in   ) :: gq0_ice(:,:)
      real(kind=kind_phys), intent(out  ) :: gq0_rain(:,:)
      real(kind=kind_phys), intent(out  ) :: gq0_snow(:,:)
      real(kind=kind_phys), intent(out  ) :: gq0_graupel(:,:)
      real(kind=kind_phys), intent(out  ) :: gq0_rain_nc(:,:)
      real(kind=kind_phys), intent(out  ) :: gq0_snow_nc(:,:)
      real(kind=kind_phys), intent(out  ) :: gq0_graupel_nc(:,:)
      real(kind=kind_phys), intent(  out) :: ice(:)
      real(kind=kind_phys), intent(  out) :: snow(:)
      real(kind=kind_phys), intent(  out) :: graupel(:)
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

      end module m_micro_post
