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
        im, levs, fprcp, mg3_as_mg2,      &
        gq0_ice, gq0_rain, gq0_snow, gq0_graupel, ten_qi, ten_qr, ten_qs, ten_qg, &
        ice, snow, graupel, dtp, errmsg, errflg)

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs, fprcp
      logical, intent(in) :: mg3_as_mg2

      real(kind=kind_phys), intent(in   ) :: gq0_ice(:,:)
      real(kind=kind_phys), intent(in   ) :: gq0_rain(:,:)
      real(kind=kind_phys), intent(in   ) :: gq0_snow(:,:)
      real(kind=kind_phys), intent(in   ) :: gq0_graupel(:,:)
      real(kind=kind_phys), intent(in   ) :: ten_qi(:,:)
      real(kind=kind_phys), intent(inout) :: ten_qr(:,:)
      real(kind=kind_phys), intent(inout) :: ten_qs(:,:)
      real(kind=kind_phys), intent(inout) :: ten_qg(:,:)
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
      real(kind=kind_phys) :: tem, new_qi, new_qr, new_qs, new_qg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      tem = dtp * con_p001 / con_day
      if (abs(fprcp) == 1 .or. mg3_as_mg2) then
        do k=1,levs
          do i=1,im
            new_qr = gq0_rain(i,k) + dtp*ten_qr(i,k)
            !zero out qr when tendencies are applied if after-application value is small or negative
            if (new_qr < qsmall) then
              ten_qr(i,k) = -gq0_rain(i,k)/dtp 
            end if
            
            new_qs = gq0_snow(i,k) + dtp*ten_qs(i,k)
            !zero out qs when tendencies are applied if after-application value is small or negative
            if (new_qs < qsmall) then
              ten_qs(i,k) = -gq0_snow(i,k)/dtp 
            end if
          enddo
        enddo
        do i=1,im
          new_qi = gq0_ice(i,1) + dtp*ten_qi(i,1)
          ice(i)  = tem * new_qi
          new_qs = gq0_snow(i,1) + dtp*ten_qs(i,1)
          snow(i) = tem * new_qs
        enddo
      elseif (fprcp > 1) then
        do k=1,levs
          do i=1,im
            new_qr = gq0_rain(i,k) + dtp*ten_qr(i,k)
            !zero out qr when tendencies are applied if after-application value is small or negative
            if (new_qr < qsmall) then
              ten_qr(i,k) = -gq0_rain(i,k)/dtp 
            end if
            
            new_qs = gq0_snow(i,k) + dtp*ten_qs(i,k)
            !zero out qs when tendencies are applied if after-application value is small or negative
            if (new_qs < qsmall) then
              ten_qs(i,k) = -gq0_snow(i,k)/dtp 
            end if
            
            
            new_qg = gq0_graupel(i,k) + dtp*ten_qg(i,k)
            !zero out qg when tendencies are applied if after-application value is small or negative
            if (new_qg < qsmall) then
              ten_qg(i,k) = -gq0_graupel(i,k)/dtp 
            end if
          enddo
        enddo
        do i=1,im
          new_qi = gq0_ice(i,1) + dtp*ten_qi(i,1)
          ice(i)     = tem * new_qi
          new_qs = gq0_snow(i,1) + dtp*ten_qs(i,1)
          snow(i) = tem * new_qs
          new_qg = gq0_graupel(i,1) + dtp*ten_qg(i,1)
          graupel(i) = tem * new_qg
        enddo

      endif

      end subroutine m_micro_post_run

      end module m_micro_post
