!> \file ozphys_2015.f
!! This file is ozone sources and sinks.


!> This module contains the CCPP-compliant Ozone 2015 photochemistry scheme.
      module ozphys_2015

      contains

!> \section arg_table_ozphys_2015_init Argument Table
!! \htmlinclude ozphys_2015_init.html
!!
      subroutine ozphys_2015_init(oz_phys_2015, errmsg, errflg)

      implicit none
      logical,          intent(in)  :: oz_phys_2015
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (.not.oz_phys_2015) then
        write (errmsg,'(*(a))') 'Logic error: oz_phys_2015 == .false.'
        errflg = 1
        return
      endif

      end subroutine ozphys_2015_init

      subroutine ozphys_2015_finalize()
      end subroutine ozphys_2015_finalize

!>\defgroup GFS_ozphys_2015 GFS Ozone Photochemistry (2015) Scheme Module
!! \brief The operational GFS currently parameterizes ozone production and
!! destruction based on monthly mean coefficients (
!! \c ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77) provided by Naval
!! Research Laboratory through CHEM2D chemistry model
!! (McCormack et al. (2006) \cite mccormack_et_al_2006).
!! \section arg_table_ozphys_2015_run Argument Table
!! \htmlinclude ozphys_2015_run.html
!!
!> \section genal_ozphys_2015 GFS ozphys_2015_run General Algorithm
!> @{
!> -  This code assumes that both prsl and po3 are from bottom to top
!!     as are all other variables.
!> -  This code is specifically for NRL parameterization and
!!     climatological T and O3 are in location 5 and 6 of prdout array
!!\author June 2015 - Shrinivas Moorthi
      subroutine ozphys_2015_run (                                      &
     &     im, levs, ko3, dt, oz, tin, po3, prsl, prdout, pl_coeff,     &
     &     delp, ldiag3d, dtend, dtidx, ntoz, index_of_process_prod_loss&
     &     , index_of_process_ozmix, index_of_process_temp,             &
     &     index_of_process_overhead_ozone, con_g, me, errmsg, errflg)
!
!
      use machine , only : kind_phys
      implicit none
!
      real(kind=kind_phys),intent(in) :: con_g
      real :: gravi
      integer, intent(in) :: im, levs, ko3, pl_coeff,me
      real(kind=kind_phys), intent(in) :: po3(ko3),                     &
     &                                    prsl(im,levs), tin(im,levs),  &
     &                                    delp(im,levs),                &
     &                                    prdout(im,ko3,pl_coeff), dt
      ! dtend may not be allocated and needs an assumed array size
      real(kind=kind_phys), intent(inout) :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:), ntoz,                          &
     &  index_of_process_prod_loss, index_of_process_ozmix,             &
     &  index_of_process_temp, index_of_process_overhead_ozone
      real(kind=kind_phys), intent(inout) :: oz(im,levs)


      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer k,kmax,kmin,l,i,j, idtend(4)
      logical              ldiag3d, flg(im), qdiag3d
      real(kind=kind_phys) pmax, pmin, tem, temp
      real(kind=kind_phys) wk1(im), wk2(im), wk3(im),prod(im,pl_coeff), &
     &                     ozib(im), colo3(im,levs+1), coloz(im,levs+1),&
     &                     ozi(im,levs)
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if(ldiag3d) then
         idtend(1) = dtidx(100+ntoz,index_of_process_prod_loss)          ! was ozp1
         idtend(2) = dtidx(100+ntoz,index_of_process_ozmix)              ! was ozp2
         idtend(3) = dtidx(100+ntoz,index_of_process_temp)               ! was ozp3
         idtend(4) = dtidx(100+ntoz,index_of_process_overhead_ozone)     ! was ozp4
      else
         idtend=0
      endif

!ccpp: save input oz in ozi
      ozi = oz
      gravi=1.0/con_g

        colo3(:,levs+1) = 0.0
        coloz(:,levs+1) = 0.0
!
      do l=levs,1,-1
        pmin =  1.0e10
        pmax = -1.0e10
!
        do i=1,im
          wk1(i) = log(prsl(i,l))
          pmin   = min(wk1(i), pmin)
          pmax   = max(wk1(i), pmax)
          prod(i,:) = 0.0
        enddo
        kmax = 1
        kmin = 1
        do k=1,ko3-1
          if (pmin < po3(k)) kmax = k
          if (pmax < po3(k)) kmin = k
        enddo
!
        do k=kmin,kmax
          temp = 1.0 / (po3(k) - po3(k+1))
          do i=1,im
            flg(i) = .false.
            if (wk1(i) < po3(k) .and. wk1(i) >= po3(k+1)) then
              flg(i) = .true.
              wk2(i) = (wk1(i) - po3(k+1)) * temp
              wk3(i) = 1.0 - wk2(i)
            endif
          enddo
          do j=1,pl_coeff
            do i=1,im
              if (flg(i)) then
                prod(i,j)  = wk2(i) * prdout(i,k,j)
     &                     + wk3(i) * prdout(i,k+1,j)
              endif
            enddo
          enddo
        enddo
!
        do j=1,pl_coeff
          do i=1,im
            if (wk1(i) < po3(ko3)) then
              prod(i,j) = prdout(i,ko3,j)
            endif
            if (wk1(i) >= po3(1)) then
              prod(i,j) = prdout(i,1,j)
            endif
          enddo
        enddo
        do i=1,im
          colo3(i,l) = colo3(i,l+1) + ozi(i,l)  * delp(i,l)*gravi
          coloz(i,l) = coloz(i,l+1) + prod(i,6) * delp(i,l)*gravi
          prod(i,2)  = min(prod(i,2), 0.0)
        enddo
!       write(1000+me,*) ' colo3=',colo3(1,l),' coloz=',coloz(1,l)
!    &,' l=',l
        do i=1,im
          ozib(i)  = ozi(i,l)            ! no filling
          tem      = prod(i,1) - prod(i,2) * prod(i,6)
     &             + prod(i,3) * (tin(i,l) - prod(i,5))
     &             + prod(i,4) * (colo3(i,l)-coloz(i,l))

!     if (me .eq. 0) print *,'ozphys_2015 tem=',tem,' prod=',prod(i,:)
!    &,' ozib=',ozib(i),' l=',l,' tin=',tin(i,l),'colo3=',colo3(i,l+1)

!ccpp            ozo(i,l) = (ozib(i)  + tem*dt) / (1.0 - prod(i,2)*dt)
          oz(i,l) = (ozib(i)  + tem*dt) / (1.0 - prod(i,2)*dt)
        enddo
        if(idtend(1)>=1) then
           dtend(:,l,idtend(1)) = dtend(:,l,idtend(1)) + ! was ozp1
     &          (prod(:,1)-prod(:,2)*prod(:,6))*dt
        endif
        if(idtend(2)>=1) then
           dtend(:,l,idtend(2)) = dtend(:,l,idtend(2)) + ! was ozp2
     &          (oz(:,l) - ozib(:))
        endif
        if(idtend(3)>=1) then
           dtend(:,l,idtend(3)) = dtend(:,l,idtend(3)) + ! was ozp3
     &          prod(:,3)*(tin(:,l)-prod(:,5))*dt
        endif
        if(idtend(4)>=1) then
           dtend(:,l,idtend(4)) = dtend(:,l,idtend(4)) + ! was ozp4
     &       prod(:,4) * (colo3(:,l)-coloz(:,l))*dt
        endif
      enddo                                ! vertical loop
!
      return
      end subroutine ozphys_2015_run

!> @}

      end module ozphys_2015
