!> \file ozphys.f
!! This file is ozone sources and sinks (previous version).


!> This module contains the CCPP-compliant Ozone photochemistry scheme.
      module ozphys

      contains

! \brief Brief description of the subroutine
!
!> \section arg_table_ozphys_init Argument Table
!! \htmlinclude ozphys_init.html
!!
      subroutine ozphys_init(oz_phys, errmsg, errflg)

      implicit none
      logical,          intent(in)  :: oz_phys
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (.not.oz_phys) then
        write (errmsg,'(*(a))') 'Logic error: oz_phys == .false.'
        errflg = 1
        return
      endif

      end subroutine ozphys_init

      subroutine ozphys_finalize()
      end subroutine ozphys_finalize

!>\defgroup GFS_ozphys GFS ozphys Main
!! \brief The operational GFS currently parameterizes ozone production and
!! destruction based on monthly mean coefficients (\c global_o3prdlos.f77) provided by Naval
!! Research Laboratory through CHEM2D chemistry model
!! (McCormack et al. (2006) \cite mccormack_et_al_2006).
!! \section arg_table_ozphys_run Argument Table
!! \htmlinclude ozphys_run.html
!!
!> \section genal_ozphys GFS ozphys_run General Algorithm
!> @{
      subroutine ozphys_run (                                           &
     &  im, levs, ko3, dt, oz, tin, po3,                                &
     &  prsl, prdout, oz_coeff, delp, ldiag3d,                          &
     &  ntoz, dtend, dtidx, index_of_process_prod_loss,                 &
     &  index_of_process_ozmix, index_of_process_temp,                  &
     &  index_of_process_overhead_ozone, con_g, me, errmsg, errflg)
!
!     this code assumes that both prsl and po3 are from bottom to top
!     as are all other variables
!
      use machine , only : kind_phys
      implicit none
!
      ! Interface variables
      integer, intent(in) :: im, levs, ko3, oz_coeff, me
      real(kind=kind_phys), intent(inout) :: oz(im,levs)
      ! The dtend array may not be allocated and needs an assumed array size
      real(kind=kind_phys), intent(inout) :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:), ntoz,                          &
     &  index_of_process_prod_loss, index_of_process_ozmix,             &
     &  index_of_process_temp, index_of_process_overhead_ozone
      real(kind=kind_phys), intent(in) ::                               &
     &                     dt, po3(ko3), prdout(im,ko3,oz_coeff),       &
     &                     prsl(im,levs), tin(im,levs), delp(im,levs),  &
     &                     con_g
      real :: gravi
      logical, intent(in) :: ldiag3d
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
      ! Local variables
      integer k,kmax,kmin,l,i,j, idtend(4)
      logical flg(im)
      real(kind=kind_phys) pmax, pmin, tem, temp
      real(kind=kind_phys) wk1(im), wk2(im), wk3(im), prod(im,oz_coeff),
     &                     ozib(im),  colo3(im,levs+1), ozi(im,levs)
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!     save input oz in ozi
      ozi = oz
      gravi=1.0/con_g


      if(ldiag3d) then
         idtend(1) = dtidx(100+ntoz,index_of_process_prod_loss)          ! was ozp1
         idtend(2) = dtidx(100+ntoz,index_of_process_ozmix)              ! was ozp2
         idtend(3) = dtidx(100+ntoz,index_of_process_temp)               ! was ozp3
         idtend(4) = dtidx(100+ntoz,index_of_process_overhead_ozone)     ! was ozp4
      else
         idtend=0
      endif

!
!> - Calculate vertical integrated column ozone values.
      if (oz_coeff > 2) then
        colo3(:,levs+1) = 0.0
        do l=levs,1,-1
          do i=1,im
            colo3(i,l) = colo3(i,l+1) + ozi(i,l) * delp(i,l) * gravi 
          enddo
        enddo
      endif
!
!> - Apply vertically linear interpolation to the ozone coefficients. 
      do l=1,levs
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
          do j=1,oz_coeff
            do i=1,im
              if (flg(i)) then
                prod(i,j)  = wk2(i) * prdout(i,k,j)
     &                     + wk3(i) * prdout(i,k+1,j)
              endif
            enddo
          enddo
        enddo
!
        do j=1,oz_coeff
          do i=1,im
            if (wk1(i) < po3(ko3)) then
              prod(i,j) = prdout(i,ko3,j)
            endif
            if (wk1(i) >= po3(1)) then
              prod(i,j) = prdout(i,1,j)
            endif
          enddo
        enddo

        if (oz_coeff == 2) then
          do i=1,im
            ozib(i)   = ozi(i,l)           ! no filling
            oz(i,l)   = (ozib(i) + prod(i,1)*dt) / (1.0 + prod(i,2)*dt)
          enddo
!
          if(idtend(1)>=1) then
             dtend(:,l,idtend(1)) = dtend(:,l,idtend(1)) +              ! was ozp1
     &            prod(:,1)*dt
          endif
          if(idtend(2)>=1) then
             dtend(:,l,idtend(2)) = dtend(:,l,idtend(2)) +              ! was ozp2
     &            (oz(:,l) - ozib(:))
          endif
        endif
!> - Calculate the 4 terms of prognostic ozone change during time \a dt:  
!!  - ozp1(:,:) - Ozone production from production/loss ratio 
!!  - ozp2(:,:) - Ozone production from ozone mixing ratio 
!!  - ozp3(:,:) - Ozone production from temperature term at model layers 
!!  - ozp4(:,:) - Ozone production from column ozone term at model layers
        if (oz_coeff == 4) then
          do i=1,im
            ozib(i)  = ozi(i,l)            ! no filling
            tem      = prod(i,1) + prod(i,3)*tin(i,l)
     &                           + prod(i,4)*colo3(i,l+1)
!     if (me .eq. 0) print *,'ozphys tem=',tem,' prod=',prod(i,:)
!    &,' ozib=',ozib(i),' l=',l,' tin=',tin(i,l),'colo3=',colo3(i,l+1)
            oz(i,l) = (ozib(i)  + tem*dt) / (1.0 + prod(i,2)*dt)
          enddo
          if(idtend(1)>=1) then
            dtend(:,l,idtend(1)) = dtend(:,l,idtend(1)) +              ! was ozp1
     &            prod(:,1)*dt
          endif
          if(idtend(2)>=1) then
            dtend(:,l,idtend(2)) = dtend(:,l,idtend(2)) +              ! was ozp2
     &            (oz(:,l)-ozib(:))
          endif
          if(idtend(3)>=1) then
            dtend(:,l,idtend(3)) = dtend(:,l,idtend(3)) +              ! was ozp3
     &            prod(:,3)*tin(:,l)*dt
          endif
          if(idtend(4)>=1) then
            dtend(:,l,idtend(4)) = dtend(:,l,idtend(4)) +              ! was ozp4
     &            prod(:,4)*colo3(:,l+1)*dt
          endif
        endif
      enddo                                ! vertical loop
!
      return
      end subroutine ozphys_run
!> @}

      end module ozphys
