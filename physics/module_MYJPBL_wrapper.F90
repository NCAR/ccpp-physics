!> \file module_myjpbl_wrapper.F90
!!  Contains all of the code related to running the MYJ PBL scheme

      MODULE myjpbl_wrapper

      USE machine, only: kfpt => kind_phys, &
                         kind_phys

      contains

      subroutine myjpbl_wrapper_init (do_myjpbl,errmsg,errflg)
      
      logical,              intent(in)  :: do_myjpbl
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

     ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

    ! Consistency checks
      if (.not. do_myjpbl) then
        write(errmsg,fmt='(*(a))') 'Logic error: do_myjpbl=.false.'        
        errflg = 1
        return
      end if
      end subroutine myjpbl_wrapper_init

      subroutine myjpbl_wrapper_finalize ()
      end subroutine myjpbl_wrapper_finalize

!!
!> \brief This scheme (1) performs pre-myjpbl work, (2) runs the myjpbl, and (3) performs post-myjpbl work
!! \section arg_table_myjpbl_wrapper_run Argument Table
!! \htmlinclude myjpbl_wrapper_run.html
!!
!###===================================================================
  SUBROUTINE myjpbl_wrapper_run(                    &
     &  restart,do_myjsfc,                          &
     &  im,levs,dt_phs,                             &
     &  kdt,ntrac,ntke,                             &
     &  ntcw,ntiw,ntrw,ntsw,ntgl,                   &
     &  ugrs, vgrs, tgrs, qgrs,                     &
     &  prsl, prsi, phii, hprime1,                  &
     &  prsik_1, prslk_1, prslki, tsfc, qsfc,       &
     &  phy_myj_qsfc, phy_myj_thz0, phy_myj_qz0,    &
     &  phy_myj_uz0, phy_myj_vz0, phy_myj_z0base,   &
     &  phy_myj_akhs, phy_myj_akms,                 &
     &  phy_myj_chkqlm, phy_myj_elflx,              &
     &  phy_myj_a1u, phy_myj_a1t, phy_myj_a1q,      &
     &  pblh, kpbl, kinver, slmsk,                  &
     &  garea, ustar, cm, ch, wind,                 &
     &  snowd, zorl, evap, hflx,                    &
     &  dudt, dvdt, dtdt, dqdt,                     &
     &  dusfc,dvsfc,dtsfc,dqsfc,                    &
     &  dkt,xkzm_m, xkzm_h,xkzm_s, gamt,gamq,       &
     &  con_cp,con_g,con_rd,                        &
     &  me, lprnt, gen_tend, ldiag3d, dtend, dtidx, &
     &  index_of_temperature, index_of_x_wind,      &
     &  index_of_y_wind, index_of_process_pbl,      &
     &  ntqv, errmsg, errflg )

!

      use MODULE_BL_MYJPBL,      only: MYJPBL_INIT,MYJPBL

!-------------------------------------------------------------------
      implicit none

!      integer,parameter:: &
!        klog=4 &                   ! logical variables
!       ,kint=4 &                   ! integer variables
!       !,kfpt=4 &                   ! floating point variables
!       ,kfpt=8 &                   ! floating point variables
!       ,kdbl=8                     ! double precision

!-------------------------------------------------------------------
!  ---  constant parameters:
!For reference
!   real    , parameter :: karman       = 0.4
!   real    , parameter :: g            = 9.81
!   real    , parameter :: r_d          = 287.
!   real    , parameter :: cp           = 7.*r_d/2.
!
!      real, parameter :: g = 9.81, r_d=287., cp= 7.*r_d/2.
!      real, parameter :: rd=r_d, rk=cp/rd
!      real, parameter :: elwv=2.501e6, eliv=2.834e6
!      real, parameter :: reliw=eliv/elwv,
      real, parameter :: xkgdx=25000.,xkzinv=0.15

!      real, parameter :: g_inv=1./con_g, cappa=con_rd/con_cp

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      real(kind=kind_phys), intent(inout), optional :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:)
      integer, intent(in) :: index_of_temperature, index_of_x_wind, &
     &                       index_of_y_wind, index_of_process_pbl, ntqv

!MYJ-1D
      integer,intent(in) :: im, levs
      integer,intent(in) :: kdt, me
      integer,intent(in) :: ntrac,ntke,ntcw,ntiw,ntrw,ntsw,ntgl
      logical,intent(in) :: restart,do_myjsfc,lprnt,ldiag3d,gen_tend
      real(kind=kind_phys),intent(in) :: con_cp, con_g, con_rd
      real(kind=kind_phys),intent(in) :: dt_phs, xkzm_m, xkzm_h, xkzm_s

!MYJ-2D
      real(kind=kind_phys),dimension(:),intent(in) ::        &
     &     prsik_1, prslk_1, prslki, slmsk, garea,           &
           snowd, evap, hflx, cm, ch, wind, hprime1
      real(kind=kind_phys),dimension(:),intent(inout) ::     &
     &     zorl, ustar, tsfc, qsfc
      real(kind=kind_phys),dimension(:),intent(inout)   ::   &
     &        phy_myj_qsfc, phy_myj_thz0, phy_myj_qz0,       &
     &        phy_myj_uz0, phy_myj_vz0, phy_myj_z0base,      &
     &        phy_myj_akhs, phy_myj_akms,                    &
     &        phy_myj_chkqlm, phy_myj_elflx,                 &
     &        phy_myj_a1u, phy_myj_a1t, phy_myj_a1q
      real(kind=kind_phys),dimension(:),intent(out) ::       &
     &     pblh,dusfc,dvsfc,dtsfc,dqsfc,gamt,gamq
      integer,dimension(:),intent(out) :: kpbl
      integer,dimension(:),intent(in) ::  kinver

!MYJ-3D
      real(kind=kind_phys),dimension(:,:),intent(in) ::      &
              phii, prsi
      real(kind=kind_phys),dimension(:,:),intent(in) ::      &
     &        ugrs, vgrs, tgrs, prsl
!      real(kind=kind_phys),dimension(:,:),intent(inout)  :: &
!             dudt, dvdt, dtdt, dkt
      real(kind=kind_phys),dimension(:,:),intent(inout)   :: &
             dudt, dvdt, dtdt
      real(kind=kind_phys),dimension(:,:),intent(out)     :: &
             dkt

!MYJ-4D
      real(kind=kind_phys),dimension(:,:,:),intent(inout) :: &
     &       qgrs,dqdt

!LOCAL
      integer :: ntsd, k, k1, i, kx1
      integer :: i_min, i_max, k_min, k_max

      logical :: lprnt1,lprnt2
      integer :: ict, ide, lm, me1
      real(kind=kfpt) :: dt_myj, tem, tem1, tem2, ptem
      integer,dimension(im) :: kpbl_myj
      real(kind=kfpt),dimension(1:levs-1):: epsl
      real(kind=kfpt),dimension(1:levs):: epsq2
      real(kind=kfpt),dimension(im) ::                     &
           xland, sice, snowd1, ht, stdh, tsk,             &
           ustar1,z0,pblh_myj,                             &
           elflx,mixht,ct
      real(kind=kfpt), dimension(im,levs) ::               &
     &        u_myj, v_myj, t_myj, q_myj, th_myj,          &
     &        cw, dz_myj, pmid, q2, exner, del
      real(kind=kfpt), dimension(im,levs+1) :: pint
      real(kind=kfpt),dimension(im,levs) ::                &
              rublten,rvblten,rthblten,rqvblten,rqcblten
      real(kind=kfpt),dimension(im,levs) :: el_myj
      real(kind=kfpt),dimension(im) ::                     &
              dusfc1,dvsfc1,dtsfc1,dqsfc1
      real(kind=kfpt),dimension(im) :: thlm,qlm
      real(kind=kfpt),dimension(im,13) :: phy_f2d_myj
      real(kind=kfpt), dimension(im,levs) :: xcofh  &
     &        ,xkzo,xkzmo
      real(kind=kind_phys) :: g, r_d, g_inv, cappa
      real(kind=kind_phys) :: thz0, qz0, a1u, a1t, a1q
      real(kind=kind_phys) :: z0m, aa1u, aa1t, z1uov, z1tox
      real(kind=kind_phys) :: tmax,tmin,t_myj1
      real(kind=kind_phys),dimension(im)       ::           &
     &        thsfc,sfcz,tsfc1,                              &
     &        sm,work3,wind1,work4                            &
     &        ,rho,qfc1,gdx,xkzm_hx,xkzm_mx,tx1, tx2
!      real(kind=kind_phys), dimension(im,levs,ntrac) ::    &
!     &        qgrs_myj
      real(kind=kind_phys),dimension(im,levs) :: dkt2
      integer :: uidx, vidx, tidx, qidx

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0


!      if (lprnt) then
!         write(0,*)"=============================================="
!         write(0,*)"in myj wrapper..."
!      endif

      ntsd = kdt - 1

      lprnt1=.false.
      lprnt2=.false.

      if (lprnt1) then
         if(me.eq.0)print*,'ntsd=', ntsd
      end if

!prep MYJ-only variables

      r_d   = con_rd
      g     = con_g
      g_inv = 1./con_g
      cappa = con_rd/con_cp

      do i=1,im
         work3(i)=prsik_1(i) / prslk_1(i)
         sice(i)=slmsk(i)*0.5
         if(sice(i) < 0.7)sice(i)=0
         sm(i)=1.; if(slmsk(i) > 0.5 ) sm(i)=0.
         z0(i)=zorl(i)*0.01
         xland(i)=sm(i)+1.
         sfcz(i)=phii(i,1)*g_inv
         work4(i)=(1.e5/prsi(i,1))**cappa
         thsfc(i)=tsfc(i)*work4(i)             ! thsfc
      enddo

      do k=1,levs
         k1=levs+1-k
         do i=1,im
            u_myj(i,k)=ugrs(i,k1)
            v_myj(i,k)=vgrs(i,k1)
            t_myj(i,k)=tgrs(i,k1)
            q_myj(i,k)=qgrs(i,k1,1)
            cw(i,k)   =qgrs(i,k1,ntcw)
!            if(ntrw.gt.0)cw(i,k) = cw(i,k) + qgrs(i,k1,ntrw)
!            if(ntiw.gt.0)cw(i,k) = cw(i,k) + qgrs(i,k1,ntiw)
!            if(ntsw.gt.0)cw(i,k) = cw(i,k) + qgrs(i,k1,ntsw)
!            if(ntgl.gt.0)cw(i,k) = cw(i,k) + qgrs(i,k1,ntgl)
            if(ntke.gt.0)then
               q2(i,k) =max(0.02,qgrs(i,k1,ntke)*2.)
            else
               q2(i,k) =0.02
            end if
!             fmid(i,k) =phil(i,k1)
            pmid(i,k) =prsl(i,k1)
            exner(i,k)=(prsl(i,k1)*1.e-5)**cappa
            th_myj(i,k)=tgrs(i,k1)/exner(i,k)
         end do
      end do
         do k=1,levs+1
           k1=levs+2-k
           do i=1,im
             pint(i,k) =prsi(i,k1)
           end do
         end do

      do i=1,im
        gdx(i) = sqrt(garea(i))
      enddo

      do i=1,im
        kx1 = 1
        tx1(i) = 1.0 / prsi(i,1)
        tx2(i) = tx1(i)
        if(gdx(i) >= xkgdx) then
          xkzm_hx(i) = xkzm_h
          xkzm_mx(i) = xkzm_m
        else
          tem  = 1. / (xkgdx - 5.)
          tem1 = (xkzm_h - 0.01) * tem
          tem2 = (xkzm_m - 0.01) * tem
          ptem = gdx(i) - 5.
          xkzm_hx(i) = 0.01 + tem1 * ptem
          xkzm_mx(i) = 0.01 + tem2 * ptem
        endif
      enddo
      xkzo  = 0.0
      xkzmo = 0.0
      do k = 1,levs-1
        do i=1,im
          if (k < kinver(i)) then
!                                  vertical background diffusivity
            ptem      = prsi(i,k+1) * tx1(i)
            tem1      = 1.0 - ptem
            tem1      = tem1 * tem1 * 10.0
            xkzo(i,k) = xkzm_hx(i) * min(1.0, exp(-tem1))
            xkzo(i,k) = min(xkzo(i,k),xkzinv)
!                                 vertical background diffusivity for momentum
            if (ptem >= xkzm_s) then
              xkzmo(i,k) = xkzm_mx(i)
              kx1     = k + 1
            else
              if (k == kx1 .and. k > 1) tx2(i) = 1.0 / prsi(i,k)
              tem1 = 1.0 - prsi(i,k+1) * tx2(i)
              tem1 = tem1 * tem1 * 5.0
              xkzmo(i,k) = xkzm_mx(i) * min(1.0, exp(-tem1))
              xkzmo(i,k) = min(xkzmo(i,k),xkzinv)
            endif
          endif
        enddo
      enddo

! change vertical coordinate
      do k=1,levs
         k1=levs+1-k
         do i=1,im
             xcofh(i,k1)=xkzo(i,k)             ! temp use xcofh
             el_myj(i,k1)=xkzmo(i,k)           ! temp use EL_MYJ
         end do
      end do

      do k=1,levs
      do i=1,im
         xkzo(i,k)=xcofh(i,k)
         xkzmo(i,k)=el_myj(i,k)
      end do
      end do

      do k=1,levs-1
         epsq2(k)=0.02
         epsl(k)=sqrt(epsq2(k)*0.5)
!          if (xkzo(i,k) .gt. 0.01) then
!             epsl(k)=1.0
!          end if
      end do
      epsq2(levs)=epsq2(levs-1)

      do k = 1, levs
         k1 = levs-k+1
         do i = 1, im
            del(i,k) = prsi(i,k1) - prsi (i,k1+1)
            dz_myj(i,k) = (phii(i,k1+1)-phii(i,k1)) * g_inv
         enddo
      enddo

      do i = 1, im
         wind1(i)=max(wind(i),1.0)
      end do

      if(.not.do_myjsfc)then
         do i=1,im
            if(sm(i).gt.0.5.and.sice(i).le.0.5) then
               z0m=max(0.018*g_inv*ustar(i)*ustar(i),1.59E-5)
               z1uov=0.35*30.*sqrt(sqrt(z0m*ustar(i)/1.5E-5))/ustar(i)
               aa1u=cm(i)*wind1(i)*z1uov
               a1u=aa1u/(1.-aa1u)
               z1tox=0.84*z1uov
               aa1t=ch(i)*wind1(i)*z1tox
               a1t=aa1t/(1.-aa1t)
!
!               a1u=0.3
!               a1t=0.25
!
               a1q=a1t
            else
               z0m=zorl(i)*0.01
               a1u=0.
               a1t=0.
               a1q=0.
            end if
            phy_myj_a1u(i)   = a1u
            phy_myj_a1t(i)   = a1t
            phy_myj_a1q(i)   = a1q
            phy_myj_akhs(i)  = ch(i)*wind1(i)*(1.+a1t)
            phy_myj_akms(i)  = cm(i)*wind1(i)*(1.+a1u)
            phy_myj_uz0(i)   = u_myj(i,levs)*a1u/(1.+a1u)
            phy_myj_vz0(i)   = v_myj(i,levs)*a1u/(1.+a1u)
            phy_myj_z0base(i)= z0m

            if(ntsd.eq.0)then
!             if(sm(i).gt.0.5)then
                 qz0=max(evap(i)/phy_myj_akhs(i)+q_myj(i,levs),1.e-9)
                 thz0=hflx(i)/phy_myj_akhs(i)+th_myj(i,levs)
!             else
!                if(sice(i).gt.0.5)then
!                    qsfc(i)=qss_ice(i)
!                else
!                    qsfc(i)=qss_land(i)
!                end if
!             endif
              phy_myj_thz0(i)  = thz0
              phy_myj_qz0(i)   = qz0
            end if
            if(cw(i,levs).gt.1.e-9)then
               phy_myj_chkqlm(i)= 0.
            else
               phy_myj_chkqlm(i)= 1.
            end if
         end do
      end if

      if(do_myjsfc)then
         do i=1,im
            phy_myj_akhs(i)=phy_myj_akhs(i)*wind1(i)/wind(i)
            phy_myj_akms(i)=phy_myj_akms(i)*wind1(i)/wind(i)
         end do
      end if

! update qsfc, thz0, qz0 and elflx after Land/Ocean model.
      do i=1,im
         phy_myj_elflx(i) = evap(i)
         qsfc(i)=max(evap(i)/(ch(i)*wind1(i))+q_myj(i,levs),1.e-9)
         tsfc1(i)=(hflx(i)/(ch(i)*wind1(i))+th_myj(i,levs))/work4(i)
         phy_myj_qsfc(i) = qsfc(i)
         thz0 = phy_myj_thz0(i)
         thlm(i)=th_myj(i,levs)
         qlm(i)=q_myj(i,levs)
!        a1t=phy_myj_a1t(i)
!        thsfc(i)=hflx(i)/phy_myj_akhs(i)+th_myj(i,levs)
!        phy_myj_thz0(i)=((a1t*thlm(i)+thsfc(i))/(a1t+1.)+thz0)*0.5           ! thz0
         phy_myj_thz0(i)=0.5*(thz0+                 &
             hflx(i)/phy_myj_akhs(i)+th_myj(i,levs))
!        a1q=phy_myj_a1q(i)
         qz0=phy_myj_qz0(i)
!        phy_myj_qz0(i) = ((a1q*q_myj(i,levs)+qsfc(i))/(a1q+1.)+qz0)*0.5
         phy_myj_qz0(i) = 0.5*(qz0+             &
             max(evap(i)/phy_myj_akhs(i)+q_myj(i,levs),1.e-9))
      enddo

      rthblten = 0.
      rqvblten = 0.
      rqcblten= 0.
      rublten = 0.
      rvblten = 0.
!      rtrblten= 0.
      xcofh=0.
      kpbl(:)=levs-1
      ict=1      !  no longer used

      if (lprnt1) then

          if (me.eq.0.and.ntsd.lt.2)then
            print*,'Qingfu test starts PBL'
            print*,'ntsd,me,im,levs,ict=',ntsd,me,im,levs,ict
            print*,'dt_phs,sfcz,dz_myj=',dt_phs,sfcz(1),dz_myj(1,5)
            print*,'pmid,pint,th_myj=',pmid(1,5),pint(1,5),th_myj(1,5)
            print*,'t_myj,exner,q_myj=',t_myj(1,5),exner(1,5),q_myj(1,5)
            print*,'cw,u_myj,v_myj=',cw(1,5),u_myj(1,5),v_myj(1,5)
            print*,'tsfc,xland,sice,snowd=',tsfc(1),xland(1),sice(1),snowd(1)
            print*,'ustar,z0,pblh,kpbl=',ustar(1),z0(1),pblh(1),kpbl(1)
            print*,'q2,xcofh=',q2(1,5),xcofh(1,5)
!            print*,'Tbd%phy_f2d_myj(1,1-5)=',(Tbd%phy_f2d_myj(1,i),i=1,5)
!            print*,'Tbd%phy_f2d_myj(1,6-10)=',(Tbd%phy_f2d_myj(1,i),i=6,10)
!            print*,'Tbd%phy_f2d_myj(1,11-13)=',(Tbd%phy_f2d_myj(1,i),i=11,13)
            print*,'thlm,thsfc=',thlm(i),thsfc(i)
          end if

         do k=1,levs
         do i=1,im
           if(t_myj(i,k).gt.390..or.t_myj(i,k).lt.110.)then
            print*,'Qingfu test starts PBL',i,k,t_myj(i,k)
            print*,'ntsd,me,im,levs,ict=',ntsd,me,im,levs,ict
            print*,'dt_phs,sfcz,dz_myj=',dt_phs,sfcz(i),dz_myj(i,k)
            print*,'pmid,pint,th_myj=',pmid(i,k),pint(i,k),th_myj(i,k)
            print*,'t_myj,exner,q_myj=',t_myj(i,k),exner(i,k),q_myj(i,k)
            print*,'cw,u_myj,v_myj=',cw(i,k),u_myj(i,k),v_myj(i,k)
            print*,'tsfc,xland,sice,snowd=',tsfc(i),xland(i),sice(i),snowd(i)
            print*,'ustar,z0,pblh,kpbl=',ustar(i),z0(i),pblh(i),kpbl(i)
            print*,'q2,xcofh=',q2(i,k),xcofh(i,k)
            end if
          end do
          end do

          tmax=-1.e-5
          tmin=1.e5
          do k=1,levs
           k1=levs+1-k
          do i=1,im
             if(tmax.lt.t_myj(i,k1))then
                tmax=t_myj(i,k1)
                i_max=i
                k_max=k
             end if
             if(tmin.gt.t_myj(i,k1))then
               tmin=t_myj(i,k1)
               i_min=i
               k_min=k
             end if
          end do
          end do
!          print*,'before i_min,k_min,i_max,k_max=',i_min,k_min,i_max,k_max
!          print*,'ntsd,me,tmin,tmax=',ntsd,me,tmin,tmax
!          if(me.eq.me1.and.tmin.lt.113.6.or.tmax.gt.350.)then
!             i=i_max
!             print*,'before bad bad tmin,tmax=',tmin,tmax,i_min,k_min,i_max,k_max
!             print*,'ntsd,me,tmin,tmax=',ntsd,me,tmin,tmax
!          end if

      end if

      ct=0.
      ide=im
      lm=levs
      dt_myj=dt_phs
      do i=1,im
         ustar1(i)=ustar(i)
         ht(i)=phii(i,1)*g_inv
         stdh(i)=hprime1(i)
         tsk(i)=tsfc(i)
         snowd1(i)=snowd(i)
         phy_f2d_myj(i,1)  = phy_myj_qsfc(i)
         phy_f2d_myj(i,2)  = phy_myj_thz0(i)
         phy_f2d_myj(i,3)  = phy_myj_qz0(i)
         phy_f2d_myj(i,4)  = phy_myj_uz0(i)
         phy_f2d_myj(i,5)  = phy_myj_vz0(i)
         phy_f2d_myj(i,6)  = phy_myj_z0base(i)
         phy_f2d_myj(i,7)  = phy_myj_akhs(i)
         phy_f2d_myj(i,8)  = phy_myj_akms(i)
         phy_f2d_myj(i,9)  = phy_myj_chkqlm(i)
         phy_f2d_myj(i,10) = phy_myj_elflx(i)
         phy_f2d_myj(i,11) = phy_myj_a1u(i)
         phy_f2d_myj(i,12) = phy_myj_a1t(i)
         phy_f2d_myj(i,13) = phy_myj_a1q(i)
!         do k=1,13
!            phy_f2d_myj(i,k)=Tbd%phy_f2d_myj(i,k)
!         end do
      end do

!      do i = 1, im
!          rho(i) = prsl(i,1)/(r_d*tgrs(i,1)      &
!                *(0.608*qgrs(i,1,1)+1.-qgrs(i,1,ntcw)))
!          if(sm(i).lt.0.5)then
!              qfc1(i)=elwv*rho(i)
!              if(snowd(i).gt.0..or.sice(i).gt.0.5)then
!                 qfc1(i)=qfc1(i)*reliw
!              end if
!          else
!              qfc1(i)=elwv*rho(i)
!          end if
!         phy_f2d_myj(i,10)=qfc1(i)*phy_f2d_myj(i,10)    ! convert units
!      end do

      if(ntsd.eq.0.or.restart)then
         if(.not.restart) xcofh=0.
            call MYJPBL_INIT(    &
                1,ide,1,1,lm,     &
                1,ide,1,1,       &
                1,ide,1,1)
      end if

      call MYJPBL(ntsd,me,dt_myj,epsl,epsq2,ht,stdh,dz_myj,del       &
                 ,pmid,pint,th_myj,t_myj,exner,q_myj,cw,u_myj,v_myj  &
                 ,tsk,phy_f2d_myj(1:im,1),phy_f2d_myj(1:im,9)        &
                 ,phy_f2d_myj(1:im,2),phy_f2d_myj(1:im,3)            &
                 ,phy_f2d_myj(1:im,4),phy_f2d_myj(1:im,5)            &
                 ,xland,sice,snowd1                                  &
                 ,q2,xcofh,ustar1,z0,el_myj,pblh_myj,kpbl_myj,ct     &
                 ,phy_f2d_myj(1:im,7),phy_f2d_myj(1:im,8)            &
                 ,phy_f2d_myj(1:im,10),mixht,thlm,qlm                &
                 ,rublten,rvblten,rthblten,rqvblten,rqcblten         &
                 ,dusfc1,dvsfc1,dtsfc1,dqsfc1,xkzo,xkzmo,ict         &
                 ,1,ide,1,1                                          &
                 ,1,ide,1,1                                          &
                 ,1,ide,1,1,lm)

      do i=1,im
         zorl(i)=z0(i)*100.
         dusfc(i)=dusfc1(i)
         dvsfc(i)=dvsfc1(i)
         dtsfc(i)=dtsfc1(i)
         dqsfc(i)=dqsfc1(i)
         pblh(i)=pblh_myj(i)
         kpbl(i)=levs-kpbl_myj(i)
!        ustar(i)=ustar1(i)
         phy_myj_qsfc(i)   = phy_f2d_myj(i,1)
         phy_myj_thz0(i)   = phy_f2d_myj(i,2)
         phy_myj_qz0(i)    = phy_f2d_myj(i,3)
         phy_myj_uz0(i)    = phy_f2d_myj(i,4)
         phy_myj_vz0(i)    = phy_f2d_myj(i,5)
         phy_myj_z0base(i) = phy_f2d_myj(i,6)
         phy_myj_akhs(i)   = phy_f2d_myj(i,7)
         phy_myj_akms(i)   = phy_f2d_myj(i,8)
         phy_myj_chkqlm(i) = phy_f2d_myj(i,9)
         phy_myj_elflx(i)  = phy_f2d_myj(i,10)
         phy_myj_a1u(i)    = phy_f2d_myj(i,11)
         phy_myj_a1t(i)    = phy_f2d_myj(i,12)
         phy_myj_a1q(i)    = phy_f2d_myj(i,13)
!         do k=1,13
!            Tbd%phy_f2d_myj(i,k)=phy_f2d_myj(i,k)
!         end do
      end do

      dkt=0.
      do k=1,levs
         k1=levs-k+1
         do i=1,im
!            dkt(i,k)=max(xcofh(i,k1),xkzo(i,k))
           dkt(i,k)=xcofh(i,k1)
         end do
      end do
      if(ntke.gt.0)then
         do k=1,levs
            k1=levs+1-k
            qgrs(:,k,ntke)=q2(:,k1)*0.5
         end do
      end if
      gamt=0.
      gamq=0.

      do k=1,levs
         k1=levs+1-k
         do i=1,im
            dudt(i,k)=dudt(i,k)+rublten(i,k1)
            dvdt(i,k)=dvdt(i,k)+rvblten(i,k1)
            dtdt(i,k)=dtdt(i,k)+rthblten(i,k1)*exner(i,k1)
            dqdt(i,k,1)=dqdt(i,k,1)+rqvblten(i,k1)
            dqdt(i,k,ntcw)=dqdt(i,k,ntcw)+rqcblten(i,k1)
         end do
      end do
      if (ldiag3d .and. .not. gen_tend) then
        uidx = dtidx(index_of_x_wind,index_of_process_pbl)
        vidx = dtidx(index_of_y_wind,index_of_process_pbl)
        tidx = dtidx(index_of_temperature,index_of_process_pbl)
        qidx = dtidx(ntqv+100,index_of_process_pbl)
        ! NOTE: The code that was here before was wrong. It replaced the
        ! cumulative value with the instantaneous value.
        do k=1,levs
           k1=levs+1-k
           if(uidx>=1) dtend(:,k,uidx)=dtend(:,k,uidx)+rublten(:,k1)*dt_phs
           if(vidx>=1) dtend(:,k,vidx)=dtend(:,k,vidx)+rvblten(:,k1)*dt_phs
           if(tidx>=1) dtend(:,k,tidx)=dtend(:,k,tidx)+rthblten(:,k1)*exner(:,k1)*dt_phs
           if(qidx>=1) dtend(:,k,qidx)=dtend(:,k,qidx)+rqvblten(:,k1)*dt_phs
        end do
      end if

      if (lprnt1) then

         do i=1,im
            if(tsfc(i).gt.350.)then
               print*,'21tsfc,tsfc1,hflx=',tsfc(i),tsfc1(i),hflx(i)
               print*,'21qsfc,evap=',qsfc(i),evap(i)
            end if
         end do
         tmax=-1.e-5
         tmin=1.e5
         do k=1,levs
           k1=levs+1-k
           do i=1,im
!             t_myj1=t_myj(i,k1)+rthblten(i,k1)*exner(i,k1)*dt_phs
             t_myj1=t_myj(i,k1)+dtdt(i,k)*dt_phs
             if(tmax.lt.t_myj1)then
                tmax=t_myj1
                i_max=i
                k_max=k
                me1=me
             end if
             if(tmin.gt.t_myj1)then
               tmin=t_myj1
               i_min=i
               k_min=k
             end if
           end do
         end do
!          print*,'2after i_min,k_min,i_max,k_max=',i_min,k_min,i_max,k_max
!          print*,'ntsd,me,tmin,tmax=',ntsd,me,tmin,tmax

         if(me.eq.me1.and.tmin.lt.113.6.or.tmax.gt.350.)then
             i=i_max
             print*,'bad bad tmin,tmax=',tmin,tmax,i_min,k_min,i_max,k_max

             do k=1,levs
                 print*,'delt,t_myj=',k,dtdt(i,k)*dt_phs,tgrs(i,k)
             end do

             print*,'ide,levs,ntsd=',ide,lm,ntsd,dt_myj
             print*,'epsl,epsq2,ht,stdh,xland,sice,snowd1=',    &
                   epsl(I),epsq2(I),ht(I),stdh(I),xland(I),sice(I),snowd1(I)
             print*,'phy_f2d_myj=',   &
                  (phy_f2d_myj(i,k),k=1,13)
             print*,'tsk(i),ustar1,z0,pblh_myj,kpbl_myj=',    &
                     tsk(i),ustar1(i),z0(i),pblh_myj(i),kpbl_myj(i)
             print*,'mixht=',mixht(i)
             do k=1,levs
                print*,'u,v,t=',k,u_myj(i,k),v_myj(i,k),   &
                   t_myj(i,k)
             end do
             do k=1,levs
                print*,'q,th,dz_myj=',k,q_myj(i,k),TH_MYJ(i,k),dz_myj(i,k)
             end do
             do k=1,levs
                print*,'del,pmid,pint,=',k,del(i,k), &
                   pmid(i,k),pint(i,K+1)
             end do
             do k=1,levs
                print*,'exner,cw,q2=',k,exner(i,k),cw(i,k),   &
                   q2(i,k)
             end do
             do k=1,levs
                print*,'xcofh,el_myj,dkt=',k,xcofh(i,k),el_myj(i,k),dkt(i,k)
             end do
         end if

      end if     ! lprnt1

      if (lprnt2) then

         tmax=-1.e-5
         tmin=1.e5
         do k=1,levs
           k1=levs+1-k
           do i=1,im
!             t_myj1=t_myj(i,k1)+rthblten(i,k1)*exner(i,k1)*dt_phs
             t_myj1=t_myj(i,k1)+dtdt(i,k)*dt_phs
             if(tmax.lt.t_myj1)then
                tmax=t_myj1
                i_max=i
                k_max=k
                me1=me
             end if
             if(tmin.gt.t_myj1)then
               tmin=t_myj1
               i_min=i
               k_min=k
             end if
           end do
         end do
         print*,'2after me i_min,k_min,i_max,k_max=',me,i_min,k_min,i_max,k_max
         print*,'ntsd,tmin,tmax=',ntsd,tmin,tmax
         print*,'dtdt(i,j)=',dtdt(i_max,k_max)*dt_phs,t_myj(i_max,k_max)

         tmax=-1.e-5
         tmin=1.e5
         do k=1,levs
           k1=levs+1-k
           do i=1,im
!             t_myj1=t_myj(i,k1)+rthblten(i,k1)*exner(i,k1)*dt_phs
             t_myj1=ugrs(i,k)+dudt(i,k)*dt_phs
!             t_myj1=dudt(i,k)*dt_phs
             if(tmax.lt.t_myj1)then
                tmax=t_myj1
                i_max=i
                k_max=k
             end if
             if(tmin.gt.t_myj1)then
               tmin=t_myj1
               i_min=i
               k_min=k
             end if
           end do
         end do
         print*,'3after i_min,k_min,i_max,k_max=',i_min,k_min,i_max,k_max
         print*,'ntsd,me,tmin,tmax=',ntsd,me,tmin,tmax
         print*,'dudt(i,k)=',dudt(i_max,k_max)*dt_phs,ugrs(i_max,k_max)

         if(tmax.gt.200.or.tmin.lt.-200)then
           print*,'bad,bad,bad=',dudt(i_max,k_max)*dt_phs,ugrs(i_max,k_max)
           do k=1,levs
              print*,'k,dudt*dt_phs,ugrs=',k,dudt(i_max,k)*dt_phs,ugrs(i_max,k)
           end do
         end if

         tmax=-1.e-5
         tmin=1.e5
         do k=1,levs
           k1=levs+1-k
           do i=1,im
!             t_myj1=t_myj(i,k1)+rthblten(i,k1)*exner(i,k1)*dt_phs
             t_myj1=vgrs(i,k)+dvdt(i,k)*dt_phs
!             t_myj1=dvdt(i,k)*dt_phs
             if(tmax.lt.t_myj1)then
                tmax=t_myj1
                i_max=i
                k_max=k
             end if
             if(tmin.gt.t_myj1)then
               tmin=t_myj1
               i_min=i
               k_min=k
             end if
           end do
         end do
         print*,'4after i_min,k_min,i_max,k_max=',i_min,k_min,i_max,k_max
         print*,'ntsd,me,tmin,tmax=',ntsd,me,tmin,tmax
         print*,'dvdt(i,k)=',dvdt(i_max,k_max)*dt_phs,vgrs(i_max,k_max)

         tmax=-1.e-5
         tmin=1.e5
         do k=1,levs
           k1=levs+1-k
           do i=1,im
!             t_myj1=q_myj(i,k1)+rthblten(i,k1)*exner(i,k1)*dt_phs
             t_myj1=q_myj(i,k1)+dqdt(i,k,1)*dt_phs
!             t_myj1=dqdt(i,k,1)*dt_phs
             if(tmax.lt.t_myj1)then
                tmax=t_myj1
                i_max=i
                k_max=k
             end if
             if(tmin.gt.t_myj1)then
               tmin=t_myj1
               i_min=i
               k_min=k
             end if
           end do
         end do
         print*,'5after i_min,k_min,i_max,k_max=',i_min,k_min,i_max,k_max
         print*,'ntsd,me,tmin,tmax=',ntsd,me,tmin,tmax
         print*,'dqdt(i,k)=',dqdt(i_max,k_max,1)*dt_phs,qgrs(i_max,k_max,1)

         tmax=-1.e-5
         tmin=1.e5
         do k=1,levs
           k1=levs+1-k
           do i=1,im
!             t_myj1=t_myj(i,k1)+rthblten(i,k1)*exner(i,k1)*dt_phs
             t_myj1=cw(i,k1)+dqdt(i,k,ntcw)*dt_phs
!             t_myj1=dqdt(i,k,ntcw)*dt_phs
             if(tmax.lt.t_myj1)then
                tmax=t_myj1
                i_max=i
                k_max=k
             end if
             if(tmin.gt.t_myj1)then
               tmin=t_myj1
               i_min=i
               k_min=k
             end if
           end do
         end do
         print*,'6after i_min,k_min,i_max,k_max=',i_min,k_min,i_max,k_max
         print*,'ntsd,me,tmin,tmax=',ntsd,me,tmin,tmax
         print*,'dqdt(i,k,ntcw)=',dqdt(i_max,k_max,ntcw)*dt_phs,qgrs(i_max,k_max,ntcw)

      end if    ! lprnt2

!       if (lprnt) then
!          print*
!          print*,"===Finished with myj_bl_driver; output:"
!          print*
!       endif

  END SUBROUTINE myjpbl_wrapper_run

!###=================================================================

END MODULE myjpbl_wrapper
