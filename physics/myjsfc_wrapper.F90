!> \file myjsfc_wrapper.F90
!!  Contains all of the code related to running the MYJ surface layer scheme

      MODULE myjsfc_wrapper

      USE machine, only: kfpt => kind_phys, &
                         kind_phys

      contains

      subroutine myjsfc_wrapper_init (do_myjsfc,    &
      & errmsg,errflg)

        logical,          intent(in)  :: do_myjsfc
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

      ! Consistency checks
        if (.not. do_myjsfc) then
          write(errmsg,fmt='(*(a))') 'Logic error: do_myjsfc = .false.'
          errflg = 1
          return
        end if
      end subroutine myjsfc_wrapper_init

!!
!> \brief This scheme (1) performs pre-myjsfc work, (20 runs the myj sfc layer scheme, and (3) performs post-myjsfc work
!! \section arg_table_myjsfc_wrapper_run Argument Table
!! \htmlinclude myjsfc_wrapper_run.html
!!
!###===================================================================
 SUBROUTINE myjsfc_wrapper_run(                    &
     &  restart,                                   &
     &  im,levs,                                   &
     &  kdt,ntrac,ntke,                            &
     &  ntcw,ntiw,ntrw,ntsw,ntgl,                  &
     &  iter,flag_iter,                            &
     &  ugrs, vgrs, tgrs, qgrs,                    &
     &  prsl, prsi, phii,                          &
     &  prsik_1, prslk_1, tsfc, qsfc,              &
     &  phy_myj_qsfc, phy_myj_thz0, phy_myj_qz0,   &
     &  phy_myj_uz0, phy_myj_vz0, phy_myj_z0base,  &
     &  phy_myj_akhs, phy_myj_akms,                &
     &  phy_myj_chkqlm, phy_myj_elflx,             &
     &  phy_myj_a1u, phy_myj_a1t, phy_myj_a1q,     &
     &  pblh, slmsk, zorl, ustar, rib,             &
     &  cm,ch,stress,ffm,ffh,fm10,fh2,             &
     &  landfrac,lakefrac,oceanfrac,fice,          &
     &  z0rl_wat,  z0rl_lnd,  z0rl_ice,            &   ! intent(inout)
     &  ustar_wat, ustar_lnd, ustar_ice,           &   ! intent(inout)
     &  cm_wat,    cm_lnd,    cm_ice,              &   ! intent(inout)
     &  ch_wat,    ch_lnd,    ch_ice,              &   ! intent(inout)
     &  rb_wat,    rb_lnd,    rb_ice,              &   ! intent(inout)
     &  stress_wat,stress_lnd,stress_ice,          &   ! intent(inout)
     &  fm_wat,    fm_lnd,    fm_ice,              &   ! intent(inout)
     &  fh_wat,    fh_lnd,    fh_ice,              &   ! intent(inout)
     &  fm10_wat,  fm10_lnd,  fm10_ice,            &   ! intent(inout)
     &  fh2_wat,   fh2_lnd,   fh2_ice,             &   ! intent(inout)
     &  wind,      con_cp,    con_g,    con_rd,    &
     &  me, lprnt, errmsg, errflg )             ! intent(inout)
!

      use MODULE_SF_JSFC, only: JSFC_INIT,JSFC

!-------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------

!      integer,parameter:: &
!         klog=4 &                   ! logical variables
!        ,kint=4 &                   ! integer variables
!        !,kfpt=4 &                   ! floating point variables
!        ,kfpt=8 &                   ! floating point variables
!        ,kdbl=8                     ! double precision
!
!  ---  constant parameters:
!      real(kind=kind_phys), parameter :: karman  = 0.4

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!For reference
!     real    , parameter :: karman       = 0.4
!     real    , parameter :: g            = 9.81
!     real    , parameter :: r_d          = 287.
!     real    , parameter :: cp           = 7.*r_d/2.
!     real    , parameter :: r_v          = 461.6
!     real    , parameter :: cpv          = 4.*r_v
!     real    , parameter :: rcp          = r_d/cp

!      real, parameter :: g_inv=1/g, cappa=r_d/cp

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

!MYJ-1D
      integer,intent(in) :: im, levs
      integer,intent(in) :: kdt, iter, me
      integer,intent(in) :: ntrac,ntke,ntcw,ntiw,ntrw,ntsw,ntgl
      logical,intent(in) :: restart, lprnt
      real(kind=kind_phys),intent(in) :: con_cp, con_g, con_rd

!MYJ-2D
      logical,dimension(:),intent(in) :: flag_iter
      real(kind=kind_phys),dimension(:),intent(in)      ::   &
     &        prsik_1, prslk_1, tsfc, qsfc, slmsk
      real(kind=kind_phys),dimension(:),intent(inout)   ::   &
     &        phy_myj_qsfc, phy_myj_thz0, phy_myj_qz0,       &
     &        phy_myj_uz0, phy_myj_vz0, phy_myj_z0base,      &
     &        phy_myj_akhs, phy_myj_akms,                    &
     &        phy_myj_chkqlm, phy_myj_elflx,                 &
     &        phy_myj_a1u, phy_myj_a1t, phy_myj_a1q
      real(kind=kind_phys),dimension(:),intent(inout)   ::   &
     &        pblh, zorl, ustar, rib
      real(kind=kind_phys),dimension(:),intent(inout)   ::   &
     &        cm, ch, stress, ffm, ffh, fm10, fh2
      real(kind=kind_phys), dimension(:), intent(inout) ::   &
     &        landfrac, lakefrac, oceanfrac, fice
      real(kind=kind_phys), dimension(:), intent(inout) ::   &
     &                    z0rl_wat,  z0rl_lnd,  z0rl_ice,    &
     &                   ustar_wat, ustar_lnd, ustar_ice,    &
     &                      cm_wat,    cm_lnd,    cm_ice,    &
     &                      ch_wat,    ch_lnd,    ch_ice,    &
     &                      rb_wat,    rb_lnd,    rb_ice,    &
     &                  stress_wat,stress_lnd,stress_ice,    &
     &                      fm_wat,    fm_lnd,    fm_ice,    &
     &                      fh_wat,    fh_lnd,    fh_ice,    &
     &                    fm10_wat,  fm10_lnd,  fm10_ice,    &
     &                     fh2_wat,   fh2_lnd,   fh2_ice,    &
     &                      wind


!MYJ-3D
      real(kind=kind_phys),dimension(:,:),intent(in) ::      &
              phii, prsi
      real(kind=kind_phys),dimension(:,:),intent(in) ::      &
     &        ugrs, vgrs, tgrs, prsl
!MYJ-4D
      real(kind=kind_phys),dimension(:,:,:),intent(in) ::    &
     &       qgrs

!LOCAL
      logical :: lprnt1, lprnt2
      integer :: ntsd, k, k1, i, n, ide, jde, kde

      real(kind=kind_phys) :: g, r_d, g_inv, cappa
      real(kind=kfpt),dimension(levs)       :: epsq2
      real(kind=kfpt),dimension(im)         ::           &
           sfcz,tsk,xland,mavail,rmol,                   &
           ustar1,z0,rib1,sm,pblh_myj
      real(kind=kfpt),dimension(im,13)      ::           &
     &     phy_f2d_myj
      real(kind=kfpt), dimension(im,levs)   ::           &
     &     u_myj, v_myj, t_myj, q_myj, th_myj,           &
     &     cw, dz_myj, pmid, q2, exner
      real(kind=kfpt), dimension(im,levs+1) :: pint
      real(kind=kfpt),dimension(im)         ::           &
     &     cm1,ch1,stress1,ffm1,ffh1,wind1,ffm10,ffh2
!      real(kind=kind_phys), dimension(im,levs,ntrac) :: &
!     &     qgrs_myj

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ntsd = kdt-1

      lprnt1 =.false.
      lprnt2 =.false.

      if (lprnt2) then
         if(me.eq.0)then
           print*,'in myj surface layer wrapper...'
           print*,'ntsd,iter=',ntsd,iter
         end if
      endif

      r_d   = con_rd
      g     = con_g
      g_inv = 1./con_g
      cappa = con_rd/con_cp

      if (ntsd==0.and.iter==1)then
        do i=1,im
           if(flag_iter(i))then
              phy_myj_qsfc(i)   = qgrs(i,1,1)     ! qsfc(:)
              phy_myj_thz0(i)   = tsfc(i)         ! thz0
              phy_myj_qz0(i)    = qgrs(i,1,1)     ! qz0(:)
              phy_myj_uz0(i)    = 0.              ! uz0(:)
              phy_myj_vz0(i)    = 0.              ! vz0(:)
              phy_myj_z0base(i) = zorl(i)*0.01    ! z0base
              phy_myj_akhs(i)   = 0.01            ! akhs(:)
              phy_myj_akms(i)   = 0.01            ! akms(:)
              phy_myj_chkqlm(i) = 0.              ! chkqlm(:)
              phy_myj_elflx(i)  = 0.              ! elflx(:)
              phy_myj_a1u(i)    = 0.              ! a1u
              phy_myj_a1t(i)    = 0.              ! a1t
              phy_myj_a1q(i)    = 0.              ! a1q
           end if
        end do
      end if

!prep MYJ-only variables
      do i=1,im
         sm(i)=1.; if(slmsk(i) > 0.5 ) sm(i)=0.
         xland(i)=sm(i)+1.
         sfcz(i)=phii(i,1)*g_inv
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
              q2(i,k) =qgrs(i,k1,ntke)*2.
            else
              q2(i,k) =0.02
            end if
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

      do k = 1, levs
         k1 = levs-k+1
         do i = 1, im
            dz_myj(i,k) = (phii(i,k1+1)-phii(i,k1)) * g_inv
         enddo
      enddo

      if (lprnt1) then
         if(me==0.and.ntsd.lt.2)then
            k=63
            k1=levs+1-k
            print*,'Qingfu starts MYJSFC'
            print*,'ntsd,iter,me,1=',ntsd,iter,me
            print*,'ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,ntke=',   &
                    ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,ntke
            print*,'im,levs,ntsd=',im,levs,ntsd
            do i=10,40,40
            print*,'Qingfu before MYJ surface kdt,i,k1=',kdt,i,k1
            print*,'sfcz,dz_myj,th_myj,tsfc,qsfc=',sfcz(i),dz_myj(i,k),   &
                    th_myj(i,k),tsfc(i),qsfc(i)
             print*,'sm,z0,xland=',                          &
                   sm(i),z0(i),xland(i)
!            print*,'phy_f2d_myj(i,1:13)=',               &
!                   (phy_f2d_myj(i,n),n=1,13)
            print*,'u_myj,v_myj=',                           &
                   u_myj(i,k),v_myj(i,k)
            print*,'t_myj,q_myj,cw,q2=',                     &
                   t_myj(i,k),q_myj(i,k),cw(i,k),q2(i,k)
            print*,'phii,pint,pmid',                         &
                    phii(i,k1),pint(i,k),pmid(i,k)
            print*,'exner,th_myj=',exner(i,k),th_myj(i,k)
            end do
         end if
      endif

!-----------------------------------------------------------------------
      ide=im+1
      jde=2
      kde=levs+1

      do i = 1, im
         epsq2(i)=0.02
         mavail(i)=1.0
         tsk(i)=tsfc(i)
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
         z0(i)=zorl(i)*0.01
         rmol(i)=0.
         rib1(I)=rib(i)
         pblh_myj(i)=pblh(i)
         ustar1(i)=ustar(i)
         cm1(i)=0.
         ch1(i)=0.
         stress1(i)=0.
         ffm1(i)=0.
         ffh1(i)=0.
         wind1(i)=0.
         ffm10(i)=0.
         ffh2(i)=0.
      end do

      if((ntsd==0.and.iter.eq.1).or.restart)then
         call JSFC_INIT(ustar1,restart                       &
     &                 ,1,ide,1,jde,1,kde                    &
     &                 ,1,im,1,1,1,levs                      &
     &                 ,1,im,1,1,1,levs)
      end if

      call JSFC(flag_iter,iter,me                            &
     &         ,ntsd,epsq2,sfcz,dz_myj                       &
     &         ,pmid,pint,th_myj,t_myj,q_myj,cw              &
     &         ,u_myj,v_myj,q2,tsk                           &
     &         ,phy_f2d_myj(1:im,1),phy_f2d_myj(1:im,2)      &
     &         ,phy_f2d_myj(1:im,3),phy_f2d_myj(1:im,4)      &
     &         ,phy_f2d_myj(1:im,5),xland                    &
     &         ,ustar1,z0,phy_f2d_myj(1:im,6)                &
     &         ,pblh_myj,mavail,rmol                         &
     &         ,phy_f2d_myj(1:im,7),phy_f2d_myj(1:im,8)      &
     &         ,phy_f2d_myj(1:im,9),phy_f2d_myj(1:im,10)     &
     &         ,rib1,cm1,ch1,stress1,ffm1,ffh1,wind1,ffm10,ffh2  &
     &         ,phy_f2d_myj(1:im,11),phy_f2d_myj(1:im,12)    &
     &         ,phy_f2d_myj(1:im,13)                         &
     &         ,1,im,1,1,1,levs                              &
     &         ,1,im,1,1,1,levs                              &
     &         ,1,im,1,1,1,levs, errmsg, errflg)

      do i = 1, im
         if(flag_iter(i))then
            zorl(i) = z0(i)*100.

            phy_myj_qsfc(i)   = phy_f2d_myj(i,1)
            phy_myj_thz0(i)   = phy_f2d_myj(i,2)
            phy_myj_qz0(i)    = phy_f2d_myj(i,3)
            phy_myj_uz0(i)    = phy_f2d_myj(i,4)
            phy_myj_vz0(i)    = phy_f2d_myj(i,5)
            phy_myj_z0base(i) = phy_f2d_myj(i,6)
            phy_myj_akhs(i)   = phy_f2d_myj(i,7)
            phy_myj_akms(i)   = phy_f2d_myj(i,8)
            phy_myj_chkqlm(i) = phy_f2d_myj(i,9)
            phy_myj_elflx(i)  = - phy_f2d_myj(i,10)    ! change flux definition
            phy_myj_a1u(i)    = phy_f2d_myj(i,11)
            phy_myj_a1t(i)    = phy_f2d_myj(i,12)
            phy_myj_a1q(i)    = phy_f2d_myj(i,13)

            rib(I)=rib1(i)
            pblh(I)=pblh_myj(i)
            cm(I)=cm1(i)
            ch(I)=ch1(i)
            stress(I)=stress1(i)
            ffm(I)=ffm1(i)
            ffh(I)=ffh1(i)
            wind(I)=wind1(i)
            fm10(I)=ffm10(i)
            fh2(I)=ffh2(i)
            ustar(i)=ustar1(i)
         end if
      end do

      if (lprnt1) then

        if(me==0.and.ntsd.lt.10)then
           print*,'ntsd,iter,me,2=',ntsd,iter,me
           do i=10,40,40
             if(flag_iter(i))then
               print*,'Qingfu after MYJ surface kdt,i,k1=',kdt,i,k1
               print*,'xland,cm,ch=',xland(i),cm(i),ch(i)
               print*,'ustar,z0,stress=',ustar(i),z0(i),stress(i)
               print*,'ffm,ffh,wind,fm10,fh2=',ffm(i),ffh(i),wind(i),fm10(i),fh2(i)
               print*,'phy_f2d_myj(9,1:13)=',    &
                     (phy_f2d_myj(i,n),n=1,13)
               print*,'u_myj,v_myj=',  &
                     u_myj(i,k),v_myj(i,k)
               print*,'t_myj,q_myj,cw,q2=',  &
                     t_myj(i,k),q_myj(i,k),cw(i,k),q2(i,k)
               print*,'phii,pint,pmid',  &
                     phii(i,k1),pint(i,k),pmid(i,k)
               print*,'exner,th_myj=',exner(i,k),th_myj(i,k)
               print*,'Qingfu finish MYJSFC'
            end if
          end do
        end if

        do k=1,levs
           k1=levs+1-k
           do i=1,im
             if(t_myj(i,k).gt.320..or.t_myj(i,k).lt.150.)then
                print*,'xland,cm,ch=',xland(i),cm(i),ch(i)
                print*,'ustar,z0,stress=',ustar(i),z0(i),stress(i)
                print*,'ffm,ffh,wind,fm10,fh2=',ffm(i),ffh(i),wind(i),fm10(i),fh2(i)
                print*,'phy_f2d_myj(9,1:13)=',    &
                      (phy_f2d_myj(i,n),n=1,13)
                print*,'u_myj,v_myj=',  &
                      u_myj(i,k),v_myj(i,k)
                print*,'t_myj,q_myj,cw,q2=',  &
                      t_myj(i,k),q_myj(i,k),cw(i,k),q2(i,k)
                print*,'phii,pint,pmid',  &
                      phii(i,k1),pint(i,k),pmid(i,k)
                print*,'exner,th_myj=',exner(i,k),th_myj(i,k)
                print*,'Qingfu finish MYJSFC'
             end if
           end do
        end do

      end if

      do i = 1, im
         if(flag_iter(i))then
                z0rl_wat(i) = zorl(i)
                  cm_wat(i) = cm(i)
                  ch_wat(i) = ch(i)
                  rb_wat(i) = rib(i)
              stress_wat(i) = stress(i)
                  fm_wat(i) = ffm(i)
                  fh_wat(i) = ffh(i)
               ustar_wat(i) = ustar(i)
                fm10_wat(i) = fm10(i)
                 fh2_wat(i) = fh2(i)

                z0rl_lnd(i) = zorl(i)
                  cm_lnd(i) = cm(i)
                  ch_lnd(i) = ch(i)
                  rb_lnd(i) = rib(i)
              stress_lnd(i) = stress(i)
                  fm_lnd(i) = ffm(i)
                  fh_lnd(i) = ffh(i)
               ustar_lnd(i) = ustar(i)
                fm10_lnd(i) = fm10(i)
                 fh2_lnd(i) = fh2(i)

                z0rl_ice(i) = zorl(i)
                  cm_ice(i) = cm(i)
                  ch_ice(i) = ch(i)
                  rb_ice(i) = rib(i)
              stress_ice(i) = stress(i)
                  fm_ice(i) = ffm(i)
                  fh_ice(i) = ffh(i)
               ustar_ice(i) = ustar(i)
                fm10_ice(i) = fm10(i)
                 fh2_ice(i) = fh2(i)
            end if
      end do

      if (lprnt2) then
        if(me==0.and.ntsd.lt.10)then
          print*,'ntsd,iter,me,3=',ntsd,iter,me
          do i=10,40,40
             if(flag_iter(i))then
               print*,'Qingfu after MYJ surface kdt,i,k1,3=',kdt,i,k1
               print*,'Qingfu test after MYJ surface kdt,i=',kdt,i,slmsk(i)
               print*,'a1u,a1t,a1q=',(phy_f2d_myj(i,k),k=11,13)
               print*,'zorl,cm,ch,rb,stress=',z0(i),    &
                       cm(i),ch(i),   &
                       rib(i),stress(i)
               print*,'ffmm,ffhh,ustar,fm10,fh2,wind=', ffm(i), &
                       ffh(i),ustar(i),fm10(i),fh2(i),wind(i)
               print*,'cm(i),ch(i)=',    &
                       (0.4/ffm(i))**2,(0.4/ffm(i)*0.4/ffh(i))
             end if
          end do
        endif
      endif


  END SUBROUTINE myjsfc_wrapper_run

!###=================================================================

END MODULE myjsfc_wrapper
