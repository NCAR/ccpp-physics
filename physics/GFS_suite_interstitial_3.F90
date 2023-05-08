!> \file GFS_suite_interstitial_3.F90
!!  Contains code to setup convectively-transported tracers, calculate critical relative humidity, and save cloud number concentrations

  module GFS_suite_interstitial_3

  contains

!> \section arg_table_GFS_suite_interstitial_3_run Argument Table
!! \htmlinclude GFS_suite_interstitial_3_run.html
!!
    subroutine GFS_suite_interstitial_3_run (otsptflag,                 &
               im, levs, nn, cscnv,imfshalcnv, imfdeepcnv,              &
               imfshalcnv_samf, imfdeepcnv_samf, imfdeepcnv_unified,    &
               imfshalcnv_unified,progsigma,                            &
               first_time_step, restart,                                &
               satmedmf, trans_trac, do_shoc, ltaerosol, ntrac, ntcw,   &
               ntiw, ntclamt, ntrw, ntsw, ntrnc, ntsnc, ntgl, ntgnc,    &
               xlon, xlat, gt0, gq0, sigmain,sigmaout,qmicro,           &
               imp_physics, imp_physics_mg,                             &
               imp_physics_zhao_carr, imp_physics_zhao_carr_pdf,        &
               imp_physics_gfdl, imp_physics_thompson, dtidx, ntlnc,    &
               imp_physics_wsm6, imp_physics_fer_hires, prsi, ntinc,    &
               imp_physics_nssl,                                        &
               prsl, prslk, rhcbot,rhcpbl, rhctop, rhcmax, islmsk,      &
               work1, work2, kpbl, kinver, ras, me, save_lnc, save_inc, &
               ldiag3d, qdiag3d, index_of_process_conv_trans,           &
               clw, rhc, save_qc, save_qi, save_tcp, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! interface variables
      logical, intent(in)     :: otsptflag(:)!  on/off switch for tracer transport (size ntrac)
      integer,              intent(in   )                   :: im, levs, nn, ntrac, ntcw, ntiw, ntclamt, ntrw, ntsw,&
        ntrnc, ntsnc, ntgl, ntgnc, imp_physics, imp_physics_mg, imp_physics_zhao_carr, imp_physics_zhao_carr_pdf,   &
        imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,imp_physics_fer_hires,  &
        imp_physics_nssl, me, index_of_process_conv_trans
      integer,              intent(in   ), dimension(:)     :: islmsk, kpbl, kinver
      logical,              intent(in   )                   :: cscnv, satmedmf, trans_trac, do_shoc, ltaerosol, ras, progsigma
      logical,              intent(in   )                   :: first_time_step, restart
      integer,              intent(in   )                   :: imfshalcnv, imfdeepcnv, imfshalcnv_samf,imfdeepcnv_samf
      integer,              intent(in   )                   :: imfshalcnv_unified,imfdeepcnv_unified
      integer,                                          intent(in) :: ntinc, ntlnc
      logical,                                          intent(in) :: ldiag3d, qdiag3d
      integer,              dimension(:,:),             intent(in) :: dtidx
      real,                 dimension(:,:),            intent(out) :: save_lnc, save_inc

      real(kind=kind_phys), intent(in   )                   :: rhcbot, rhcmax, rhcpbl, rhctop
      real(kind=kind_phys), intent(in   ), dimension(:)     :: work1, work2
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: prsl, prslk
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: prsi
      real(kind=kind_phys), intent(in   ), dimension(:)     :: xlon, xlat
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: gt0
      real(kind=kind_phys), intent(in   ), dimension(:,:,:) :: gq0

      real(kind=kind_phys), intent(inout   ), dimension(:,:)   :: sigmain
      real(kind=kind_phys), intent(inout   ), dimension(:,:)   :: sigmaout,qmicro
      real(kind=kind_phys), intent(inout), dimension(:,:)   :: rhc, save_qc
      ! save_qi is not allocated for Zhao-Carr MP
      real(kind=kind_phys), intent(inout), dimension(:,:)   :: save_qi
      real(kind=kind_phys), intent(inout), dimension(:,:)   :: save_tcp
      real(kind=kind_phys), intent(inout), dimension(:,:,:) :: clw

      character(len=*),     intent(  out)                   :: errmsg
      integer,              intent(  out)                   :: errflg

      ! local variables
      integer :: i,k,n,tracers,kk
      real(kind=kind_phys) :: tem, tem1, tem2
      real(kind=kind_phys), dimension(im) :: tx1, tx2, tx3, tx4

      !real(kind=kind_phys),parameter :: slope_mg = 0.02, slope_upmg = 0.04,  &
      !                   turnrhcrit = 0.900, turnrhcrit_upper = 0.150
      ! in the following inverse of slope_mg and slope_upmg are specified
      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys
      real(kind=kind_phys), parameter :: slope_mg   = 50.0_kind_phys,   &
                                         slope_upmg = 25.0_kind_phys

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! In case of using prognostic updraf area fraction, initialize area fraction here
      ! since progsigma_calc is called from both deep and shallow schemes.
      if(((imfshalcnv == imfshalcnv_samf) .or. (imfdeepcnv == imfdeepcnv_samf) &
          .or. (imfshalcnv == imfshalcnv_unified) .or. (imfdeepcnv == imfdeepcnv_unified)) &
          .and. progsigma)then
         if(first_time_step .and. .not. restart)then
            do k=1,levs
               do i=1,im
                  sigmain(i,k)=0.0
                  sigmaout(i,k)=0.0
                  qmicro(i,k)=0.0
               enddo
            enddo
         endif
         do k=1,levs
            do i=1,im
               sigmaout(i,k)=0.0
            enddo
         enddo
      endif


      if (cscnv .or. satmedmf .or. trans_trac .or. ras) then
        tracers = 2
        do n=2,ntrac
!          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
!               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
!               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc) then
            IF ( otsptflag(n) ) THEN
            tracers = tracers + 1
            do k=1,levs
              do i=1,im
                clw(i,k,tracers) = gq0(i,k,n)
              enddo
            enddo
          endif
        enddo
      endif ! end if_ras or cfscnv or samf

      if (ntcw > 0) then
        if (imp_physics == imp_physics_mg .and. rhcpbl < 0.5_kind_phys) then ! compute rhc for GMAO macro physics cloud pdf
          do i=1,im
            tx1(i) = one / prsi(i,1)
            tx2(i) = one - rhcmax*work1(i)-rhcbot*work2(i)

            kk     = min(kinver(i), max(2,kpbl(i)))
            tx3(i) = prsi(i,kk)*tx1(i)
            tx4(i) = rhcpbl - rhctop*abs(cos(xlat(i)))
          enddo
          do k = 1, levs
            do i = 1, im
              tem  = prsl(i,k) * tx1(i)
              tem1 = min(max((tem-tx3(i))*slope_mg, -20.0_kind_phys), 20.0_kind_phys)
              ! Using rhcpbl and rhctop from the namelist instead of 0.3 and 0.2
              ! and rhcbot represents pbl top critical relative humidity
              tem2 = min(max((tx4(i)-tem)*slope_upmg, -20.0_kind_phys), 20.0_kind_phys) ! Anning
              if (islmsk(i) > 0) then
                tem1 = one / (one+exp(tem1+tem1))
              else
                tem1 = 2.0_kind_phys / (one+exp(tem1+tem1))
              endif
              tem2 = one / (one+exp(tem2))

              rhc(i,k) = min(rhcmax, max(0.7_kind_phys, one-tx2(i)*tem1*tem2))
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              kk = max(10,kpbl(i))
#ifdef SINGLE_PREC
              if (k < kk) then
                tem    = rhcbot - (rhcbot-rhcpbl) * (one-prslk(i,k)) / max(one-prslk(i,kk),1e-7)
              else
                tem    = rhcpbl - (rhcpbl-rhctop) * (prslk(i,kk)-prslk(i,k)) / max(prslk(i,kk),1e-7)
              endif
#else
              if (k < kk) then
                tem    = rhcbot - (rhcbot-rhcpbl) * (one-prslk(i,k)) / (one-prslk(i,kk))
              else
                tem    = rhcpbl - (rhcpbl-rhctop) * (prslk(i,kk)-prslk(i,k)) / prslk(i,kk)
              endif
#endif
              tem      = rhcmax * work1(i) + tem * work2(i)
              rhc(i,k) = max(zero, min(one,tem))
            enddo
          enddo
        endif
      else
        rhc(:,:) = 1.0
      endif

      if (imp_physics == imp_physics_zhao_carr .or. imp_physics == imp_physics_zhao_carr_pdf) then   ! zhao-carr microphysics
        !GF* move to GFS_MP_generic_pre (from gscond/precpd)
        ! do i=1,im
        !   psautco_l(i) = Model%psautco(1)*work1(i) + Model%psautco(2)*work2(i)
        !   prautco_l(i) = Model%prautco(1)*work1(i) + Model%prautco(2)*work2(i)
        ! enddo
        !*GF
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntcw)
          enddo
        enddo
      elseif (imp_physics == imp_physics_gfdl) then
        clw(1:im,:,1) = gq0(1:im,:,ntcw)
      elseif (imp_physics == imp_physics_thompson) then
        do k=1,levs
          do i=1,im
            clw(i,k,1)    = gq0(i,k,ntiw)                    ! ice
            clw(i,k,2)    = gq0(i,k,ntcw)                    ! water
            save_tcp(i,k) = gt0(i,k)
          enddo
        enddo
        if(ltaerosol) then
          save_qi(:,:) = clw(:,:,1)
          save_qc(:,:) = clw(:,:,2)
        else
          save_qi(:,:) = clw(:,:,1)
        endif
      else if (imp_physics == imp_physics_nssl ) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntiw)                    ! cloud ice
            clw(i,k,2) = gq0(i,k,ntcw)                    ! cloud droplets
          enddo
        enddo
          save_qi(:,:) = clw(:,:,1)
          save_qc(:,:) = clw(:,:,2)
      elseif (imp_physics == imp_physics_wsm6 .or. imp_physics == imp_physics_mg .or. imp_physics == imp_physics_fer_hires) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
            clw(i,k,2) = gq0(i,k,ntcw)                    ! water
          enddo
        enddo
      endif

      if(imp_physics == imp_physics_thompson .and. ldiag3d .and. qdiag3d) then
         if(dtidx(100+ntlnc,index_of_process_conv_trans)>0) then
            save_lnc = gq0(:,:,ntlnc)
         endif
         if(dtidx(100+ntinc,index_of_process_conv_trans)>0) then
            save_inc = gq0(:,:,ntinc)
         endif
      endif

    end subroutine GFS_suite_interstitial_3_run

  end module GFS_suite_interstitial_3
