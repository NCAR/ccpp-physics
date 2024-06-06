!> \file GFS_suite_interstitial_4.F90
!!  Contains code to calculate tendencies of tracers due to convective transport, updates tracers after convective transport, and updates cloud condensation nuclei.

  module GFS_suite_interstitial_4

  contains

!> \section arg_table_GFS_suite_interstitial_4_run Argument Table
!! \htmlinclude GFS_suite_interstitial_4_run.html
!!
    subroutine GFS_suite_interstitial_4_run (im, levs, ltaerosol, tracers_total, ntrac, ntcw, ntiw, ntclamt, &
      ntrw, ntsw, ntrnc, ntsnc, ntgl, ntgnc, ntlnc, ntinc, ntccn, nn, imp_physics, imp_physics_gfdl, imp_physics_thompson,  &
      imp_physics_nssl, nssl_invertccn, nssl_ccn_on,                                                  &
      imp_physics_zhao_carr, imp_physics_zhao_carr_pdf, convert_dry_rho, dtf, save_qc, save_qi, con_pi, dtidx, dtend,&
      index_of_process_conv_trans, gq0, clw, prsl, save_tcp, con_rd, con_eps, nssl_cccn, nwfa, spechum, ldiag3d,     &
      qdiag3d, save_lnc, save_inc, ntk, ntke, otsptflag, errmsg, errflg)

      use machine,               only: kind_phys
      use module_mp_thompson_make_number_concentrations, only: make_IceNumber, make_DropletNumber

      implicit none

      ! interface variables

      logical, intent(in)     :: otsptflag(:)! on/off switch for tracer transport by updraft and
      integer,              intent(in   )                   :: im, levs, tracers_total, ntrac, ntcw, ntiw, ntclamt, ntrw, &
        ntsw, ntrnc, ntsnc, ntgl, ntgnc, ntlnc, ntinc, ntccn, nn, imp_physics, imp_physics_gfdl, imp_physics_thompson,    &
        imp_physics_zhao_carr, imp_physics_zhao_carr_pdf, imp_physics_nssl

      logical,                                  intent(in) :: ltaerosol, convert_dry_rho
      logical,                                  intent(in) :: nssl_ccn_on, nssl_invertccn

      real(kind=kind_phys), intent(in   )                   :: con_pi, dtf
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: save_qc
      ! save_qi is not allocated for Zhao-Carr MP
      real(kind=kind_phys), intent(in   ), dimension(:,:)   :: save_qi, save_lnc, save_inc

      ! dtend and dtidx are only allocated if ldiag3d
      logical, intent(in)                                   :: ldiag3d, qdiag3d
      real(kind=kind_phys), dimension(:,:,:), intent(inout), optional :: dtend
      integer,              dimension(:,:),   intent(in)    :: dtidx
      integer,                                intent(in)    :: index_of_process_conv_trans,ntk,ntke

      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: gq0
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: clw
      real(kind=kind_phys), dimension(:,:),   intent(in) :: prsl
      real(kind=kind_phys),                   intent(in) :: con_rd, con_eps, nssl_cccn
      real(kind=kind_phys), dimension(:,:),   intent(in), optional :: nwfa
      real(kind=kind_phys), dimension(:,:),   intent(in) :: save_tcp
      real(kind=kind_phys), dimension(:,:),   intent(in) :: spechum

      character(len=*),     intent(  out)                   :: errmsg
      integer,              intent(  out)                   :: errflg

      ! local variables
      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys
      integer :: i,k,n,tracers,idtend
      real(kind=kind_phys) :: liqm, icem, xccn, xcwmas, xccw, xcimas, qccn

      real(kind=kind_phys) :: rho, orho
      real(kind=kind_phys), dimension(im,levs) :: qv_mp !< kg kg-1 (dry mixing ratio)
      real(kind=kind_phys), dimension(im,levs) :: qc_mp !< kg kg-1 (dry mixing ratio)
      real(kind=kind_phys), dimension(im,levs) :: qi_mp !< kg kg-1 (dry mixing ratio)
      real(kind=kind_phys), dimension(im,levs) :: nc_mp !< kg-1 (dry mixing ratio)
      real(kind=kind_phys), dimension(im,levs) :: ni_mp !< kg-1 (dry mixing ratio)

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! This code was previously in GFS_SCNV_generic_post, but it really belongs
      ! here, because it fixes the convective transportable_tracers mess for Zhao-Carr
      ! and GFDL MP from GFS_suite_interstitial_3. This whole code around clw(:,:,2)
      ! being set to -999 for Zhao-Carr MP (which doesn't have cloud ice) and GFDL-MP
      ! (which does have cloud ice, but for some reason it was decided to code it up
      ! in the same way as for Zhao-Carr, nowadays unnecessary and confusing) needs
      ! to be cleaned up. The convection schemes doing something different internally
      ! based on clw(i,k,2) being -999.0 or not is not a good idea.
      do k=1,levs
        do i=1,im
          if (clw(i,k,2) <= -999.0) clw(i,k,2) = 0.0
        enddo
      enddo

      if(ldiag3d) then
         if(ntk>0 .and. ntk<=size(clw,3)) then
            idtend=dtidx(100+ntke,index_of_process_conv_trans)
            if(idtend>=1) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,ntk)-gq0(:,:,ntk)
            endif
         endif
         if(ntcw>0) then
            if (imp_physics == imp_physics_zhao_carr     .or. &
                 imp_physics == imp_physics_zhao_carr_pdf .or. &
                 imp_physics == imp_physics_gfdl) then
               idtend=dtidx(100+ntcw,index_of_process_conv_trans)
               if(idtend>=1) then
                  dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,1)+clw(:,:,2) - gq0(:,:,ntcw)
               endif
            else if(ntiw>0) then
               idtend=dtidx(100+ntiw,index_of_process_conv_trans)
               if(idtend>=1) then
                  dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,1)-gq0(:,:,ntiw)
               endif
               idtend=dtidx(100+ntcw,index_of_process_conv_trans)
               if(idtend>=1) then
                  dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,2)-gq0(:,:,ntcw)
               endif
            else
               idtend=dtidx(100+ntcw,index_of_process_conv_trans)
               if(idtend>=1) then
                  dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,1)+clw(:,:,2) - gq0(:,:,ntcw)
               endif
            endif
         endif
      endif

!  --- update the tracers due to deep & shallow cumulus convective transport
!           (except for suspended water and ice)

      if (tracers_total > 0) then
        tracers = 2
        do n=2,ntrac
!         if ( n /= ntcw .and. n /= ntiw .and. n /= ntclamt) then
!          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
!               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
!               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc  &
!               .and. &
!             n /= nthl  .and. n /= nthnc .and. n /= ntgv    .and. &
!             n /= nthv .and. n /= ntccn  &
!                                                               ) then
           IF ( otsptflag(n) ) THEN                                                    
              tracers = tracers + 1
            if(n/=ntk .and. n/=ntlnc .and. n/=ntinc .and. n /= ntcw .and. n /= ntiw) then
               idtend=dtidx(100+n,index_of_process_conv_trans)
               if(idtend>=1) then
                  dtend(:,:,idtend) = dtend(:,:,idtend) + clw(:,:,tracers)-gq0(:,:,n)
               endif
            endif
            do k=1,levs
              do i=1,im
                gq0(i,k,n) = clw(i,k,tracers)
              enddo
            enddo
          endif
        enddo
      endif

      if (ntcw > 0) then

!  for microphysics
        if (imp_physics == imp_physics_zhao_carr     .or. &
            imp_physics == imp_physics_zhao_carr_pdf .or. &
            imp_physics == imp_physics_gfdl) then
           gq0(1:im,:,ntcw) = clw(1:im,:,1) + clw(1:im,:,2)

        elseif (ntiw > 0) then
          do k=1,levs
            do i=1,im
              gq0(i,k,ntiw) = clw(i,k,1)                     ! ice
              gq0(i,k,ntcw) = clw(i,k,2)                     ! water
            enddo
          enddo

          if ( imp_physics == imp_physics_nssl ) then
              liqm =  con_pi/6.*1.e3*(18.e-6)**3  ! 4./3.*con_pi*1.e-12
              icem =  con_pi/6.*1.e3*(120.e-6)**3 ! 4./3.*con_pi*3.2768*1.e-14*890.
              qccn = nssl_cccn/1.225 !1.225 is a reference air density and should match what is used in module_mp_nssl_2mom.F90 (rho00)
              do k=1,levs
                do i=1,im
                   ! check number of available ccn
                   IF ( nssl_ccn_on ) THEN
                     IF ( nssl_invertccn ) THEN
                       xccn = qccn - gq0(i,k,ntccn)
                     ELSE
                       xccn = gq0(i,k,ntccn)
                     ENDIF
                   ELSE
                     xccn = Max(0.0, qccn - gq0(i,k,ntlnc))
                   ENDIF
                   
                   IF ( gq0(i,k,ntlnc) > 0.0 .and. save_qc(i,k) > 0.0 ) THEN
                      xcwmas = Max( liqm, clw(i,k,2)/gq0(i,k,ntlnc) )
                   ELSE
                      xcwmas = liqm
                   ENDIF

                   IF ( gq0(i,k,ntinc) > 0.0 .and. save_qi(i,k) > 0.0 ) THEN
                      xcimas = Max( liqm, clw(i,k,1)/gq0(i,k,ntinc) )
                   ELSE
                      xcimas = icem
                   ENDIF
                   
                  IF ( xccn > 0.0 ) THEN
                  xccw = Min( xccn, max(0.0, (clw(i,k,2)-save_qc(i,k))) / xcwmas )
                  gq0(i,k,ntlnc) = gq0(i,k,ntlnc) + xccw 
                  IF ( nssl_ccn_on ) THEN
                     IF ( nssl_invertccn ) THEN
                       ! ccn are activated CCN, so add
                       gq0(i,k,ntccn) = gq0(i,k,ntccn) + xccw
                     ELSE
                       ! ccn are unactivated CCN, so subtract
                       gq0(i,k,ntccn) = gq0(i,k,ntccn) - xccw
                     ENDIF
                  ENDIF
                  ENDIF

                  gq0(i,k,ntinc) = gq0(i,k,ntinc)  &
                           +  max(0.0, (clw(i,k,1)-save_qi(i,k))) / xcimas
                enddo
              enddo
          endif

          if (imp_physics == imp_physics_thompson .and. (ntlnc>0 .or. ntinc>0)) then
            if_convert_dry_rho: if (convert_dry_rho) then
              do k=1,levs
                do i=1,im
                  !> - Convert specific humidity to dry mixing ratio
                  qv_mp(i,k) = spechum(i,k) / (one-spechum(i,k))
                  !> - Density of air in kg m-3 and inverse density
                  rho = con_eps*prsl(i,k) / (con_rd*save_tcp(i,k)*(qv_mp(i,k)+con_eps))
                  orho = one/rho
                  if (ntlnc>0) then
                    !> - Convert moist mixing ratio to dry mixing ratio
                    qc_mp(i,k) = (clw(i,k,2)-save_qc(i,k)) / (one-spechum(i,k))
                    !> - Convert number concentration from moist to dry
                    nc_mp(i,k) = gq0(i,k,ntlnc) / (one-spechum(i,k))
                    nc_mp(i,k) = max(zero, nc_mp(i,k) + make_DropletNumber(qc_mp(i,k) * rho, nwfa(i,k)*rho) * orho)
                    !> - Convert number concentrations from dry to moist
                    gq0(i,k,ntlnc) = nc_mp(i,k) / (one+qv_mp(i,k))
                  endif
                  if (ntinc>0) then
                    !> - Convert moist mixing ratio to dry mixing ratio
                    qi_mp(i,k) = (clw(i,k,1)-save_qi(i,k)) / (one-spechum(i,k))
                    !> - Convert number concentration from moist to dry
                    ni_mp(i,k) = gq0(i,k,ntinc) / (one-spechum(i,k)) 
                    ni_mp(i,k) = max(zero, ni_mp(i,k) + make_IceNumber(qi_mp(i,k) * rho, save_tcp(i,k)) * orho)
                    !> - Convert number concentrations from dry to moist
                    gq0(i,k,ntinc) = ni_mp(i,k) / (one+qv_mp(i,k))
                  endif
                enddo
              enddo
            else
              do k=1,levs
                do i=1,im
                  !> - Density of air in kg m-3 and inverse density
                  rho = con_eps*prsl(i,k) / (con_rd*save_tcp(i,k)*(spechum(i,k)+con_eps))
                  orho = one/rho
                  if (ntlnc>0) then
                    !> - Update cloud water mixing ratio
                    qc_mp(i,k) = (clw(i,k,2)-save_qc(i,k))
                    !> - Update cloud water number concentration
                    gq0(i,k,ntlnc) = max(zero, gq0(i,k,ntlnc) + make_DropletNumber(qc_mp(i,k) * rho, nwfa(i,k)*rho) * orho)
                  endif
                  if (ntinc>0) then
                    !> - Update cloud ice mixing ratio
                    qi_mp(i,k) = (clw(i,k,1)-save_qi(i,k))
                    !> - Update cloud ice number concentration
                    gq0(i,k,ntinc) = max(zero, gq0(i,k,ntinc) + make_IceNumber(qi_mp(i,k) * rho, save_tcp(i,k)) * orho)
                  endif
                enddo
              enddo
            end if if_convert_dry_rho
            if(ldiag3d .and. qdiag3d) then
              idtend = dtidx(100+ntlnc,index_of_process_conv_trans)
              if(idtend>0) then
                dtend(:,:,idtend) = dtend(:,:,idtend) + gq0(:,:,ntlnc) - save_lnc
              endif
              idtend = dtidx(100+ntinc,index_of_process_conv_trans)
              if(idtend>0) then
                dtend(:,:,idtend) = dtend(:,:,idtend) + gq0(:,:,ntinc) - save_inc
              endif
            endif
          endif

        else
          do k=1,levs
            do i=1,im
              gq0(i,k,ntcw) = clw(i,k,1) + clw(i,k,2)
            enddo
          enddo
        endif   ! end if_ntiw

      else
        do k=1,levs
          do i=1,im
            clw(i,k,1) = clw(i,k,1) + clw(i,k,2)
          enddo
        enddo
      endif   ! end if_ntcw

    end subroutine GFS_suite_interstitial_4_run

  end module GFS_suite_interstitial_4
