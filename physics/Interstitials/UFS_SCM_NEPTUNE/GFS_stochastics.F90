!> \file GFS_stochastics.F90
!! This file contains code previously in GFS_stochastics_driver.

    module GFS_stochastics

      contains

!>\defgroup gfs_stoch_mod GFS Stochastics Physics Module
!> @{
!! This is the GFS stochastics physics driver module.
!!
!> \section arg_table_GFS_stochastics_init Argument Table
!! \htmlinclude GFS_stochastics_init.html
!!
!>\section gfs_stochyini_general GFS_stochastics_init General Algorithm
!! This is the GFS stochastic physics initialization.
!! -# define vertical tapering for CA global
      subroutine GFS_stochastics_init (si,vfact_ca,km,do_ca,ca_global, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none
      real(kind_phys), dimension(:),         intent(in)    :: si
      real(kind_phys), dimension(:),         intent(inout) :: vfact_ca
      integer,                               intent(in)    :: km
      logical,                               intent(in)    :: do_ca
      logical,                               intent(in)    :: ca_global
      character(len=*),                      intent(out)   :: errmsg
      integer,                               intent(out)   :: errflg
      integer :: k,nz

      errmsg = ''
      errflg = 0
      if (do_ca .and. ca_global) then
         nz=min(km,size(vfact_ca))
         vfact_ca(:)=0.0
         do k=1,nz
            if (si(k) .lt. 0.1 .and. si(k) .gt. 0.025) then
               vfact_ca(k) = (si(k)-0.025)/(0.1-0.025)
            else if (si(k) .lt. 0.025) then
               vfact_ca(k) = 0.0
            else
               vfact_ca(k) = 1.0
            endif
         enddo
         vfact_ca(2)=vfact_ca(3)*0.5
         vfact_ca(1)=0.0
      endif
      end subroutine GFS_stochastics_init

!> \section arg_table_GFS_stochastics_run Argument Table
!! \htmlinclude GFS_stochastics_run.html
!!
!>\section gfs_stochy_general GFS_stochastics_run General Algorithm
!! This is the GFS stochastic physics driver.
!! Routines are called prior to radiation and physics steps to handle:
!! -# sets up various time/date variables
!! -# sets up various triggers
!! -# defines random seed indices for radiation (in a reproducible way)
!! -# interpolates coefficients for prognostic ozone calculation
!! -# performs surface data cycling via the GFS gcycle routine
      subroutine GFS_stochastics_run (im, km, kdt, delt, do_sppt, pert_mp, use_zmtnblck, &
                                      do_shum ,do_skeb, do_ca,ca_global,ca1,vfact_ca,    &
                                      zmtnblck, sppt_wts, skebu_wts, skebv_wts, shum_wts,&
                                      diss_est, ugrs, vgrs, tgrs, qgrs_wv,               &
                                      qgrs_cw, qgrs_rw, qgrs_sw, qgrs_iw, qgrs_gl,       &
                                      gu0, gv0, gt0, gq0_wv, dtdtnp,                     &
                                      gq0_cw, gq0_rw, gq0_sw, gq0_iw, gq0_gl,            &
                                      rain, rainc, tprcp, totprcp, cnvprcp,              &
                                      totprcpb, cnvprcpb, cplflx, cpllnd,                &
                                      rain_cpl, snow_cpl, drain_cpl, dsnow_cpl,          &
                                      ntcw,ntrw,ntsw,ntiw,ntgl,                          &
                                      errmsg, errflg)

         use machine,               only: kind_phys

         implicit none

         integer,                               intent(in)    :: im
         integer,                               intent(in)    :: km
         integer,                               intent(in)    :: kdt
         real(kind_phys),                       intent(in)    :: delt     
         logical,                               intent(in)    :: do_sppt
         logical,                               intent(in)    :: pert_mp
         logical,                               intent(in)    :: do_ca
         logical,                               intent(in)    :: ca_global
         logical,                               intent(in)    :: use_zmtnblck
         logical,                               intent(in)    :: do_shum
         logical,                               intent(in)    :: do_skeb
         real(kind_phys), dimension(:),         intent(in)    :: zmtnblck
         ! sppt_wts only allocated if do_sppt == .true.
         real(kind_phys), dimension(:,:),       intent(inout), optional :: sppt_wts
         ! skebu_wts, skebv_wts only allocated if do_skeb == .true.
         real(kind_phys), dimension(:,:),       intent(in), optional    :: skebu_wts
         real(kind_phys), dimension(:,:),       intent(in), optional    :: skebv_wts
         ! shum_wts only allocated if do_shum == .true.
         real(kind_phys), dimension(:,:),       intent(in), optional    :: shum_wts
         real(kind_phys), dimension(:,:),       intent(in)    :: diss_est
         real(kind_phys), dimension(:,:),       intent(in)    :: ugrs
         real(kind_phys), dimension(:,:),       intent(in)    :: vgrs
         real(kind_phys), dimension(:,:),       intent(in)    :: tgrs
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_wv
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_cw
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_rw
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_sw
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_iw
         real(kind_phys), dimension(:,:),       intent(in)    :: qgrs_gl
         real(kind_phys), dimension(:,:),       intent(inout) :: gu0
         real(kind_phys), dimension(:,:),       intent(inout) :: gv0
         real(kind_phys), dimension(:,:),       intent(inout) :: gt0
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_wv
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_cw
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_rw
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_sw
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_iw
         real(kind_phys), dimension(:,:),       intent(inout) :: gq0_gl
         integer, intent(in) ::      ntcw
         integer, intent(in) ::      ntrw
         integer, intent(in) ::      ntsw
         integer, intent(in) ::      ntiw
         integer, intent(in) ::      ntgl
         real(kind_phys), dimension(:,:),       intent(inout), optional :: dtdtnp
         real(kind_phys), dimension(:),         intent(in)    :: rain
         real(kind_phys), dimension(:),         intent(in)    :: rainc
         real(kind_phys), dimension(:),         intent(inout) :: tprcp
         real(kind_phys), dimension(:),         intent(inout) :: totprcp
         real(kind_phys), dimension(:),         intent(inout) :: cnvprcp
         real(kind_phys), dimension(:),         intent(inout) :: totprcpb
         real(kind_phys), dimension(:),         intent(inout) :: cnvprcpb
         logical,                               intent(in)    :: cplflx
         logical,                               intent(in)    :: cpllnd
         ! rain_cpl only allocated if cplflx == .true. or cplchm == .true. or cpllnd == .true.
         real(kind_phys), dimension(:),         intent(inout), optional :: rain_cpl
         ! snow_cpl only allocated if cplflx == .true. or cplchm == .true.
         real(kind_phys), dimension(:),         intent(inout), optional :: snow_cpl
         ! drain_cpl, dsnow_cpl only allocated if cplflx == .true. or cplchm == .true.
         real(kind_phys), dimension(:),         intent(in), optional    :: drain_cpl
         real(kind_phys), dimension(:),         intent(in), optional    :: dsnow_cpl
         real(kind_phys), dimension(:),         intent(in)    :: vfact_ca
         real(kind_phys), dimension(:),         intent(in), optional :: ca1
         character(len=*),                      intent(out)   :: errmsg
         integer,                               intent(out)   :: errflg

         !--- local variables
         integer :: k, i
         real(kind=kind_phys) :: upert, vpert, tpert, qpert, qnew, sppt_vwt
         real(kind=kind_phys), dimension(1:im,1:km) :: ca

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (do_sppt) then
           do k=1,km
             do i=1,im
               sppt_vwt=1.0
               if (zmtnblck(i).EQ.0.0) then
                  sppt_vwt=1.0
               else 
                  if (k.GT.zmtnblck(i)+2) then
                     sppt_vwt=1.0
                  endif
                  if (k.LE.zmtnblck(i)) then
                     sppt_vwt=0.0
                  endif
                  if (k.EQ.zmtnblck(i)+1) then
                     sppt_vwt=0.333333
                  endif
                  if (k.EQ.zmtnblck(i)+2) then
                     sppt_vwt=0.666667
                  endif
               endif
               if (use_zmtnblck)then
                  sppt_wts(i,k)=(sppt_wts(i,k)-1)*sppt_vwt+1.0
               endif

               upert = (gu0(i,k) - ugrs(i,k))   * sppt_wts(i,k)
               vpert = (gv0(i,k) - vgrs(i,k))   * sppt_wts(i,k)
               tpert = (gt0(i,k) - tgrs(i,k) - (delt*dtdtnp(i,k))) * sppt_wts(i,k)
               qpert = (gq0_wv(i,k) - qgrs_wv(i,k)) * sppt_wts(i,k)

               gu0(i,k)  = ugrs(i,k)+upert
               gv0(i,k)  = vgrs(i,k)+vpert

               !negative humidity check
               qnew = qgrs_wv(i,k)+qpert
               if (qnew >= 1.0e-10) then
                  gq0_wv(i,k) = qnew
                  gt0(i,k) = tgrs(i,k) + tpert + (delt*dtdtnp(i,k))
               endif
               if (pert_mp) then
                  if (ntcw>0) then
                     qpert = (gq0_cw(i,k) - qgrs_cw(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_cw(i,k)+qpert
                     gq0_cw(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_cw(i,k) = 0.0
                     endif
                  endif
                  if (ntrw>0) then
                     qpert = (gq0_rw(i,k) - qgrs_rw(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_rw(i,k)+qpert
                     gq0_rw(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_rw(i,k) = 0.0
                     endif
                  endif
                  if (ntsw>0) then
                     qpert = (gq0_sw(i,k) - qgrs_sw(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_sw(i,k)+qpert
                     gq0_sw(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_sw(i,k) = 0.0
                     endif
                  endif
                  if (ntiw>0) then
                     qpert = (gq0_iw(i,k) - qgrs_iw(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_iw(i,k)+qpert
                     gq0_iw(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_iw(i,k) = 0.0
                     endif
                  endif
                  if (ntgl>0) then
                     qpert = (gq0_gl(i,k) - qgrs_gl(i,k)) * sppt_wts(i,k)
                     qnew = qgrs_gl(i,k)+qpert
                     gq0_gl(i,k) = qnew
                     if (qnew < 0.0) then
                        gq0_gl(i,k) = 0.0
                     endif
                  endif
               endif
             enddo
           enddo

           ! instantaneous precip rate going into land model at the next time step
           tprcp(:) = sppt_wts(:,15)*tprcp(:)
           totprcp(:) = totprcp(:) + (sppt_wts(:,15) - 1 )*rain(:)
           ! acccumulated total and convective preciptiation
           cnvprcp(:) = cnvprcp(:) + (sppt_wts(:,15) - 1 )*rainc(:)
           ! bucket precipitation adjustment due to sppt
           totprcpb(:) = totprcpb(:) + (sppt_wts(:,15) - 1 )*rain(:)
           cnvprcpb(:) = cnvprcpb(:) + (sppt_wts(:,15) - 1 )*rainc(:)

           if (cplflx .or. cpllnd) then
               rain_cpl(:) = rain_cpl(:) + (sppt_wts(:,15) - 1.0)*drain_cpl(:)
           endif
           if (cplflx) then
               snow_cpl(:) = snow_cpl(:) + (sppt_wts(:,15) - 1.0)*dsnow_cpl(:)
           endif
           !zero out radiative heating tendency for next physics step
           dtdtnp(:,:)=0.0

         endif

         if (do_ca .and. ca_global) then

          !if(kdt == 1)then
          !endif
   
            do k = 1,km
               do i = 1,im
                  sppt_vwt=1.0
                  if (zmtnblck(i).EQ.0.0) then
                     sppt_vwt=1.0
                  else
                     if (k.GT.zmtnblck(i)+2) then
                        sppt_vwt=1.0
                     endif
                     if (k.LE.zmtnblck(i)) then
                        sppt_vwt=0.0
                     endif
                     if (k.EQ.zmtnblck(i)+1) then
                        sppt_vwt=0.333333
                     endif
                     if (k.EQ.zmtnblck(i)+2) then
                        sppt_vwt=0.666667
                     endif
                  endif

                  ca(i,k)=((ca1(i)-1.)*sppt_vwt*vfact_ca(k))+1.0

                  upert = (gu0(i,k)   - ugrs(i,k))   * ca(i,k)
                  vpert = (gv0(i,k)   - vgrs(i,k))   * ca(i,k)
                  tpert = (gt0(i,k)   - tgrs(i,k) - (delt*dtdtnp(i,k))) * ca(i,k)
                  qpert = (gq0_wv(i,k)   - qgrs_wv(i,k)) * ca(i,k)
                  gu0(i,k)  = ugrs(i,k)+upert
                  gv0(i,k)  = vgrs(i,k)+vpert
                  !negative humidity check                                                                                                                                                                                                                     
                  qnew = qgrs_wv(i,k)+qpert
                  if (qnew >= 1.0e-10) then
                     gq0_wv(i,k) = qnew
                     gt0(i,k)   = tgrs(i,k) + tpert + (delt*dtdtnp(i,k))
                  endif
                  if (pert_mp) then
                     if (ntcw>0) then
                        qpert = (gq0_cw(i,k) - qgrs_cw(i,k)) * ca(i,k)
                        qnew = qgrs_cw(i,k)+qpert
                        gq0_cw(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_cw(i,k) = 0.0
                        endif
                     endif
                     if (ntrw>0) then
                        qpert = (gq0_rw(i,k) - qgrs_rw(i,k)) * ca(i,k)
                        qnew = qgrs_rw(i,k)+qpert
                        gq0_rw(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_rw(i,k) = 0.0
                        endif
                     endif
                     if (ntsw>0) then
                        qpert = (gq0_sw(i,k) - qgrs_sw(i,k)) * ca(i,k)
                        qnew = qgrs_sw(i,k)+qpert
                        gq0_sw(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_sw(i,k) = 0.0
                        endif
                     endif
                     if (ntiw>0) then
                        qpert = (gq0_iw(i,k) - qgrs_iw(i,k)) * ca(i,k)
                        qnew = qgrs_iw(i,k)+qpert
                        gq0_iw(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_iw(i,k) = 0.0
                        endif
                     endif
                     if (ntgl>0) then
                        qpert = (gq0_gl(i,k) - qgrs_gl(i,k)) * ca(i,k)
                        qnew = qgrs_gl(i,k)+qpert
                        gq0_gl(i,k) = qnew
                        if (qnew < 0.0) then
                           gq0_gl(i,k) = 0.0
                        endif
                     endif
                  endif
               enddo
            enddo
       
            ! instantaneous precip rate going into land model at the next time step                                                                                                                                                                         
            tprcp(:) = ca(:,15)*tprcp(:)
            totprcp(:) = totprcp(:) + (ca(:,15) - 1 )*rain(:)
            ! acccumulated total and convective preciptiation                                                                                                                                                                                               
            cnvprcp(:) = cnvprcp(:)      + (ca(:,15) - 1 )*rainc(:)
            ! bucket precipitation adjustment due to sppt                                                                                                                                                                                                   
            totprcpb(:)      = totprcpb(:)      + (ca(:,15) - 1 )*rain(:)
            cnvprcpb(:)      = cnvprcpb(:)      + (ca(:,15) - 1 )*rainc(:)
            
            if (cplflx .or. cpllnd) then
               rain_cpl(:) = rain_cpl(:) + (ca(:,15) - 1.0)*drain_cpl(:)
            endif
            if (cplflx) then
               snow_cpl(:) = snow_cpl(:) + (ca(:,15) - 1.0)*dsnow_cpl(:)
            endif
            !zero out radiative heating tendency for next physics step
            dtdtnp(:,:)=0.0


         endif

         if (do_shum) then
           do k=1,km
             gq0_wv(:,k) = gq0_wv(:,k)*(1.0 + shum_wts(:,k))
           end do
         endif
         
         if (do_skeb) then
           do k=1,km
             gu0(:,k) = gu0(:,k)+skebu_wts(:,k)*(diss_est(:,k))
             gv0(:,k) = gv0(:,k)+skebv_wts(:,k)*(diss_est(:,k))
           enddo
         endif

      end subroutine GFS_stochastics_run
!> @}
    end module GFS_stochastics
