!> \file GFS_stochastics.f90
!! This file contains code previously in GFS_stochastics_driver.

    module GFS_stochastics

      contains

      subroutine GFS_stochastics_init ()
      end subroutine GFS_stochastics_init

      subroutine GFS_stochastics_finalize()
      end subroutine GFS_stochastics_finalize


!>\defgroup gfs_stoch GFS Stochastics Physics Module
!! This module
!> @{
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
      subroutine GFS_stochastics_run (im, km, kdt, do_sppt, use_zmtnblck, do_shum,       &
                                      do_skeb, do_ca,ca_global,ca1,si,vfact_ca,          &
                                      zmtnblck, sppt_wts, skebu_wts, skebv_wts, shum_wts,&
                                      sppt_wts_inv, skebu_wts_inv, skebv_wts_inv,        &
                                      shum_wts_inv, diss_est,                            &
                                      ugrs, vgrs, tgrs, qgrs, gu0, gv0, gt0, gq0, dtdtr, &
                                      rain, rainc, tprcp, totprcp, cnvprcp,              &
                                      totprcpb, cnvprcpb, cplflx,                        &
                                      rain_cpl, snow_cpl, drain_cpl, dsnow_cpl,          &
                                      errmsg, errflg)

         use machine,               only: kind_phys

         implicit none

         integer,                               intent(in)    :: im
         integer,                               intent(in)    :: km
         integer,                               intent(in)    :: kdt
         logical,                               intent(in)    :: do_sppt
         logical,                               intent(in)    :: do_ca
         logical,                               intent(in)    :: ca_global
         logical,                               intent(in)    :: use_zmtnblck
         logical,                               intent(in)    :: do_shum
         logical,                               intent(in)    :: do_skeb
         real(kind_phys), dimension(1:im),      intent(in)    :: zmtnblck
         ! sppt_wts only allocated if do_sppt == .true.
         real(kind_phys), dimension(:,:),       intent(inout) :: sppt_wts
         ! skebu_wts, skebv_wts only allocated if do_skeb == .true.
         real(kind_phys), dimension(:,:),       intent(in)    :: skebu_wts
         real(kind_phys), dimension(:,:),       intent(in)    :: skebv_wts
         ! shum_wts only allocated if do_shum == .true.
         real(kind_phys), dimension(:,:),       intent(in)    :: shum_wts
         ! inverse/flipped weights are always allocated
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: sppt_wts_inv
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: skebu_wts_inv
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: skebv_wts_inv
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: shum_wts_inv
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: diss_est
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: ugrs
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: vgrs
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: tgrs
         real(kind_phys), dimension(1:im,1:km), intent(in)    :: qgrs
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: gu0
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: gv0
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: gt0
         real(kind_phys), dimension(1:im,1:km), intent(inout) :: gq0
         ! dtdtr only allocated if do_sppt == .true.
         real(kind_phys), dimension(:,:),       intent(in)    :: dtdtr
         real(kind_phys), dimension(1:im),      intent(in)    :: rain
         real(kind_phys), dimension(1:im),      intent(in)    :: rainc
         real(kind_phys), dimension(1:im),      intent(inout) :: tprcp
         real(kind_phys), dimension(1:im),      intent(inout) :: totprcp
         real(kind_phys), dimension(1:im),      intent(inout) :: cnvprcp
         real(kind_phys), dimension(1:im),      intent(inout) :: totprcpb
         real(kind_phys), dimension(1:im),      intent(inout) :: cnvprcpb
         logical,                               intent(in)    :: cplflx
         ! rain_cpl, snow_cpl only allocated if cplflx == .true. or cplchm == .true.
         real(kind_phys), dimension(:),         intent(inout) :: rain_cpl
         real(kind_phys), dimension(:),         intent(inout) :: snow_cpl
         ! drain_cpl, dsnow_cpl only allocated if cplflx == .true. or cplchm == .true.
         real(kind_phys), dimension(:),         intent(in)    :: drain_cpl
         real(kind_phys), dimension(:),         intent(in)    :: dsnow_cpl
         real(kind_phys), dimension(1:km),      intent(in)    :: si
         real(kind_phys), dimension(1:km),      intent(inout) :: vfact_ca
         real(kind_phys), dimension(1:im),      intent(in)    :: ca1
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
               sppt_wts_inv(i,k)=sppt_wts(i,k)

               upert = (gu0(i,k) - ugrs(i,k))   * sppt_wts(i,k)
               vpert = (gv0(i,k) - vgrs(i,k))   * sppt_wts(i,k)
               tpert = (gt0(i,k) - tgrs(i,k) - dtdtr(i,k)) * sppt_wts(i,k)
               qpert = (gq0(i,k) - qgrs(i,k)) * sppt_wts(i,k)

               gu0(i,k)  = ugrs(i,k)+upert
               gv0(i,k)  = vgrs(i,k)+vpert

               !negative humidity check
               qnew = qgrs(i,k)+qpert
               if (qnew >= 1.0e-10) then
                  gq0(i,k) = qnew
                  gt0(i,k) = tgrs(i,k) + tpert + dtdtr(i,k)
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

            if (cplflx) then
               rain_cpl(:) = rain_cpl(:) + (sppt_wts(:,15) - 1.0)*drain_cpl(:)
               snow_cpl(:) = snow_cpl(:) + (sppt_wts(:,15) - 1.0)*dsnow_cpl(:)
            endif

         endif

         if (do_ca .and. ca_global) then

          if(kdt == 1)then
            do k=1,km
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
                  tpert = (gt0(i,k)   - tgrs(i,k) - dtdtr(i,k)) * ca(i,k)
                  qpert = (gq0(i,k)   - qgrs(i,k)) * ca(i,k)
                  gu0(i,k)  = ugrs(i,k)+upert
                  gv0(i,k)  = vgrs(i,k)+vpert
                  !negative humidity check                                                                                                                                                                                                                     
                  qnew = qgrs(i,k)+qpert
                  if (qnew >= 1.0e-10) then
                     gq0(i,k) = qnew
                     gt0(i,k)   = tgrs(i,k) + tpert + dtdtr(i,k)
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
            
            if (cplflx) then
               rain_cpl(:) = rain_cpl(:) + (ca(:,15) - 1.0)*drain_cpl(:)
               snow_cpl(:) = snow_cpl(:) + (ca(:,15) - 1.0)*dsnow_cpl(:)
            endif


         endif

         if (do_shum) then
           do k=1,km
             gq0(:,k) = gq0(:,k)*(1.0 + shum_wts(:,k))
             shum_wts_inv(:,k) = shum_wts(:,k)
           end do
         endif
         
         if (do_skeb) then
           do k=1,km
             gu0(:,k) = gu0(:,k)+skebu_wts(:,k)*(diss_est(:,k))
             gv0(:,k) = gv0(:,k)+skebv_wts(:,k)*(diss_est(:,k))
             skebu_wts_inv(:,k) = skebu_wts(:,k)
             skebv_wts_inv(:,k) = skebv_wts(:,k)
           enddo
         endif

      end subroutine GFS_stochastics_run

    end module GFS_stochastics
!> @}
