!> \file GFS_DCNV_generic_pre.F90
!!  Contains code related to deep convective schemes to be used within the GFS physics suite.

      module GFS_DCNV_generic_pre

      contains

!> \brief Interstitial scheme called prior to any deep convective scheme to save state variables for calculating tendencies after the deep convective scheme is executed
!! \section arg_table_GFS_DCNV_generic_pre_run Argument Table
!! \htmlinclude GFS_DCNV_generic_pre_run.html
!!
    subroutine GFS_DCNV_generic_pre_run (im, levs, ldiag3d, qdiag3d, do_cnvgwd, cplchm,  &
                                         gu0, gv0, gt0, gq0, nsamftrac, ntqv,            &
                                         save_u, save_v, save_t, save_q, clw,            &
                                         ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,   &
                                         ntgnc, nthl, nthnc, nthv, ntgv,                 &
                                         ntrz, ntgz, nthz, ntsigma,                      &
                                         cscnv, satmedmf, trans_trac, ras, ntrac,        &
                                         dtidx, index_of_process_dcnv, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      integer, intent(in) :: im, levs, nsamftrac, ntqv, index_of_process_dcnv, dtidx(:,:), &
           ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,ntrac,ntgnc,nthl,nthnc,nthv,ntgv,  &
           ntrz, ntgz, nthz, ntsigma
      logical, intent(in) :: ldiag3d, qdiag3d, do_cnvgwd, cplchm
      real(kind=kind_phys), dimension(:,:),   intent(in)    :: gu0
      real(kind=kind_phys), dimension(:,:),   intent(in)    :: gv0
      real(kind=kind_phys), dimension(:,:),   intent(in)    :: gt0
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: gq0
      real(kind=kind_phys), dimension(:,:),   intent(inout) :: save_u
      real(kind=kind_phys), dimension(:,:),   intent(inout) :: save_v
      real(kind=kind_phys), dimension(:,:),   intent(inout) :: save_t
      real(kind=kind_phys), dimension(:,:,:), intent(inout) :: save_q
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg
      logical, intent(in) :: cscnv, satmedmf, trans_trac, ras
      real(kind=kind_phys), parameter :: zero    = 0.0d0
      real(kind=kind_phys), dimension(:,:,:), intent(in) :: clw

      integer :: i, k, n, tracers

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (ldiag3d) then
        do k=1,levs
          do i=1,im
            save_t(i,k) = gt0(i,k)
            save_u(i,k) = gu0(i,k)
            save_v(i,k) = gv0(i,k)
          enddo
        enddo
      elseif (do_cnvgwd) then
        do k=1,levs
          do i=1,im
            save_t(i,k) = gt0(i,k)
          enddo
        enddo
      endif

      if ((ldiag3d.and.qdiag3d) .or. cplchm) then
         if (cscnv .or. satmedmf .or. trans_trac .or. ras) then
            tracers = 2
            do n=2,ntrac
               if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
                    n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
                    n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc   .and. &
                    n /= nthl  .and. n /= nthnc .and. n /= nthv    .and. &
                    n /= ntrz  .and. n /= ntgz  .and. n /= nthz    .and. &
                    n /= ntgv  .and. n/= ntsigma) then
                  tracers = tracers + 1
                  if(dtidx(100+n,index_of_process_dcnv)>0) then
                     save_q(:,:,n) = clw(:,:,tracers)
                  endif
               endif
            enddo
         else
            do n=2,ntrac
               if(dtidx(100+n,index_of_process_dcnv)>0) then
                  save_q(:,:,n) = gq0(:,:,n)
               endif
            enddo
         endif ! end if_ras or cfscnv or samf
         save_q(:,:,ntqv) = gq0(:,:,ntqv)
      endif

    end subroutine GFS_DCNV_generic_pre_run

    end module GFS_DCNV_generic_pre
