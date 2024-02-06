!> \file GFS_SCNV_generic_pre.F90
!!  Contains code related to shallow convective schemes to be run prior to shallow convection for GFS-based physics suites.

      module GFS_SCNV_generic_pre

      contains

!> \section arg_table_GFS_SCNV_generic_pre_run Argument Table
!! \htmlinclude GFS_SCNV_generic_pre_run.html
!!
      subroutine GFS_SCNV_generic_pre_run (im, levs, ldiag3d, qdiag3d, gu0, gv0, gt0, gq0, &
        save_u, save_v, save_t, save_q, ntqv, nsamftrac, flag_for_scnv_generic_tend,       &
        dtidx, index_of_process_scnv, ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,ntgnc,  &
        ntsigma, cscnv, satmedmf, trans_trac, ras, ntrac, clw, errmsg, errflg)

        use machine,               only: kind_phys

        implicit none

        integer, intent(in) :: im, levs, ntqv, nsamftrac, index_of_process_scnv, dtidx(:,:)
        integer, intent(in) :: ntcw,ntiw,ntclamt,ntrw,ntsw,ntrnc,ntsnc,ntgl,ntgnc, ntsigma,ntrac
        logical, intent(in) :: ldiag3d, qdiag3d, flag_for_scnv_generic_tend
        real(kind=kind_phys), dimension(:,:),   intent(in) :: gu0, gv0, gt0
        real(kind=kind_phys), dimension(:,:,:), intent(in) :: gq0
        real(kind=kind_phys), dimension(:,:,:), intent(inout) :: save_q
        real(kind=kind_phys), dimension(:,:),   intent(inout) :: save_u, save_v, save_t
        character(len=*),                 intent(out) :: errmsg
        integer,                          intent(out) :: errflg
        logical, intent(in) :: cscnv, satmedmf, trans_trac, ras
        real(kind=kind_phys), dimension(:,:,:), intent(in) :: clw

        integer :: i, k, n, tracers

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (ldiag3d .and. flag_for_scnv_generic_tend) then
          do k=1,levs
            do i=1,im
              save_u(i,k)   = gu0(i,k)
              save_v(i,k)   = gv0(i,k)
              save_t(i,k)   = gt0(i,k)
            enddo
          enddo
          if (qdiag3d) then
             if (cscnv .or. satmedmf .or. trans_trac .or. ras) then
                tracers = 2
                do n=2,ntrac
                   if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
                        n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
                        n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc   .and. n /= ntsigma) then
                      tracers = tracers + 1
                      if(dtidx(100+n,index_of_process_scnv)>0) then
                         save_q(:,:,n) = clw(:,:,tracers)
                      endif
                   endif
                enddo
             else
                do n=2,ntrac
                   if(dtidx(100+n,index_of_process_scnv)>0) then
                      save_q(:,:,n) = gq0(:,:,n)
                   endif
                enddo
             endif ! end if_ras or cfscnv or samf
             save_q(:,:,ntqv) = gq0(:,:,ntqv)
          endif
        endif

    end subroutine GFS_SCNV_generic_pre_run


    end module GFS_SCNV_generic_pre
