module phys_tend

   use machine, only: kind_phys

   implicit none

   private

   public phys_tend_init, phys_tend_run, phys_tend_finalize

contains

   subroutine phys_tend_init()
   end subroutine phys_tend_init

   subroutine phys_tend_finalize()
   end subroutine phys_tend_finalize

!> \section arg_table_phys_tend_run Argument Table
!! \htmlinclude phys_tend_run.html
!!
   subroutine phys_tend_run(ldiag3d, dtend, dtidx, ntracp100, &
       index_for_cause_physics, index_for_cause_non_physics, &
       ncause, errmsg, errflg)

       ! Interface variables
       logical, intent(in) :: ldiag3d
       real(kind=kind_phys), optional, intent(inout) :: dtend
       integer, intent(in) :: dtidx(:,:), index_for_cause_physics, index_for_cause_non_physics, ntracp100
       character(len=*), intent(out) :: errmsg
       integer, intent(out)          :: errflg

       integer :: itrac, iphys, icause, idtend
       logical :: first

       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if(.not.ldiag3d) then
          return
       endif

       do itrac=2,ntracp100
          first=.true.
          iphys = dtidx(itrac,index_for_cause_physics)
          if(iphys<2) then
             cycle ! No physics tendency requested for this tracer
          endif
          do icause=1,ncause
             if(icause==index_for_cause_physics .or. &
                  icuase==index_for_cause_non_physics) then
                cycle ! Don't sum up the sums.
             endif
             idtend = dtidx(itrac,icause)
             if(idtend>1) then
                ! This tendency was calculated for this tracer, so
                ! accumulate it into the total physics tendency.
                if(first) then
                   dtend(:,:,iphys) = dtend(:,:,idtend)
                   first=.false.
                else
                   dtend(:,:,iphys) = dtend(:,:,iphys) + dtend(:,:,idtend)
                endif
             endif
          enddo
          if(first) then
             ! No physics tendencies were calculated for this tracer,
             ! so total physics tendency is 0.
             dtend(:,:,iphys) = 0
          endif
       enddo

   end subroutine phys_tend_run

end module phys_tend
