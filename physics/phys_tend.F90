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
       index_of_process_physics, index_of_process_photochem,  &
       nprocess, nprocess_summed, is_photochem, ntoz, errmsg, errflg)

       ! Interface variables
       logical, intent(in) :: ldiag3d, is_photochem(:)
       real(kind=kind_phys), optional, intent(inout) :: dtend(:,:,:)
       integer, intent(in) :: dtidx(:,:), index_of_process_physics, ntoz, &
         ntracp100, nprocess, nprocess_summed, index_of_process_photochem
       character(len=*), intent(out) :: errmsg
       integer, intent(out)          :: errflg

       integer :: ichem, iphys, itrac
       logical :: all_true(nprocess)

       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if(.not.ldiag3d) then
          return
       endif

       all_true = .true.

       ! Total photochemical tendencies
       itrac=ntoz+100
       ichem = dtidx(itrac,index_of_process_photochem)
       if(ichem>=1) then
          call sum_it(ichem,itrac,is_photochem)
       endif


       do itrac=2,ntracp100
          ! Total physics tendencies
          iphys = dtidx(itrac,index_of_process_physics)
          if(iphys>=1) then
             call sum_it(iphys,itrac,all_true)
          endif
       enddo

     contains
       
       subroutine sum_it(isum,itrac,sum_me)
         implicit none
         integer, intent(in) :: isum ! third index of dtend of summary process
         integer, intent(in) :: itrac ! tracer or state variable being summed
         logical, intent(in) :: sum_me(nprocess) ! false = skip this process
         logical :: first
         integer :: idtend, iprocess

         first=.true.
         do iprocess=1,nprocess
            if(iprocess>nprocess_summed) then
               exit ! Don't sum up the sums.
            else if(.not.sum_me(iprocess)) then
               cycle ! We were asked to skip this one.
            endif
            idtend = dtidx(itrac,iprocess)
            if(idtend>=1) then
               ! This tendency was calculated for this tracer, so
               ! accumulate it into the total tendency.
               if(first) then
                  dtend(:,:,isum) = dtend(:,:,idtend)
                  first=.false.
               else
                  dtend(:,:,isum) = dtend(:,:,isum) + dtend(:,:,idtend)
               endif
            endif
         enddo
         if(first) then
            ! No tendencies were calculated, so sum is 0:
            dtend(:,:,isum) = 0
         endif
       end subroutine sum_it
       
   end subroutine phys_tend_run

end module phys_tend
