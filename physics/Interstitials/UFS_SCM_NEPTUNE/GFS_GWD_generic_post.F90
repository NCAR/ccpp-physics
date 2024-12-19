!> \file GFS_gwd_generic_post.F90 
!! This file contains the CCPP-compliant orographic gravity wave drag post
!! interstitial codes.
module GFS_GWD_generic_post

contains

!> \section arg_table_GFS_GWD_generic_post_run Argument Table
!! \htmlinclude GFS_GWD_generic_post_run.html
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!>  @{
      subroutine GFS_GWD_generic_post_run(lssav, ldiag3d, dtf, dusfcg, dvsfcg, dudt, dvdt, dtdt,          &
      &  dugwd, dvgwd, flag_for_gwd_generic_tend, dtend, dtidx, index_of_temperature, index_of_x_wind,  &
      &  index_of_y_wind, index_of_process_orographic_gwd, errmsg, errflg)

      use machine, only : kind_phys
      implicit none
      
      logical, intent(in) :: lssav, ldiag3d, flag_for_gwd_generic_tend
      
      real(kind=kind_phys), intent(in) :: dusfcg(:), dvsfcg(:)
      real(kind=kind_phys), intent(in) :: dudt(:,:), dvdt(:,:), dtdt(:,:)
      real(kind=kind_phys), intent(in) :: dtf
      
      real(kind=kind_phys), intent(inout) :: dugwd(:), dvgwd(:)

      ! dtend only allocated only if ldiag3d is .true.
      real(kind=kind_phys), intent(inout), optional :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:), index_of_temperature,          &
     &  index_of_x_wind, index_of_y_wind, index_of_process_orographic_gwd
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: idtend

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lssav) then
        dugwd(:) = dugwd(:) + dusfcg(:)*dtf
        dvgwd(:) = dvgwd(:) + dvsfcg(:)*dtf

        if (ldiag3d .and. flag_for_gwd_generic_tend) then
          idtend = dtidx(index_of_temperature, index_of_process_orographic_gwd)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + dtdt*dtf
          endif

          idtend = dtidx(index_of_x_wind, index_of_process_orographic_gwd)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + dudt*dtf
          endif

          idtend = dtidx(index_of_y_wind, index_of_process_orographic_gwd)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + dvdt*dtf
          endif
        endif
      endif

    end subroutine GFS_GWD_generic_post_run
!> @}

end module GFS_GWD_generic_post
