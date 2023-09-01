!> \file GFS_suite_stateout_update.f90
!!  Contains code to update the state variables due to process-split physics from accumulated tendencies during that phase.
!!  Also, set bounds on the mass-weighted rime factor when using Ferrier-Aligo microphysics.

  module GFS_suite_stateout_update

  contains

!> \section arg_table_GFS_suite_stateout_update_run Argument Table
!! \htmlinclude GFS_suite_stateout_update_run.html
!!
    subroutine GFS_suite_stateout_update_run (im, levs, ntrac, dtp,  &
                     tgrs, ugrs, vgrs, qgrs, dudt, dvdt, dtdt, dqdt, &
                     gt0, gu0, gv0, gq0, ntiw, nqrimef, imp_physics, &
                     imp_physics_fer_hires, epsq, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! Interface variables
      integer,              intent(in )                   :: im
      integer,              intent(in )                   :: levs
      integer,              intent(in )                   :: ntrac
      integer,              intent(in )                   :: imp_physics,imp_physics_fer_hires
      integer,              intent(in )                   :: ntiw, nqrimef
      real(kind=kind_phys), intent(in )                   :: dtp, epsq

      real(kind=kind_phys), intent(in ), dimension(:,:)   :: tgrs, ugrs, vgrs
      real(kind=kind_phys), intent(in ), dimension(:,:,:) :: qgrs
      real(kind=kind_phys), intent(in ), dimension(:,:)   :: dudt, dvdt, dtdt
      real(kind=kind_phys), intent(in ), dimension(:,:,:) :: dqdt
      real(kind=kind_phys), intent(out), dimension(:,:)   :: gt0, gu0, gv0
      real(kind=kind_phys), intent(out), dimension(:,:,:) :: gq0

      character(len=*),     intent(out)                   :: errmsg
      integer,              intent(out)                   :: errflg

      integer                       :: i, k
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      gt0(:,:)   = tgrs(:,:)   + dtdt(:,:)   * dtp
      gu0(:,:)   = ugrs(:,:)   + dudt(:,:)   * dtp
      gv0(:,:)   = vgrs(:,:)   + dvdt(:,:)   * dtp
      gq0(:,:,:) = qgrs(:,:,:) + dqdt(:,:,:) * dtp
      
      if (imp_physics == imp_physics_fer_hires) then
       do k=1,levs
         do i=1,im
           if(gq0(i,k,ntiw) > epsq) then
             gq0(i,k,nqrimef) = max(1., gq0(i,k,nqrimef)/gq0(i,k,ntiw))
           else
             gq0(i,k,nqrimef) = 1.
           end if
         end do
       end do
      end if

    end subroutine GFS_suite_stateout_update_run

  end module GFS_suite_stateout_update