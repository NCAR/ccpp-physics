! ###########################################################################################
!> \file GFS_physics_post.F90
!!
!! This module contains GFS specific calculations (e.g. diagnostics) and suite specific
!! code (e.g Saving fields for subsequent physics timesteps).  For interoperability across a 
!! wide range of hosts, CCPP compliant schemes should avoid including such calculations. This 
!! module/scheme is intended for such "host-specific" computations.
!!
! ###########################################################################################
module GFS_physics_post
  use machine, only : kind_phys, kind_dbl_prec, kind_sngl_prec
  implicit none
  public GFS_physics_post_run
contains

! ###########################################################################################
! SUBROUTINE GFS_physics_post_run
! ###########################################################################################
!! \section arg_table_GFS_physics_post_run Argument Table
!! \htmlinclude GFS_physics_post_run.html
!!
  subroutine GFS_physics_post_run(nCol, nLev, ntoz, ntracp100, nprocess, nprocess_summed,   &
       dtidx, is_photochem, ldiag3d, ip_physics, ip_photochem, ip_prod_loss, ip_ozmix,      &
       ip_temp, ip_overhead_ozone, do3_dt_prd, do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz,       &
       dtend, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: &
         nCol,           & ! Horizontal dimension
         nLev,           & ! Number of vertical layers
         ntoz,           & ! Index for ozone mixing ratio
         ntracp100,      & ! Number of tracers plus 100
         nprocess,       & ! Number of processes that cause changes in state variables 
         nprocess_summed,& ! Number of causes in dtidx per tracer summed for total physics tendency
         ip_physics,     & ! Index for process in diagnostic tendency output
         ip_photochem,   & ! Index for process in diagnostic tendency output
         ip_prod_loss,   & ! Index for process in diagnostic tendency output
         ip_ozmix,       & ! Index for process in diagnostic tendency output
         ip_temp,        & ! Index for process in diagnostic tendency output
         ip_overhead_ozone ! Index for process in diagnostic tendency output    
    integer, intent(in), dimension(:,:) :: &
         dtidx             ! Bookkeeping indices for GFS diagnostic tendencies
    logical, intent(in) :: &
         ldiag3d           ! Flag for 3d diagnostic fields
    logical, intent(in), dimension(:) :: &
         is_photochem      ! Flags for photochemistry processes to sum

    ! Inputs (optional)
    real(kind=kind_phys), intent(in), dimension(:,:), pointer, optional :: &
         do3_dt_prd,     & ! Physics tendency: production and loss effect
         do3_dt_ozmx,    & ! Physics tendency: ozone mixing ratio effect
         do3_dt_temp,    & ! Physics tendency: temperature effect
         do3_dt_ohoz       ! Physics tendency: overhead ozone effect

    ! Outputs
    real(kind=kind_phys), intent(inout), dimension(:,:,:), optional :: &
         dtend             ! Diagnostic tendencies for state variables
    character(len=*), intent(out) :: &
         errmsg            ! CCPP error message
    integer, intent(out) :: &
         errflg            ! CCPP error flag

    ! Locals
    integer :: idtend, ichem, iphys, itrac
    logical :: all_true(nprocess)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if(.not.ldiag3d) then
       return
    endif

    ! #######################################################################################
    !
    ! Ozone physics diagnostics
    !
    ! #######################################################################################
    idtend = dtidx(100+ntoz,ip_prod_loss)
    if (idtend >= 1 .and. associated(do3_dt_prd)) then  
       dtend(:,:,idtend) = dtend(:,:,idtend) + do3_dt_prd
    endif
    !
    idtend = dtidx(100+ntoz,ip_ozmix)
    if (idtend >= 1 .and. associated(do3_dt_ozmx)) then
       dtend(:,:,idtend) = dtend(:,:,idtend) + do3_dt_ozmx
    endif
    !
    idtend = dtidx(100+ntoz,ip_temp)
    if (idtend >= 1 .and. associated(do3_dt_temp)) then
       dtend(:,:,idtend) = dtend(:,:,idtend) + do3_dt_temp
    endif
    !
    idtend = dtidx(100+ntoz,ip_overhead_ozone)
    if (idtend >= 1 .and. associated(do3_dt_ohoz)) then
       dtend(:,:,idtend) = dtend(:,:,idtend) + do3_dt_ohoz
    endif

    ! #######################################################################################
    !
    ! Total (photochemical) tendencies.
    !
    ! #######################################################################################
    itrac = ntoz+100
    ichem = dtidx(itrac, ip_photochem)
    if(ichem >= 1) then
       call sum_it(ichem, itrac, is_photochem)
    endif

    ! #######################################################################################
    !
    ! Total (physics) tendencies
    !
    ! #######################################################################################
    all_true = .true.
    do itrac = 2,ntracp100
       iphys = dtidx(itrac,ip_physics)
       if(iphys >= 1) then
          call sum_it(iphys, itrac, all_true)
       endif
    enddo

  contains

    subroutine sum_it(isum,itrac,sum_me)
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
  end subroutine GFS_physics_post_run
end module GFS_physics_post
