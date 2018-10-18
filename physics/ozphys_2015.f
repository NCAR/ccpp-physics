!> \file ozphys_2015.f
!! This file is ozone sources and sinks.


!> This module contains the CCPP-compliant Ozone 2015 photochemistry scheme.
      module ozphys_2015

      contains

!> \section arg_table_ozphys_2015_init Argument Table
!!
      subroutine ozphys_2015_init()
      end subroutine ozphys_2015_init

! \brief Brief description of the subroutine
!
!> \section arg_table_ozphys_2015_finalize Argument Table
!!
      subroutine ozphys_2015_finalize()
      end subroutine ozphys_2015_finalize


!>\defgroup GFS_ozphys_2015 GFS ozphys_2015 Main
!! @{
!! \brief The operational GFS currently parameterizes ozone production
!and
!! destruction based on monthly mean coefficients (
!! \c ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77) provided by Naval
!! Research Laboratory through CHEM2D chemistry model
!! (McCormack et al. (2006) \cite mccormack_et_al_2006).
!! \section arg_table_ozphys_2015_run Argument Table
!! | local_name     | standard_name                                                            | long_name                                                                  | units   | rank | type      | kind      | intent | optional |
!! |----------------|--------------------------------------------------------------------------|----------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | ix             | horizontal_dimension                                                     | horizontal dimension                                                       | count   |    0 | integer   |           | in     | F        |
!! | im             | horizontal_loop_extent                                                   | horizontal loop extent                                                     | count   |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                                       | number of vertical layers                                                  | count   |    0 | integer   |           | in     | F        |
!! | ko3            | vertical_dimension_of_ozone_forcing_data                                 | number of vertical layers in ozone forcing data                            | count   |    0 | integer   |           | in     | F        |
!! | dt             | time_step_for_physics                                                    | physics time step                                                          | s       |    0 | real      | kind_phys | in     | F        |
!! | oz             | ozone_concentration_updated_by_physics                                   | ozone concentration updated by physics                                     | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | tin            | air_temperature_updated_by_physics                                       | updated air temperature                                                    | K       |    2 | real      | kind_phys | in     | F        |
!! | po3            | natural_log_of_ozone_forcing_data_pressure_levels                        | natural log of ozone forcing data pressure levels                          | log(Pa) |    1 | real      | kind_phys | in     | F        |
!! | prsl           | air_pressure                                                             | mid-layer pressure                                                         | Pa      |    2 | real      | kind_phys | in     | F        |
!! | prdout         | ozone_forcing                                                            | ozone forcing data                                                         | various |    3 | real      | kind_phys | in     | F        |
!! | pl_coeff       | number_of_coefficients_in_ozone_forcing_data                             | number of coefficients in ozone forcing data                               | index   |    0 | integer   |           | in     | F        |
!! | delp           | air_pressure_difference_between_midlayers                                | difference between mid-layer pressures                                     | Pa      |    2 | real      | kind_phys | in     | F        |
!! | ldiag3d        | flag_diagnostics_3D                                                      | flag for calculating 3-D diagnostic fields                                 | flag    |    0 | logical   |           | in     | F        |
!! | ozp1           | cumulative_change_in_ozone_concentration_due_to_production_and_loss_rate | cumulative change in ozone concentration due to production and loss rate   | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | ozp2           | cumulative_change_in_ozone_concentration_due_to_ozone_mixing_ratio       | cumulative change in ozone concentration due to ozone mixing ratio         | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | ozp3           | cumulative_change_in_ozone_concentration_due_to_temperature              | cumulative change in ozone concentration due to temperature                | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | ozp4           | cumulative_change_in_ozone_concentration_due_to_overhead_ozone_column    | cumulative change in ozone concentration due to overhead ozone column      | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | con_g          | gravitational_acceleration                                               | gravitational acceleration                                                 | m s-2   |    0 | real      | kind_phys | in     | F        |
!! | me             | mpi_rank                                                                 | rank of the current MPI task                                               | index   |    0 | integer   |           | in     | F        |
!! | errmsg         | ccpp_error_message                                                       | error message for error handling in CCPP                                   | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                          | error flag for error handling in CCPP                                      | flag    |    0 | integer   |           | out    | F        |
!!
!> \section genal_ozphys_2015 GFS ozphys_2015_run General Algorithm
!! @{
      subroutine ozphys_2015_run (                                      &
     &                        ix, im, levs, ko3, dt, oz, tin, po3,      &
     &                        prsl, prdout, pl_coeff, delp, ldiag3d,    &
     &                        ozp1,ozp2,ozp3,ozp4,con_g,                &
     &                        me, errmsg, errflg)
!
!     this code assumes that both prsl and po3 are from bottom to top
!     as are all other variables
!     This code is specifically for NRL parameterization and
!     climatological T and O3 are in location 5 and 6 of prdout array
! June 2015 - Shrinivas Moorthi
!
      use machine , only : kind_phys
      implicit none
!
      real(kind=kind_phys),intent(in) :: con_g
      real :: gravi
      integer, intent(in) :: im, ix, levs, ko3, pl_coeff,me
      real(kind=kind_phys), intent(in) :: po3(ko3),                     &
     &                                    prsl(ix,levs), tin(ix,levs),  &
     &                                    delp(ix,levs),                &
     &                                    prdout(ix,ko3,pl_coeff), dt
      real(kind=kind_phys), intent(inout) :: ozp1(ix,levs),             &
     &                  ozp2(ix,levs), ozp3(ix,levs),ozp4(ix,levs)
      real(kind=kind_phys), intent(inout) :: oz(ix,levs)


      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer k,kmax,kmin,l,i,j
      logical              ldiag3d, flg(im)
      real(kind=kind_phys) pmax, pmin, tem, temp
      real(kind=kind_phys) wk1(im), wk2(im), wk3(im),prod(im,pl_coeff), &
     &                     ozib(im), colo3(im,levs+1), coloz(im,levs+1),&
     &                     ozi(ix,levs)
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!ccpp: save input oz in ozi
      ozi = oz
      gravi=1.0/con_g

        colo3(:,levs+1) = 0.0
        coloz(:,levs+1) = 0.0
!
      do l=levs,1,-1
        pmin =  1.0e10
        pmax = -1.0e10
!
        do i=1,im
          wk1(i) = log(prsl(i,l))
          pmin   = min(wk1(i), pmin)
          pmax   = max(wk1(i), pmax)
          prod(i,:) = 0.0
        enddo
        kmax = 1
        kmin = 1
        do k=1,ko3-1
          if (pmin < po3(k)) kmax = k
          if (pmax < po3(k)) kmin = k
        enddo
!
        do k=kmin,kmax
          temp = 1.0 / (po3(k) - po3(k+1))
          do i=1,im
            flg(i) = .false.
            if (wk1(i) < po3(k) .and. wk1(i) >= po3(k+1)) then
              flg(i) = .true.
              wk2(i) = (wk1(i) - po3(k+1)) * temp
              wk3(i) = 1.0 - wk2(i)
            endif
          enddo
          do j=1,pl_coeff
            do i=1,im
              if (flg(i)) then
                prod(i,j)  = wk2(i) * prdout(i,k,j)
     &                     + wk3(i) * prdout(i,k+1,j)
              endif
            enddo
          enddo
        enddo
!
        do j=1,pl_coeff
          do i=1,im
            if (wk1(i) < po3(ko3)) then
              prod(i,j) = prdout(i,ko3,j)
            endif
            if (wk1(i) >= po3(1)) then
              prod(i,j) = prdout(i,1,j)
            endif
          enddo
        enddo
        do i=1,im
          colo3(i,l) = colo3(i,l+1) + ozi(i,l)  * delp(i,l)*gravi
          coloz(i,l) = coloz(i,l+1) + prod(i,6) * delp(i,l)*gravi
          prod(i,2)  = min(prod(i,2), 0.0)
        enddo
!       write(1000+me,*) ' colo3=',colo3(1,l),' coloz=',coloz(1,l)
!    &,' l=',l
        do i=1,im
          ozib(i)  = ozi(i,l)            ! no filling
          tem      = prod(i,1) - prod(i,2) * prod(i,6)
     &             + prod(i,3) * (tin(i,l) - prod(i,5))
     &             + prod(i,4) * (colo3(i,l)-coloz(i,l))

!     if (me .eq. 0) print *,'ozphys_2015 tem=',tem,' prod=',prod(i,:)
!    &,' ozib=',ozib(i),' l=',l,' tin=',tin(i,l),'colo3=',colo3(i,l+1)

!ccpp            ozo(i,l) = (ozib(i)  + tem*dt) / (1.0 - prod(i,2)*dt)
          oz(i,l) = (ozib(i)  + tem*dt) / (1.0 - prod(i,2)*dt)
        enddo
        if (ldiag3d) then     !     ozone change diagnostics
          do i=1,im
            ozp1(i,l) = ozp1(i,l) + (prod(i,1)-prod(i,2)*prod(i,6))*dt
!ccpp            ozp(i,l,2) = ozp(i,l,2) + (ozo(i,l) - ozib(i))
            ozp2(i,l) = ozp2(i,l) + (oz(i,l) - ozib(i))
            ozp3(i,l) = ozp3(i,l) + prod(i,3)*(tin(i,l)-prod(i,5))*dt
            ozp4(i,l) = ozp4(i,l) + prod(i,4)
     &                              * (colo3(i,l)-coloz(i,l))*dt
          enddo
        endif
      enddo                                ! vertical loop
!
      return
      end subroutine ozphys_2015_run

!! @}
!! @}

      end module ozphys_2015
