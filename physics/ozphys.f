!> \file ozphys.f
!! This file is ozone sources and sinks.

!> \defgroup GFS_ozn GFS Ozone Sources and Sinks
!! The operational GFS currently parameterizes ozone production and
!! destruction based on monthly mean coefficients provided by Naval
!! Research Laboratory through CHEM2D chemistry model
!! (McCormack et al. 2006 \cite mccormack_et_al_2006).
!! Monthly and zonal mean ozone production rate and ozone destruction
!! rate per unit ozone mixing ratio were provided by NRL based on
!! CHEM2D model.
!! Original version of these terms were provided by NASA/DAO based on
!! NASA 2D Chemistry model - GSM is capable of running both versions
!!
!! \section intra_oz Intraphysics Communication
!! - Routine OZPHYS is called from GBPHYS after call to RAYLEIGH_DAMP
!! @{



      module ozphys_pre

      contains

!! \section arg_table_ozphys_pre_init Argument Table
!!
      subroutine ozphys_pre_init()
      end subroutine ozphys_pre_init


!! \section arg_table_ozphys_pre_run Argument Table
!!
      subroutine ozphys_pre_run()
      end subroutine ozphys_pre_run


!! \section arg_table_ozphys_pre_finalize Argument Table
!!
      subroutine ozphys_pre_finalize()
      end subroutine ozphys_pre_finalize


      end module ozphys_pre




      module ozphys

      contains

!> \ingroup GFS_ozphys
!! \brief Brief description of the subroutine
!!
!! \section arg_table_ozphys_init Argument Table
!!
      subroutine ozphys_init()
      end subroutine ozphys_init

!>
!! \section arg_table_ozphys_run Argument Table
!! | local var name | longname                                          | description                                       | units   | rank | type    | kind      | intent | optional |
!! |----------------|---------------------------------------------------|---------------------------------------------------|---------|------|---------|-----------|--------|----------|
!! | ix             | horizontal_dimension                              | horizontal dimension                              | index   | 0    | integer | default   | in     | F        |
!! | im             | horizontal_loop_extent                            | horizontal loop extent                            | index   | 0    | integer | default   | in     | F        |
!! | levs           | vertical_dimension                                | number of vertical layers                         | index   | 0    | integer | default   | in     | F        |
!! | ko3            | vertical_dimension_of_ozone_forcing_data          | number of vertical layers in ozone forcing data   | index   | 0    | integer | default   | in     | F        |
!! | dt             | time_step_for_physics                             | physics time step                                 | s       | 0    | real    | kind_phys | in     | F        |
!! | ozi            | ozone_concentration_updated_by_physics            | ozone concentration                               | kg kg-1 | 2    | real    | kind_phys | in     | F        |
!! | ozo            | ozone_concentration_updated_by_physics            | ozone concentration updated by physics            | kg kg-1 | 2    | real    | kind_phys | out    | F        |
!! | tin            | air_temperature_updated_by_physics                | updated air temperature                           | K       | 2    | real    | kind_phys | in     | F        |
!! | po3            | natural_log_of_ozone_forcing_data_pressure_levels | natural log of ozone forcing data pressure levels | log(Pa) | 1    | real    | kind_phys | in     | F        |
!! | prsl           | air_pressure                                      | mid-layer pressure                                | Pa      | 2    | real    | kind_phys | in     | F        |
!! | prdout         | ozone_forcing                                     | ozone forcing data                                | various | 3    | real    | kind_phys | in     | F        |
!! | pl_coeff       | number_of_coefficients_in_ozone_forcing_data      | number of coefficients in ozone forcing data      | index   | 0    | integer | default   | in     | F        |
!! | delp           | air_pressure_difference_between_midlayers         | difference between mid-layer pressures            | Pa      | 2    | real    | kind_phys | in     | F        |
!! | ldiag3d        | flag_diagnostics_3D                               | flag for calculating 3-D diagnostic fields        | flag    | 0    | logical | default   | in     | F        |
!! | ozp            | change_in_ozone_concentration                     | change in ozone concentration                     | kg kg-1 | 3    | real    | kind_phys | inout  | F        |
!! | me             | mpi_rank                                          | rank of the current MPI task                      | index   | 0    | integer | default   | in     | F        |
!!
!! \param[in] ix,im     integer, horizontal dimension and num of used pts
!! \param[in] levs      integer, vertical layer dimension
!! \param[in] ko3       integer, number of layers for ozone data
!! \param[in] dt        real, physics time step in seconds
!! \param[in] ozi       real, updated ozone
!! \param[out] ozo      real, updated ozone
!! \param[in] tin       real, updated temperature
!! \param[in] po3       real, (ko3), ozone forcing data level pressure
!!                      (ln(Pa))
!! \param[in] prsl      real, (ix,levs),mean layer pressure
!! \param[in] prdout    real, (ix,ko3,pl_coeff),ozone forcing data
!! \param[in] pl_coeff  integer, number coefficients in ozone forcing
!! \param[in] delp      real, (ix,levs)
!! \param[in] ldiag3d   logical, flag for 3d diagnostic fields
!! \param[out] ozp      real, ozone change due to physics
!! \param[in] me        integer, pe number - used for debug prints
!! \section gen_al General Algorithm
!> @{
      subroutine ozphys_run (                                           &
     &  ix, im, levs, ko3, dt, ozi, ozo, tin, po3,                      &
     &  prsl, prdout, pl_coeff, delp, ldiag3d,                          &
     &  ozp, me)
!
!     this code assumes that both prsl and po3 are from bottom to top
!     as are all other variables
!
      use machine , only : kind_phys
      use physcons, only : grav => con_g
      implicit none
!
      real, parameter :: gravi=1.0/grav
      integer im, ix, levs, ko3, pl_coeff, me
      real(kind=kind_phys) ozi(ix,levs),   ozo(ix,levs), po3(ko3),      &
     &                     prsl(ix,levs),  tin(ix,levs), delp(ix,levs), &
     &                     prdout(ix,ko3,pl_coeff),                     &
     &                     ozp(ix,levs,pl_coeff), dt
!
      integer k,kmax,kmin,l,i,j
      logical              ldiag3d, flg(im)
      real(kind=kind_phys) pmax, pmin, tem, temp
      real(kind=kind_phys) wk1(im), wk2(im), wk3(im), prod(im,pl_coeff),
     &                     ozib(im),  colo3(im,levs+1)
!
      if (pl_coeff > 2) then
        colo3(:,levs+1) = 0.0
        do l=levs,1,-1
          do i=1,im
            colo3(i,l) = colo3(i,l+1) + ozi(i,l) * delp(i,l) * gravi
          enddo
        enddo
      endif
!
      do l=1,levs
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

        if (pl_coeff == 2) then
          do i=1,im
            ozib(i)   = ozi(i,l)           ! no filling
            ozo(i,l)  = (ozib(i) + prod(i,1)*dt) / (1.0 + prod(i,2)*dt)
          enddo
!
          if (ldiag3d) then     !     ozone change diagnostics
            do i=1,im
              ozp(i,l,1) = ozp(i,l,1) + prod(i,1)*dt
              ozp(i,l,2) = ozp(i,l,2) + (ozo(i,l) - ozib(i))
            enddo
          endif
        endif

        if (pl_coeff == 4) then
          do i=1,im
            ozib(i)  = ozi(i,l)            ! no filling
            tem      = prod(i,1) + prod(i,3)*tin(i,l)
     &                           + prod(i,4)*colo3(i,l+1)
!     if (me .eq. 0) print *,'ozphys tem=',tem,' prod=',prod(i,:)
!    &,' ozib=',ozib(i),' l=',l,' tin=',tin(i,l),'colo3=',colo3(i,l+1)
            ozo(i,l) = (ozib(i)  + tem*dt) / (1.0 + prod(i,2)*dt)
          enddo
          if (ldiag3d) then     !     ozone change diagnostics
            do i=1,im
              ozp(i,l,1) = ozp(i,l,1) + prod(i,1)*dt
              ozp(i,l,2) = ozp(i,l,2) + (ozo(i,l) - ozib(i))
              ozp(i,l,3) = ozp(i,l,3) + prod(i,3)*tin(i,l)*dt
              ozp(i,l,4) = ozp(i,l,4) + prod(i,4)*colo3(i,l+1)*dt
            enddo
          endif
        endif

      enddo                                ! vertical loop
!
      return
      end subroutine ozphys_run
!! @}
!> @}

!> \ingroup GFS_ozphys
!! \brief Brief description of the subroutine
!!
!! \section arg_table_ozphys_finalize Argument Table
!!
      subroutine ozphys_finalize()
      end subroutine ozphys_finalize


      end module ozphys




      module ozphys_post

      contains

!! \section arg_table_ozphys_post_init Argument Table
!!
      subroutine ozphys_post_init()
      end subroutine ozphys_post_init


!! \section arg_table_ozphys_post_run Argument Table
!! | local var name | longname                                     | description                                  | units   | rank | type                       | kind      | intent | optional |
!! |----------------|----------------------------------------------|----------------------------------------------|---------|------|----------------------------|-----------|--------|----------|
!! | ix             | horizontal_dimension                         | horizontal dimension                         | index   | 0    | integer                    | default   | in     | F        |
!! | levs           | vertical_dimension                           | number of vertical layers                    | index   | 0    | integer                    | default   | in     | F        |
!! | pl_coeff       | number_of_coefficients_in_ozone_forcing_data | number of coefficients in ozone forcing data | index   | 0    | integer                    | default   | in     | F        |
!! | ozp            | change_in_ozone_concentration                | change in ozone concentration                | kg kg-1 | 3    | real                       | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                            | GFS diagnostics derived data type variable   | DDT     | 0    | GFS_typedefs%GFS_diag_type |           | inout  | F        |
!!
      subroutine ozphys_post_run(ix, levs, pl_coeff, ozp, Diag)

      use GFS_typedefs, only: GFS_diag_type
      use machine,      only: kind_phys

      implicit none

      integer, intent(in) :: ix, levs, pl_coeff
      real(kind=kind_phys), intent(in) :: ozp(ix,levs,pl_coeff)
      type(GFS_diag_type), intent(inout) :: Diag

      Diag%dq3dt(:,:,6) = ozp(:,:,1)
      Diag%dq3dt(:,:,7) = ozp(:,:,2)
      Diag%dq3dt(:,:,8) = ozp(:,:,3)
      Diag%dq3dt(:,:,9) = ozp(:,:,4)

      return

      end subroutine ozphys_post_run


!! \section arg_table_ozphys_post_finalize Argument Table
!!
      subroutine ozphys_post_finalize()
      end subroutine ozphys_post_finalize


      end module ozphys_post

