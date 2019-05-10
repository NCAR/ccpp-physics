!>\file rrtmgp_lw_pre.f90
!! This file contains a call to module_radiation_surface::setemis() to
!! setup surface emissivity for LW radiation.
      module rrtmgp_lw_pre
      contains

!>\defgroup rrtmgp_lw_pre GFS RRTMGP scheme pre
!! @{
!> \section arg_table_rrtmgp_lw_pre_init Argument Table
!!
      subroutine rrtmgp_lw_pre_init ()
      end subroutine rrtmgp_lw_pre_init 

!> \section arg_table_rrtmgp_lw_pre_run Argument Table
!! | local_name       | standard_name                             | long_name                                                          | units    | rank |  type                 |   kind    | intent | optional |
!! |------------------|-------------------------------------------|--------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model            | GFS_control_type_instance                 | Fortran DDT containing FV3-GFS model control parameters            | DDT      |    0 | GFS_control_type      |           | in     | F        |
!! | Grid             | GFS_grid_type_instance                    | Fortran DDT containing FV3-GFS grid and interpolation related data | DDT      |    0 | GFS_grid_type         |           | in     | F        |
!! | Sfcprop          | GFS_sfcprop_type_instance                 | Fortran DDT containing FV3-GFS surface fields                      | DDT      |    0 | GFS_sfcprop_type      |           | in     | F        |
!! | Radtend          | GFS_radtend_type_instance                 | Fortran DDT containing FV3-GFS radiation tendencies                | DDT      |    0 | GFS_radtend_type      |           | inout  | F        |
!! | im               | horizontal_loop_extent                    | horizontal loop extent                                             | count    |    0 | integer               |           | in     | F        |
!! | tsfg             | surface_ground_temperature_for_radiation  | surface ground temperature for radiation                           | K        |    1 | real                  | kind_phys | in     | F        |
!! | tsfa             | surface_air_temperature_for_radiation     | lowest model layer air temperature for radiation                   | K        |    1 | real                  | kind_phys | in     | F        |
!! | errmsg           | ccpp_error_message                        | error message for error handling in CCPP                           | none     |    0 | character             | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                           | error flag for error handling in CCPP                              | flag     |    0 | integer               |           | out    | F        |
!! | kdist_lw         | K_distribution_file_for_RRTMGP_LW_scheme  | DDT containing spectral information for RRTMGP LW radiation scheme | DDT      |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | sfc_emiss_byband | surface_longwave_emissivity_in_each_band  | surface lw emissivity in fraction in each LW band                  | frac     |    2 | real                  | kind_phys | out    | F        |
!!

      subroutine rrtmgp_lw_pre_run (Model, Grid, Sfcprop, Radtend, im, tsfg, tsfa, kdist_lw, sfc_emiss_byband, errmsg, errflg)
    
      use machine,                   only: kind_phys

      use GFS_typedefs,              only: GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type,           &
                                           GFS_sfcprop_type         
      use module_radiation_surface,  only: setemis
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp

      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_sfcprop_type),         intent(in)    :: Sfcprop
      type(GFS_grid_type),            intent(in)    :: Grid
      integer, intent(in)                           :: im
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) ::  tsfa, tsfg
      type(ty_gas_optics_rrtmgp),intent(in) :: &
           kdist_lw        ! DDT containing LW spectral information
      real(kind_phyd),dimension(kdist_lw%get_nband(),im),intent(out) :: sfc_emiss_byband
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      integer :: ij

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lslwr) then
!>  - Call module_radiation_surface::setemis(),to setup surface
!! emissivity for LW radiation.
        call setemis (Grid%xlon, Grid%xlat, Sfcprop%slmsk,        &        !  ---  inputs
                     Sfcprop%snowd, Sfcprop%sncovr, Sfcprop%zorl, &
                     tsfg, tsfa, Sfcprop%hprim, IM,               &
                      Radtend%semis)                              !  ---  outputs
        do ij=1,kdist_lw%get_nband()
           sfc_emiss_byband(ij,:) = Radtend%semis
        enddo
      endif

      end subroutine rrtmgp_lw_pre_run

!> \section arg_table_rrtmgp_lw_pre_finalize Argument Table
!!
       subroutine rrtmgp_lw_pre_finalize ()
       end subroutine rrtmgp_lw_pre_finalize
!! @}
       end module rrtmgp_lw_pre
