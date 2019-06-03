!>\file rrtmgp_sw_pre.f90
!! This file contains a subroutine to module_radiation_surface::setalb() to
!! setup surface albedo for SW radiation.
module rrtmgp_sw_pre
  use machine,                   only: kind_phys
  use GFS_typedefs,              only: GFS_control_type,           &
                                       GFS_grid_type,              &
                                       GFS_radtend_type,           &
                                       GFS_sfcprop_type
  use module_radiation_surface,  only: NF_ALBD, setalb
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  implicit none
contains
  
  subroutine rrtmgp_sw_pre_init ()
  end subroutine rrtmgp_sw_pre_init

!> \section arg_table_rrtmgp_sw_pre_run Argument Table
!! | local_name            | standard_name                                               | long_name                                                          | units    | rank |  type                |   kind    | intent | optional |
!! |-----------------------|-------------------------------------------------------------|--------------------------------------------------------------------|----------|------|----------------------|-----------|--------|----------|
!! | Model                 | GFS_control_type_instance                                   | Fortran DDT containing FV3-GFS model control parameters            | DDT      |    0 | GFS_control_type     |           | in     | F        |
!! | Grid                  | GFS_grid_type_instance                                      | Fortran DDT containing FV3-GFS grid and interpolation related data | DDT      |    0 | GFS_grid_type        |           | in     | F        |
!! | Sfcprop               | GFS_sfcprop_type_instance                                   | Fortran DDT containing FV3-GFS surface fields                      | DDT      |    0 | GFS_sfcprop_type     |           | in     | F        |
!! | Radtend               | GFS_radtend_type_instance                                   | Fortran DDT containing FV3-GFS radiation tendencies                | DDT      |    0 | GFS_radtend_type     |           | inout  | F        |
!! | im                    | horizontal_loop_extent                                      | horizontal loop extent                                             | count    |    0 | integer              |           | in     | F        |
!! | nday                  | daytime_points_dimension                                    | daytime points dimension                                           | count    |    0 | integer              |           | out    | F        |
!! | idxday                | daytime_points                                              | daytime points                                                     | index    |    1 | integer              |           | out    | F        |
!! | tsfg                  | surface_ground_temperature_for_radiation                    | surface ground temperature for radiation                           | K        |    1 | real                 | kind_phys | in     | F        |
!! | tsfa                  | surface_air_temperature_for_radiation                       | lowest model layer air temperature for radiation                   | K        |    1 | real                 | kind_phys | in     | F        |
!! | sw_gas_props          | coefficients_for_sw_gas_optics                              | DDT containing spectral information for RRTMGP SW radiation scheme | DDT      |    0 | ty_gas_optics_rrtmgp |           | in     | F        |
!! | sfc_alb_nir_dir       | surface_shortwave_albedo_near_infrared_direct_in_each_band  | surface sw near-infrared direct albedo in each SW band             | frac     |    2 | real                 | kind_phys | out    | F        |
!! | sfc_alb_nir_dif       | surface_shortwave_albedo_near_infrared_diffuse_in_each_band | surface sw near-infrared diffuse albedo in each SW band            | frac     |    2 | real                 | kind_phys | out    | F        |
!! | sfc_alb_uvvis_dir     | surface_shortwave_albedo_uv_visible_direct_in_each_band     | surface sw uv-visible direct albedo in each SW band                | frac     |    2 | real                 | kind_phys | out    | F        |
!! | sfc_alb_uvvis_dif     | surface_shortwave_albedo_uv_visible_diffuse_in_each_band    | surface sw uv-visible diffuse albedo in each SW band               | frac     |    2 | real                 | kind_phys | out    | F        |
!! | alb1d                 | surface_albedo_perturbation                                 | surface albedo perturbation                                        | frac     |    1 | real                 | kind_phys | in     | F        |
!! | errmsg                | ccpp_error_message                                          | error message for error handling in CCPP                           | none     |    0 | character            | len=*     | out    | F        |
!! | errflg                | ccpp_error_flag                                             | error flag for error handling in CCPP                              | flag     |    0 | integer              |           | out    | F        |
!!
  subroutine rrtmgp_sw_pre_run (Model, Grid, Sfcprop, Radtend, im, sw_gas_props, &
       nday, idxday, tsfg, tsfa, sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, &
       sfc_alb_uvvis_dif, alb1d, errmsg, errflg)

    ! Inputs
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props    ! RRTMGP DDT containing spectral information for SW calculation
    type(GFS_control_type),         intent(in)    :: Model
    type(GFS_radtend_type),         intent(inout) :: Radtend
    type(GFS_sfcprop_type),         intent(in)    :: Sfcprop
    type(GFS_grid_type),            intent(in)    :: Grid
    integer,                        intent(in)    :: im
    integer,                        intent(out)   :: nday
    integer, dimension(size(Grid%xlon,1)), intent(out) :: idxday
    real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in)  ::  tsfa, tsfg
    real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in)  :: alb1d

    ! Outputs
    real(kind_phys),dimension(sw_gas_props%get_nband(),IM),intent(out) :: &
         sfc_alb_nir_dir,   & ! Shortwave surface albedo (nIR-direct) 
         sfc_alb_nir_dif,   & ! Shortwave surface albedo (nIR-diffuse)
         sfc_alb_uvvis_dir, & ! Shortwave surface albedo (uvvis-direct)
         sfc_alb_uvvis_dif    ! Shortwave surface albedo (uvvis-diffuse)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    integer :: i, iBand
    real(kind=kind_phys), dimension(size(Grid%xlon,1),NF_ALBD) :: sfcalb
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (Model%lsswr) then
       ! Check for daytime points for SW radiation.
       nday = 0
       idxday = 0
       do i = 1, IM
          if (Radtend%coszen(i) >= 0.0001) then
             nday = nday + 1
             idxday(nday) = i
          endif
       enddo
       
       ! Call module_radiation_surface::setalb() to setup surface albedo.
       call setalb (Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr,&    !  ---  inputs:
                    Sfcprop%snoalb, Sfcprop%zorl, Radtend%coszen,&
                    tsfg, tsfa, Sfcprop%hprim, Sfcprop%alvsf,    &
                    Sfcprop%alnsf, Sfcprop%alvwf, Sfcprop%alnwf, &
                    Sfcprop%facsf, Sfcprop%facwf, Sfcprop%fice,  &
                    Sfcprop%tisfc, IM,                           &
                    alb1d, Model%pertalb,                        &    !  mg, sfc-perts
                    sfcalb)                                           !  ---  outputs
       
       ! Approximate mean surface albedo from vis- and nir-  diffuse values.
       Radtend%sfalb(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
    else
       nday   = 0
       idxday = 0
       sfcalb = 0.0
    endif
      
    ! Spread across all SW bands
    do iBand=1,sw_gas_props%get_nband()
       sfc_alb_nir_dir(iBand,1:IM)   = sfcalb(1:IM,1)
       sfc_alb_nir_dif(iBand,1:IM)   = sfcalb(1:IM,2)
       sfc_alb_uvvis_dir(iBand,1:IM) = sfcalb(1:IM,3)
       sfc_alb_uvvis_dif(iBand,1:IM) = sfcalb(1:IM,4)
    enddo

  end subroutine rrtmgp_sw_pre_run

  subroutine rrtmgp_sw_pre_finalize ()
  end subroutine rrtmgp_sw_pre_finalize

end module rrtmgp_sw_pre
