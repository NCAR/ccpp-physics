!> \file GFS_surface_generic.f90
!!  Contains code related to all GFS surface schemes.

      module GFS_surface_generic_pre

      contains

      subroutine GFS_surface_generic_pre_init ()
      end subroutine GFS_surface_generic_pre_init

      subroutine GFS_surface_generic_pre_finalize()
      end subroutine GFS_surface_generic_pre_finalize

!> \section arg_table_GFS_surface_generic_pre_run Argument Table
!! | local_name     | standard_name                                                                | long_name                                                        | units      | rank | type      |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                           | count      |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                                           | number of vertical levels                                        | count      |    0 | integer   |           | in     | F        |
!! | vfrac          | vegetation_area_fraction                                                     | areal fractional cover of green vegetation                       | frac       |    1 | real      | kind_phys | in     | F        |
!! | islmsk         | sea_land_ice_mask                                                            | landmask: sea/land/ice=0/1/2                                     | flag       |    1 | integer   |           | in     | F        |
!! | isot           | soil_type_dataset_choice                                                     | soil type dataset choice                                         | index      |    0 | integer   |           | in     | F        |
!! | ivegsrc        | vegetation_type_dataset_choice                                               | land use dataset choice                                          | index      |    0 | integer   |           | in     | F        |
!! | stype          | soil_type_classification_real                                                | soil type for lsm                                                | index      |    1 | real      | kind_phys | in     | F        |
!! | vtype          | vegetation_type_classification_real                                          | vegetation type for lsm                                          | index      |    1 | real      | kind_phys | in     | F        |
!! | slope          | surface_slope_classification_real                                            | sfc slope type for lsm                                           | index      |    1 | real      | kind_phys | in     | F        |
!! | prsik_1        | dimensionless_exner_function_at_lowest_model_interface                       | dimensionless Exner function at lowest model interface           | none       |    1 | real      | kind_phys | in     | F        |
!! | prslk_1        | dimensionless_exner_function_at_lowest_model_layer                           | dimensionless Exner function at lowest model layer               | none       |    1 | real      | kind_phys | in     | F        |
!! | semis          | surface_longwave_emissivity                                                  | surface lw emissivity in fraction                                | frac       |    1 | real      | kind_phys | in     | F        |
!! | adjsfcdlw      | surface_downwelling_longwave_flux                                            | surface downwelling longwave flux at current time                | W m-2      |    1 | real      | kind_phys | in     | F        |
!! | tsfc           | surface_skin_temperature                                                     | surface skin temperature                                         | K          |    1 | real      | kind_phys | in     | F        |
!! | phil           | geopotential                                                                 | geopotential at model layer centers                              | m2 s-2     |    2 | real      | kind_phys | in     | F        |
!! | con_g          | gravitational_acceleration                                                   | gravitational acceleration                                       | m s-2      |    0 | real      | kind_phys | in     | F        |
!! | sigmaf         | bounded_vegetation_area_fraction                                             | areal fractional cover of green vegetation bounded on the bottom | frac       |    1 | real      | kind_phys | inout  | F        |
!! | soiltyp        | soil_type_classification                                                     | soil type at each grid cell                                      | index      |    1 | integer   |           | inout  | F        |
!! | vegtype        | vegetation_type_classification                                               | vegetation type at each grid cell                                | index      |    1 | integer   |           | inout  | F        |
!! | slopetyp       | surface_slope_classification                                                 | surface slope type at each grid cell                             | index      |    1 | integer   |           | inout  | F        |
!! | work3          | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer      | ratio      |    1 | real      | kind_phys | inout  | F        |
!! | gabsbdlw       | surface_downwelling_longwave_flux_absorbed_by_ground                         | total sky surface downward longwave flux absorbed by the ground  | W m-2      |    1 | real      | kind_phys | inout  | F        |
!! | tsurf          | surface_skin_temperature_after_iteration                                     | surface skin temperature after iteration                         | K          |    1 | real      | kind_phys | inout  | F        |
!! | zlvl           | height_above_ground_at_lowest_model_layer                                    | layer 1 height above ground (not MSL)                            | m          |    1 | real      | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                                           | error message for error handling in CCPP                         | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                              | error flag for error handling in CCPP                            | flag       |    0 | integer   |           | out    | F        |
!!
      subroutine GFS_surface_generic_pre_run (im, levs, vfrac, islmsk, isot, ivegsrc, stype, vtype, slope, &
        prsik_1, prslk_1, semis, adjsfcdlw, tsfc, phil, con_g, sigmaf, soiltyp, vegtype, slopetyp, work3, gabsbdlw, &
        tsurf, zlvl, errmsg, errflg)

        use machine,               only: kind_phys

        implicit none

        integer, intent(in) :: im, levs, isot, ivegsrc
        integer, dimension(im), intent(in) :: islmsk
        integer, dimension(im), intent(inout) :: soiltyp, vegtype, slopetyp

        real(kind=kind_phys), intent(in) :: con_g
        real(kind=kind_phys), dimension(im), intent(in) :: vfrac, stype, vtype, slope, prsik_1, prslk_1, &
          semis, adjsfcdlw, tsfc
        real(kind=kind_phys), dimension(im,levs), intent(in) :: phil

        real(kind=kind_phys), dimension(im), intent(inout) :: sigmaf, work3, gabsbdlw, tsurf, zlvl

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        integer :: i
        real(kind=kind_phys) :: onebg

        onebg  = 1.0/con_g

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        do i=1,im
          sigmaf(i) = max(vfrac(i),0.01 )
          if (islmsk(i) == 2) then
            if (isot == 1) then
              soiltyp(i)  = 16
            else
              soiltyp(i)  = 9
            endif
            if (ivegsrc == 1) then
              vegtype(i)  = 15
            elseif(ivegsrc == 2) then
              vegtype(i)  = 13
            endif
            slopetyp(i) = 9
          else
            soiltyp(i)  = int( stype(i)+0.5 )
            vegtype(i)  = int( vtype(i)+0.5 )
            slopetyp(i) = int( slope(i)+0.5 )    !! clu: slope -> slopetyp
          endif

          work3(i)   = prsik_1(i) / prslk_1(i)
        end do

        !  ---  convert lw fluxes for land/ocean/sea-ice models
        !  note: for sw: adjsfcdsw and adjsfcnsw are zenith angle adjusted downward/net fluxes.
        !        for lw: adjsfcdlw is (sfc temp adjusted) downward fluxe with no emiss effect.
        !                adjsfculw is (sfc temp adjusted) upward fluxe including emiss effect.
        !        one needs to be aware that that the absorbed downward lw flux (used by land/ocean
        !        models as downward flux) is not the same as adjsfcdlw but a value reduced by
        !        the factor of emissivity.  however, the net effects are the same when seeing
        !        it either above the surface interface or below.
        !
        !   - flux above the interface used by atmosphere model:
        !        down: adjsfcdlw;    up: adjsfculw = sfcemis*sigma*T**4 + (1-sfcemis)*adjsfcdlw
        !        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)
        !   - flux below the interface used by lnd/oc/ice models:
        !        down: sfcemis*adjsfcdlw;  up: sfcemis*sigma*T**4
        !        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)

        !  --- ...  define the downward lw flux absorbed by ground
        gabsbdlw(:) = semis(:) * adjsfcdlw(:)

        do i=1,im
          tsurf(i)        = tsfc(i)
          zlvl(i)    = phil(i,1) * onebg
        end do

      end subroutine GFS_surface_generic_pre_run

      end module GFS_surface_generic_pre

      module GFS_surface_generic_post

      contains

      subroutine GFS_surface_generic_post_init ()
      end subroutine GFS_surface_generic_post_init

      subroutine GFS_surface_generic_post_finalize()
      end subroutine GFS_surface_generic_post_finalize

!> \section arg_table_GFS_surface_generic_post_run Argument Table
!! | local_name     | standard_name                                       | long_name                                                            | units      | rank | type                  |    kind   | intent | optional |
!! |----------------|-----------------------------------------------------|----------------------------------------------------------------------|------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                | Fortran DDT containing FV3-GFS model control parameters              | DDT        |    0 | GFS_control_type      |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                   | Fortran DDT containing FV3-GFS grid and interpolation related data   | DDT        |    0 | GFS_grid_type         |           | in     | F        |
!! | ep1d           | surface_upward_potential_latent_heat_flux           | surface upward potential latent heat flux                            | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | gflx           | upward_heat_flux_in_soil                            | upward soil heat flux                                                | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | evbs           | soil_upward_latent_heat_flux                        | soil upward latent heat flux                                         | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | evcw           | canopy_upward_latent_heat_flux                      | canopy upward latent heat flux                                       | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | trans          | transpiration_flux                                  | total plant transpiration rate                                       | kg m-2 s-1 |    1 | real                  | kind_phys | in     | F        |
!! | sbsno          | snow_deposition_sublimation_upward_latent_heat_flux | latent heat flux from snow depo/subl                                 | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | snowc          | surface_snow_area_fraction                          | surface snow area fraction                                           | frac       |    1 | real                  | kind_phys | in     | F        |
!! | snohf          | snow_freezing_rain_upward_latent_heat_flux          | latent heat flux due to snow and frz rain                            | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                   | Fortran DDT containing FV3-GFS fields targeted for diagnostic output | DDT        |    0 | GFS_diag_type         |           | inout  | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                | Fortran DDT containing FV3-GFS surface fields                        | DDT        |    0 | GFS_sfcprop_type      |           | inout  | F        |
!! | errmsg         | ccpp_error_message                                  | error message for error handling in CCPP                             | none       |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                     | error flag for error handling in CCPP                                | flag       |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_surface_generic_post_run (Model, Grid, ep1d, gflx, evbs, evcw, trans, sbsno, snowc, snohf, Diag, &
                                               Sfcprop, errmsg, errflg)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_sfcprop_type, GFS_diag_type

        implicit none

        type(GFS_control_type),           intent(in) :: Model
        type(GFS_grid_type),              intent(in) :: Grid
        type(GFS_sfcprop_type),           intent(inout) :: Sfcprop
        type(GFS_diag_type),              intent(inout) :: Diag

        real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in)  :: ep1d, gflx, evbs, evcw, trans, sbsno, snowc, snohf

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        Diag%epi(:)     = ep1d(:)
        Diag%gfluxi(:)  = gflx(:)

        if (Model%lssav) then
          Diag%gflux(:)   = Diag%gflux(:)  + gflx(:)  * Model%dtf
          Diag%evbsa(:)   = Diag%evbsa(:)  + evbs(:)  * Model%dtf
          Diag%evcwa(:)   = Diag%evcwa(:)  + evcw(:)  * Model%dtf
          Diag%transa(:)  = Diag%transa(:) + trans(:) * Model%dtf
          Diag%sbsnoa(:)  = Diag%sbsnoa(:) + sbsno(:) * Model%dtf
          Diag%snowca(:)  = Diag%snowca(:) + snowc(:) * Model%dtf
          Diag%snohfa(:)  = Diag%snohfa(:) + snohf(:) * Model%dtf
          Diag%ep(:)      = Diag%ep(:)     + ep1d(:)  * Model%dtf

          Diag%tmpmax(:)  = max(Diag%tmpmax(:),Sfcprop%t2m(:))
          Diag%tmpmin(:)  = min(Diag%tmpmin(:),Sfcprop%t2m(:))

          Diag%spfhmax(:) = max(Diag%spfhmax(:),Sfcprop%q2m(:))
          Diag%spfhmin(:) = min(Diag%spfhmin(:),Sfcprop%q2m(:))
        endif

      end subroutine GFS_surface_generic_post_run

      end module GFS_surface_generic_post
