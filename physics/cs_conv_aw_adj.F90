!>  \file cs_conv_aw_adj.F90
!! This file contains a subroutine to adjusts surface rainrate for conservation for CSAW.  

!>\defgroup mod_cs_conv_aw_adj CSAW adjustment Module
!! This module adjusts surface rainrate for conservation.
!> @{
module cs_conv_aw_adj

   implicit none

   private

   public :: cs_conv_aw_adj_init, cs_conv_aw_adj_run, cs_conv_aw_adj_finalize

   contains

   subroutine cs_conv_aw_adj_init()
   end subroutine cs_conv_aw_adj_init

   subroutine cs_conv_aw_adj_finalize()
   end subroutine cs_conv_aw_adj_finalize

!>\ingroup cs_scheme
!> This subroutine adjusts surface rainrate for conservation.
!> \section arg_table_cs_conv_aw_adj_run Argument Table
!! | local_name      | standard_name                                                 | long_name                                                                        | units   | rank | type      |    kind   | intent | optional |
!! |-----------------|---------------------------------------------------------------|----------------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | im              | horizontal_dimension                                          | horizontal dimension                                                             | count   |    0 | integer   |           | in     | F        |
!! | levs            | vertical_dimension                                            | number of veritcal levels                                                        | count   |    0 | integer   |           | in     | F        |
!! | do_cscnv        | flag_for_Chikira_Sugiyama_deep_convection                     | flag for Chikira-Sugiyama convection                                             | flag    |    0 | logical   |           | in     | F        |
!! | do_aw           | flag_for_Arakawa_Wu_adjustment                                | flag for Arakawa Wu scale-aware adjustment                                       | flag    |    0 | logical   |           | in     | F        |
!! | do_shoc         | flag_for_shoc                                                 | flag for SHOC                                                                    | flag    |    0 | logical   |           | in     | F        |
!! | ntrac           | number_of_tracers                                             | number of tracers                                                                | count   |    0 | integer   |           | in     | F        |
!! | ncld            | number_of_hydrometeors                                        | number of hydrometeors                                                           | count   |    0 | integer   |           | in     | F        |
!! | ntcw            | index_for_liquid_cloud_condensate                             | tracer index for cloud condensate (or liquid water)                              | index   |    0 | integer   |           | in     | F        |
!! | ntclamt         | index_for_cloud_amount                                        | tracer index for cloud amount integer                                            | index   |    0 | integer   |           | in     | F        |
!! | nncl            | number_of_tracers_for_cloud_condensate                        | number of tracers for cloud condensate                                           | count   |    0 | integer   |           | in     | F        |
!! | con_g           | gravitational_acceleration                                    | gravitational acceleration                                                       | m s-2   |    0 | real      | kind_phys | in     | F        |
!! | sigmafrac       | convective_updraft_area_fraction                              | convective updraft area fraction                                                 | frac    |    2 | real      | kind_phys | in     | F        |
!! | gt0             | air_temperature_updated_by_physics                            | temperature updated by physics                                                   | K       |    2 | real      | kind_phys | inout  | F        |
!! | gq0             | tracer_concentration_updated_by_physics                       | tracer concentration updated by physics                                          | kg kg-1 |    3 | real      | kind_phys | inout  | F        |
!! | save_t          | air_temperature_save                                          | air temperature before entering a physics scheme                                 | K       |    2 | real      | kind_phys | in     | F        |
!! | save_q          | tracer_concentration_save                                     | tracer concentration before entering a physics scheme                            | kg kg-1 |    3 | real      | kind_phys | in     | F        |
!! | prsi            | air_pressure_at_interface                                     | air pressure at model layer interfaces                                           | Pa      |    2 | real      | kind_phys | in     | F        |
!! | cldfrac         | cloud_fraction_for_MG                                         | cloud fraction used by Morrison-Gettelman MP                                     | frac    |    2 | real      | kind_phys | inout  | F        |
!! | subcldfrac      | subgrid_scale_cloud_fraction_from_shoc                        | subgrid-scale cloud fraction from the SHOC scheme                                | frac    |    2 | real      | kind_phys | inout  | F        |
!! | prcp            | lwe_thickness_of_explicit_precipitation_amount                | explicit precipitation (rain, ice, snow, graupel, ...) on physics timestep       | m       |    1 | real      | kind_phys | inout  | F        |
!! | imp_physics     | flag_for_microphysics_scheme                                  | choice of microphysics scheme                                                    | flag    |    0 | integer   |           | in     | F        |
!! | imp_physics_mg  | flag_for_morrison_gettelman_microphysics_scheme               | choice of Morrison-Gettelman microphysics scheme                                 | flag    |    0 | integer   |           | in     | F        |
!! | errmsg          | ccpp_error_message                                            | error message for error handling in CCPP                                         | none    |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                               | error flag for error handling in CCPP                                            | flag    |    0 | integer   |           | out    | F        |
!!
!\section gen_cs_conv_aw_adj_run CPT cs_conv_aw_adj_run General Algorithm
   subroutine cs_conv_aw_adj_run(im, levs, do_cscnv, do_aw, do_shoc, &
                ntrac, ncld, ntcw, ntclamt, nncl, con_g, sigmafrac,  &
                gt0, gq0, save_t, save_q, prsi, cldfrac, subcldfrac, &
                prcp, imp_physics, imp_physics_mg, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

! --- interface variables
      integer,                                    intent(in)    :: im, levs
      logical,                                    intent(in)    :: do_cscnv, do_aw, do_shoc
      integer,                                    intent(in)    :: ntrac, ncld, ntcw, ntclamt, nncl
      real(kind_phys),                            intent(in)    :: con_g
      real(kind_phys),  dimension(im,levs),       intent(inout) :: sigmafrac
      real(kind_phys),  dimension(im,levs),       intent(inout) :: gt0
      real(kind_phys),  dimension(im,levs,ntrac), intent(inout) :: gq0
      real(kind_phys),  dimension(im,levs),       intent(in)    :: save_t
      real(kind_phys),  dimension(im,levs,ntrac), intent(in)    :: save_q
      real(kind_phys),  dimension(im,levs+1),     intent(in)    :: prsi
      real(kind_phys),  dimension(im,levs),       intent(inout) :: cldfrac
      real(kind_phys),  dimension(im,levs),       intent(inout) :: subcldfrac
      real(kind_phys),  dimension(im),            intent(inout) :: prcp
      integer,                                    intent(in   ) :: imp_physics, imp_physics_mg
      character(len=*),                           intent(  out) :: errmsg
      integer,                                    intent(  out) :: errflg

! --- local variables
      real(kind=kind_phys), dimension(im) :: temrain1
      real(kind=kind_phys)                :: tem1, tem2
      real(kind=kind_phys)                :: onebg
      integer :: i, k, n

! --- initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !if (do_cscnv .and. do_aw) then

      onebg = 1.0_kind_phys/con_g

!  Arakawa-Wu adjustment of large-scale microphysics tendencies:
!  reduce by factor of (1-sigma)
!  these are microphysics increments. We want to keep (1-sigma) of the increment,
!  we will remove sigma*increment from final values
!         fsigma = 0.  ! don't apply any AW correction, in addition comment next line
!         fsigma = sigmafrac

!  adjust sfc rainrate for conservation
!  vertically integrate reduction of water increments, reduce precip by that amount

      temrain1(:) = 0.0
      do k = 1,levs
        do i = 1,im
          tem1        = sigmafrac(i,k)
          gt0(i,k)    = gt0(i,k) - tem1 * (gt0(i,k)-save_t(i,k))
          tem2        = tem1 * (gq0(i,k,1)-save_q(i,k,1))
          gq0(i,k,1)  = gq0(i,k,1) - tem2
          temrain1(i) = temrain1(i) - (prsi(i,k)-prsi(i,k+1)) * tem2 * onebg
        enddo
      enddo
! add convective clouds if shoc is true and not MG microphysics
      if (do_shoc .and. imp_physics /= imp_physics_mg) then
        do k = 1,levs
          do i = 1,im
            subcldfrac(i,k) = min(1.0, subcldfrac(i,k) + sigmafrac(i,k))
          enddo
        enddo
      endif
      !
      do n=ntcw,ntcw+nncl-1
        do k = 1,levs
          do i = 1,im
            tem1        = sigmafrac(i,k) * (gq0(i,k,n)-save_q(i,k,n))
            gq0(i,k,n)  = gq0(i,k,n) - tem1
            temrain1(i) = temrain1(i) - (prsi(i,k)-prsi(i,k+1)) * tem1 * onebg
          enddo
        enddo
      enddo
      !
      do i = 1,im
        prcp(i) = max(prcp(i) - temrain1(i)*0.001, 0.0_kind_phys)
      enddo

      !endif

      return

   end subroutine cs_conv_aw_adj_run


end module cs_conv_aw_adj
!> @}
