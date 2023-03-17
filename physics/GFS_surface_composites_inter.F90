!> \file GFS_surface_composites_inter.F90
!!  Contains code related to generating composites for all GFS surface schemes.

module GFS_surface_composites_inter

   use machine, only: kind_phys

   implicit none

   private

   public GFS_surface_composites_inter_run

contains

!> \section arg_table_GFS_surface_composites_inter_run Argument Table
!! \htmlinclude GFS_surface_composites_inter_run.html
!!
   subroutine GFS_surface_composites_inter_run (im, dry, icy, wet, semis_wat, semis_lnd, semis_ice, &
                                                adjsfcdlw, gabsbdlw_lnd, gabsbdlw_ice, gabsbdlw_wat,&
                                                adjsfcusw, adjsfcdsw, adjsfcnsw, use_lake_model, errmsg, errflg)

      implicit none

      ! Interface variables
      integer,                            intent(in   ) :: im
      logical,              dimension(:), intent(in   ) :: dry, icy
      logical,              dimension(:), intent(inout) :: wet
      real(kind=kind_phys), dimension(:), intent(in   ) :: semis_wat, semis_lnd, semis_ice,  &
                                                           adjsfcdlw, adjsfcdsw, adjsfcnsw
      real(kind=kind_phys), dimension(:), intent(inout) :: gabsbdlw_lnd, gabsbdlw_ice, gabsbdlw_wat
      real(kind=kind_phys), dimension(:), intent(out)   :: adjsfcusw
      integer, dimension(:), intent(in) :: use_lake_model

      ! CCPP error handling
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !  ---  convert lw fluxes for land/ocean/sea-ice models - requires dcyc2t3 to set adjsfcdlw
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
      ! surface upwelling shortwave flux at current time is in adjsfcusw

      !  --- ...  define the downward lw flux absorbed by ground
      do i=1,im
        if (dry(i)) gabsbdlw_lnd(i) = semis_lnd(i) * adjsfcdlw(i)
        if (icy(i)) gabsbdlw_ice(i) = semis_ice(i) * adjsfcdlw(i)
        if (wet(i)) gabsbdlw_wat(i) = semis_wat(i) * adjsfcdlw(i)
        adjsfcusw(i) = adjsfcdsw(i) - adjsfcnsw(i)
      enddo

   end subroutine GFS_surface_composites_inter_run

end module GFS_surface_composites_inter
