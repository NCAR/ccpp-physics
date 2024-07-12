!> \file m_micro_pre.F90
!! This file contains subroutines that prepare data for the Morrison-Gettelman microphysics scheme
!! as part of the GFS physics suite.
      module m_micro_pre

      implicit none

      contains

!! \section arg_table_m_micro_pre_run Argument Table
!! \htmlinclude m_micro_pre_run.html
!!
      subroutine m_micro_pre_run (im, levs, do_shoc, skip_macro, fprcp, mg3_as_mg2, gq0_ice, gq0_water, gq0_rain,  &
        gq0_snow, gq0_graupel, gq0_rain_nc, gq0_snow_nc, gq0_graupel_nc, cld_shoc, cnvc, cnvw, tcr, tcrf, gt0,     &
        qrn, qsnw, qgl, ncpr, ncps, ncgl, cld_frc_MG, clw_water, clw_ice, clcn, errmsg, errflg )

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs, fprcp
      logical, intent(in) :: do_shoc, mg3_as_mg2
      logical, intent(inout) :: skip_macro
      real(kind=kind_phys), intent(in) :: tcr, tcrf

      real(kind=kind_phys), intent(in) ::                               &
          gq0_ice(:,:), gq0_water(:,:), gq0_rain(:,:), gq0_snow(:,:),   &
          gq0_graupel(:,:), gq0_rain_nc(:,:), gq0_snow_nc(:,:),         &
          gq0_graupel_nc(:,:), cnvc(:,:), cnvw(:,:), gt0(:,:)
      real(kind=kind_phys), intent(in), optional :: cld_shoc(:,:)
      real(kind=kind_phys), intent(inout), optional ::                    &
          qrn(:,:), qsnw(:,:), qgl(:,:), ncpr(:,:), ncps(:,:), ncgl(:,:), &
          cld_frc_MG(:,:)

      real(kind=kind_phys), intent(out) :: clw_ice(:,:), clw_water(:,:)

      real(kind=kind_phys), intent(in), optional :: clcn(:,:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k
      real(kind=kind_phys) :: tem

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !       Acheng used clw here for other code to run smoothly and minimum change
      !       to make the code work. However, the nc and clw should be treated
      !       in other procceses too.  August 28/2015; Hope that can be done next
      !       year. I believe this will make the physical interaction more reasonable
      !       Anning 12/5/2015 changed ntcw hold liquid only
      skip_macro = do_shoc
      if (do_shoc) then
        if (fprcp == 0) then
          do k=1,levs
            do i=1,im
              clw_ice(i,k)    = gq0_ice(i,k)
              clw_water(i,k)  = gq0_water(i,k)
              cld_frc_MG(i,k) = cld_shoc(i,k)
            enddo
          enddo
        else if ((abs(fprcp) == 1) .or. mg3_as_mg2) then
          do k=1,levs
            do i=1,im
              clw_ice(i,k)    = gq0_ice(i,k)
              clw_water(i,k)  = gq0_water(i,k)
              qrn(i,k)        = gq0_rain(i,k)
              qsnw(i,k)       = gq0_snow(i,k)
              ncpr(i,k)       = gq0_rain_nc(i,k)
              ncps(i,k)       = gq0_snow_nc(i,k)
              cld_frc_MG(i,k) = cld_shoc(i,k)
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              clw_ice(i,k)    = gq0_ice(i,k)
              clw_water(i,k)  = gq0_water(i,k)
              qrn(i,k)        = gq0_rain(i,k)
              qsnw(i,k)       = gq0_snow(i,k)
              qgl(i,k)        = gq0_graupel(i,k)
              ncpr(i,k)       = gq0_rain_nc(i,k)
              ncps(i,k)       = gq0_snow_nc(i,k)
              ncgl(i,k)       = gq0_graupel_nc(i,k)
              cld_frc_MG(i,k) = cld_shoc(i,k)
            enddo
          enddo
        end if
      else
        if (fprcp == 0 ) then
          do k=1,levs
            do i=1,im
              clw_ice(i,k)   = gq0_ice(i,k)
              clw_water(i,k) = gq0_water(i,k)
            enddo
          enddo
        elseif (abs(fprcp) == 1 .or. mg3_as_mg2) then
          do k=1,levs
            do i=1,im
              clw_ice(i,k)   = gq0_ice(i,k)
              clw_water(i,k) = gq0_water(i,k)
              qrn(i,k)       = gq0_rain(i,k)
              qsnw(i,k)      = gq0_snow(i,k)
              ncpr(i,k)      = gq0_rain_nc(i,k)
              ncps(i,k)      = gq0_snow_nc(i,k)
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              clw_ice(i,k)   = gq0_ice(i,k)
              clw_water(i,k) = gq0_water(i,k)
              qrn(i,k)       = gq0_rain(i,k)
              qsnw(i,k)      = gq0_snow(i,k)
              qgl(i,k)       = gq0_graupel(i,k)
              ncpr(i,k)      = gq0_rain_nc(i,k)
              ncps(i,k)      = gq0_snow_nc(i,k)
              ncgl(i,k)      = gq0_graupel_nc(i,k)
            enddo
          enddo
        endif
      end if

      ! add convective cloud fraction
      do k = 1,levs
        do i = 1,im
          cld_frc_MG(i,k) = min(1.0, cld_frc_MG(i,k) + clcn(i,k))
        enddo
      enddo

      end subroutine m_micro_pre_run

      end module m_micro_pre
