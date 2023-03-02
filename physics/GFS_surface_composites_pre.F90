!> \file GFS_surface_composites_pre.F90
!!  Contains code related to generating composites for all GFS surface schemes.

module GFS_surface_composites_pre

   use machine, only: kind_phys

   implicit none

   private

   public GFS_surface_composites_pre_run

   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys, epsln = 1.0e-10_kind_phys

!  real(kind=kind_phys), parameter :: huge      = 9.9692099683868690E36 ! NetCDF float FillValue

contains

!> \section arg_table_GFS_surface_composites_pre_run Argument Table
!! \htmlinclude GFS_surface_composites_pre_run.html
!!
   subroutine GFS_surface_composites_pre_run (im, lkm, frac_grid,                                                         &
                                 flag_cice, cplflx, cplice, cplwav2atm, lsm, lsm_ruc,                                     &
                                 landfrac, lakefrac, lakedepth, oceanfrac, frland,                                        &
                                 dry, icy, lake, use_flake, wet, hice, cice, zorlo, zorll, zorli,                         &
                                 snowd,            snowd_lnd, snowd_ice, tprcp, tprcp_wat,                                &
                                 tprcp_lnd, tprcp_ice, uustar, uustar_wat, uustar_lnd, uustar_ice,                        &
                                 weasd,            weasd_lnd, weasd_ice, ep1d_ice, tsfc, tsfco, tsfcl, tsfc_wat,          &
                                           tisfc, tsurf_wat, tsurf_lnd, tsurf_ice,                                        &
                                 gflx_ice, tgice, islmsk, islmsk_cice, slmsk, qss, qss_wat, qss_lnd, qss_ice,             &
                                 min_lakeice, min_seaice, kdt, huge, errmsg, errflg)

      implicit none

      ! Interface variables
      integer,                             intent(in   ) :: im, lkm, kdt, lsm, lsm_ruc
      logical,                             intent(in   ) :: cplflx, cplice, cplwav2atm, frac_grid
      logical, dimension(:),              intent(inout)  :: flag_cice
      logical,              dimension(:), intent(inout)  :: dry, icy, lake, use_flake, wet
      real(kind=kind_phys), dimension(:), intent(in   )  :: landfrac, lakefrac, lakedepth, oceanfrac
      real(kind=kind_phys), dimension(:), intent(inout)  :: cice, hice
      real(kind=kind_phys), dimension(:), intent(  out)  :: frland
      real(kind=kind_phys), dimension(:), intent(in   )  :: snowd, tprcp, uustar, weasd, qss

      real(kind=kind_phys), dimension(:), intent(inout)  :: tsfc, tsfco, tsfcl, tisfc
      real(kind=kind_phys), dimension(:), intent(inout)  :: snowd_lnd, snowd_ice, tprcp_wat,            &
                    tprcp_lnd, tprcp_ice, tsfc_wat, tsurf_wat,tsurf_lnd, tsurf_ice,                     &
                    uustar_wat, uustar_lnd, uustar_ice, weasd_lnd, weasd_ice,                           &
                    qss_wat, qss_lnd, qss_ice, ep1d_ice, gflx_ice
      real(kind=kind_phys),                intent(in   ) :: tgice
      integer,              dimension(:), intent(inout)  :: islmsk, islmsk_cice
      real(kind=kind_phys), dimension(:), intent(inout)  :: slmsk
      real(kind=kind_phys),               intent(in   )  :: min_lakeice, min_seaice, huge
      !
      real(kind=kind_phys), dimension(:), intent(inout)  :: zorlo, zorll, zorli
      !
      real(kind=kind_phys), parameter :: timin = 173.0_kind_phys  ! minimum temperature allowed for snow/ice

      real(kind=kind_phys) :: tem

      ! CCPP error handling
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (frac_grid) then  ! cice is ice fraction wrt water area
        do i=1,im
          frland(i) = landfrac(i)
          if (frland(i) > zero) dry(i) = .true.
          if (frland(i) < one) then
            if (oceanfrac(i) > zero) then
              if (cice(i) >= min_seaice) then
                icy(i)  = .true.
                tisfc(i) = max(timin, min(tisfc(i), tgice))
                if (cplflx)  then
                  islmsk_cice(i) = 4
                  flag_cice(i)   = .true.
                else
                  islmsk_cice(i) = 2
                  flag_cice(i)   = .false.
                endif
                islmsk(i) = 2
              else
                cice(i)        = zero
                hice(i)        = zero
                flag_cice(i)   = .false.
                islmsk_cice(i) = 0
                islmsk(i)      = 0
                icy(i)         = .false.
              endif
              if (cice(i) < one) then
                wet(i) = .true. ! some open ocean
                if (.not. cplflx .and. icy(i)) tsfco(i) = max(tisfc(i), tgice)
              endif
            else
              if (cice(i) >= min_lakeice) then
                icy(i)    = .true.
                islmsk(i) = 2
                tisfc(i)  = max(timin, min(tisfc(i), tgice))
              else
                cice(i)   = zero
                hice(i)   = zero
                islmsk(i) = 0
                icy(i)    = .false.
              endif
              islmsk_cice(i) = islmsk(i)
              flag_cice(i)   = .false.
              if (cice(i) < one) then
                wet(i) = .true. ! some open lake
                if (icy(i)) tsfco(i) = max(tisfc(i), tgice)
              endif
            endif
          else            ! all land
            cice(i)        = zero
            hice(i)        = zero
            islmsk_cice(i) = 1
            islmsk(i)      = 1
            wet(i)         = .false.
            icy(i)         = .false.
            flag_cice(i)   = .false.
          endif
        enddo

      else

        do i = 1, IM
          if (islmsk(i) == 1) then
!           tsfcl(i)  = tsfc(i)
            dry(i)    = .true.
            frland(i) = one
            cice(i)   = zero
            hice(i)   = zero
            icy(i)    = .false.
          else
            frland(i) = zero
            if (oceanfrac(i) > zero) then
              if (cice(i) >= min_seaice) then
                icy(i)   = .true.
                tisfc(i) = max(timin, min(tisfc(i), tgice))
                ! This cplice namelist option was added to deal with the
                ! situation of the FV3ATM-HYCOM coupling without an active sea
                ! ice (e.g., CICE6) component. By default, the cplice is true
                ! when cplflx is .true. (e.g., for the S2S application).
                ! Whereas, for the HAFS FV3ATM-HYCOM coupling, cplice is set as
                ! .false.. In the future HAFS FV3ATM-MOM6 coupling, the cplflx
                ! could be .true., while cplice being .false..
                if (cplice .and. cplflx)  then
                  islmsk_cice(i) = 4
                  flag_cice(i)   = .true.
                else
                  islmsk_cice(i) = 2
                  flag_cice(i)   = .false.
                endif
                islmsk(i) = 2
              else
                cice(i)        = zero
                hice(i)        = zero
                flag_cice(i)   = .false.
                islmsk(i)      = 0
                islmsk_cice(i) = 0
                icy(i)         = .false.
              endif
              if (cice(i) < one) then
                wet(i) = .true. ! some open ocean
                if (cplice) then
                  if (.not. cplflx .and. icy(i)) tsfco(i) = max(tisfc(i), tgice)
                else
                  if (icy(i)) tsfco(i) = max(tisfc(i), tgice)
                endif
              endif
            else
              if (cice(i) >= min_lakeice) then
                icy(i) = .true.
                tisfc(i) = max(timin, min(tisfc(i), tgice))
                islmsk(i) = 2
              else
                cice(i)   = zero
                hice(i)   = zero
                islmsk(i) = 0
                icy(i)    = .false.
              endif
              islmsk_cice(i) = islmsk(i)
              flag_cice(i)   = .false.
              if (cice(i) < one) then
                wet(i) = .true. ! some open lake
                if (icy(i)) tsfco(i) = max(tisfc(i), tgice)
              endif
            endif
          endif
        enddo
      endif ! frac_grid

      do i=1,im
        tprcp_wat(i) = tprcp(i)
        tprcp_lnd(i) = tprcp(i)
        tprcp_ice(i) = tprcp(i)

        if (wet(i)) then                   ! Water
          uustar_wat(i) = uustar(i)
            tsfc_wat(i) = tsfco(i)
           tsurf_wat(i) = tsfco(i)
               zorlo(i) = max(1.0e-5, min(one, zorlo(i)))
        ! DH*
        else
          zorlo(i) = huge
        ! *DH
        endif
        if (dry(i)) then                   ! Land
          uustar_lnd(i) = uustar(i)
           if(lsm /= lsm_ruc) weasd_lnd(i) = weasd(i)
           tsurf_lnd(i) = tsfcl(i)
        ! DH*
        else
          zorll(i) = huge
        ! *DH
        !mjz
          tsfcl(i) = huge
        endif
        if (icy(i)) then                   ! Ice
          uustar_ice(i) = uustar(i)
           if(lsm /= lsm_ruc) weasd_ice(i) = weasd(i)
           tsurf_ice(i) = tisfc(i)
            ep1d_ice(i) = zero
            gflx_ice(i) = zero
               zorli(i) = max(1.0e-5, min(one, zorli(i)))
        ! DH*
        else
          zorli(i) = huge
        ! *DH
        endif
        if (nint(slmsk(i)) /= 1) slmsk(i)  = islmsk(i)
      enddo

! to prepare to separate lake from ocean under water category
      do i = 1, im
        if ((wet(i) .or. icy(i)) .and. lakefrac(i) > zero) then
          lake(i) = .true.
          if (lkm == 1 .and. lakefrac(i) >= 0.15 .and. lakedepth(i) > one) then
            use_flake(i) = .true.
          else
            use_flake(i) = .false.
          endif
        else
          lake(i) = .false.
          use_flake(i) = .false.
        endif
      enddo
!
      if (frac_grid) then
        do i=1,im
          if (dry(i)) then
            if (icy(i)) then
              if (kdt == 1 .or. (.not. cplflx .or. lakefrac(i) > zero)) then 
                tem = one / (cice(i)*(one-frland(i)))
                snowd_ice(i) = max(zero, (snowd(i) - snowd_lnd(i)*frland(i)) * tem)
                weasd_ice(i) = max(zero, (weasd(i) - weasd_lnd(i)*frland(i)) * tem)
              endif
            endif
          elseif (icy(i)) then
            if (kdt == 1 .or. (.not. cplflx .or. lakefrac(i) > zero)) then 
              tem = one / cice(i)
              snowd_lnd(i) = zero
              snowd_ice(i) = snowd(i) * tem
              weasd_lnd(i) = zero
              weasd_ice(i) = weasd(i) * tem
            endif
          endif
        enddo
      elseif(lsm /= lsm_ruc) then ! do not do snow initialization  with RUC lsm
        do i=1,im
          if (icy(i)) then
            if (kdt == 1 .or. (.not. cplflx .or. lakefrac(i) > zero)) then
              snowd_lnd(i) = zero
              weasd_lnd(i) = zero
              tem = one / cice(i)
              snowd_ice(i) = snowd(i) * tem
              weasd_ice(i) = weasd(i) * tem
            endif
          endif
        enddo
      endif

!     write(0,*)' minmax of ice snow=',minval(snowd_ice),maxval(snowd_ice)

   end subroutine GFS_surface_composites_pre_run

end module GFS_surface_composites_pre
