!>\file GFS_radiation_surface.F90
!! This file contains calls to module_radiation_surface::setemis() to set up
!! surface emissivity for LW radiation and to module_radiation_surface::setalb()
!! to set up surface albedo for SW radiation.

      module GFS_radiation_surface

      use machine,                   only: kind_phys

      contains

!>\defgroup GFS_radiation_surface_mod GFS Radiation Surface Module
!! This module contains calls to module_radiation_surface::setemis() to set up
!! surface emissivity for LW radiation and to module_radiation_surface::setalb()
!! to set up surface albedo for SW radiation.
!> @{
!> \section arg_table_GFS_radiation_surface_init Argument Table
!! \htmlinclude GFS_radiation_surface_init.html
!!
      subroutine GFS_radiation_surface_init (me, ialb, iems, semis_file, con_pi, errmsg, errflg)

      use module_radiation_surface, only: sfc_init

      implicit none

      integer,                              intent(in)  :: me, ialb, iems
      character(len=26),                    intent(in)  :: semis_file
      real(kind_phys),                      intent(in)  :: con_pi
      character(len=*),                     intent(out) :: errmsg
      integer,                              intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if ( me == 0 ) then
        print *,'In GFS_radiation_surface_init, before calling sfc_init'
        print *,'ialb=',ialb,' iems=',iems
      end if

      ! Call surface initialization routine
      call sfc_init ( me, ialb, iems, semis_file, con_pi, errmsg, errflg )

      end subroutine GFS_radiation_surface_init


!> \section arg_table_GFS_radiation_surface_run Argument Table
!! \htmlinclude GFS_radiation_surface_run.html
!!
      subroutine GFS_radiation_surface_run (                            &
        ialb, im, nf_albd, frac_grid, lslwr, lsswr, lsm, lsm_noahmp,    &
        lsm_ruc, xlat, xlon, slmsk, lndp_type, n_var_lndp, sfc_alb_pert,&
        lndp_var_list, lndp_prt_list, landfrac, snodl, snodi, sncovr,   &
        sncovr_ice, fice, zorl, hprime, tsfg, tsfa, tisfc, coszen,      &
        cplice, min_seaice, min_lakeice, lakefrac, use_flake,           &
        alvsf, alnsf, alvwf, alnwf, facsf, facwf,                       &
        semis_lnd, semis_ice, semis_wat, snoalb, use_cice_alb, con_ttp, &
        albdvis_lnd, albdnir_lnd, albivis_lnd, albinir_lnd,             &
        albdvis_ice, albdnir_ice, albivis_ice, albinir_ice,             &
        semisbase, semis, sfcalb, sfc_alb_dif, errmsg, errflg)

      use module_radiation_surface,  only: f_zero, f_one,  &
                                           epsln,          &
                                           setemis, setalb

      implicit none

      integer,               intent(in) :: im, nf_albd, ialb
      logical,               intent(in) :: frac_grid, lslwr, lsswr, use_cice_alb, cplice
      integer,               intent(in) :: lsm, lsm_noahmp, lsm_ruc, lndp_type, n_var_lndp
      real(kind=kind_phys),  intent(in) :: min_seaice, min_lakeice, con_ttp
      logical, dimension(:), intent(in) :: use_flake

      real(kind=kind_phys), dimension(:),   intent(in)  :: xlat, xlon, slmsk,           &
                                                           sfc_alb_pert, lndp_prt_list, &
                                                           landfrac, lakefrac,          &
                                                           snodl, snodi, sncovr,        &
                                                           sncovr_ice, fice, zorl,      &
                                                           hprime, tsfg, tsfa, tisfc,   &
                                                           coszen, alvsf, alnsf, alvwf, &
                                                           alnwf, facsf, facwf, snoalb
      character(len=3)    , dimension(:),   intent(in)  :: lndp_var_list
      real(kind=kind_phys), dimension(:),   intent(in)  :: albdvis_ice, albdnir_ice,    &
                                                           albivis_ice, albinir_ice

      real(kind=kind_phys), dimension(:),   intent(inout) :: albdvis_lnd, albdnir_lnd,  &
                                                             albivis_lnd, albinir_lnd,  &
                                                             semis_lnd,   semis_ice, semis_wat
      real(kind=kind_phys), dimension(:),   intent(inout) :: semisbase, semis
      real(kind=kind_phys), dimension(:,:), intent(inout) :: sfcalb
      real(kind=kind_phys), dimension(:),   intent(inout) :: sfc_alb_dif

      character(len=*),                     intent(out) :: errmsg
      integer,                              intent(out) :: errflg

      ! Local variables
      integer                             :: i
      real(kind=kind_phys)                :: lndp_alb
      real(kind=kind_phys), dimension(im) :: cimin, fracl, fraci, fraco
      logical,              dimension(im) :: icy

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Return immediately if neither shortwave nor longwave radiation are called
      if (.not. lsswr .and. .not. lslwr) return

      do i=1,im
        if (lakefrac(i) > f_zero) then
          cimin(i) = min_lakeice
        else
          cimin(i) = min_seaice
        endif
      enddo

      ! Set up land/ice/ocean fractions for emissivity and albedo calculations
      if (.not. frac_grid) then
        do i=1,im
          if (slmsk(i) == 1) then
            fracl(i) = f_one
            fraci(i) = f_zero
            fraco(i) = f_zero
            icy(i)   = .false.
          else
            fracl(i) = f_zero
            fraco(i) = f_one
            if(fice(i) < cimin(i)) then
              fraci(i) = f_zero
              icy(i)   = .false.
            else
              fraci(i) = fraco(i) * fice(i)
              icy(i)   = .true.
            endif
            fraco(i) = max(f_zero, fraco(i)-fraci(i))
          endif
        enddo
      else
        do i=1,im
          fracl(i) = landfrac(i)
          fraco(i) = max(f_zero, f_one - fracl(i))
          if(fice(i) < cimin(i)) then
            fraci(i) = f_zero
            icy(i)   = .false.
          else
            fraci(i) = fraco(i) * fice(i)
            icy(i)   = .true.
          endif
          fraco(i) = max(f_zero, fraco(i)-fraci(i))
        enddo
      endif

      if (lslwr) then
!>  - Call module_radiation_surface::setemis(),to set up surface
!! emissivity for LW radiation.
        call setemis (lsm, lsm_noahmp, lsm_ruc, frac_grid, cplice,  &
                      use_flake, lakefrac, xlon, xlat, slmsk,       &
!                     frac_grid, min_seaice, xlon, xlat, slmsk,     &
                      snodl, snodi, sncovr, sncovr_ice, zorl, tsfg, &
                      tsfa, hprime, semis_lnd, semis_ice, semis_wat,&
                      im, fracl, fraco, fraci, icy,                 & !  ---  inputs
                      semisbase, semis)                               !  ---  outputs
      endif

      if (lsswr) then
!>  - Set surface albedo perturbation, if requested
        lndp_alb = -999.
        if (lndp_type==1) then
          do i =1,n_var_lndp
            if (lndp_var_list(i) == 'alb') then
                lndp_alb = lndp_prt_list(i)
            endif
          enddo
        endif

!>  - Call module_radiation_surface::setalb(),to set up surface
!! albedor for SW radiation.

        call setalb (slmsk, lsm, lsm_noahmp, lsm_ruc, use_cice_alb, snodi, sncovr, sncovr_ice, &
                     snoalb, zorl, coszen, tsfg, tsfa, hprime, frac_grid, lakefrac,            &
!                    snoalb, zorl, coszen, tsfg, tsfa, hprime, frac_grid, min_seaice,          &
                     alvsf, alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc,                    &
                     albdvis_lnd, albdnir_lnd, albivis_lnd, albinir_lnd,                       &
                     albdvis_ice, albdnir_ice, albivis_ice, albinir_ice,                       &
                     im, nf_albd, sfc_alb_pert, lndp_alb, fracl, fraco, fraci, icy, ialb,      &
                     con_ttp,                                                                  & !  ---  inputs
                     sfcalb )                                                                    !  ---  outputs

!> -# Approximate mean surface albedo from vis- and nir- diffuse values.
        sfc_alb_dif(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
      endif

      end subroutine GFS_radiation_surface_run

!> @}
       end module GFS_radiation_surface
