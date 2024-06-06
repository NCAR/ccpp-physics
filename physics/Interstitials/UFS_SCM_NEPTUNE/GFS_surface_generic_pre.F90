!> \file GFS_surface_generic_pre.F90
!!  Contains code related to running prior to all GFS surface schemes.

      module GFS_surface_generic_pre

      use machine, only: kind_phys

      implicit none

      private

      public GFS_surface_generic_pre_init, GFS_surface_generic_pre_run

      real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

      contains

!>\defgroup mod_GFS_surface_generic_pre GFS surface_generic_pre module
!! This module contains code related to running prior to all GFS surface schemes.
!> @{
!> \section arg_table_GFS_surface_generic_pre_init Argument Table
!! \htmlinclude GFS_surface_generic_pre_init.html
!!
      subroutine GFS_surface_generic_pre_init (nthreads, im, slmsk, isot, ivegsrc, stype,scolor, vtype, slope, &
                                               vtype_save, stype_save,scolor_save, slope_save, errmsg, errflg)

        implicit none

        ! Interface variables
        integer,                       intent(in)    :: nthreads, im, isot, ivegsrc
        real(kind_phys), dimension(:), intent(in)    :: slmsk
        integer,         dimension(:), intent(inout) :: vtype, stype, scolor,slope
        integer,         dimension(:), intent(out)   :: vtype_save, stype_save,scolor_save, slope_save

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        integer, dimension(1:im) :: islmsk
        integer :: i

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        islmsk = nint(slmsk)

        ! Save current values of vegetation, soil and slope type
        vtype_save(:) = vtype(:)
        stype_save(:) = stype(:)
        scolor_save(:) = scolor(:)
        slope_save(:) = slope(:)

        call update_vegetation_soil_slope_type(nthreads, im, isot, ivegsrc, islmsk, vtype, stype,scolor, slope)

      end subroutine GFS_surface_generic_pre_init

!> \section arg_table_GFS_surface_generic_pre_run Argument Table
!! \htmlinclude GFS_surface_generic_pre_run.html
!!
      subroutine GFS_surface_generic_pre_run (nthreads, im, levs, vfrac, islmsk, isot, ivegsrc, stype, scolor,vtype, slope, &
                          prsik_1, prslk_1, tsfc, phil, con_g, sigmaf, work3, zlvl,                        &
                          lndp_type, n_var_lndp, sfc_wts, lndp_var_list, lndp_prt_list,                    &
                          z01d, zt1d, bexp1d, xlai1d, vegf1d, lndp_vgf,                                    &
                          cplflx, flag_cice, islmsk_cice, slimskin_cpl,                                    &
                          wind, u1, v1, cnvwind, smcwlt2, smcref2, vtype_save, stype_save,scolor_save, slope_save,     &
                          errmsg, errflg)

        use surface_perturbation,  only: cdfnor

        implicit none

        ! Interface variables
        integer, intent(in) :: nthreads, im, levs, isot, ivegsrc
        integer, dimension(:), intent(in) :: islmsk

        real(kind=kind_phys), intent(in) :: con_g
        real(kind=kind_phys), dimension(:), intent(in) :: vfrac, prsik_1, prslk_1
        integer, dimension(:), intent(inout) :: vtype, stype,scolor, slope
        integer, dimension(:), intent(out)   :: vtype_save(:), stype_save(:),scolor_save(:), slope_save(:)

        real(kind=kind_phys), dimension(:), intent(inout) :: tsfc
        real(kind=kind_phys), dimension(:,:), intent(in) :: phil

        real(kind=kind_phys), dimension(:), intent(inout) :: sigmaf, work3, zlvl

        ! Stochastic physics / surface perturbations
        integer,                              intent(in)  :: lndp_type, n_var_lndp
        character(len=3),     dimension(:),   intent(in), optional  :: lndp_var_list
        real(kind=kind_phys), dimension(:),   intent(in), optional  :: lndp_prt_list
        real(kind=kind_phys), dimension(:,:), intent(in), optional  :: sfc_wts
        real(kind=kind_phys), dimension(:),   intent(out) :: z01d
        real(kind=kind_phys), dimension(:),   intent(out) :: zt1d
        real(kind=kind_phys), dimension(:),   intent(out) :: bexp1d
        real(kind=kind_phys), dimension(:),   intent(out) :: xlai1d
        real(kind=kind_phys), dimension(:),   intent(out) :: vegf1d
        real(kind=kind_phys),                 intent(out) :: lndp_vgf

        logical,                              intent(in)    :: cplflx
        real(kind=kind_phys), dimension(:),   intent(in), optional    :: slimskin_cpl
        logical,              dimension(:),   intent(inout) :: flag_cice
        integer,              dimension(:),   intent(out)   :: islmsk_cice

        real(kind=kind_phys), dimension(:),   intent(out) :: wind
        real(kind=kind_phys), dimension(:),   intent(in ) :: u1, v1
        ! surface wind enhancement due to convection
        real(kind=kind_phys), dimension(:),   intent(inout ), optional :: cnvwind
        !
        real(kind=kind_phys), dimension(:),   intent(out)    :: smcwlt2, smcref2

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        integer              :: i, k
        real(kind=kind_phys) :: onebg, cdfz

        ! Set constants
        onebg  = 1.0/con_g

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Scale random patterns for surface perturbations with perturbation size
        ! Turn vegetation fraction pattern into percentile pattern
        lndp_vgf=-999.

        if (lndp_type==1) then
          do k =1,n_var_lndp
            select case(lndp_var_list(k))
            case ('rz0')
                z01d(:) = lndp_prt_list(k)* sfc_wts(:,k)
            case ('rzt')
                 zt1d(:) = lndp_prt_list(k)* sfc_wts(:,k)
            case ('shc')
                 bexp1d(:) = lndp_prt_list(k) * sfc_wts(:,k)
            case ('lai')
                xlai1d(:) = lndp_prt_list(k)* sfc_wts(:,k)
            case ('vgf')
        ! note that the pertrubed vegfrac is being used in sfc_drv, but not sfc_diff
              do i=1,im
                call cdfnor(sfc_wts(i,k),cdfz)
                vegf1d(i) = cdfz
              enddo
              lndp_vgf = lndp_prt_list(k)
            end select
          enddo
        endif

        ! End of stochastic physics / surface perturbation

        ! Save current values of vegetation, soil and slope type
        vtype_save(:) = vtype(:)
        stype_save(:) = stype(:)
        scolor_save(:) = scolor(:)
        slope_save(:) = slope(:)

        call update_vegetation_soil_slope_type(nthreads, im, isot, ivegsrc, islmsk, vtype, stype,scolor, slope)

        do i=1,im
          sigmaf(i) = max(vfrac(i), 0.01_kind_phys)
          islmsk_cice(i) = islmsk(i)

          work3(i)   = prsik_1(i) / prslk_1(i)

          zlvl(i)    = phil(i,1) * onebg
          smcwlt2(i) = zero
          smcref2(i) = zero

          wind(i)  = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))   &
                         + max(zero, min(cnvwind(i), 30.0_kind_phys)), one)
         !wind(i)  = max(sqrt(Statein%ugrs(i,1)*Statein%ugrs(i,1) + &
         !                         Statein%vgrs(i,1)*Statein%vgrs(i,1))  &
         !              + max(zero, min(Tbd%phy_f2d(i,Model%num_p2d), 30.0)), one)
          cnvwind(i) = zero

        enddo

      if (cplflx) then
        do i=1,im
          islmsk_cice(i) = nint(slimskin_cpl(i))
          flag_cice(i)   = (islmsk_cice(i) == 4)
        enddo
      endif

      end subroutine GFS_surface_generic_pre_run

      subroutine update_vegetation_soil_slope_type(nthreads, im, isot, ivegsrc, islmsk, vtype, stype,scolor, slope)

        implicit none

        integer, intent(in)    :: nthreads, im, isot, ivegsrc, islmsk(:)
        integer, intent(inout) :: vtype(:), stype(:),scolor(:), slope(:)
        integer :: i

!$OMP  parallel do num_threads(nthreads) default(none) private(i) &
!$OMP      shared(im, isot, ivegsrc, islmsk, vtype, stype,scolor, slope)

! scolor is a place holder now, how to update soil color based on the mask/veg/sot src

        do i=1,im
          if (islmsk(i) == 2) then
            if (isot == 1) then
              stype(i) = 16
            else
              stype(i) = 9
            endif
            if (ivegsrc == 0 .or. ivegsrc == 4) then
              vtype(i) = 24
            elseif (ivegsrc == 1) then
              vtype(i) = 15
            elseif (ivegsrc == 2) then
              vtype(i) = 13
            elseif (ivegsrc == 3 .or. ivegsrc == 5) then
              vtype(i) = 15
            endif
            slope(i)  = 9
          else
            if (vtype(i)  < 1) vtype(i)  = 17
            if (slope(i) < 1) slope(i) = 1
          endif
        enddo
!$OMP end parallel do

      end subroutine update_vegetation_soil_slope_type
!> @}

      end module GFS_surface_generic_pre
