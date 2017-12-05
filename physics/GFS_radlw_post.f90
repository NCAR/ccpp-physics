      subroutine Post_radiation (Radtend, tsfa, lm, kd, htlwc, htlw0, &
          Model, Coupling, Grid, htswc, htsw0, scmpsw, sfcalb, Diag,  &
          nday, Statein, im, kt, kb, raddt, cldsa, mtopa, mbota,      &
          clouds, aerodp)

        implicit none

        integer, intent(in) :: lm, kd, im, kt, kb
        type(GFS_grid_type),    intent(in)     :: Grid
        type(GFS_control_type), intent(in)     :: Model
        type(GFS_statein_type), intent(in)     :: Statein
        type(GFS_radtend_type), intent(inout)  :: Radtend
        type(GFS_coupling_type), intent(inout) :: Coupling
        type(GFS_diag_type), intent(inout)     :: Diag

        real(kind = kind_phys), dimension(Size (Grid%xlon, 1), Model%levr + &
          LTP), intent(in) :: htlw0, htlwc, htswc, htsw0
        real(kind = kind_phys), dimension(Size (Grid%xlon, 1)), intent(in) :: tsfa
        type(cmpfsw_type), dimension(size(Grid%xlon,1)), intent(inout) :: scmpsw
        real(kind = kind_phys), dimension(Size (Grid%xlon, 1), NF_ALBD), intent(in) :: sfcalb
        integer, intent(in) :: nday
        real(kind = kind_phys), intent(in) :: raddt
        real(kind = kind_phys), dimension(Size (Grid%xlon, 1), 5), intent(in) :: cldsa
        integer, dimension(size(Grid%xlon, 1), 3), intent(in) :: mbota, mtopa
        real(kind = kind_phys), dimension(Size (Grid%xlon, 1), Model%levr + &
            LTP, NF_CLDS), intent(in) :: clouds
        real(kind = kind_phys), dimension(Size (Grid%xlon, 1), NSPC1), intent(in) :: aerodp


        !pedro Save LW results
        !pedro call Post_lw (Radtend, tsfa, lm, kd, htlwc, htlw0, Model, Coupling, Grid)

!>  - Save surface air temp for diurnal adjustment at model t-steps

      if (Model%lslwr) then
        Radtend%tsflw (:) = tsfa(:)

        do k = 1, LM
          k1 = k + kd
            Radtend%htrlw(:,k) = htlwc(:,k1)
        enddo
        ! --- repopulate the points above levr
        if (Model%levr < Model%levs) then
          do k = LM,Model%levs
            Radtend%htrlw (:,k) = Radtend%htrlw (:,LM)
          enddo
        endif

        if (Model%lwhtr) then
          do k = 1, lm
            k1 = k + kd
            Radtend%lwhc(:,k) = htlw0(:,k1)
          enddo
          ! --- repopulate the points above levr
          if (Model%levr < Model%levs) then
            do k = LM,Model%levs
              Radtend%lwhc(:,k) = Radtend%lwhc(:,LM)
            enddo
          endif
        endif

! --- radiation fluxes for other physics processes
        Coupling%sfcdlw(:) = Radtend%sfcflw(:)%dnfxc

      endif                                ! end_if_lslwr


            ! post SW
        call Save_sw_heating_rate (Radtend, Model, Grid, htswc, lm, kd, &
            Model%lsswr)

        call Save_sw_heating_rate_csk (Radtend, Model, Grid, htsw0, lm, &
            kd, Model%lsswr)

            ! Surface down and up spectral component fluxes
            ! Save two spectral bands' surface downward and upward fluxes for output.
        call Save_sw_fluxes (Coupling, scmpsw, Grid, sfcalb, Model%lsswr)

            ! Night time: set SW heating rates and fluxes to zero
        call Zero_out_heatrate_flux (Radtend, Diag, scmpsw, Coupling, &
            Grid, Model, nday, Model%lsswr)

        call Save_more_sw_fluxes (Radtend, Coupling, Model%lsswr)


          ! Collect the fluxr data for wrtsfc
        call Organize_output (Diag, Model, Grid, Radtend, Statein, &
            Coupling, im, kd, kt, kb, lm, scmpsw, raddt, cldsa,    &
            mtopa, mbota, clouds, aerodp)

      end subroutine Post_radiation

