!>\file gcycle.F90
!! This file repopulates specific time-varying surface properties for
!! atmospheric forecast runs.

!>\ingroup mod_GFS_phys_time_vary
!! This subroutine repopulates specific time-varying surface properties for
!! atmospheric forecast runs.
  SUBROUTINE GCYCLE (nblks, nthrds, Model, Grid, Sfcprop, Cldprop)
!
!
    USE MACHINE,      only: kind_phys
    USE PHYSCONS,     only: PI => con_PI
    USE GFS_typedefs, only: GFS_control_type, GFS_grid_type, &
                            GFS_sfcprop_type, GFS_cldprop_type
    implicit none

    integer,                  intent(in)    :: nblks, nthrds
    type(GFS_control_type),   intent(in)    :: Model
    type(GFS_grid_type),      intent(in)    :: Grid(nblks)
    type(GFS_sfcprop_type),   intent(inout) :: Sfcprop(nblks)
    type(GFS_cldprop_type),   intent(inout) :: Cldprop(nblks)

!
!     Local variables
!     ---------------
    integer              ::                     &
           I_INDEX(Model%nx*Model%ny),          &
           J_INDEX(Model%nx*Model%ny)

    real(kind=kind_phys) ::                     &
           RLA (Model%nx*Model%ny),             &
           RLO (Model%nx*Model%ny),             &
        SLMASK (Model%nx*Model%ny),             &
          OROG (Model%nx*Model%ny),             &
       OROG_UF (Model%nx*Model%ny),             &
        SLIFCS (Model%nx*Model%ny),             &
        TSFFCS (Model%nx*Model%ny),             &
        SNOFCS (Model%nx*Model%ny),             &
        ZORFCS (Model%nx*Model%ny),             &
        TG3FCS (Model%nx*Model%ny),             &
        CNPFCS (Model%nx*Model%ny),             &
        AISFCS (Model%nx*Model%ny),             &
!       F10MFCS(Model%nx*Model%ny),             &
        VEGFCS (Model%nx*Model%ny),             &
        VETFCS (Model%nx*Model%ny),             &
        SOTFCS (Model%nx*Model%ny),             &
         CVFCS (Model%nx*Model%ny),             &
        CVBFCS (Model%nx*Model%ny),             &
        CVTFCS (Model%nx*Model%ny),             &
        SWDFCS (Model%nx*Model%ny),             &
        SIHFCS (Model%nx*Model%ny),             &
        SICFCS (Model%nx*Model%ny),             &
        SITFCS (Model%nx*Model%ny),             &
        VMNFCS (Model%nx*Model%ny),             &
        VMXFCS (Model%nx*Model%ny),             &
        SLPFCS (Model%nx*Model%ny),             &
        ABSFCS (Model%nx*Model%ny),             &
        ALFFC1 (Model%nx*Model%ny*2),           &
        ALBFC1 (Model%nx*Model%ny*4),           &
        SMCFC1 (Model%nx*Model%ny*Model%lsoil), &
        STCFC1 (Model%nx*Model%ny*Model%lsoil), &
        SLCFC1 (Model%nx*Model%ny*Model%lsoil)

    logical          :: lake(Model%nx*Model%ny)

    character(len=6) :: tile_num_ch
    real(kind=kind_phys), parameter :: pifac=180.0/pi
    real(kind=kind_phys)            :: sig1t, dt_warm
    integer :: npts, len, nb, ix, jx, ls, ios, ll
    logical :: exists
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     if (Model%me .eq. 0) print *,' nlats=',nlats,' lonsinpe='
!    *,lonsinpe(0,1)

      tile_num_ch = "      "
      if (Model%tile_num < 10) then
        write(tile_num_ch, "(a4,i1)") "tile", Model%tile_num
      else
        write(tile_num_ch, "(a4,i2)") "tile", Model%tile_num
      endif

      len = 0
      do jx = Model%jsc, (Model%jsc+Model%ny-1)
        do ix = Model%isc, (Model%isc+Model%nx-1)
          len = len + 1
          i_index(len) = ix
          j_index(len) = jx
        enddo
      enddo

      sig1t = 0.0_kind_phys
      npts  = Model%nx*Model%ny
!
      len = 0
      do nb = 1,nblks
        do ix = 1,size(Grid(nb)%xlat,1)
          len = len + 1
          RLA     (len)          = Grid(nb)%xlat      (ix) * pifac
          RLO     (len)          = Grid(nb)%xlon      (ix) * pifac
          OROG    (len)          = Sfcprop(nb)%oro    (ix)
          OROG_UF (len)          = Sfcprop(nb)%oro_uf (ix)
          SLIFCS  (len)          = Sfcprop(nb)%slmsk  (ix)
          if ( Model%nstf_name(1) > 0 ) then
            TSFFCS(len)          = Sfcprop(nb)%tref   (ix)
          else
            TSFFCS(len)          = Sfcprop(nb)%tsfc   (ix)
          endif
          SNOFCS  (len)          = Sfcprop(nb)%weasd  (ix)
          ZORFCS  (len)          = Sfcprop(nb)%zorll  (ix)
          if (SLIFCS(len) > 1.9_kind_phys .and. .not. Model%frac_grid) then
            ZORFCS  (len)        = Sfcprop(nb)%zorli  (ix)
          elseif (SLIFCS(len) < 0.1_kind_phys .and. .not. Model%frac_grid) then
            ZORFCS  (len)        = Sfcprop(nb)%zorlo  (ix)
          endif
          TG3FCS  (len)          = Sfcprop(nb)%tg3    (ix)
          CNPFCS  (len)          = Sfcprop(nb)%canopy (ix)
!         F10MFCS (len)          = Sfcprop(nb)%f10m   (ix)
          VEGFCS  (len)          = Sfcprop(nb)%vfrac  (ix)
          VETFCS  (len)          = Sfcprop(nb)%vtype  (ix)
          SOTFCS  (len)          = Sfcprop(nb)%stype  (ix)
          CVFCS   (len)          = Cldprop(nb)%cv     (ix)
          CVBFCS  (len)          = Cldprop(nb)%cvb    (ix)
          CVTFCS  (len)          = Cldprop(nb)%cvt    (ix)
          SWDFCS  (len)          = Sfcprop(nb)%snowd  (ix)
          SIHFCS  (len)          = Sfcprop(nb)%hice   (ix)
          SICFCS  (len)          = Sfcprop(nb)%fice   (ix)
          SITFCS  (len)          = Sfcprop(nb)%tisfc  (ix)
          VMNFCS  (len)          = Sfcprop(nb)%shdmin (ix)
          VMXFCS  (len)          = Sfcprop(nb)%shdmax (ix)
          SLPFCS  (len)          = Sfcprop(nb)%slope  (ix)
          ABSFCS  (len)          = Sfcprop(nb)%snoalb (ix)

          ALFFC1  (len       )   = Sfcprop(nb)%facsf  (ix)
          ALFFC1  (len + npts)   = Sfcprop(nb)%facwf  (ix)

          ALBFC1  (len         ) = Sfcprop(nb)%alvsf  (ix)
          ALBFC1  (len + npts  ) = Sfcprop(nb)%alvwf  (ix)
          ALBFC1  (len + npts*2) = Sfcprop(nb)%alnsf  (ix)
          ALBFC1  (len + npts*3) = Sfcprop(nb)%alnwf  (ix)

          do ls = 1,Model%lsoil
            SMCFC1 (len + (ls-1)*npts) = Sfcprop(nb)%smc (ix,ls)
            STCFC1 (len + (ls-1)*npts) = Sfcprop(nb)%stc (ix,ls)
            SLCFC1 (len + (ls-1)*npts) = Sfcprop(nb)%slc (ix,ls)
          enddo

          IF (SLIFCS(len) < 0.1_kind_phys .OR. SLIFCS(len) > 1.5_kind_phys) THEN
             SLMASK(len) = 0.0_kind_phys
          ELSE
             SLMASK(len) = 1.0_kind_phys
          ENDIF

          IF (SLIFCS(len) > 1.99_kind_phys) THEN
            AISFCS(len) = 1.0_kind_phys
          ELSE
            AISFCS(len) = 0.0_kind_phys
          ENDIF
          if (Sfcprop(nb)%lakefrac(ix) > 0.0_kind_phys) then
            lake(len) = .true.
          else
            lake(len) = .false.
          endif

!     if (Model%me .eq. 0)
!    &   print *,' len=',len,' rla=',rla(len),' rlo=',rlo(len)
        ENDDO                 !-----END BLOCK SIZE LOOP------------------------------
      ENDDO                   !-----END BLOCK LOOP-------------------------------

! check
!     call mymaxmin(slifcs,len,len,1,'slifcs')
!     call mymaxmin(slmask,len,len,1,'slmsk')
!
#ifndef INTERNAL_FILE_NML
      inquire (file=trim(Model%fn_nml),exist=exists)
      if (.not. exists) then
        write(6,*) 'gcycle:: namelist file: ',trim(Model%fn_nml),' does not exist'
        stop
      else
        open (unit=Model%nlunit, file=trim(Model%fn_nml), action='READ', status='OLD', iostat=ios)
        rewind (Model%nlunit)
      endif
#endif
      CALL SFCCYCLE (9998, npts, Model%lsoil, SIG1T, Model%fhcyc, &
                     Model%idate(4), Model%idate(2),              &
                     Model%idate(3), Model%idate(1),              &
                     Model%phour, RLA, RLO, SLMASK,               &
!                    Model%fhour, RLA, RLO, SLMASK,               &
                     OROG, OROG_UF, Model%USE_UFO, Model%nst_anl, &
                     SIHFCS, SICFCS, SITFCS, SWDFCS, SLCFC1,      &
                     VMNFCS, VMXFCS, SLPFCS, ABSFCS, TSFFCS,      &
                     SNOFCS, ZORFCS, ALBFC1, TG3FCS, CNPFCS,      &
                     SMCFC1, STCFC1, SLIFCS, AISFCS,              &
                     VEGFCS, VETFCS, SOTFCS, ALFFC1, CVFCS,       &
                     CVBFCS, CVTFCS, Model%me, nthrds,            &
                     Model%nlunit, size(Model%input_nml_file),    &
                     Model%input_nml_file,                        &
                     lake, Model%min_lakeice, Model%min_seaice,   &
                     Model%ialb, Model%isot, Model%ivegsrc,       &
                     trim(tile_num_ch), i_index, j_index)
#ifndef INTERNAL_FILE_NML
      close (Model%nlunit)
#endif

      len = 0
      do nb = 1,nblks
        do ix = 1,size(Grid(nb)%xlat,1)
          len = len + 1
          Sfcprop(nb)%slmsk  (ix) = SLIFCS  (len)
          if ( Model%nstf_name(1) > 0 ) then
             Sfcprop(nb)%tref(ix) = TSFFCS  (len)
!           if ( Model%nstf_name(2) == 0 ) then
!             dt_warm = (Sfcprop(nb)%xt(ix) + Sfcprop(nb)%xt(ix) ) &
!                     / Sfcprop(nb)%xz(ix)
!             Sfcprop(nb)%tsfco(ix) = Sfcprop(nb)%tref(ix)         &
!                                   + dt_warm - Sfcprop(nb)%dt_cool(ix)
!           endif
          else
             Sfcprop(nb)%tsfc(ix)  = TSFFCS  (len)
             Sfcprop(nb)%tsfco(ix) = TSFFCS  (len)
          endif
          Sfcprop(nb)%weasd  (ix) = SNOFCS  (len)
          Sfcprop(nb)%zorll  (ix) = ZORFCS  (len)
          if (SLIFCS(len) > 1.9_kind_phys .and. .not. Model%frac_grid) then
            Sfcprop(nb)%zorli(ix) = ZORFCS  (len)
          elseif (SLIFCS(len) < 0.1_kind_phys .and. .not. Model%frac_grid) then
            Sfcprop(nb)%zorlo(ix) = ZORFCS  (len)
          endif
          Sfcprop(nb)%tg3    (ix) = TG3FCS  (len)
          Sfcprop(nb)%canopy (ix) = CNPFCS  (len)
!         Sfcprop(nb)%f10m   (ix) = F10MFCS (len)
          Sfcprop(nb)%vfrac  (ix) = VEGFCS  (len)
          Sfcprop(nb)%vtype  (ix) = VETFCS  (len)
          Sfcprop(nb)%stype  (ix) = SOTFCS  (len)
          Cldprop(nb)%cv     (ix) = CVFCS   (len)
          Cldprop(nb)%cvb    (ix) = CVBFCS  (len)
          Cldprop(nb)%cvt    (ix) = CVTFCS  (len)
          Sfcprop(nb)%snowd  (ix) = SWDFCS  (len)
          Sfcprop(nb)%hice   (ix) = SIHFCS  (len)
          Sfcprop(nb)%fice   (ix) = SICFCS  (len)
          Sfcprop(nb)%tisfc  (ix) = SITFCS  (len)
          Sfcprop(nb)%shdmin (ix) = VMNFCS  (len)
          Sfcprop(nb)%shdmax (ix) = VMXFCS  (len)
          Sfcprop(nb)%slope  (ix) = SLPFCS  (len)
          Sfcprop(nb)%snoalb (ix) = ABSFCS  (len)

          Sfcprop(nb)%facsf  (ix) = ALFFC1  (len       )
          Sfcprop(nb)%facwf  (ix) = ALFFC1  (len + npts)

          Sfcprop(nb)%alvsf  (ix) = ALBFC1  (len         )
          Sfcprop(nb)%alvwf  (ix) = ALBFC1  (len + npts  )
          Sfcprop(nb)%alnsf  (ix) = ALBFC1  (len + npts*2)
          Sfcprop(nb)%alnwf  (ix) = ALBFC1  (len + npts*3)
          do ls = 1,Model%lsoil
            ll = len + (ls-1)*npts
            Sfcprop(nb)%smc (ix,ls) = SMCFC1 (ll)
            Sfcprop(nb)%stc (ix,ls) = STCFC1 (ll)
            Sfcprop(nb)%slc (ix,ls) = SLCFC1 (ll)
            if (ls<=Model%kice) Sfcprop(nb)%tiice (ix,ls) = STCFC1 (ll)
          enddo
        ENDDO                 !-----END BLOCK SIZE LOOP--------------------------
      ENDDO                   !-----END BLOCK LOOP-------------------------------

! check
!     call mymaxmin(slifcs,len,len,1,'slifcs')
!
!     if (Model%me .eq. 0) print*,'executed gcycle during hour=',fhour

      RETURN
      END
