!>\file gcycle.F90
!! This file repopulates specific time-varying surface properties for
!! atmospheric forecast runs.

module gcycle_mod

  implicit none

  private

  public gcycle

contains

!>\ingroup mod_GFS_phys_time_vary
!! This subroutine repopulates specific time-varying surface properties for
!! atmospheric forecast runs.
  subroutine gcycle (me, nthrds, nx, ny, isc, jsc, nsst, tile_num, nlunit,        &
      input_nml_file, lsoil, lsoil_lsm, kice, idate, ialb, isot, ivegsrc,         &
      use_ufo, nst_anl, fhcyc, phour, landfrac, lakefrac, min_seaice, min_lakeice,&
      frac_grid, smc, slc, stc, smois, sh2o, tslb, tiice, tg3, tref, tsfc,        &
      tsfco, tisfc, hice, fice, facsf, facwf, alvsf, alvwf, alnsf, alnwf,         &
      zorli, zorll, zorlo, weasd, slope, snoalb, canopy, vfrac, vtype,            &
      stype, shdmin, shdmax, snowd, cv, cvb, cvt, oro, oro_uf,                    &
      xlat_d, xlon_d, slmsk, imap, jmap)
!
!
    use machine,      only: kind_phys, kind_io8
    implicit none

    integer,              intent(in)    :: me, nthrds, nx, ny, isc, jsc, nsst,    &
                                           tile_num, nlunit, lsoil, lsoil_lsm, kice
    integer,              intent(in)    :: idate(:), ialb, isot, ivegsrc
    character(len=*),     intent(in)    :: input_nml_file(:)
    logical,              intent(in)    :: use_ufo, nst_anl, frac_grid
    real(kind=kind_phys), intent(in)    :: fhcyc, phour, landfrac(:), lakefrac(:),&
                                           min_seaice, min_lakeice,               &
                                           xlat_d(:), xlon_d(:)
    real(kind=kind_phys), intent(inout) :: smc(:,:),   &
                                           slc(:,:),   &
                                           stc(:,:),   &
                                           smois(:,:), &
                                           sh2o(:,:),  &
                                           tslb(:,:),  &
                                           tiice(:,:), &
                                           tg3(:),     &
                                           tref(:),    &
                                           tsfc(:),    &
                                           tsfco(:),   &
                                           tisfc(:),   &
                                           hice(:),    &
                                           fice(:),    &
                                           facsf(:),   &
                                           facwf(:),   &
                                           alvsf(:),   &
                                           alvwf(:),   &
                                           alnsf(:),   &
                                           alnwf(:),   &
                                           zorli(:),   &
                                           zorll(:),   &
                                           zorlo(:),   &
                                           weasd(:),   &
                                           slope(:),   &
                                           snoalb(:),  &
                                           canopy(:),  &
                                           vfrac(:),   &
                                           vtype(:),   &
                                           stype(:),   &
                                           shdmin(:),  &
                                           shdmax(:),  &
                                           snowd(:),   &
                                           cv(:),      &
                                           cvb(:),     &
                                           cvt(:),     &
                                           oro(:),     &
                                           oro_uf(:),  &
                                           slmsk(:)

    integer,              intent(in)    :: imap(:), jmap(:)
!
!     Local variables
!     ---------------
    real(kind=kind_phys) ::                  &
        slmskl (nx*ny),                      &
        slmskw (nx*ny),                      &
        TSFFCS (nx*ny),                      &
        ZORFCS (nx*ny),                      &
        AISFCS (nx*ny),                      &
        ALFFC1 (nx*ny*2),                    &
        ALBFC1 (nx*ny*4),                    &
        SMCFC1 (nx*ny*max(lsoil,lsoil_lsm)), &
        STCFC1 (nx*ny*max(lsoil,lsoil_lsm)), &
        SLCFC1 (nx*ny*max(lsoil,lsoil_lsm))


    real (kind=kind_io8) :: min_ice(nx*ny)
    character(len=6)     :: tile_num_ch
    real(kind=kind_phys) :: sig1t, dt_warm
    integer              :: npts, nb, ix, jx, ls, ios, ll
    logical              :: exists
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     if (Model%me .eq. 0) print *,' nlats=',nlats,' lonsinpe='
!    *,lonsinpe(0,1)
!
      tile_num_ch = "      "
      if (tile_num < 10) then
        write(tile_num_ch, "(a4,i1)") "tile", tile_num
      else
        write(tile_num_ch, "(a4,i2)") "tile", tile_num
      endif
!
      sig1t = 0.0_kind_phys
      npts  = nx*ny
!
      if ( nsst > 0 ) then
        TSFFCS = tref
      else
        TSFFCS = tsfco
      endif
!
      do ix=1,npts
        if (landfrac(ix) > -1.0e-6_kind_phys) then
          slmskl(ix) = ceiling(landfrac(ix))
          slmskw(ix) = floor(landfrac(ix)+1.0e-6_kind_phys)
        endif

        if (lakefrac(ix) > 0.0_kind_phys) then
          min_ice(ix) = min_lakeice
        else
          min_ice(ix) = min_seaice
        endif

        zorfcs(ix) = zorll (ix)
        if (nint(slmskl(ix)) /= 1 ) then
          if (fice(ix) >= min_ice(ix)) then
            zorfcs(ix) = zorli(ix)
          else
            zorfcs(ix) = zorlo(ix)
          endif
        endif

        IF (fice(ix) >= min_ice(ix)) THEN
          AISFCS(ix) = 1.0_kind_phys
        ELSE
          AISFCS(ix) = 0.0_kind_phys
        ENDIF
!
        ALFFC1(ix         ) = facsf(ix)
        ALFFC1(ix + npts  ) = facwf(ix)
!
        ALBFC1(ix         ) = alvsf(ix)
        ALBFC1(ix + npts  ) = alvwf(ix)
        ALBFC1(ix + npts*2) = alnsf(ix)
        ALBFC1(ix + npts*3) = alnwf(ix)
!
        do ls = 1,max(lsoil,lsoil_lsm)
          ll = ix + (ls-1)*npts
          if (lsoil == lsoil_lsm) then
            SMCFC1(ll) = smc(ix,ls)
            STCFC1(ll) = stc(ix,ls)
            SLCFC1(ll) = slc(ix,ls)
          else
            SMCFC1(ll) = smois(ix,ls)
            STCFC1(ll) = tslb(ix,ls)
            SLCFC1(ll) = sh2o(ix,ls)
          endif
        enddo
      enddo
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
      CALL SFCCYCLE (9998, npts, max(lsoil,lsoil_lsm), sig1t, fhcyc, &
                     idate(4), idate(2), idate(3), idate(1),         &
                     phour, xlat_d, xlon_d, slmskl, slmskw,          &
                     oro, oro_uf, use_ufo, nst_anl,                  &
                     hice, fice, tisfc, snowd, slcfc1,               &
                     shdmin, shdmax, slope, snoalb, tsffcs,          &
                     weasd, zorfcs, albfc1, tg3, canopy,             &
                     smcfc1, stcfc1, slmsk, aisfcs,                  &
                     vfrac, vtype, stype, alffc1, cv,                &
                     cvb, cvt, me, nthrds,                           &
                     nlunit, size(input_nml_file), input_nml_file,   &
                     min_ice, ialb, isot, ivegsrc,                   &
                     trim(tile_num_ch), imap, jmap)
#ifndef INTERNAL_FILE_NML
      close (Model%nlunit)
#endif
!
      if ( nsst > 0 ) then
        tref = TSFFCS
      else
        tsfco = TSFFCS
      endif
!
      do ix=1,npts
        zorll(ix) = ZORFCS(ix)
        if (.not. frac_grid) then
          if (slmsk(ix) > 1.9_kind_phys) then
            zorli(ix) = ZORFCS(ix)
          elseif (slmsk(ix) < 0.1_kind_phys) then
            zorlo(ix) = ZORFCS(ix)
          endif
        else
          if (nint(slmskw(ix))  == 0 .and. nint(slmskl(ix)) /= 1) then
             if (fice(ix) >= min_ice(ix)) then
              zorli(ix) = ZORFCS(ix)
            else
              zorlo(ix) = ZORFCS(ix)
            endif
          endif
        endif
!
        facsf(ix) = ALFFC1(ix         )
        facwf(ix) = ALFFC1(ix + npts  )
!
        alvsf(ix) = ALBFC1(ix         )
        alvwf(ix) = ALBFC1(ix + npts  )
        alnsf(ix) = ALBFC1(ix + npts*2)
        alnwf(ix) = ALBFC1(ix + npts*3)
!
        do ls = 1,max(lsoil,lsoil_lsm)
          ll = ix + (ls-1)*npts
          if(lsoil == lsoil_lsm) then
            smc(ix,ls) = SMCFC1(ll)
            stc(ix,ls) = STCFC1(ll)
            slc(ix,ls) = SLCFC1(ll)
          else
            smois(ix,ls) = SMCFC1(ll)
            tslb(ix,ls)  = STCFC1(ll)
            sh2o(ix,ls)  = SLCFC1(ll)
          endif
!         if (ls <= kice) tiice(ix,ls) = STCFC1(ll)
        enddo
      enddo
!
!     if (Model%me .eq. 0) print*,'executed gcycle during hour=',fhour
!
      RETURN
      END

end module gcycle_mod
