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
  subroutine gcycle (me, nthrds, nx, ny, isc, jsc, nsst, tile_num, nlunit, &
      input_nml_file, lsoil, lsoil_lsm, kice, idate, ialb, isot, ivegsrc,  &
      use_ufo, nst_anl, fhcyc, phour, lakefrac, min_seaice, min_lakeice,   &
      frac_grid, smc, slc, stc, smois, sh2o, tslb, tiice, tg3, tref, tsfc, &
      tsfco, tisfc, hice, fice, facsf, facwf, alvsf, alvwf, alnsf, alnwf,  &
      zorli, zorll, zorlo, weasd, slope, snoalb, canopy, vfrac, vtype,     &
      stype, shdmin, shdmax, snowd, cv, cvb, cvt, oro, oro_uf,             &
      xlat_d, xlon_d, slmsk, imap, jmap)
!
!
    use machine,      only: kind_phys
    implicit none

    integer,              intent(in)    :: me, nthrds, nx, ny, isc, jsc, nsst, &
                                           tile_num, nlunit, lsoil, lsoil_lsm, kice
    integer,              intent(in)    :: idate(:), ialb, isot, ivegsrc
    character(len=*),     intent(in)    :: input_nml_file(:)
    logical,              intent(in)    :: use_ufo, nst_anl, frac_grid
    real(kind=kind_phys), intent(in)    :: fhcyc, phour, lakefrac(:), &
                                           min_seaice, min_lakeice,   &
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
        SLMASK (nx*ny),                      &
        TSFFCS (nx*ny),                      &
        ZORFCS (nx*ny),                      &
        AISFCS (nx*ny),                      &
        ALFFC1 (nx*ny*2),                    &
        ALBFC1 (nx*ny*4),                    &
        SMCFC1 (nx*ny*max(lsoil,lsoil_lsm)), &
        STCFC1 (nx*ny*max(lsoil,lsoil_lsm)), &
        SLCFC1 (nx*ny*max(lsoil,lsoil_lsm))


    logical              :: lake(nx*ny)
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
        TSFFCS = tsfc
      end if
!
      do ix=1,npts
        ZORFCS(ix) = zorll (ix)
        if (slmsk(ix) > 1.9_kind_phys .and. .not. frac_grid) then
          ZORFCS(ix) = zorli  (ix)
        elseif (slmsk(ix) < 0.1_kind_phys .and. .not. frac_grid) then
          ZORFCS(ix) = zorlo  (ix)
        endif
        ! DH* Why not 1.9 as for ZORFCS?
        IF (slmsk(ix) > 1.99_kind_phys) THEN
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
        !
        IF (slmsk(ix) < 0.1_kind_phys .OR. slmsk(ix) > 1.5_kind_phys) THEN
           SLMASK(ix) = 0.0_kind_phys
        ELSE
           SLMASK(ix) = 1.0_kind_phys
        ENDIF
        !
        if (lakefrac(ix) > 0.0_kind_phys) then
          lake(ix) = .true.
        else
          lake(ix) = .false.
        endif
      end do
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
                     phour, xlat_d, xlon_d, slmask,                  &
                     oro, oro_uf, use_ufo, nst_anl,                  &
                     hice, fice, tisfc, snowd, slcfc1,               &
                     shdmin, shdmax, slope, snoalb, tsffcs,          &
                     weasd, zorfcs, albfc1, tg3, canopy,             &
                     smcfc1, stcfc1, slmsk, aisfcs,                  &
                     vfrac, vtype, stype, alffc1, cv,                &
                     cvb, cvt, me, nthrds,                           &
                     nlunit, size(input_nml_file), input_nml_file,   &
                     lake, min_lakeice, min_seaice,                  &
                     ialb, isot, ivegsrc,                            &
                     trim(tile_num_ch), imap, jmap)
#ifndef INTERNAL_FILE_NML
      close (Model%nlunit)
#endif
!
      if ( nsst > 0 ) then
        tref = TSFFCS
      else
        tsfc  = TSFFCS
        tsfco = TSFFCS
      end if
!
      do ix=1,npts
        zorll(ix) = ZORFCS(ix)
        if (slmsk(ix) > 1.9_kind_phys .and. .not. frac_grid) then
          zorli(ix) = ZORFCS(ix)
        elseif (slmsk(ix) < 0.1_kind_phys .and. .not. frac_grid) then
          zorlo(ix) = ZORFCS(ix)
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
            tslb(ix,ls) = STCFC1(ll)
            sh2o(ix,ls) = SLCFC1(ll)
          endif
          if (ls<=kice) tiice(ix,ls) = STCFC1(ll)
        enddo
      enddo
!
!     if (Model%me .eq. 0) print*,'executed gcycle during hour=',fhour
!
      RETURN
      END

end module gcycle_mod
