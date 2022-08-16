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
  subroutine gcycle (me, nthrds, nx, ny, isc, jsc, nsst, tile_num, nlunit, fn_nml, &
      input_nml_file, lsoil, lsoil_lsm, kice, idate, ialb, isot, ivegsrc,          &
      use_ufo, nst_anl, fhcyc, phour, landfrac, lakefrac, min_seaice, min_lakeice, &
      frac_grid, smc, slc, stc, smois, sh2o, tslb, tiice, tg3, tref, tsfc,         &
      tsfco, tisfc, hice, fice, facsf, facwf, alvsf, alvwf, alnsf, alnwf,          &
      zorli, zorll, zorlo, weasd, slope, snoalb, canopy, vfrac, vtype,             &
      stype, shdmin, shdmax, snowd, cv, cvb, cvt, oro, oro_uf,                     &
      xlat_d, xlon_d, slmsk, imap, jmap, errmsg, errflg)
!
!
    use machine,      only: kind_phys, kind_io8
    implicit none

    integer,              intent(in)    :: me, nthrds, nx, ny, isc, jsc, nsst, &
                                           tile_num, nlunit, lsoil, lsoil_lsm, kice
    integer,              intent(in)    :: idate(:), ialb, isot, ivegsrc
    character(len = 64), intent(in)     :: fn_nml
    character(len=*),     intent(in)    :: input_nml_file(:)
    logical,              intent(in)    :: use_ufo, nst_anl, frac_grid
    real(kind=kind_phys), intent(in)    :: fhcyc, phour, landfrac(:), lakefrac(:), &
                                           min_seaice, min_lakeice,                &
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
                                           snoalb(:),  &
                                           canopy(:),  &
                                           vfrac(:),   &
                                           shdmin(:),  &
                                           shdmax(:),  &
                                           snowd(:),   &
                                           cv(:),      &
                                           cvb(:),     &
                                           cvt(:),     &
                                           oro(:),     &
                                           oro_uf(:),  &
                                           slmsk(:)
    integer,              intent(inout) :: vtype(:),   &
                                           stype(:),   &
                                           slope(:)

    integer,              intent(in)    :: imap(:), jmap(:)
    character(len=*),     intent(out)   :: errmsg
    integer,              intent(out)   :: errflg

!
!     Local variables
!     ---------------
!   real(kind=kind_phys) ::                  &
    real(kind=kind_io8) ::                   &
        slmskl (nx*ny),                      &
        slmskw (nx*ny),                      &
        slpfcs (nx*ny),                      &
        vegfcs (nx*ny),                      &
        sltfcs (nx*ny),                      &
        TSFFCS (nx*ny),                      &
        ZORFCS (nx*ny),                      &
        AISFCS (nx*ny),                      &
        ALFFC1 (nx*ny*2),                    &
        ALBFC1 (nx*ny*4),                    &
        SMCFC1 (nx*ny*max(lsoil,lsoil_lsm)), &
        STCFC1 (nx*ny*max(lsoil,lsoil_lsm)), &
        SLCFC1 (nx*ny*max(lsoil,lsoil_lsm))


    real (kind=kind_io8) :: min_ice(nx*ny)
    integer              :: i_indx(nx*ny), j_indx(nx*ny)
    character(len=6)     :: tile_num_ch
    real(kind=kind_phys) :: sig1t
    integer              :: npts, nb, ix, jx, ls, ios, ll
    logical              :: exists

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

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
      end if
! integer to real/double precision
      slpfcs = real(slope)
      vegfcs = real(vtype)
      sltfcs = real(stype)
!
      if (frac_grid) then
        do ix=1,npts
!         if (landfrac(ix) > -1.0e-8_kind_phys) then
          if (landfrac(ix) > 0.0_kind_phys) then
            slmskl(ix) = ceiling(landfrac(ix)-1.0e-6_kind_phys)
            slmskw(ix) =   floor(landfrac(ix)+1.0e-6_kind_phys)
          else
            if (nint(slmsk(ix)) == 1) then
              slmskl(ix) = 1.0_kind_phys
              slmskw(ix) = 1.0_kind_phys
            else
              slmskl(ix) = 0.0_kind_phys
              slmskw(ix) = 0.0_kind_phys
            endif
          endif
          ZORFCS(ix) = zorll(ix)
          if (nint(slmskl(ix)) == 0) then
            if (slmsk(ix) > 1.99_kind_phys) then
              ZORFCS(ix) = zorli(ix)
            else
              ZORFCS(ix) = zorlo(ix)
            endif
          endif
        enddo
      else
        do ix=1,npts
          if (nint(slmsk(ix)) == 1) then
            slmskl(ix) = 1.0_kind_phys
            slmskw(ix) = 1.0_kind_phys
          else
            slmskl(ix) = 0.0_kind_phys
            slmskw(ix) = 0.0_kind_phys
          endif
          ZORFCS(ix) = zorll(ix)
          if (slmsk(ix) > 1.99_kind_phys) then
            ZORFCS(ix) = zorli(ix)
          elseif (slmsk(ix) < 0.1_kind_phys) then
            ZORFCS(ix) = zorlo(ix)
          endif
        enddo
      endif
      do ix=1,npts
        i_indx(ix) = imap(ix) + isc - 1
        j_indx(ix) = jmap(ix) + jsc - 1

        if (lakefrac(ix) > 0.0_kind_phys) then
          min_ice(ix) = min_lakeice
        else
          min_ice(ix) = min_seaice
        endif

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
      enddo
!
#ifndef INTERNAL_FILE_NML
      inquire (file=trim(fn_nml),exist=exists)
      if (.not. exists) then
        write(6,*) 'gcycle:: namelist file: ',trim(fn_nml),' does not exist'
        errflg = 1
        errmsg = 'ERROR(gcycle): namelist file: ',trim(fn_nml),' does not exist.'
        return
      else
        open (unit=nlunit, file=trim(fn_nml), action='READ', status='OLD', iostat=ios)
        rewind (nlunit)
      endif
#endif
      CALL SFCCYCLE (9998, npts, max(lsoil,lsoil_lsm), sig1t, fhcyc, &
                     idate(4), idate(2), idate(3), idate(1),         &
                     phour, xlat_d, xlon_d, slmskl, slmskw,          &
                     oro, oro_uf, use_ufo, nst_anl,                  &
                     hice, fice, tisfc, snowd, slcfc1,               &
                     shdmin, shdmax, slpfcs, snoalb, tsffcs,         &
                     weasd, zorfcs, albfc1, tg3, canopy,             &
                     smcfc1, stcfc1, slmsk, aisfcs,                  &
                     vfrac, vegfcs, sltfcs, alffc1, cv,              &
                     cvb, cvt, me, nthrds,                           &
                     nlunit, size(input_nml_file), input_nml_file,   &
                     min_ice, ialb, isot, ivegsrc,                   &
                     trim(tile_num_ch), i_indx, j_indx)
#ifndef INTERNAL_FILE_NML
      close (nlunit)
#endif
!
      if ( nsst > 0 ) then
        tref = TSFFCS
      else
        tsfc  = TSFFCS
        tsfco = TSFFCS
      endif
!
! real/double precision to integer
      slope = int(slpfcs)
      vtype = int(vegfcs)
      stype = int(sltfcs)
!
      do ix=1,npts
        zorll(ix) = ZORFCS(ix)
        if (nint(slmskl(ix)) == 0) then
          if (slmsk(ix) > 1.99_kind_phys) then
            zorli(ix) = ZORFCS(ix)
          elseif (slmsk(ix) < 0.1_kind_phys) then
            zorlo(ix) = ZORFCS(ix)
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
            tslb(ix,ls) = STCFC1(ll)
            sh2o(ix,ls) = SLCFC1(ll)
          endif
!         if (ls<=kice) tiice(ix,ls) = STCFC1(ll)
        enddo
      enddo
!
!     if (Model%me .eq. 0) print*,'executed gcycle during hour=',fhour
!
      RETURN
      END

end module gcycle_mod
