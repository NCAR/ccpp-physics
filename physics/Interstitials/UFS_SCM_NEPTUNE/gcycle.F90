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
      stype, scolor, shdmin, shdmax, snowd, cv, cvb, cvt, oro, oro_uf,             &
      cplflx, oceanfrac,                                                           &
      xlat_d, xlon_d, slmsk, imap, jmap, errmsg, errflg)
!
!
    use machine,      only: kind_phys, kind_io8
    use sfccyc_module, only: sfccycle
    implicit none

    integer,              intent(in)    :: me, nthrds, nx, ny, isc, jsc, nsst, &
                                           tile_num, nlunit, lsoil, lsoil_lsm, kice
    integer,              intent(in)    :: idate(:), ialb, isot, ivegsrc
    character(len = 64), intent(in)     :: fn_nml
    character(len=*),     intent(in)    :: input_nml_file(:)
    logical,              intent(in)    :: use_ufo, nst_anl, frac_grid, cplflx
    real(kind=kind_phys), intent(in)    :: fhcyc, phour, landfrac(:), lakefrac(:), &
                                           min_seaice, min_lakeice,oceanfrac(:), &
                                           xlat_d(:), xlon_d(:)
    real(kind=kind_phys), intent(inout), optional ::   &
                                           smois(:,:), &
                                           sh2o(:,:),  &
                                           tslb(:,:),  &
                                           tref(:)
    real(kind=kind_phys), intent(inout) :: smc(:,:),   &
                                           slc(:,:),   &
                                           stc(:,:),   &
                                           tiice(:,:), &
                                           tg3(:),     &
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
                                           scolor(:),  &
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
        slcfcs (nx*ny),                      &               !soil color
        TSFFCS (nx*ny),                      &
        ZORFCS (nx*ny),                      &
        AISFCS (nx*ny),                      &
        ALFFC1 (nx*ny*2),                    &
        ALBFC1 (nx*ny*4),                    &
        SMCFC1 (nx*ny*max(lsoil,lsoil_lsm)), &
        STCFC1 (nx*ny*max(lsoil,lsoil_lsm)), &
        SLCFC1 (nx*ny*max(lsoil,lsoil_lsm))

!
! declare the variables (arrays) for cplflx, surface type dependent gcycle changes
!
    real(kind=kind_io8) ::                   &
        hice_save    (nx*ny),                &        ! sea or lake ice thickness
        fice_save    (nx*ny),                &        ! sea or lake ice fraction
        snowd_save   (nx*ny),                &        ! water equivalent snow depth
        snoalb_save  (nx*ny),                &        ! maximum snow albedo
        tisfc_save   (nx*ny),                &        ! surface skin temperature over (sea or lake) ice
        weasd_save   (nx*ny)                          ! water equiv of acc snow depth over land and (sea or lake)  ice


    real (kind=kind_io8) :: min_ice(nx*ny)
    integer              :: i_indx(nx*ny), j_indx(nx*ny)
    character(len=6)     :: tile_num_ch
    real(kind=kind_phys) :: sig1t(nx*ny)
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
! Some surface variables need to be updated by gcycle with coupled mode, and nsst mode dependent. A few variables are saved 
! in order to be able to update them over the specific surface types only after call sfccycle
!
      if ( cplflx ) then
        hice_save   = hice 
        fice_save   = fice
        snowd_save  = snowd
        snoalb_save = snoalb
        tisfc_save  = tisfc
        weasd_save  = weasd
      endif

      if ( nsst > 0 ) then
        TSFFCS = tref
      else
        TSFFCS = tsfco
      endif


! integer to real/double precision
      slpfcs = real(slope)
      vegfcs = real(vtype)
      sltfcs = real(stype)
      slcfcs = real(scolor)         !soil color
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
                     vfrac, vegfcs, sltfcs, slcfcs,alffc1, cv,       &   !slcfcs: soil color
                     cvb, cvt, me, nthrds,                           &
                     nlunit, size(input_nml_file), input_nml_file,   &
                     min_ice, ialb, isot, ivegsrc,                   &
                     trim(tile_num_ch), i_indx, j_indx)
#ifndef INTERNAL_FILE_NML
      close (nlunit)
#endif
!
! The gcycle resulted change is applied to some variables in the way of the coupled mode dependent, water surface type (ocean or lake)
! dependent and nsst mode dependent
!
      if ( cplflx ) then
!       In coupled mode, keep these variables the same as is (before sfccycle is called) over ocean
        do ix=1,npts
          if ( oceanfrac(ix) > 0.0_kind_phys ) then      
            hice(ix)   = hice_save(ix) 
            fice(ix)   = fice_save(ix)
            snowd(ix)  = snowd_save(ix)
            snoalb(ix) = snoalb_save(ix)
            tisfc(ix)  = tisfc_save(ix)
            weasd(ix)  = weasd_save(ix)
          endif
        enddo
!       In the coupled mode and when NSST is on, update tref over non-ocean 
        if ( nsst > 0 ) then       
          do ix=1,npts
            if ( oceanfrac(ix) == 0.0_kind_phys ) then 
              tref(ix)  = TSFFCS(ix) 
            endif
          enddo
!       In the coupled mode and when NSST is off, update tsfc and tsfco over non-ocean
        else             
          do ix=1,npts
            if ( oceanfrac(ix) == 0.0_kind_phys ) then 
              tsfc(ix)  = TSFFCS(ix)
              tsfco(ix) = TSFFCS(ix)
            endif
          enddo
        endif
!     The same as before (this modification) in uncoupled mode
      else
        if ( nsst > 0 ) then
          tref = TSFFCS
        else
          tsfc  = TSFFCS
          tsfco = TSFFCS
        endif
      endif
!
! real/double precision to integer
      slope = int(slpfcs)
      vtype = int(vegfcs)
      stype = int(sltfcs)
      scolor = int(slcfcs)  !soil color
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
