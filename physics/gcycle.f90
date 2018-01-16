  SUBROUTINE GCYCLE (blocksize, Model, Grid, Sfcprop, Cldprop, Sfccycle)
!
!
    USE MACHINE,      only: kind_phys
    USE PHYSCONS,     only: PI => con_PI
    USE GFS_typedefs, only: GFS_control_type, GFS_grid_type, &
                            GFS_sfcprop_type, GFS_cldprop_type, &
                            GFS_sfccycle_type
    implicit none

    integer,                  intent(in)    :: blocksize
    type(GFS_control_type),   intent(in)    :: Model
    type(GFS_grid_type),      intent(in)    :: Grid
    type(GFS_sfcprop_type),   intent(inout) :: Sfcprop
    type(GFS_cldprop_type),   intent(inout) :: Cldprop
    type(GFS_sfccycle_type),  intent(inout) :: Sfccycle

!
!     Local variables
!     ---------------
    real(kind=kind_phys) ::             &
           RLA (blocksize),             &
           RLO (blocksize),             &
        SLMASK (blocksize),             &
          OROG (blocksize),             &
       OROG_UF (blocksize),             &
        SLIFCS (blocksize),             &
        TSFFCS (blocksize),             &
        SNOFCS (blocksize),             &
        ZORFCS (blocksize),             &
        TG3FCS (blocksize),             &
        CNPFCS (blocksize),             &
        AISFCS (blocksize),             &
        F10MFCS(blocksize),             &
        VEGFCS (blocksize),             &
        VETFCS (blocksize),             &
        SOTFCS (blocksize),             &
         CVFCS (blocksize),             &
        CVBFCS (blocksize),             &
        CVTFCS (blocksize),             &
        SWDFCS (blocksize),             &
        SIHFCS (blocksize),             &
        SICFCS (blocksize),             &
        SITFCS (blocksize),             &
        VMNFCS (blocksize),             &
        VMXFCS (blocksize),             &
        SLPFCS (blocksize),             &
        ABSFCS (blocksize),             &
        ALFFC1 (blocksize*2),           &
        ALBFC1 (blocksize*4),           &
        SMCFC1 (blocksize*Model%lsoil), &
        STCFC1 (blocksize*Model%lsoil), &
        SLCFC1 (blocksize*Model%lsoil), &
        GLACIR (blocksize),             &
        AMXICE (blocksize),             &
        TSFCL0 (blocksize)

    real(kind=kind_phys)    :: sig1t, pifac
    integer :: npts, ix, ls, ios
    logical :: exists
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     if (Model%me .eq. 0) print *,' nlats=',nlats,' lonsinpe='
!    *,lonsinpe(0,1)

      sig1t = 0.0
      npts  = blocksize
!
      pifac = 180.0 / pi

!$OMP parallel do default (shared) private (ix,ls)
      do ix = 1,size(Grid%xlat,1)
        RLA     (ix)          =    Grid%xlat   (ix) * pifac
        RLO     (ix)          =    Grid%xlon   (ix) * pifac
        OROG    (ix)          = Sfcprop%oro    (ix)
        OROG_UF (ix)          = Sfcprop%oro_uf (ix)
        SLIFCS  (ix)          = Sfcprop%slmsk  (ix)
        TSFFCS  (ix)          = Sfcprop%tsfc   (ix)
        SNOFCS  (ix)          = Sfcprop%weasd  (ix)
        ZORFCS  (ix)          = Sfcprop%zorl   (ix)
        TG3FCS  (ix)          = Sfcprop%tg3    (ix)
        CNPFCS  (ix)          = Sfcprop%canopy (ix)
        F10MFCS (ix)          = Sfcprop%f10m   (ix)
        VEGFCS  (ix)          = Sfcprop%vfrac  (ix)
        VETFCS  (ix)          = Sfcprop%vtype  (ix)
        SOTFCS  (ix)          = Sfcprop%stype  (ix)
        CVFCS   (ix)          = Cldprop%cv     (ix)
        CVBFCS  (ix)          = Cldprop%cvb    (ix)
        CVTFCS  (ix)          = Cldprop%cvt    (ix)
        SWDFCS  (ix)          = Sfcprop%snowd  (ix)
        SIHFCS  (ix)          = Sfcprop%hice   (ix)
        SICFCS  (ix)          = Sfcprop%fice   (ix)
        SITFCS  (ix)          = Sfcprop%tisfc  (ix)
        VMNFCS  (ix)          = Sfcprop%shdmin (ix)
        VMXFCS  (ix)          = Sfcprop%shdmax (ix)
        SLPFCS  (ix)          = Sfcprop%slope  (ix)
        ABSFCS  (ix)          = Sfcprop%snoalb (ix)

        ALFFC1  (ix       )   = Sfcprop%facsf  (ix)
        ALFFC1  (ix + npts)   = Sfcprop%facwf  (ix)

        ALBFC1  (ix         ) = Sfcprop%alvsf  (ix)
        ALBFC1  (ix + npts  ) = Sfcprop%alvwf  (ix)
        ALBFC1  (ix + npts*2) = Sfcprop%alnsf  (ix)
        ALBFC1  (ix + npts*3) = Sfcprop%alnwf  (ix)

        do ls = 1,Model%lsoil
          SMCFC1 (ix + (ls-1)*npts) = Sfcprop%smc (ix,ls)
          STCFC1 (ix + (ls-1)*npts) = Sfcprop%stc (ix,ls)
          SLCFC1 (ix + (ls-1)*npts) = Sfcprop%slc (ix,ls)
        enddo

        IF (SLIFCS(ix) .LT. 0.1 .OR. SLIFCS(ix) .GT. 1.5) THEN
           SLMASK(ix) = 0
        ELSE
           SLMASK(ix) = 1
        ENDIF

        IF (SLIFCS(ix) .EQ. 2) THEN
          AISFCS(ix) = 1.
        ELSE
          AISFCS(ix) = 0.
        ENDIF

        GLACIR  (ix)          = Sfccycle%glacir (ix)
        AMXICE  (ix)          = Sfccycle%amxice (ix)
        TSFCL0  (ix)          = Sfccycle%tsfcl0 (ix)

  !     if (Model%me .eq. 0)
  !    &   print *,' ix=',ix,' rla=',rla(ix),' rlo=',rlo(ix)
      enddo
!$OMP end parallel do

! check
!     call mymaxmin(slifcs,ix,ix,1,'slifcs')
!     call mymaxmin(slmask,ix,ix,1,'slmsk')
!
      inquire (file=trim(Model%fn_nml),exist=exists)
      if (.not. exists) then
        write(6,*) 'gcycle:: namelist file: ',trim(Model%fn_nml),' does not exist'
        stop
      else
        open (unit=Model%nlunit, file=trim(Model%fn_nml), action='READ', status='OLD', iostat=ios)
      endif
      CALL SFCCYCLE_SUB (9998, npts, Model%lsoil, SIG1T, Model%fhcyc, &
                     Model%idate(4), Model%idate(2),                  &
                     Model%idate(3), Model%idate(1),                  &
                     Model%fhour, RLA, RLO, SLMASK,                   &
                     OROG, OROG_UF, Model%USE_UFO, Model%nst_anl,     &
                     SIHFCS, SICFCS, SITFCS, SWDFCS, SLCFC1,          &
                     VMNFCS, VMXFCS, SLPFCS, ABSFCS, TSFFCS,          &
                     SNOFCS, ZORFCS, ALBFC1, TG3FCS, CNPFCS,          &
                     SMCFC1, STCFC1, SLIFCS, AISFCS, F10MFCS,         &
                     VEGFCS, VETFCS, SOTFCS, ALFFC1, CVFCS,           &
                     CVBFCS, CVTFCS,                                  &
                     GLACIR, AMXICE, TSFCL0,                          &
                     Model%me, Model%nlunit,                          &
                     Model%ialb,                                      &
                     Sfccycle%ifp, Sfccycle%clima,                    &
                     Model%isot, Model%ivegsrc)
      close (Model%nlunit)

!$OMP parallel do default (shared) private (ix,ls)
      do ix = 1,size(Grid%xlat,1)
        Sfcprop%slmsk  (ix) = SLIFCS  (ix)
        Sfcprop%tsfc   (ix) = TSFFCS  (ix)
        Sfcprop%weasd  (ix) = SNOFCS  (ix)
        Sfcprop%zorl   (ix) = ZORFCS  (ix)
        Sfcprop%tg3    (ix) = TG3FCS  (ix)
        Sfcprop%canopy (ix) = CNPFCS  (ix)
        Sfcprop%f10m   (ix) = F10MFCS (ix)
        Sfcprop%vfrac  (ix) = VEGFCS  (ix)
        Sfcprop%vtype  (ix) = VETFCS  (ix)
        Sfcprop%stype  (ix) = SOTFCS  (ix)
        Cldprop%cv     (ix) = CVFCS   (ix)
        Cldprop%cvb    (ix) = CVBFCS  (ix)
        Cldprop%cvt    (ix) = CVTFCS  (ix)
        Sfcprop%snowd  (ix) = SWDFCS  (ix)
        Sfcprop%hice   (ix) = SIHFCS  (ix)
        Sfcprop%fice   (ix) = SICFCS  (ix)
        Sfcprop%tisfc  (ix) = SITFCS  (ix)
        Sfcprop%shdmin (ix) = VMNFCS  (ix)
        Sfcprop%shdmax (ix) = VMXFCS  (ix)
        Sfcprop%slope  (ix) = SLPFCS  (ix)
        Sfcprop%snoalb (ix) = ABSFCS  (ix)

        Sfcprop%facsf  (ix) = ALFFC1  (ix       )
        Sfcprop%facwf  (ix) = ALFFC1  (ix + npts)

        Sfcprop%alvsf  (ix) = ALBFC1  (ix         )
        Sfcprop%alvwf  (ix) = ALBFC1  (ix + npts  )
        Sfcprop%alnsf  (ix) = ALBFC1  (ix + npts*2)
        Sfcprop%alnwf  (ix) = ALBFC1  (ix + npts*3)
        do ls = 1,Model%lsoil
          Sfcprop%smc (ix,ls) = SMCFC1 (ix + (ls-1)*npts)
          Sfcprop%stc (ix,ls) = STCFC1 (ix + (ls-1)*npts)
          Sfcprop%slc (ix,ls) = SLCFC1 (ix + (ls-1)*npts)
        enddo

        Sfccycle%glacir (ix) = GLACIR  (ix)
        Sfccycle%amxice (ix) = AMXICE  (ix)
        Sfccycle%tsfcl0 (ix) = TSFCL0  (ix)

      enddo
!$OMP end parallel do

! check
!     call mymaxmin(slifcs,ix,ix,1,'slifcs')
!
!     if (Model%me .eq. 0) print*,'executed gcycle during hour=',fhour
      
      RETURN
      END
