!>\file aerinterp.F90
!! This file contains subroutines of reading and interpolating
!! aerosol data for MG microphysics.

!>\ingroup mod_GFS_phys_time_vary
!! This module contain subroutines of reading and interpolating
!! aerosol data for MG microphysics.
module aerinterp

    implicit none

    private read_netfaer

    public :: read_aerdata, setindxaer, aerinterpol,read_aerdataf

contains

      logical function netcdf_check(status, errmsg, errflg, why)
      use netcdf
      implicit none
      character(len=*), intent(inout) :: errmsg
      integer, intent(out) :: errflg
      integer, intent(in) :: status
      character(len=*), intent(in) :: why

      netcdf_check = (status == NF90_NOERR)

      if(netcdf_check) then
        errflg = 0
        errmsg = ' '
      else
        errflg = 1
        errmsg = trim(why) // ': ' // trim(nf90_strerror(status))
      endif

      END function netcdf_check

      SUBROUTINE read_aerdata (me, master, iflip, idate, errmsg, errflg)
      use machine, only: kind_phys, kind_io4, kind_io8
      use aerclm_def
      use netcdf

!--- in/out
      integer, intent(in) :: me, master, iflip, idate(4)
      character(len=*), intent(inout) :: errmsg
      integer, intent(inout) :: errflg

!--- locals
      integer      :: ncid, varid, ndims, hmx
      integer      :: i, j, k, n, ii, imon, klev
      character    :: fname*50, mn*2, vname*10
      logical      :: file_exist
      integer :: dimids(NF90_MAX_VAR_DIMS)
      integer :: dimlen(NF90_MAX_VAR_DIMS)

      errflg = 0
      errmsg = ' '

!
!! ===================================================================
      if (me == master) then
         if ( iflip == 0 )  then             ! data from toa to sfc
          print *, "GFS is top-down"
         else
          print *, "GFS is bottom-up"
         endif
      endif
!
!! ===================================================================
!! check if all necessary files exist
!! ===================================================================
      do imon = 1, 12
         write(mn,'(i2.2)') imon
         fname=trim("aeroclim.m"//mn//".nc")
         inquire (file = fname, exist = file_exist)
         if (.not. file_exist) then
            errmsg = 'Error in read_aerdata: file ' // trim(fname) // ' not found'
            errflg = 1
            return
         endif
      enddo
!
!! ===================================================================
!! fetch dim spec and lat/lon from m01 data set
!! ===================================================================
      fname=trim("aeroclim.m"//'01'//".nc")
      ncid = -1
      if(.not.netcdf_check(nf90_open(fname , nf90_NOWRITE, ncid), &
           errmsg, errflg, 'open '//trim(fname))) then
        return
      endif

      vname =  trim(specname(1))
      varid = -1
      if(.not.netcdf_check(nf90_inq_varid(ncid, vname, varid), &
           errmsg, errflg, 'find id of '//trim(vname)//' var')) then
        return
      endif
      ndims = 0
      if(.not.netcdf_check(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), &
           errmsg, errflg, 'inquire details about '//trim(vname)//' var')) then
        return
      endif
      do i=1,ndims
        if(.not.netcdf_check(nf90_inquire_dimension(ncid, dimids(i), len=dimlen(i)), &
             errmsg, errflg, 'inquire details about dimension')) then
          return
        endif
      enddo

! specify latsaer, lonsaer, hmx
      lonsaer = dimlen(1)
      latsaer = dimlen(2)
      levsw = dimlen(3)

      if(me==master) then
         print *, 'MERRA2 dim: ',dimlen(1:ndims)
      endif

! allocate arrays

      if (.not. allocated(aer_lat)) then
        allocate(aer_lat(latsaer))
        allocate(aer_lon(lonsaer))
      endif

! construct lat/lon array
      varid = -1
      if(.not.netcdf_check(nf90_inq_varid(ncid, 'lat', varid), &
           errmsg, errflg, 'find id of lat var')) then
        return
      endif
      aer_lat = 0
      if(.not.netcdf_check(nf90_get_var(ncid, varid, aer_lat, (/ 1, 1, 1 /), (/latsaer, 1, 1/)), &
           errmsg, errflg, 'read lat var')) then
        return
      endif
      varid = -1
      if(.not.netcdf_check(nf90_inq_varid(ncid, 'lon', varid), &
           errmsg, errflg, 'find id of lon var')) then
        return
      endif
      aer_lon = 0
      if(.not.netcdf_check(nf90_get_var(ncid, varid, aer_lon, (/ 1, 1, 1 /), (/lonsaer, 1, 1/)), &
           errmsg, errflg, 'read lon var')) then
        return
      endif
      if(.not.netcdf_check(nf90_close(ncid), errmsg, errflg, 'close '//trim(fname))) then
        return
      endif
      END SUBROUTINE read_aerdata
!
!**********************************************************************
      SUBROUTINE read_aerdataf ( me, master, iflip, idate, FHOUR, errmsg, errflg)
      use machine, only: kind_phys, kind_dbl_prec
      use aerclm_def

!--- in/out
      integer, intent(in) :: me, master, iflip, idate(4)
      character(len=*), intent(inout) :: errmsg
      integer, intent(inout) :: errflg
      real(kind=kind_phys), intent(in) :: fhour
!--- locals
      integer      :: i, j, k, n, ii, imon, klev, n1, n2
      logical      :: file_exist
      integer  IDAT(8),JDAT(8)
      real(kind=kind_phys) rjday
      real(kind=kind_dbl_prec) rinc(5)
      integer jdow, jdoy, jday

      integer, allocatable  :: invardims(:)
!
      if (.not. allocated(aerin)) then
        allocate(aerin(iamin:iamax,jamin:jamax,levsaer,ntrcaerm,timeaer))
        allocate(aer_pres(iamin:iamax,jamin:jamax,levsaer,timeaer))
      endif

! allocate local working arrays
!!  found interpolation months
      IDAT = 0
      IDAT(1) = IDATE(4)
      IDAT(2) = IDATE(2)
      IDAT(3) = IDATE(3)
      IDAT(5) = IDATE(1)
      RINC = 0.
      RINC(2) = FHOUR
      CALL W3MOVDAT(RINC,IDAT,JDAT)
!
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow,jdoy,jday)
      rjday = jdoy + jdat(5) / 24.
      IF (RJDAY < aer_time(1)) RJDAY = RJDAY+365.
!
      n2 = 13
      do j=2, 12
       if (rjday < aer_time(j)) then
          n2 = j
          exit
       endif
      enddo
      n1 = n2 - 1
      if (n2 > 12) n2 = n2 -12
!! ===================================================================
      call read_netfaer(n1, iflip, 1, errmsg, errflg)
      if(errflg/=0) return
      call read_netfaer(n2, iflip, 2, errmsg, errflg)
      if(errflg/=0) return
!! ===================================================================
      n1sv=n1
      n2sv=n2
!---
      END SUBROUTINE read_aerdataf
!
      SUBROUTINE setindxaer(npts,dlat,jindx1,jindx2,ddy,dlon,           &
                            iindx1,iindx2,ddx,me,master)
!
      USE MACHINE,  ONLY: kind_phys
      use aerclm_def, only: aer_lat, jaero=>latsaer,                    &
                            aer_lon, iaero=>lonsaer
!
      implicit none
!
      integer me, master
      integer npts, JINDX1(npts),JINDX2(npts),IINDX1(npts),IINDX2(npts)
      real(kind=kind_phys) dlat(npts),DDY(npts),dlon(npts),DDX(npts)
!
      integer i,j

      DO J=1,npts
        jindx2(j) = jaero + 1
        do i=1,jaero
          if (dlat(j) < aer_lat(i)) then
            jindx2(j) = i
            exit
          endif
        enddo
        jindx1(j) = max(jindx2(j)-1,1)
        jindx2(j) = min(jindx2(j),jaero)
        if (jindx2(j) .ne. jindx1(j)) then
          DDY(j) = (dlat(j)            - aer_lat(jindx1(j))) &
                 / (aer_lat(jindx2(j)) - aer_lat(jindx1(j)))
        else
          ddy(j) = 1.0
        endif

      ENDDO

      DO J=1,npts
        iindx2(j) = iaero + 1
        do i=1,iaero
          if (dlon(j) < aer_lon(i)) then
            iindx2(j) = i
            exit
          endif
        enddo
        iindx1(j) = max(iindx2(j)-1,1)
        iindx2(j) = min(iindx2(j),iaero)
        if (iindx2(j) .ne. iindx1(j)) then
          ddx(j) = (dlon(j)            - aer_lon(iindx1(j))) &
                 / (aer_lon(iindx2(j)) - aer_lon(iindx1(j)))
        else
          ddx(j) = 1.0
        endif
      ENDDO

      RETURN
      END SUBROUTINE setindxaer
!
!**********************************************************************
!**********************************************************************
!
      SUBROUTINE aerinterpol( me,master,nthrds,npts,IDATE,FHOUR,iflip, jindx1,jindx2, &
                             ddy,iindx1,iindx2,ddx,lev,prsl,aerout, errmsg,errflg)
!
      use machine, only: kind_phys, kind_dbl_prec
      use aerclm_def

      implicit none
      integer, intent(inout) :: errflg
      character(*), intent(inout) :: errmsg
      integer, intent(in) :: iflip
      integer   i1,i2, iday,j,j1,j2,l,npts,nc,n1,n2,lev,k,i,ii, klev
      real(kind=kind_phys) fhour,temj, tx1, tx2,temi, tem
      real(kind=kind_phys), dimension(npts) :: temij,temiy,temjx,ddxy
      
!

      integer  JINDX1(npts), JINDX2(npts), iINDX1(npts), iINDX2(npts)
      integer  me,idate(4), master, nthrds
      integer  IDAT(8),JDAT(8)
!
      real(kind=kind_phys) DDY(npts), ddx(npts),ttt
      real(kind=kind_phys) aerout(npts,lev,ntrcaer)
      real(kind=kind_phys) aerpm(npts,levsaer,ntrcaer)
      real(kind=kind_phys) prsl(npts,lev), aerpres(npts,levsaer)
      real(kind=kind_phys) rjday
      real(kind=kind_dbl_prec) rinc(5)
      integer jdow, jdoy, jday
!
      errflg = 0
      errmsg = ' '
      IDAT = 0
      IDAT(1) = IDATE(4)
      IDAT(2) = IDATE(2)
      IDAT(3) = IDATE(3)
      IDAT(5) = IDATE(1)
      RINC = 0.
      RINC(2) = FHOUR
      CALL W3MOVDAT(RINC,IDAT,JDAT)
!
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow,jdoy,jday)
      rjday = jdoy + jdat(5) / 24.
      IF (RJDAY < aer_time(1)) RJDAY = RJDAY+365.
!
      n2 = 13
      do j=2, 12
       if (rjday < aer_time(j)) then
          n2 = j
          exit
       endif
      enddo
      n1 = n2 - 1
      if (n2 > 12) n2 = n2 -12
!     need to read a new month 
      if (n1.ne.n1sv) then
#ifdef DEBUG
        if (me == master) write(*,*)"read in a new month MERRA2", n2
#endif
        DO ii = 1, ntrcaerm
          do j = jamin, jamax
            do k = 1, levsaer
              do i = iamin, iamax
                aerin(i,j,k,ii,1) = aerin(i,j,k,ii,2)
              enddo   !i-loop (lon)
            enddo     !k-loop (lev)
          enddo       !j-loop (lat)
        ENDDO         ! ii-loop (ntracaerm)
!! ===================================================================
        call read_netfaer(n2, iflip, 2, errmsg, errflg)
        n1sv=n1
        n2sv=n2
      end if
!
      tx1 = (aer_time(n2) - rjday) / (aer_time(n2) - aer_time(n1))
      tx2 = 1.0 - tx1
      if (n2 > 12) n2 = n2 -12

      do j=1,npts
         TEMJ     = 1.0 - DDY(J)
         TEMI     = 1.0 - DDX(J)
         temij(j) = TEMI*TEMJ
         temiy(j) = TEMI*DDY(j)
         temjx(j) = TEMJ*DDX(j)
         ddxy(j)  = DDX(j)*DDY(J)
      enddo

#ifndef __GFORTRAN__
!$OMP parallel num_threads(nthrds) default(none)             &
!$OMP          shared(npts,ntrcaer,aerin,aer_pres,prsl)      &
!$OMP          shared(ddx,ddy,jindx1,jindx2,iindx1,iindx2)   &
!$OMP          shared(aerpm,aerpres,aerout,lev,nthrds) &
!$OMP          shared(temij,temiy,temjx,ddxy)                &
!$OMP          private(l,j,k,ii,i1,i2,j1,j2,tem)             &
!$OMP          copyin(tx1,tx2) firstprivate(tx1,tx2)

!$OMP do
#endif
      DO L=1,levsaer
        DO J=1,npts
          J1    = JINDX1(J)
          J2    = JINDX2(J)
          I1    = IINDX1(J)
          I2    = IINDX2(J)
          DO ii=1,ntrcaer
           aerpm(j,L,ii) =                                                  &
           tx1*(TEMIJ(j)*aerin(I1,J1,L,ii,1)+DDXY(j)*aerin(I2,J2,L,ii,1)  &
               +TEMIY(j)*aerin(I1,J2,L,ii,1)+temjx(j)*aerin(I2,J1,L,ii,1))&
          +tx2*(TEMIJ(j)*aerin(I1,J1,L,ii,2)+DDXY(j)*aerin(I2,J2,L,ii,2)  &
               +TEMIY(j)*aerin(I1,J2,L,ii,2)+temjx(j)*aerin(I2,J1,L,ii,2))
          ENDDO

          aerpres(j,L) =                                                    &
           tx1*(TEMIJ(j)*aer_pres(I1,J1,L,1)+DDXY(j)*aer_pres(I2,J2,L,1)  &
               +TEMIY(j)*aer_pres(I1,J2,L,1)+temjx(j)*aer_pres(I2,J1,L,1))&
          +tx2*(TEMIJ(j)*aer_pres(I1,J1,L,2)+DDXY(j)*aer_pres(I2,J2,L,2)  &
               +TEMIY(j)*aer_pres(I1,J2,L,2)+temjx(j)*aer_pres(I2,J1,L,2))
        ENDDO
      ENDDO
#ifndef __GFORTRAN__
!$OMP end do

! don't flip, input is the same direction as GFS  (bottom-up)
!$OMP do
#endif
      DO J=1,npts
        DO L=1,lev
           if(prsl(j,L) >= aerpres(j,1)) then
              DO ii=1, ntrcaer
               aerout(j,L,ii) = aerpm(j,1,ii)        !! sfc level
              ENDDO
           else if(prsl(j,L) <= aerpres(j,levsaer)) then
              DO ii=1, ntrcaer
               aerout(j,L,ii) = aerpm(j,levsaer,ii)  !! toa top
              ENDDO
           else
             DO  k=1, levsaer-1      !! from sfc to toa
              IF(prsl(j,L) <= aerpres(j,k) .and. prsl(j,L)>aerpres(j,k+1)) then
                 i1 = k
                 i2 = min(k+1,levsaer)
                 exit
              ENDIF
             ENDDO
             tem  = 1.0 / (aerpres(j,i1) - aerpres(j,i2))
             tx1  = (prsl(j,L) - aerpres(j,i2)) * tem
             tx2  = (aerpres(j,i1) - prsl(j,L)) * tem
             DO ii = 1, ntrcaer
               aerout(j,L,ii) = aerpm(j,i1,ii)*tx1 + aerpm(j,i2,ii)*tx2
             ENDDO
           endif
        ENDDO   !L-loop
      ENDDO     !J-loop
#ifndef __GFORTRAN__
!$OMP end do

!$OMP end parallel
#endif

      RETURN
      END SUBROUTINE aerinterpol

      subroutine read_netfaer(nf, iflip,nt, errmsg,errflg)
      use machine, only: kind_phys, kind_io4, kind_io8
      use aerclm_def
      use netcdf
      integer, intent(in) :: iflip, nf, nt
      integer, intent(inout) :: errflg
      character(*), intent(inout) :: errmsg
      integer      :: ncid, varid, i,j,k,ii,klev
      character    :: fname*50, mn*2, vname*10
      real(kind=kind_io4),allocatable,dimension(:,:,:) :: buff
      real(kind=kind_io4),allocatable,dimension(:,:,:,:):: buffx
      real(kind=kind_io4),allocatable,dimension(:,:)   :: pres_tmp
      
!! ===================================================================
      allocate (buff(lonsaer, latsaer, levsw))
      allocate (pres_tmp(lonsaer, levsw))
      allocate (buffx(lonsaer, latsaer, levsw, 1))

      errflg = 0
      errmsg = ' '
      buff = 0
      pres_tmp = 0
      buffx = 0

      write(mn,'(i2.2)') nf 
      fname=trim("aeroclim.m"//mn//".nc")
      ncid = -1
      if(.not.netcdf_check(nf90_open(fname , nf90_NOWRITE, ncid), &
           errmsg, errflg, 'open '//trim(fname))) then
        return
      endif

! ====> construct 3-d pressure array (Pa)
      varid = -1
      if(.not.netcdf_check(nf90_inq_varid(ncid, "DELP", varid), &
           errmsg, errflg, 'find id of DELP var')) then
        return
      endif
      if(.not.netcdf_check(nf90_get_var(ncid, varid, buff), &
           errmsg, errflg, 'read DELP var')) then
        return
      endif

      do j = jamin, jamax
        do i = iamin, iamax
! constract pres_tmp (top-down), note input is top-down
          pres_tmp(i,1) = 0.
          do k=2, levsw
            pres_tmp(i,k) = pres_tmp(i,k-1)+buff(i,j,k)
          enddo    !k-loop
        enddo     !i-loop (lon)

! extract pres_tmp to fill aer_pres (in  Pa)
        do k = 1, levsaer
          if ( iflip == 0 )  then             ! data from toa to sfc
            klev = k
          else                                ! data from sfc to top
            klev = ( levsw - k ) + 1
          endif
          do i = iamin, iamax
            aer_pres(i,j,k,nt)    = 1.d0*pres_tmp(i,klev)
          enddo     !i-loop (lon)
        enddo     !k-loop (lev)
      enddo     !j-loop (lat)

! ====> construct 4-d aerosol array (kg/kg)
! merra2 data is top down
! for GFS, iflip 0: toa to sfc; 1: sfc to toa
      DO ii = 1, ntrcaerm
        vname=trim(specname(ii))
        varid = -1
        if(.not.netcdf_check(nf90_inq_varid(ncid, vname, varid), &
             errmsg, errflg, 'get id of '//trim(vname)//' var')) then
          return
        endif
        if(.not.netcdf_check(nf90_get_var(ncid, varid, buffx), &
             errmsg, errflg, 'read '//trim(vname)//' var')) then
          return
        endif

        do j = jamin, jamax
          do k = 1, levsaer
! input is from toa to sfc
            if ( iflip == 0 )  then             ! data from toa to sfc
              klev = k
            else                                ! data from sfc to top
              klev = ( levsw - k ) + 1
            endif
            do i = iamin, iamax
              aerin(i,j,k,ii,nt) = 1.d0*buffx(i,j,klev,1)
              if(aerin(i,j,k,ii,nt) < 0 .or. aerin(i,j,k,ii,nt) > 1.)  then
                aerin(i,j,k,ii,nt) = 1.e-15
              endif
            enddo   !i-loop (lon)
          enddo     !k-loop (lev)
        enddo       !j-loop (lat)

      ENDDO         ! ii-loop (ntracaerm)

! close the file
      if(.not.netcdf_check(nf90_close(ncid), errmsg, errflg, 'close '//trim(fname))) then
        return
      endif
      deallocate (buff, pres_tmp)
      deallocate (buffx)
      END SUBROUTINE read_netfaer

end module aerinterp

