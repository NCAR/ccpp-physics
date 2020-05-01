!>\file aerinterp.F90
!! This file contains subroutines of reading and interpolating
!! aerosol data for MG microphysics.

!>\ingroup mod_GFS_phys_time_vary
!! This module contain subroutines of reading and interpolating
!! aerosol data for MG microphysics.
module aerinterp

    implicit none

    private

    public :: read_aerdata, setindxaer, aerinterpol

contains

      SUBROUTINE read_aerdata (me, master, iflip, idate, errmsg, errflg)
      use machine, only: kind_phys, kind_io4, kind_io8
      use aerclm_def
      use netcdf

!--- in/out
      integer, intent(in) :: me, master, iflip, idate(4)
      character(len=*), intent(inout) :: errmsg
      integer, intent(inout) :: errflg

!--- locals
      integer      :: ncid, varid, ndims, dim1, dim2, dim3, hmx
      integer      :: i, j, k, n, ii, imon, klev
      character    :: fname*50, mn*2, vname*10
      logical      :: file_exist

      integer, allocatable  :: invardims(:)
      real(kind=kind_io4),allocatable,dimension(:,:,:) :: buff
      real(kind=kind_io4),allocatable,dimension(:,:,:,:):: buffx
      real(kind=kind_io4),allocatable,dimension(:,:)   :: pres_tmp
      real(kind=kind_io8),allocatable,dimension(:)     :: aer_lati
      real(kind=kind_io8),allocatable,dimension(:)     :: aer_loni
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
!! fetch dim spec and lat/lon from m01 data set
!! ===================================================================
      fname=trim("aeroclim.m"//'01'//".nc")
      inquire (file = fname, exist = file_exist)
      if (.not. file_exist) then
         errmsg = 'Error in read_aerdata: file ' // trim(fname) // ' not found'
         errflg = 1
         return
      endif
      call nf_open(fname , nf90_NOWRITE, ncid)

      vname =  trim(specname(1))
      call nf_inq_varid(ncid, vname, varid)
      call nf_inq_varndims(ncid, varid, ndims)

      if(.not. allocated(invardims)) allocate(invardims(3))
      call nf_inq_vardimid(ncid,varid,invardims)
      call nf_inq_dimlen(ncid, invardims(1), dim1)
      call nf_inq_dimlen(ncid, invardims(2), dim2)
      call nf_inq_dimlen(ncid, invardims(3), dim3)

! specify latsaer, lonsaer, hmx
      lonsaer = dim1
      latsaer = dim2
      hmx = int(dim1/2)       ! to swap long from W-E to E-W

      if(me==master) then
         print *, 'MERRA2 dim: ',dim1, dim2, dim3
      endif

! allocate arrays
      if (.not. allocated(aer_loni)) then
        allocate (aer_loni(lonsaer))
        allocate (aer_lati(latsaer))
      endif

      if (.not. allocated(aer_lat)) then
        allocate(aer_lat(latsaer))
        allocate(aer_lon(lonsaer))
        allocate(aerin(lonsaer,latsaer,levsaer,ntrcaerm,timeaer))
        allocate(aer_pres(lonsaer,latsaer,levsaer,timeaer))
      endif

! construct lat/lon array
      call nf_inq_varid(ncid, 'lat', varid)
      call nf_get_var(ncid, varid, aer_lati)
      call nf_inq_varid(ncid, 'lon', varid)
      call nf_get_var(ncid, varid, aer_loni)

      do i = 1, hmx     ! flip from (-180,180) to (0,360)
        if(aer_loni(i)<0.)  aer_loni(i)=aer_loni(i)+360.
        aer_lon(i+hmx) = aer_loni(i)
        aer_lon(i)     = aer_loni(i+hmx)
      enddo

      do i = 1, latsaer
        aer_lat(i)     = aer_lati(i)
      enddo

      call nf_close(ncid)

! allocate local working arrays
      if (.not. allocated(buff)) then
        allocate (buff(lonsaer, latsaer, dim3))
        allocate (pres_tmp(lonsaer,dim3))
      endif
      if (.not. allocated(buffx)) then
        allocate (buffx(lonsaer, latsaer, dim3,1))
      endif

!! ===================================================================
!! loop thru m01 - m12 for aer/pres array
!! ===================================================================
      do imon = 1, timeaer
       write(mn,'(i2.2)') imon
       fname=trim("aeroclim.m"//mn//".nc")
       inquire (file = fname, exist = file_exist)
       if (.not. file_exist) then
         errmsg = 'Error in read_aerdata: file ' // trim(fname) // ' not found'
         errflg = 1
         return
       endif

       call nf_open(fname , nf90_NOWRITE, ncid)

! ====> construct 3-d pressure array (Pa)
       call nf_inq_varid(ncid, "DELP", varid)
       call nf_get_var(ncid, varid, buff)

       do j = 1, latsaer
        do i = 1, lonsaer
! constract pres_tmp (top-down), note input is top-down
         pres_tmp(i,1) = 0.
         do k=2, dim3
          pres_tmp(i,k) = pres_tmp(i,k-1)+buff(i,j,k)
         enddo    !k-loop
        enddo     !i-loop (lon)

! extract pres_tmp to fill aer_pres (in  Pa)
        do k = 1, levsaer
         if ( iflip == 0 )  then             ! data from toa to sfc
           klev = k
         else                                ! data from sfc to top
           klev = ( dim3 - k ) + 1
         endif
         do i = 1, hmx
         aer_pres(i+hmx,j,k,imon)= 1.d0*pres_tmp(i,klev)
         aer_pres(i,j,k,imon)    = 1.d0*pres_tmp(i+hmx,klev)
         enddo     !i-loop (lon)
        enddo     !k-loop (lev)
       enddo     !j-loop (lat)

! ====> construct 4-d aerosol array (kg/kg)
! merra2 data is top down
! for GFS, iflip 0: toa to sfc; 1: sfc to toa
       DO ii = 1, ntrcaerm
         vname=trim(specname(ii))
         call nf_inq_varid(ncid, vname, varid)
         call nf_get_var(ncid, varid, buffx)

         do j = 1, latsaer
         do k = 1, levsaer
! input is from toa to sfc
          if ( iflip == 0 )  then             ! data from toa to sfc
            klev = k
          else                                ! data from sfc to top
            klev = ( dim3 - k ) + 1
          endif
          do i = 1, hmx
          aerin(i+hmx,j,k,ii,imon) = 1.d0*buffx(i,j,klev,1)
          if(aerin(i+hmx,j,k,ii,imon)<0.or.aerin(i+hmx,j,k,ii,imon)>1.)  then
            aerin(i+hmx,j,k,ii,imon) = 0.
          end if
          aerin(i,j,k,ii,imon) = 1.d0*buffx(i+hmx,j,klev,1)
          if(aerin(i,j,k,ii,imon)<0.or.aerin(i,j,k,ii,imon)>1.)  then
            aerin(i,j,k,ii,imon) = 0.
          end if
          enddo    !i-loop (lon)
         enddo     !k-loop (lev)
         enddo     !j-loop (lat)

       ENDDO           ! ii-loop (ntracaerm)

! close the file
       call nf_close(ncid)
      enddo      !imon-loop
!---
      deallocate (aer_loni, aer_lati)
      deallocate (buff, pres_tmp)
      deallocate (buffx)

      END SUBROUTINE read_aerdata
!
!**********************************************************************
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
      SUBROUTINE aerinterpol(me,master,npts,IDATE,FHOUR,jindx1,jindx2, &
                 ddy,iindx1,iindx2,ddx,lev,prsl,aerout)
!
      USE MACHINE,  ONLY : kind_phys
      use aerclm_def
      implicit none
      integer   i1,i2, iday,j,j1,j2,l,npts,nc,n1,n2,lev,k,i,ii
      real(kind=kind_phys) fhour,temj, tx1, tx2,temi
!

      integer  JINDX1(npts), JINDX2(npts),iINDX1(npts),iINDX2(npts)
      integer  me,idate(4), master
      integer  IDAT(8),JDAT(8)
!
      real(kind=kind_phys) DDY(npts), ddx(npts),ttt
      real(kind=kind_phys) aerout(npts,lev,ntrcaer)
      real(kind=kind_phys) aerpm(npts,levsaer,ntrcaer)
      real(kind=kind_phys) prsl(npts,lev), aerpres(npts,levsaer)
      real(kind=kind_phys) RINC(5), rjday
      integer jdow, jdoy, jday
      real(4) rinc4(5)
      integer w3kindreal,w3kindint
!
      IDAT=0
      IDAT(1)=IDATE(4)
      IDAT(2)=IDATE(2)
      IDAT(3)=IDATE(3)
      IDAT(5)=IDATE(1)
      RINC=0.
      RINC(2)=FHOUR
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal==4) then
        rinc4=rinc
        CALL W3MOVDAT(RINC4,IDAT,JDAT)
      else
        CALL W3MOVDAT(RINC,IDAT,JDAT)
      endif
!
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow,jdoy,jday)
      rjday = jdoy + jdat(5) / 24.
      IF (RJDAY .LT. aer_time(1)) RJDAY = RJDAY+365.
!
      n2 = 13
      do j=2, 12
       if (rjday .lt. aer_time(j)) then
          n2 = j
          exit
       endif
      enddo
      n1 = n2 - 1
!
      tx1 = (aer_time(n2) - rjday) / (aer_time(n2) - aer_time(n1))
      tx2 = 1.0 - tx1
      if (n2 > 12) n2 = n2 -12

!
      DO L=1,levsaer
        DO J=1,npts
          J1  = JINDX1(J)
          J2  = JINDX2(J)
          TEMJ = 1.0 - DDY(J)
          I1  = IINDX1(J)
          I2  = IINDX2(J)
          TEMI = 1.0 - DDX(J)
          DO ii=1,ntrcaer
           aerpm(j,L,ii) =                                                        &
           tx1*(TEMI*TEMJ*aerin(I1,J1,L,ii,n1)+DDX(j)*DDY(J)*aerin(I2,J2,L,ii,n1)&
               +TEMI*DDY(j)*aerin(I1,J2,L,ii,n1)+DDX(j)*TEMJ*aerin(I2,J1,L,ii,n1))&
          +tx2*(TEMI*TEMJ*aerin(I1,J1,L,ii,n2)+DDX(j)*DDY(J)*aerin(I2,J2,L,ii,n2) &
               +TEMI*DDY(j)*aerin(I1,J2,L,ii,n2)+DDX(j)*TEMJ*aerin(I2,J1,L,ii,n2))
          ENDDO

          aerpres(j,L) =                                                         &
           tx1*(TEMI*TEMJ*aer_pres(I1,J1,L,n1)+DDX(j)*DDY(J)*aer_pres(I2,J2,L,n1)&
               +TEMI*DDY(j)*aer_pres(I1,J2,L,n1)+DDX(j)*TEMJ*aer_pres(I2,J1,L,n1))&
          +tx2*(TEMI*TEMJ*aer_pres(I1,J1,L,n2)+DDX(j)*DDY(J)*aer_pres(I2,J2,L,n2) &
               +TEMI*DDY(j)*aer_pres(I1,J2,L,n2)+DDX(j)*TEMJ*aer_pres(I2,J1,L,n2))

        ENDDO
      ENDDO

! don't flip, input is the same direction as GFS  (bottom-up)
      DO J=1,npts
        DO L=1,lev
           if(prsl(j,L).ge.aerpres(j,1)) then
              DO ii=1, ntrcaer
               aerout(j,L,ii)=aerpm(j,1,ii)        !! sfc level
              ENDDO
           else if(prsl(j,L).le.aerpres(j,levsaer)) then
              DO ii=1, ntrcaer
               aerout(j,L,ii)=aerpm(j,levsaer,ii)  !! toa top
              ENDDO
           else
             DO  k=1, levsaer-1      !! from sfc to toa
              IF(prsl(j,L)<aerpres(j,k).and.prsl(j,L)>aerpres(j,k+1)) then
                 i1=k
                 i2=min(k+1,levsaer)
                 exit
              ENDIF
             ENDDO
             temi = prsl(j,L)-aerpres(j,i2)
             temj = aerpres(j,i1) - prsl(j,L)
             tx1 = temi/(aerpres(j,i1) - aerpres(j,i2))
             tx2 = temj/(aerpres(j,i1) - aerpres(j,i2))
             DO ii = 1, ntrcaer
           aerout(j,L,ii)= aerpm(j,i1,ii)*tx1 + aerpm(j,i2,ii)*tx2
             ENDDO
           endif
        ENDDO   !L-loop
      ENDDO     !J-loop
!
      RETURN
      END SUBROUTINE aerinterpol

end module aerinterp

