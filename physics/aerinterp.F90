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

      SUBROUTINE read_aerdata (me, master, iflip, idate )

      use machine,  only: kind_phys
      use aerclm_def
      use netcdf

!--- in/out
      integer, intent(in) :: me, master, iflip, idate(4)

!--- locals
      integer      :: ncid, varid
      integer      :: i, j, k, n, ii, ijk, imon, klev
      character    :: fname*50, mn*2, fldname*10
      logical      :: file_exist
      real(kind=4), allocatable, dimension(:,:,:)   :: ps_clm
      real(kind=4), allocatable, dimension(:,:,:,:) :: delp_clm
      real(kind=4), allocatable, dimension(:,:,:,:) :: aer_clm
      real(kind=4), allocatable, dimension(:,:,:,:) :: airden_clm
      real(kind=4), allocatable, dimension(:)       :: pres_tmp

      allocate (delp_clm(lonsaer,latsaer,lmerra,1))
      allocate (aer_clm(lonsaer,latsaer,lmerra,1))
      allocate (airden_clm(lonsaer,latsaer,lmerra,1))
      allocate (ps_clm(lonsaer,latsaer,1))
      allocate (pres_tmp(lmerra))

! allocate aerclm_def arrays: aerin and aer_pres
      allocate (aerin(lonsaer,latsaer,levsaer,ntrcaer,timeaer))
      allocate (aer_pres(lonsaer,latsaer,levsaer,timeaer))

      if (me == master) then
         if ( iflip == 0 )  then             ! data from toa to sfc
          print *, "EJ, GFS is top-down"
         else
          print *, "EJ, GFS is bottom-up"
         endif
      endif

      do imon = 1, timeaer
       !ijk = imon + idate(2)+int(idate(3)/16)-2
       !if ( ijk .le. 0 )  ijk = 12
       !if ( ijk .eq. 13 ) ijk = 1
       !if ( ijk .eq. 14 ) ijk = 2
       write(mn,'(i2.2)') imon
       fname=trim("merra2C.aerclim.2003-2014.m"//mn//".nc")
       if (me == master) print *, "EJ,aerosol climo:", fname, &
                         "for imon:",imon,idate

       inquire (file = fname, exist = file_exist)
       if ( file_exist ) then
        if (me == master) print *,     &
                         "EJ, aerosol climo found; proceed the run" 
       else
        print *,"EJ, Error! aerosol climo not found; abort the run"
        stop 555
       endif

       call nf_open(fname, NF90_NOWRITE, ncid)

! merra2 data is top down
! for GFS, iflip 0: toa to sfc; 1: sfc to toa

! read aerosol mixing ratio arrays (kg/kg)
! construct 4-d aerosol mass concentration (kg/m3)
       call nf_inq_varid(ncid, 'AIRDENS', varid)
       call nf_get_var(ncid, varid, airden_clm)
!      if(me==master)  print *, "EJ, read airdens", airden_clm(1,1,:,1)

       do ii = 1, ntrcaer
        fldname=specname(ii)
        call nf_inq_varid(ncid, fldname, varid)
        call nf_get_var(ncid, varid, aer_clm)
!       if(me==master)  print *, "EJ, read ", fldname, aer_clm(1,1,:,1)
        do i = 1, lonsaer
        do j = 1, latsaer
        do k = 1, levsaer 
! input is from toa to sfc
         if ( iflip == 0 )  then             ! data from toa to sfc
           klev = k
         else                                ! data from sfc to top
           klev = ( lmerra - k ) + 1
         endif 
         aerin(i,j,k,ii,imon) = aer_clm(i,j,klev,1)*airden_clm(i,j,klev,1) 
        enddo     !k-loop (lev)
        enddo     !j-loop (lat)
        enddo     !i-loop (lon)
       enddo      !ii-loop (ntrac)

! aer_clm is top-down (following MERRA2)
! aerin is bottom-up (following GFS)

!      if ( imon == 1 .and.  me == master ) then
!        print *, 'EJ, du1(1,1) :', aerin(1,1,:,1,imon)
!      endif

! construct 3-d pressure array (Pa)
       call nf_inq_varid(ncid, "PS", varid)
       call nf_get_var(ncid, varid, ps_clm)
       call nf_inq_varid(ncid, "DELP", varid)
       call nf_get_var(ncid, varid, delp_clm)

!      if ( imon == 1 .and.  me == master ) then
!        print *, 'EJ, ps_clm:', ps_clm(1,1,1)
!        print *, 'EJ, delp_clm:', delp_clm(1,1,:,1)
!      endif

       do i = 1, lonsaer
       do j = 1, latsaer

! constract pres_tmp (top-down)
        pres_tmp(1) = 0.
        do k=2, lmerra
         pres_tmp(k) = pres_tmp(k-1) + delp_clm(i,j,k,1)
        enddo
!       if (imon==1 .and.  me==master .and. i==1 .and. j==1 ) then
!        print *, 'EJ, pres_tmp:', pres_tmp(:)
!       endif

! extract pres_tmp to fill aer_pres
        do k = 1, levsaer
         if ( iflip == 0 )  then             ! data from toa to sfc
           klev = k
         else                                ! data from sfc to top
           klev = ( lmerra - k ) + 1
         endif 
         aer_pres(i,j,k,imon)=  pres_tmp(klev)
        enddo     !k-loop (lev)
!       if (imon==1 .and.  me==master .and. i==1 .and. j==1 ) then
!        print *, 'EJ, aer_pres:', aer_pres(i,j,:,imon) 
!       endif

       enddo     !j-loop (lat)
       enddo     !i-loop (lon)

!      if (imon==1 .and.  me==master ) then
!        print *, 'EJx, aer_pres_i1:',(aer_pres(1,1:180,levsaer,imon) )
!      endif

! construct lat/lon array
       if (imon == 1 ) then
        call nf_inq_varid(ncid, "lat", varid)
        call nf_get_var(ncid, varid, aer_lat)
        call nf_inq_varid(ncid, "lon", varid)
        call nf_get_var(ncid, varid, aer_lon)
        do i = 1, lonsaer
          if(aer_lon(i) < 0.)  aer_lon(i) = aer_lon(i) + 360.
        enddo
!       if (imon==1 .and. me == master) then
! print *, "EJ, lat:", aer_lat(:)
! print *, "EJ, lon:", aer_lon(:)
!       endif
       endif

! close the file
       call nf_close(ncid)
      enddo      !imon-loop

!---
      deallocate (ps_clm, delp_clm, pres_tmp, aer_clm, airden_clm )
      if (me == master) then
        write(*,*) 'Reading in GOCART aerosols data'
      endif

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
 
!       if (me == master .and. j<= 3) then
!       print *,'EJj,',j,' dlat=',dlat(j),' jindx12=',jindx1(j),&
!         jindx2(j),' aer_lat=',aer_lat(jindx1(j)),              &
!         aer_lat(jindx2(j)),' ddy=',ddy(j)
!       endif
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
!       if (me == master .and. j<= 3) then
!       print *,'EJi,',j,' dlon=',dlon(j),' iindx12=',iindx1(j),&
!        iindx2(j),' aer_lon=',aer_lon(iindx1(j)),              &
!        aer_lon(iindx2(j)),' ddx=',ddx(j)
!       endif
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
      real(kind=kind_phys) aerout(npts,lev,ntrcaer),aerpm(npts,levsaer,ntrcaer)
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
!     if(me==master) print *,'EJ, IDAT ',IDAT(1:3), IDAT(5)
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
!     if(me==master)print *,'EJ,rjday=',rjday, ';aer_time,tx1,tx='   &
!    ,        aer_time(n1),aer_time(n2),tx1,tx2,n1,n2
!    
!     if(me==master) then
!      DO L=1,levsaer
!       print *,'EJ,aerin(n1,n2)=',L,aerin(1,1,L,1,n1),aerin(1,1,L,1,n2)
!      ENDDO
!     endif

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

!         IF(me==master .and. j==1) THEN
!          print *, 'EJ,aer/ps:',L,aerpm(j,L,1),aerpres(j,L)
!          if(L==1) then
!             print *, 'EJ, wgt:',TEMI*TEMJ,DDX(j)*DDY(J),TEMI*DDY(j),DDX(j)*TEMJ
!             print *, 'EJ, aerx:',aerin(I1,J1,L,ii,n1), &
!             aerin(I2,J2,L,ii,n1), aerin(I1,J2,L,ii,n1), aerin(I2,J1,L,ii,n1)
!             print *, 'EJ, aery:',aerin(I1,J1,L,ii,n2), &
!             aerin(I2,J2,L,ii,n2), aerin(I1,J2,L,ii,n2), aerin(I2,J1,L,ii,n2)
!          endif
!         ENDIF 
        ENDDO
      ENDDO

! note: input is set to be same as GFS 
      DO J=1,npts
        DO L=1,lev
           if(prsl(j,l).ge.aerpres(j,levsaer)) then
              DO ii=1, ntrcaer
               aerout(j,l,ii)=aerpm(j,levsaer,ii)
              ENDDO
           else if(prsl(j,l).le.aerpres(j,1)) then
              DO ii=1, ntrcaer
               aerout(j,l,ii)=aerpm(j,1,ii)
              ENDDO
           else
             DO  k=levsaer-1,1,-1
               IF(prsl(j,l)>aerpres(j,k)) then
                 i1=k
                 i2=min(k+1,levsaer)
                 exit
               end if
             end do
             DO ii = 1, ntrcaer
             aerout(j,l,ii)=aerpm(j,i1,ii)+(aerpm(j,i2,ii)-aerpm(j,i1,ii))&
                 /(aerpres(j,i2)-aerpres(j,i1))*(prsl(j,l)-aerpres(j,i1))
!            IF(me==master .and. j==1 .and. ii==1) then
!             print *, 'EJ, aerout:',aerout(j,l,ii), aerpm(j,i1,ii), &
!                aerpm(j,i2,ii), aerpres(j,i2), aerpres(j,i1), prsl(j,l)
!            ENDIF
             ENDDO
           endif
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE aerinterpol

end module aerinterp
