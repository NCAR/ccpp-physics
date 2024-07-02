!>\file iccninterp.F90
!! This file contains subrouines of reading and interplating
!! IN and CCN data.

!>\ingroup mod_GFS_phys_time_vary
!! This module contains subroutines of reading and interplating
!! IN and CCN data.
module iccninterp

    implicit none

    private

    public :: read_cidata, setindxci, ciinterpol

contains

      SUBROUTINE read_cidata (me, master)
      use machine,  only: kind_phys
      use iccn_def
      use netcdf
!--- in/out
      integer, intent(in) :: me
      integer, intent(in) :: master
!--- locals
      integer :: ncerr
      integer :: i, n, k, ncid, varid,j,it
      real(kind=kind_phys), allocatable, dimension(:) :: hyam,hybm
      real(kind=4), allocatable, dimension(:,:,:) :: ci_ps

      allocate (hyam(kcipl), hybm(kcipl), ci_ps(lonscip,latscip,timeci))
      allocate (ciplin(lonscip,latscip,kcipl,timeci))
      allocate (ccnin(lonscip,latscip,kcipl,timeci))
      allocate (ci_pres(lonscip,latscip,kcipl,timeci))
      ncerr = nf90_open("cam5_4_143_NAAI_monclimo2.nc", NF90_NOWRITE, ncid)
      ncerr = nf90_inq_varid(ncid, "lat", varid)
      ncerr = nf90_get_var(ncid, varid, ci_lat)
      ncerr = nf90_inq_varid(ncid, "lon", varid)
      ncerr = nf90_get_var(ncid, varid, ci_lon)
      ncerr = nf90_inq_varid(ncid, "PS", varid)
      ncerr = nf90_get_var(ncid, varid, ci_ps)
      ncerr = nf90_inq_varid(ncid, "hyam", varid)
      ncerr = nf90_get_var(ncid, varid, hyam)
      ncerr = nf90_inq_varid(ncid, "hybm", varid)
      ncerr = nf90_get_var(ncid, varid, hybm)
      ncerr = nf90_inq_varid(ncid, "NAAI", varid)
      ncerr = nf90_get_var(ncid, varid, ciplin)
      do it = 1,timeci
        do k=1, kcipl
          ci_pres(:,:,k,it)=hyam(k)*1.e5+hybm(k)*ci_ps(:,:,it)
        end do
      end do
      ncerr = nf90_close(ncid)
      ncerr = nf90_open("cam5_4_143_NPCCN_monclimo2.nc", NF90_NOWRITE, ncid)
      ncerr = nf90_inq_varid(ncid, "NPCCN", varid)
      ncerr = nf90_get_var(ncid, varid, ccnin)
      ncerr = nf90_close(ncid)
!---
      deallocate (hyam, hybm, ci_ps)
      if (me == master) then
        write(*,*) 'Reading in ICCN data',ci_time
      endif

      END SUBROUTINE read_cidata
!
!**********************************************************************
!
      SUBROUTINE setindxci(npts,dlat,jindx1,jindx2,ddy,dlon,                &
                 iindx1,iindx2,ddx)
!
      USE MACHINE,  ONLY: kind_phys
      USE iccn_def, ONLY: jci => latscip, ci_lat,ici=>lonscip, ci_lon
!
      implicit none
!
      integer npts, JINDX1(npts),JINDX2(npts),iINDX1(npts),iINDX2(npts)
      real(kind=kind_phys) dlat(npts),DDY(npts),dlon(npts),DDX(npts)
!
      integer i,j
!
      DO J=1,npts
        jindx2(j) = jci + 1
        do i=1,jci
          if (dlat(j) < ci_lat(i)) then
            jindx2(j) = i
            exit
          endif
        enddo
        jindx1(j) = max(jindx2(j)-1,1)
        jindx2(j) = min(jindx2(j),jci)
        if (jindx2(j) .ne. jindx1(j)) then
          DDY(j) = (dlat(j)           - ci_lat(jindx1(j))) &
                 / (ci_lat(jindx2(j)) - ci_lat(jindx1(j)))
        else
          ddy(j) = 1.0
        endif
        !print *,' j=',j,' dlat=',dlat(j),' jindx12=',jindx1(j), &
        ! jindx2(j),' ci_lat=',ci_lat(jindx1(j)),               &
        ! ci_lat(jindx2(j)),' ddy=',ddy(j)
      ENDDO

      DO J=1,npts
        iindx2(j) = ici + 1
        do i=1,ici
          if (dlon(j) < ci_lon(i)) then
            iindx2(j) = i
            exit
          endif
        enddo
        iindx1(j) = max(iindx2(j)-1,1)
        iindx2(j) = min(iindx2(j),ici)
        if (iindx2(j) .ne. iindx1(j)) then
          ddx(j) = (dlon(j)           - ci_lon(iindx1(j))) &
                 / (ci_lon(iindx2(j)) - ci_lon(iindx1(j)))
        else
          ddx(j) = 1.0
        endif
        !print *,' j=',j,' dlon=',dlon(j),' iindx12=',iindx1(j), &
        ! iindx2(j),' ci_lon=',ci_lon(iindx1(j)),               &
        ! ci_lon(iindx2(j)),' ddx=',ddx(j)
      ENDDO
 
      RETURN
      END SUBROUTINE setindxci
!
!**********************************************************************
!**********************************************************************
!
      SUBROUTINE ciinterpol(me,npts,IDATE,FHOUR,jindx1,jindx2,ddy, &
                 iindx1,iindx2,ddx,lev, prsl, ciplout,ccnout)
!
      USE MACHINE,  ONLY : kind_phys, kind_dbl_prec
      use iccn_def
      implicit none
      integer   i1,i2, iday,j,j1,j2,l,npts,nc,n1,n2,lev,k,i
      real(kind=kind_phys) fhour,temj, tx1, tx2,temi
!
 
      integer  JINDX1(npts), JINDX2(npts),iINDX1(npts),iINDX2(npts)
      integer  me,idate(4)
      integer  IDAT(8),JDAT(8)
!
      real(kind=kind_phys) DDY(npts), ddx(npts),ttt
      real(kind=kind_phys) ciplout(npts,lev),cipm(npts,kcipl)
      real(kind=kind_phys) ccnout(npts,lev),ccnpm(npts,kcipl)
      real(kind=kind_phys) cipres(npts,kcipl), prsl(npts,lev)
      real(kind=kind_phys) rjday
      real(kind=kind_dbl_prec) rinc(5)
      integer jdow, jdoy, jday
!
      IDAT=0
      IDAT(1)=IDATE(4)
      IDAT(2)=IDATE(2)
      IDAT(3)=IDATE(3)
      IDAT(5)=IDATE(1)
      RINC=0.
      RINC(2)=FHOUR
      CALL W3MOVDAT(RINC,IDAT,JDAT)
!
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow,jdoy,jday)
      rjday = jdoy + jdat(5) / 24.
      IF (RJDAY .LT. ci_time(1)) RJDAY = RJDAY+365.
!
      n2 = timeci + 1
      do j=2,timeci
        if (rjday .lt. ci_time(j)) then
          n2 = j
          exit
        endif
      enddo
      n1 = n2 - 1

!
!
      tx1 = (ci_time(n2) - rjday) / (ci_time(n2) - ci_time(n1))
      if (n2 > timeci) n2 = n2 - timeci
!     if (me .eq. 0) print *,' n1=',n1,' n2=',n2,' rjday=',rjday   &
!     ,'ci_time=',ci_time(n1),ci_time(n2), ci_time(timeci+1),tx1
      tx2 = 1.0 - tx1
!
      DO L=1,kcipl
        DO J=1,npts
          J1  = JINDX1(J)
          J2  = JINDX2(J)
          TEMJ = 1.0 - DDY(J)
          I1  = IINDX1(J)
          I2  = IINDX2(J)
          TEMI = 1.0 - DDX(J)
          cipm(j,L) =                                                           & 
            tx1*(TEMI*TEMJ*ciplin(I1,J1,L,n1)+DDX(j)*DDY(J)*ciplin(I2,J2,L,n1)  &
                +TEMI*DDY(j)*ciplin(I1,J2,L,n1)+DDX(j)*TEMJ*ciplin(I2,J1,L,n1)) & 
          + tx2*(TEMI*TEMJ*ciplin(I1,J1,L,n2)+DDX(j)*DDY(J)*ciplin(I2,J2,L,n2)  &
                +TEMI*DDY(j)*ciplin(I1,J2,L,n2)+DDX(j)*TEMJ*ciplin(I2,J1,L,n2)) 

          ccnpm(j,L) =                                                           & 
            tx1*(TEMI*TEMJ*ccnin(I1,J1,L,n1)+DDX(j)*DDY(J)*ccnin(I2,J2,L,n1)  &
                +TEMI*DDY(j)*ccnin(I1,J2,L,n1)+DDX(j)*TEMJ*ccnin(I2,J1,L,n1)) & 
          + tx2*(TEMI*TEMJ*ccnin(I1,J1,L,n2)+DDX(j)*DDY(J)*ccnin(I2,J2,L,n2)  &
                +TEMI*DDY(j)*ccnin(I1,J2,L,n2)+DDX(j)*TEMJ*ccnin(I2,J1,L,n2)) 

          cipres(j,L) =                                                          & 
            tx1*(TEMI*TEMJ*ci_pres(I1,J1,L,n1)+DDX(j)*DDY(J)*ci_pres(I2,J2,L,n1)  &
                +TEMI*DDY(j)*ci_pres(I1,J2,L,n1)+DDX(j)*TEMJ*ci_pres(I2,J1,L,n1)) & 
          + tx2*(TEMI*TEMJ*ci_pres(I1,J1,L,n2)+DDX(j)*DDY(J)*ci_pres(I2,J2,L,n2)  &
                +TEMI*DDY(j)*ci_pres(I1,J2,L,n2)+DDX(j)*TEMJ*ci_pres(I2,J1,L,n2)) 
        ENDDO
      ENDDO

      DO J=1,npts
        DO L=1,lev
           ! noticed input is from top to bottom
           if(prsl(j,l).ge.cipres(j,kcipl)) then
              ciplout(j,l)=cipm(j,kcipl)
              ccnout(j,l)=ccnpm(j,kcipl)
           else if(prsl(j,l).le.cipres(j,1)) then
               ciplout(j,l)=cipm(j,1)
               ccnout(j,l)=ccnpm(j,1)
           else
             DO  k=kcipl-1,1,-1
               ! DH* There is no backstop if this condition isn't met,
               ! i.e. i1 and i2 will have values determined by the
               ! previous code (line 178) - this leads to crashes in
               ! debug mode (out of bounds), for example for regression
               ! test fv3_stretched_nest_debug. For the time being,
               ! this is 'solved' by simply switching off ICCN
               ! if MG2/3 are not used (these are the only microphysics
               ! schemes that use the ICCN data); however, this doesn't
               ! mean that the code is correct for MG2/3, it just doesn't
               ! abort if the below condition isn't met, because the code
               ! is not tested in DEBUG mode. *DH
               IF(prsl(j,l)>cipres(j,k)) then
                 i1=k
                 i2=min(k+1,kcipl)
                 exit
               end if
             end do
             ciplout(j,l)=cipm(j,i1)+(cipm(j,i2)-cipm(j,i1))     &
                 /(cipres(j,i2)-cipres(j,i1))*(prsl(j,l)-cipres(j,i1))
             ccnout(j,l)=ccnpm(j,i1)+(ccnpm(j,i2)-ccnpm(j,i1))     &
                 /(cipres(j,i2)-cipres(j,i1))*(prsl(j,l)-cipres(j,i1))
            end if
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE ciinterpol

end module iccninterp
