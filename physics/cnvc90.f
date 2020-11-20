!>\file cnvc90.f
!! This file contains the calculation of fraction of convective cloud,
!! pressure at bottom of convective cloud and at top of convective
!! cloud.
      module cnvc90

      contains

      subroutine cnvc90_init()
      end subroutine cnvc90_init

!>\defgroup GFS_cnvc90 GFS Convective Cloud Diagnostics Module
!> @{
!! This module contains the calculation of fraction of convective cloud,
!! pressure at bottom of convective cloud and at top of convective
!! cloud.
!> \section arg_table_cnvc90_run Argument Table
!! \htmlinclude cnvc90_run.html
!!
! \section gen_cnvc_run GFS cnvc90_run General Algorithm
      SUBROUTINE cnvc90_run(CLSTP,IM,RN,KBOT,KTOP,KM,PRSI,              &
     &                      ACV,ACVB,ACVT,CV,CVB,CVT,errmsg,errflg)

      USE MACHINE, ONLY :kind_phys
      implicit none

      ! Interface variables
      real(kind=kind_phys), intent(in) :: clstp
      integer,              intent(in) :: im, km
      real(kind=kind_phys), intent(in) :: RN(IM)
      integer,              intent(in) :: KBOT(IM)
      integer,              intent(in) :: KTOP(IM)
      real(kind=kind_phys), intent(in) :: prsi(IM,km+1)
      real(kind=kind_phys), intent(inout) :: ACV(IM)
      real(kind=kind_phys), intent(inout) :: ACVB(IM)
      real(kind=kind_phys), intent(inout) :: ACVT(IM)
      real(kind=kind_phys), intent(inout) :: CV(IM)
      real(kind=kind_phys), intent(inout) :: CVB(IM)
      real(kind=kind_phys), intent(inout) :: CVT(IM)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer              :: i,ibot,itop,lc,lz,n,ncc
      real(kind=kind_phys) :: ah,cc1,cc2,cvb0,p1,p2,rkbot,rktop,val
      integer              :: NMD(IM)
      real(kind=kind_phys) :: PMD(IM)
!
      real (kind=kind_phys), parameter :: cons_100=100.0
      real(kind=kind_phys) :: R_KBOT_I, R_KTOP_I
!
      PARAMETER(NCC=9)
      real(kind=kind_phys) :: CC(NCC),P(NCC)
      DATA CC/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/
      DATA P/.14,.31,.70,1.6,3.4,7.7,17.,38.,85./
      DATA CVB0/100./
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
      LZ=0
      LC=0
      IF(CLSTP.GE.1000.) LZ=1
      IF(CLSTP.GE.1100..OR.(CLSTP.LT.1000..AND.CLSTP.GE.100.)) LC=1
      AH=MOD(CLSTP,cons_100)
      IF(LZ.NE.0) THEN
        DO I=1,IM
          ACV(I)  = 0.
          ACVB(I) = CVB0
          ACVT(I) = 0.
        ENDDO
      ENDIF
      IF(LC.NE.0) THEN
        DO I=1,IM
          IF(RN(I).GT.0.) THEN
            ACV(I)  = ACV(I)+RN(I)
            R_KBOT_I= KBOT(I)
            ACVB(I) = MIN(ACVB(I),R_KBOT_I)
            R_KTOP_I= KTOP(I)
            ACVT(I) = MAX(ACVT(I),R_KTOP_I)
          ENDIF
        ENDDO
      ENDIF
      IF(AH.GT.0.01.AND.AH.LT.99.99) THEN
        DO I=1,IM
          IF(ACV(I).GT.0.) THEN
!           CVB(I) = ACVB(I)
!           CVT(I) = ACVT(I)
!       convert cvt and cvb to pressures
            ITOP   = NINT(ACVT(I))
            CVT(I) = PRSI(i,ITOP+1)
            IBOT   = NINT(ACVB(I))
            CVB(I) = PRSI(i,IBOT)
          ELSE
!           CVB(I) = CVB0
            CVB(I) = 0.
            CVT(I) = 0.
          ENDIF
          PMD(I)   = ACV(I)*(24.E+3/AH)
          NMD(I)   = 0
        ENDDO
        DO N=1,NCC
          DO I=1,IM
            IF(PMD(I).GT.P(N)) NMD(I) = N
          ENDDO
        ENDDO
        DO I=1,IM
          IF(NMD(I).EQ.0) THEN
            CV(I)  = 0.
!           CVB(I) = CVB0
            CVB(I) = 0.
            CVT(I) = 0.
          ELSEIF(NMD(I).EQ.NCC) THEN
            CV(I)  = CC(NCC)
          ELSE
            CC1    = CC(NMD(I))
            CC2    = CC(NMD(I)+1)
            P1     = P(NMD(I))
            P2     = P(NMD(I)+1)
            CV(I)  = CC1 + (CC2-CC1)*(PMD(I)-P1)/(P2-P1)
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE cnvc90_run
!> @}

      subroutine cnvc90_finalize()
      end subroutine cnvc90_finalize

      end module cnvc90

