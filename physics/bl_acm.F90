!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \file bl_acm.F90
!! This file contains the CCPP-compliant asymmetric convective model (ACM) boundary layer scheme.
!! IT COMPUTES VERTICAL MIXING IN AND ABOVE THE PBL ACCORDING TO THE ASYMMETRICAL 
!! CONVECTIVE MODEL VERSION 2  (ACM2), WHICH IS A COMBINED LOCAL NON-LOCAL CLOSURE SCHEME 
!! BASED ON THE ORIGINAL ACM (PLEIM AND CHANG 1992)
!!  REFERENCES:
!!   Pleim (2007) A combined local and non-local closure model for the atmospheric
!!                boundary layer.  Part1: Model description and testing.
!!                JAMC, 46, 1383-1395
!!   Pleim (2007) A combined local and non-local closure model for the atmospheric
!!                boundary layer.  Part2: Application and evaluation in a mesoscale
!!                meteorology model.  JAMC, 46, 1396-1409
!!
!!  REVISION HISTORY:
!!     AX        3/2005   - developed WRF version based on the MM5 PX LSM
!!     RG and JP 7/2006   - Finished WRF adaptation
!!     JP 12/2011 12/2011 - ACM2 modified so it's not dependent on first layer thickness.
!!     JP        3/2013   - WRFChem version. Mixing of chemical species are added
!!     HF        5/2016   - MPAS version
!!     JP        8/2017   - Z-coord version from MPAS version for hybrid coords
!!
!----------------------------------------------------------------------
module bl_acm

    use machine, only : kind_phys

    implicit none

    public bl_acm_init, bl_acm_run, bl_acm_finalize

    private
    
    real, parameter      :: ric    = 0.25                ! critical richardson number
    real, parameter      :: crankp = 0.5                 ! crank-nic parameter
    real, parameter      :: karman = 0.4

    contains
! -----------------------------------------------------------------------
! CCPP entry points for ACM
! -----------------------------------------------------------------------

    subroutine bl_acm_init ()
    end subroutine bl_acm_init

    subroutine bl_acm_finalize ()
    end subroutine bl_acm_finalize
      
!> \defgroup ACM bl_acm Main
!! \brief This subroutine contains all of the logic for the
!! ACM PBL scheme.
!!
!> \section arg_table_bl_acm_run Argument Table
!! \htmlinclude bl_acm_run.html
!!
! -------------------------------------------------------------------------------------------
    subroutine bl_acm_run (im, km, nvdiff, ntqvx, ntcwx, ntiwx, ntozx, xtime, dt, cpd, ep1, g, rd, qxs,&
                           tt, us, vs, ust, hfx, qfx, exner,    &
                           phii, phil, pbl, klpbl, utnp, vtnp, ttnp, qtnp,     &
                           lssav, ldiag3d, qdiag3d,                                    &
                           flag_for_pbl_generic_tend, du3dt_PBL, dv3dt_PBL,             &
                           dt3dt_PBL, dq3dt_PBL, do3dt_PBL, eddyz, eddyzm, errmsg, errflg)
      
      implicit none
      
      integer, intent(in) :: im, km, nvdiff, ntqvx, ntcwx, ntiwx, ntozx, xtime
      logical, intent(in) ::  lssav, ldiag3d, qdiag3d, flag_for_pbl_generic_tend
      real(kind=kind_phys), intent(in) :: dt, cpd, ep1, g, rd
      real(kind=kind_phys), dimension(im), intent(in) :: ust, hfx, qfx
      real(kind=kind_phys), dimension(im,km), intent(in) :: tt, exner, phil, us, vs
      real(kind=kind_phys), dimension(im,km+1), intent(in) :: phii
      real(kind=kind_phys), dimension(im,km,nvdiff), intent(in) :: qxs

      real(kind=kind_phys), dimension(im), intent(inout) :: pbl
      real(kind=kind_phys), dimension(im,im), intent(inout) :: utnp, vtnp, ttnp
      real(kind=kind_phys), dimension(im,im,nvdiff), intent(inout) :: qtnp
      real(kind=kind_phys), dimension(im,km), intent(inout) :: du3dt_PBL, dv3dt_PBL, dt3dt_PBL, dq3dt_PBL, do3dt_PBL

      integer, dimension(im), intent(out) :: klpbl
      real(kind=kind_phys), dimension(im,km), intent(out) :: eddyz, eddyzm
      
      ! error messages
      character(len=*),                               intent(out)     ::  errmsg
      integer,                                        intent(out)     ::  errflg
      
      !local variables
      integer :: i, k, l, ksrc, kmix, kpblht
      logical :: pbl_found
      
      integer, dimension(im) :: kpblh, noconv
      
      real(kind=kind_phys) :: tmpvtcon, thv1, gravi, th1, zh1, uh1, vh1, wss, tconv, dtmp, fintt, zmix, umix, vmix, tog, wssq, rdt
      real(kind=kind_phys), dimension(im) :: cpair, tst, qst, ustm, tstv, mol, rmol, wst, fint
      real(kind=kind_phys), dimension(im,km) :: theta, thetax, thetav, ux, vx, za, dzh, dzhi, dzfi, rib, cld_water
      real(kind=kind_phys), dimension(im,km+1) :: zf
      real(kind=kind_phys), dimension(im,km,nvdiff) :: qxx
      
      character*512 :: message
      
      real(kind=kind_phys), parameter :: tstv_min = 1.0e-6
      
      gravi = 1.0/g
      
      do k=1, km
        do i=1,im
          theta(i,k) = tt(i,k)/exner(i,k)
          tmpvtcon  = 1.0 + ep1 * qxs(i,k,ntqvx)
          thetav(i,k) = theta(i,k) * tmpvtcon
        end do
      end do
      
      do i = 1, im
        cpair(i)  = cpd * (1.0 + 0.84 * qxs(i,1,ntqvx))                    ! J/(K KG)
        tmpvtcon  = 1.0 + ep1 * qxs(i,1,ntqvx)                             ! COnversion factor for virtual temperature
        !pass in ustar rather than calc
        tst(i)    = -hfx(i) / ust(i)
        qst(i)    = -qfx(i) / ust(i)
        ustm(i)   = ust(i) !GJF: WRF scheme multiplies and divides by lowest model level wind speed for ustm (necessary?)
        thv1      = tmpvtcon * theta(i,1)
        tstv(i)   = tst(i)*tmpvtcon + thv1*ep1*qst(i)
        if(abs(tstv(i)) < tstv_min) then
          tstv(i) = sign(tstv_min,tstv(i))
        endif
        mol(i)    = thv1 * ust(i)**2/(karman*g*tstv(i))
        rmol(i)   = 1.0/mol(i)
        wst(i)    = ust(i) * (pbl(i)/(karman*abs(mol(i)))) ** 0.333333
      end do
      
!... Compute PBL height

!... compute the height of full- and half-sigma level above ground level
      do i = 1, im
        zf(i,1)    = 0.0
        klpbl(i)   = 1
      enddo
      
      do k = 1, km
        do i = 1, im
          !zf(i,k) = dzf(i,k) + zf(i,k-1)
          !za(i,k) = 0.5 * (zf(i,k) + zf(i,k-1))
          zf(i,k+1) = phii(i,k+1)*gravi
          za(i,k) = phil(i,k)*gravi
          dzh(i,k) = zf(i,k+1) - zf(i,k)
          dzhi(i,k)= 1./dzh(i,k)
        enddo
      enddo

      do k = 1, km-1
        do i = 1, im
          dzfi(i,k) = 1./(za(i,k+1)-za(i,k))
        enddo
      enddo
      do i = 1,im
        dzfi(i,km) = dzfi(i,km-1)
      enddo
      
!...  COMPUTE PBL WHERE RICHARDSON NUMBER = RIC (0.25) HOLTSLAG ET AL 1990  
      do i = 1, im
        
        do k = 1, km
          ksrc = k
          if (zf(i,k) .gt. 30.0) exit
        enddo

        th1 = 0.0
        zh1 = 0.0
        uh1 = 0.0
        vh1 = 0.0
        do k = 1,ksrc
          th1 = th1 + thetav(i,k)  
          zh1 = zh1 + za(i,k)
          uh1 = uh1 + us(i,k)
          vh1 = vh1 + vs(i,k)
        enddo  
        th1 = th1/ksrc
        zh1 = zh1/ksrc
        uh1 = uh1/ksrc
        vh1 = vh1/ksrc
        
        if(mol(i) < 0.0 .and. xtime > 1) then
          wss   = (ust(i) ** 3 + 0.6 * wst(i) ** 3) ** 0.33333
          tconv = -8.5 * ust(i) * tstv(i) / wss
          th1   = th1 + tconv
        endif
 
        kmix = ksrc
        do k = ksrc,km
          dtmp   = thetav(i,k) - th1
          if (dtmp < 0.0) kmix = k
        enddo
        
        if(kmix.gt.ksrc) then
          fintt = (th1 - thetav(i,kmix)) / (thetav(i,kmix+1)               &
                - thetav(i,kmix))
          zmix = fintt * (za(i,kmix+1)-za(i,kmix)) + za(i,kmix)
          umix = fintt * (us(i,kmix+1)-us(i,kmix)) + us(i,kmix)
          vmix = fintt * (vs(i,kmix+1)-vs(i,kmix)) + vs(i,kmix)
        else
          zmix = zh1
          umix = uh1
          vmix = vh1
        endif
        
        pbl_found = .false.
        do k = kmix,km
          dtmp   = thetav(i,k) - th1
          tog = 0.5 * (thetav(i,k) + th1) / g
          wssq = (us(i,k)-umix)**2                                     &
               + (vs(i,k)-vmix)**2
          if (kmix == ksrc) wssq = wssq + 100.*ust(i)*ust(i) 
          wssq = max( wssq, 0.1 )
          rib(i,k) = abs(za(i,k)-zmix) * dtmp / (tog * wssq)
          if (rib(i,k) .ge. ric) then
            pbl_found = .true.
            kpblh(i) = k
            exit
          end if
        end do
  
        if (.not. pbl_found) then
           write (message,*)' RIBX never exceeds RIC, RIB(i,kte) = ',rib(i,km),        &
                   ' THETAV(i,1) = ',thetav(i,1),' MOL=',mol(i),            &
                   ' TCONV = ',TCONV,' WST = ',WST(I),                      &
                   ' KMIX = ',kmix,' UST = ',UST(I),                       &
                   ' TST = ',TST(I),' U,V = ',US(I,1),VS(I,1),              &
                   ' I=',I
           errflg = 1
           errmsg = trim(message)
           return
        end if
      end do

      do i = 1,im
        if (kpblh(i) .gt. ksrc) then
!---------interpolate between levels -- jp 7/93
          fint(i) = (ric - rib(i,kpblh(i)-1)) / (rib(i,kpblh(i)) -       &
                     rib(i,kpblh(i)-1))
          if (fint(i) .gt. 0.5) then
            kpblht  = kpblh(i)
            fint(i) = fint(i) - 0.5
          else
            kpblht  = kpblh(i) - 1
            fint(i) = fint(i) + 0.5
          end if
          pbl(i)  = fint(i) * (zf(i,kpblht+1) - zf(i,kpblht)) +          &
                      zf(i,kpblht)
          klpbl(i) = kpblht
        else
          klpbl(i) = ksrc
          pbl(i)    = za(i,ksrc)                                                  
        end if
      end do
     
      do i = 1,im       
        noconv(i) = 0
       
! check for cbl and identify conv. vs. non-conv cells
        if (pbl(i) / mol(i) < -0.02 .and. klpbl(i) > 3        &
           .and. thetav(i,1) > thetav(i,2) .and. xtime > 1) then
          noconv(i)   = 1
          !GJF: regime isn't used internally and not needed by other CCPP schemes
!          regime(i) = 4.0                     ! free convective - acm
        endif
      enddo
      
      if (ntiwx > 0) then
        cld_water = qxs(:,:,ntcwx) + qxs(:,:,ntiwx)
      else
        cld_water = qxs(:,:,ntcwx)
      end if
      
      
      
!... Calculate Kz
      call EDDYX(zf, za, mol, pbl, ust, us, vs, tt, thetav, qxs(:,:,ntqvx), cld_water, g, &
                 rd, cpair, eddyz, eddyzm, 1, im, 1, km, 1, im, 1, km)
      !write(*,*) 'ACM:',qxs(:,:,1)
      call ACM(dt, noconv, zf, dzh, dzhi, klpbl, pbl, dzfi, mol, ust,&
                   tst, qst, eddyz, nvdiff, ntqvx, theta, qxs, thetax, qxx,  &
                   1, im, 1, km, 1, im, 1, km, 1, im, 1, km)
      !STOP
      !thetax = theta
      !qxx = qxs
      !qxx(:,:,ntqvx) = qxs(:,:,ntqvx)
      !write(*,*) 'here'
      !STOP
      !write(*,*) 'ACM:',qxx(:,:,1)
      !qxx = qxs
      !write(*,*) 'ACM (thetax): ',xtime, thetax - theta
      !write(*,*) 'ACM (qxx1): ',xtime, qxx(:,:,1) - qxs(:,:,1)
      
      call ACMM(dt, noconv, zf, dzh, dzhi, klpbl, pbl, dzfi, mol,    &
                  ustm,  eddyzm, us, vs, ux, vx,                     &
                  1, im, 1, km, 1, im, 1, km, 1, im, 1, km)

!... Calculate tendency due to PBL parameterization
      rdt = 1.0 / dt
      do k = 1, km
        do i = 1, im
          utnp(i,k)  = utnp(i,k) + (ux(i,k) - us(i,k)) * rdt
          vtnp(i,k)  = vtnp(i,k) + (vx(i,k) - vs(i,k)) * rdt
          ttnp(i,k)  = ttnp(i,k) + (thetax(i,k) - theta(i,k)) * exner(i,k) * rdt
          ! qtnp(i,k,1) = qtnp(i,k,1) + (qxx(i,k,1) - qxs(i,k,1)) * rdt
          ! qtnp(i,k,2) = qtnp(i,k,2) + (qxx(i,k,2) - qxs(i,k,2)) * rdt
          ! qtnp(i,k,3) = qtnp(i,k,3) + (qxx(i,k,3) - qxs(i,k,3)) * rdt
          
          !do l=1, nvdiff
          do l=1, 7
            qtnp(i,k,l) = qtnp(i,k,l) + (qxx(i,k,l) - qxs(i,k,l)) * rdt
          end do
          do l=8, nvdiff-1
            write(*,*) k,l, qxx(i,k,l), qxs(i,k,l)
            qtnp(i,k,l) = qtnp(i,k,l) + (qxx(i,k,l) - qxs(i,k,l)) * rdt
          end do
!          write(*,*) i,k,utnp(i,k),vtnp(i,k),ttnp(i,k),qtnp(i,k,1)
          if(lssav .and. ldiag3d .and. .not. flag_for_pbl_generic_tend) then
            dt3dt_pbl(i,k) = dt3dt_pbl(i,k) + (thetax(i,k) - theta(i,k)) * exner(i,k)
            du3dt_pbl(i,k) = du3dt_pbl(i,k) + (ux(i,k) - us(i,k))
            dv3dt_pbl(i,k) = dv3dt_pbl(i,k) + (vx(i,k) - vs(i,k))
            if (qdiag3d) then
              dq3dt_pbl(i,k) = dq3dt_pbl(i,k) + (qxx(i,k,ntqvx) - qxs(i,k,ntqvx))
              if (ntozx > 0) then
                do3dt_pbl(i,k) = do3dt_pbl(i,k) + (qxx(i,k,ntozx) - qxs(i,k,ntozx))
              end if
            end if
          end if
        end do
      end do
     
    end subroutine bl_acm_run
   
    SUBROUTINE EDDYX(ZF,  ZA,     MOL, PBL,  UST,               &
                     US,    VS,  TT,  THETAV,            &
                     QVS,   QCS, G, RD, CPAIR,                    &
                     EDDYZ, EDDYZM, its,ite, kts,kte,ims,ime,kms,kme )


!**********************************************************************
!   Two methods for computing Kz:
!   1.  Boundary scaling similar to Holtslag and Boville (1993)
!   2.  Local Kz computed as function of local Richardson # and vertical 
!       wind shear, similar to LIU & CARROLL (1996)
!
!**********************************************************************
!
!-- ZF              height of full level
!-- ZA              height of half level
!-- MOL             Monin-Obukhov length in 1D form
!-- PBL             PBL height in 1D form
!-- UST             friction velocity U* in 1D form (m/s)
!-- US              U wind 
!-- VS              V wind
!-- TT              temperature
!-- THETAV          potential virtual temperature
!-- QVS             water vapor mixing ratio (Kg/Kg)
!-- QCS             cloud mixing ratio (water + ice)(Kg/Kg)
!-- G               gravity
!-- RD              gas constant for dry air (j/kg/k)
!-- CPAIR           specific heat of moist air (M^2 S^-2 K^-1)
!-- EDDYZ           eddy diffusivity for heat KZ
!-- EDDYZM          eddy diffusivity for momentum KM
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      IMPLICIT NONE

!.......Arguments
  
!... Integer
      INTEGER,  INTENT(IN)   ::    its,ite, kts,kte,ims,ime,kms,kme
!... Real
      REAL , DIMENSION( ims:ime ),          INTENT(IN)  :: PBL, UST
      REAL ,                                INTENT(IN)  :: G, RD
      REAL , DIMENSION( its:ite ),          INTENT(IN)  :: MOL, CPAIR

      REAL , DIMENSION( ims:ime, kms:kme ), INTENT(IN)  :: US,VS, TT,   &
                                                           QVS, QCS
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(IN) :: ZA, THETAV
      REAL, DIMENSION( its:ite, 0:kte )  , INTENT(IN) :: ZF
      
      REAL , DIMENSION( its:ite, kts:kte ), INTENT(OUT) :: EDDYZ,EDDYZM

!.......Local variables

!... Integer
      INTEGER  :: ILX, KL, KLM, K, I

!... Real
      REAL     :: ZOVL, PHIH, WT, ZSOL, ZFUNC, DZF, SS, GOTH, EDYZ
      REAL     :: RI, QMEAN, TMEAN, XLV, ALPH, CHI, ZK, SQL, DENSF, KZO
      REAL     :: FH, FM
      REAL     :: WM, EDYZM, PHIM
!... Parameters
      REAL, PARAMETER :: RV     = 461.5
      REAL, PARAMETER :: RC     = 0.25
      REAL, PARAMETER :: RLAM   = 80.0
      REAL, PARAMETER :: GAMH   = 16.0 !Dyer74 !15.0  !  Holtslag and Boville (1993)
      REAL, PARAMETER :: GAMM   = 16.0 !Dyer74
      REAL, PARAMETER :: BETAH  = 5.0   !  Holtslag and Boville (1993)  BETAM = BETAH
      REAL, PARAMETER :: P      = 2.0   ! ZFUNC exponent
      REAL, PARAMETER :: EDYZ0  = 0.01  ! New Min Kz
      REAL, PARAMETER :: PR     = 0.8   ! Prandtl #
!      REAL, PARAMETER :: EDYZ0  = 0.1
!--   IMVDIF      imvdif=1 for moist adiabat vertical diffusion
      INTEGER, PARAMETER :: imvdif = 1
!
      ILX = ite 
      KL  = kte
      KLM = kte - 1
      
      DO K = kts,KLM
        DO I = its,ILX
          EDYZ = 0.0
          ZOVL = 0.0
          DZF  = ZA(I,K+1) - ZA(I,K)
          KZO = EDYZ0
!--------------------------------------------------------------------------
          IF (ZF(I,K) .LT. PBL(I)) THEN
            ZOVL = ZF(I,K) / MOL(I)
            IF (ZOVL .LT. 0.0) THEN
              IF (ZF(I,K) .LT. 0.1 * PBL(I)) THEN
                PHIH = 1.0 / SQRT(1.0 - GAMH * ZOVL)
                PHIM = (1.0 - GAMM * ZOVL)**(-0.25)
                WT   = UST(I) / PHIH
                WM   = UST(I) / PHIM
              ELSE
                ZSOL = 0.1 * PBL(I) / MOL(I)
                PHIH = 1.0 / SQRT(1.0 - GAMH * ZSOL)
                PHIM = (1.0 - GAMM * ZSOL)**(-0.25)
                WT   = UST(I) / PHIH
                WM   = UST(I) / PHIM
              ENDIF
            ELSE IF (ZOVL .LT. 1.0) THEN
              PHIH = 1.0 + BETAH * ZOVL
              WT   = UST(I) / PHIH
              WM   = WT
            ELSE
              PHIH = BETAH + ZOVL
              WT   = UST(I) / PHIH
              WM   = WT
            ENDIF
            ZFUNC      = ZF(I,K) * (1.0 - ZF(I,K) / PBL(I)) ** P
            EDYZ = KARMAN * WT * ZFUNC
            EDYZM = KARMAN * WM * ZFUNC
          ENDIF
!--------------------------------------------------------------------------!--------------------------------------------------------------------------
          SS   = ((US(I,K+1) - US(I,K)) ** 2 + (VS(I,K+1) - VS(I,K)) ** 2)   &
                  / (DZF * DZF) + 1.0E-9
          GOTH = 2.0 * G / (THETAV(I,K+1) + THETAV(I,K))
          RI   = GOTH * (THETAV(I,K+1) - THETAV(I,K)) / (DZF * SS)
!--------------------------------------------------------------------------
!         Adjustment to vert diff in Moist air
          IF(imvdif.eq.1)then
            IF (QCS(I,K) > 0.01E-3 .OR. QCS(I,K+1) > 0.01E-3) THEN
              QMEAN = 0.5 * (QVS(I,K) + QVS(I,K+1))
              TMEAN = 0.5 * (TT(I,K) + TT(I,K+1))
              XLV   = (2.501 - 0.00237 * (TMEAN - 273.15)) * 1.E6
              ALPH  =  XLV * QMEAN / RD / TMEAN
              CHI   =  XLV * XLV * QMEAN / CPAIR(I) / RV / TMEAN / TMEAN
              RI    = (1.0 + ALPH) * (RI -G * G / SS / TMEAN / CPAIR(I) *       &
                      ((CHI - ALPH) / (1.0 + CHI)))
            ENDIF
          ENDIF
!--------------------------------------------------------------------------
            
          ZK  = 0.4 * ZF(I,K)
          SQL = (ZK * RLAM / (RLAM + ZK)) ** 2
            
          IF (RI .GE. 0.0) THEN
!	          IF (ZF(I,K).LT.PBL(I).AND.ZOVL.GT.0.0) THEN
!	            FH = MAX((1.-ZF(I,K)/PBL(I))**2,0.01) * PHIH **(-2)
!                  SQL = ZK ** 2
!	          ELSE
!	            FH = (MAX(1.-RI/RC,0.01))**2
!	          ENDIF
            FH=1./(1.+10.*RI+50.*RI**2+5000.*RI**4)+0.0012  !pleim5
            FM= PR*FH + 0.00104

            EDDYZ(I,K) = KZO + SQRT(SS) * FH * SQL
            EDDYZM(I,K) = KZO + SQRT(SS) * FM * SQL
          ELSE
            EDDYZ(I,K) = KZO + SQRT(SS * (1.0 - 25.0 * RI)) * SQL
            EDDYZM(I,K) = EDDYZ(I,K) * PR
          ENDIF
  
          IF(EDYZ.GT.EDDYZ(I,K)) THEN
            EDDYZ(I,K) = EDYZ
            EDDYZM(I,K) = MIN(EDYZM,EDYZ*0.8)  !PR
          ENDIF

          EDDYZ(I,K) = MIN(1000.0,EDDYZ(I,K))
          EDDYZ(I,K) = MAX(KZO,EDDYZ(I,K))
          EDDYZM(I,K) = MIN(1000.0,EDDYZM(I,K))
          EDDYZM(I,K) = MAX(KZO,EDDYZM(I,K))

        ENDDO             ! for I loop
      ENDDO               ! for k loop
!
      DO I = its,ILX
        EDDYZ(I,KL) = 0.0 ! EDDYZ(I,KLM) -- changed jp 3/08
        EDDYZM(I,KL) = 0.0
      ENDDO

    END SUBROUTINE EDDYX
    
    SUBROUTINE ACM (DTPBL, NOCONV, ZF, DZH, DZHI,       &
                    KLPBL, PBL,   DZFI, MOL,  UST,                  &
                    TST, QST,   EDDYZ, NVDIFF, NTQVX,              &
                    THETA,  QXS,     &
                    THETAX, QXX,     &
                    ids,ide, kds,kde,                      &
                    ims,ime, kms,kme,                      &
                    its,ite, kts,kte)
 !**********************************************************************
 !   PBL model called the Asymmetric Convective Model, Version 2 (ACM2) 
 !   -- See top of module for summary and references
 !
 !---- REVISION HISTORY:
 !   AX     3/2005 - developed WRF version based on ACM2 in the MM5 PX LSM
 !   JP and RG 8/2006 - updates
 !   JP     3/2013 - Chem additions
 !
 !**********************************************************************
 !  ARGUMENTS:
 !-- DTPBL           PBL time step
 !-- NOCONV          If free convection =0, no; =1, yes
 !-- ZF              Height for full layer
 !-- DZH             Layer thickness --> ZF(K) - ZF(K-1)
 !-- DZHI            Inverse of layer thickness
 !-- KLPBL           PBL level at K index
 !-- PBL             PBL height in m
 !-- DZFI            Inverse layer thickness --> 1/(Z(K+1)-Z(K))
 !-- MOL             Monin-Obukhov length in 1D form
 !-- UST             U* in 1D form
 !-- TST             Theta* in 1D form
 !-- QST             Q* in 1D form
 !-- EDDYZ           eddy diffusivity KZ
 !-- US              U wind 
 !-- VS              V wind
 !-- THETA           potential temperature
 !-- QVS             water vapor mixing ratio (Kg/Kg)
 !-- QCS             cloud mixing ratio (Kg/Kg)
 !-- QIS             ice mixing ratio (Kg/Kg)
 !-- UX              new U wind 
 !-- VX              new V wind
 !-- THETAX          new potential temperature
 !-- QVX             new water vapor mixing ratio (Kg/Kg)
 !-- QCX             new cloud mixing ratio (Kg/Kg)
 !-- QIX             new ice mixing ratio (Kg/Kg)
 !-- CHEM            Chemical species mixing ratios (ppm)  WRFChem   
 !-- VD              Dry deposition velocity (m/s)         WRFChem
 !-- NCHEM           Number of chemical species            WRFChem
 !-----------------------------------------------------------------------
 !-----------------------------------------------------------------------
 
       IMPLICIT NONE
 
 !.......Arguments
 
 !... Integer
       INTEGER,  INTENT(IN)      ::      nvdiff, ntqvx,    &
                                         ids,ide, kds,kde, &
                                         ims,ime, kms,kme, &
                                         its,ite, kts,kte
       INTEGER,  DIMENSION( its:ite ), INTENT(IN)  :: NOCONV
       INTEGER,  DIMENSION( ims:ime ), INTENT(IN)  :: KLPBL
 
 !... Real
       REAL , DIMENSION( ims:ime ),          INTENT(IN)  :: PBL, UST
       REAL ,                                INTENT(IN)  :: DTPBL
       REAL , DIMENSION( its:ite ),          INTENT(IN)  :: MOL, TST, &
                                                            QST
       REAL , DIMENSION( its:ite, kts:kte ), INTENT(IN)  :: DZHI, DZH, DZFI
       REAL , DIMENSION( its:ite, 0:kte ),   INTENT(IN)  :: ZF
       REAL , DIMENSION( its:ite, kts:kte ), INTENT(INOUT)  :: EDDYZ
       REAL , DIMENSION( ims:ime, kms:kme ), INTENT(IN)  :: THETA
       REAL , DIMENSION( its:ite, kts:kte ), INTENT(OUT) :: THETAX
       REAL , DIMENSION( its:ite, kms:kme, nvdiff), INTENT(IN) :: QXS
       REAL , DIMENSION( its:ite, kts:kte, nvdiff), INTENT(OUT) :: QXX
 
 !.......Local variables
 
 !... Parameters
       REAL,    PARAMETER :: XX    = 0.5          ! FACTOR APPLIED TO CONV MIXING TIME STEP
        
 !... Integer
       INTEGER :: ILX, KL, KLM, I, K, NSP, NSPX, NLP, NL, JJ, L,LL
       INTEGER :: KCBLMX
       INTEGER, DIMENSION( its:ite ) :: KCBL
 
 !... Real
       REAL                               :: MBMAX, HOVL, MEDDY, MBAR
       REAL                               :: EKZ, RZ, FM, WSPD, DTS, DTRAT, F1
       REAL, DIMENSION( its:ite )         :: FSACM, DTLIM
       REAL, DIMENSION( kts:kte, its:ite) :: MBARKS, MDWN
       REAL, DIMENSION( kts:kte )         :: XPLUS, XMINUS
       REAL  DELC
       REAL, DIMENSION( kts:kte )                :: AI, BI, CI, EI !, Y
       REAL, ALLOCATABLE, DIMENSION( : , : )     :: DI, UI    
       REAL, ALLOCATABLE, DIMENSION( : , : )     :: FS
       REAL, ALLOCATABLE, DIMENSION( : , : , : ) :: VCI
 
       CHARACTER*80 :: message
 
 !
 !--Start Executable ----
 
       ILX = ite
       KL  = kte
       KLM = kte - 1
       NSP = nvdiff + 1
       !write(*,*) 'NSP:',NSP
       !write(*,*) 'QXS(ntqvx):',QXS(:,:,ntqvx)
       NSPX = NSP
 
       KCBLMX = 0
       MBMAX  = 0.0
       
 !...Allocate species variables
       ALLOCATE (DI( 1:NSPX,kts:kte ))       
       ALLOCATE (UI( 1:NSPX,kts:kte ))  
       ALLOCATE (FS( 1:NSPX, its:ite )) 
       ALLOCATE (VCI( 1:NSPX,its:ite,kts:kte  ))
 
 !---COMPUTE ACM MIXING RATE
       DO I = its, ILX
         DTLIM(I)  = DTPBL
         KCBL(I)   = 1
         FSACM(I)  = 0.0
 
         IF (NOCONV(I) .EQ. 1) THEN
           KCBL(I) = KLPBL(I)
 
 !-------MBARKS IS UPWARD MIXING RATE; MDWN IS DOWNWARD MIXING RATE
 !--New couple ACM & EDDY-------------------------------------------------------------
           HOVL     = -PBL(I) / MOL(I)
           FSACM(I) = 1./(1.+((KARMAN/(HOVL))**0.3333)/(0.72*KARMAN))
           MEDDY    = EDDYZ(I,1) * DZFI(i,1) / (PBL(I) - ZF(i,1))
           MBAR     = MEDDY * FSACM(I)
           DO K = kts,KCBL(I)-1
             EDDYZ(I,K) = EDDYZ(I,K) * (1.0 - FSACM(I))
           ENDDO
 
           MBMAX = AMAX1(MBMAX,MBAR)
           DO K = kts+1,KCBL(I)
             MBARKS(K,I) = MBAR
             MDWN(K,I)   = MBAR * (PBL(I) - ZF(i,K-1)) * DZHI(i,K)
           ENDDO
           MBARKS(1,I) = MBAR
           MBARKS(KCBL(I),I) = MDWN(KCBL(I),I)
           MDWN(KCBL(I)+1,I) = 0.0
         ENDIF
       ENDDO                              ! end of I loop
       
       DO K = kts,KLM
         DO I = its,ILX
           EKZ = EDDYZ(I,K) * DZFI(i,K) * DZHI(i,K)
           DTLIM(I) = AMIN1(0.75 / EKZ,DTLIM(I))
         ENDDO
       ENDDO
        
       DO I = its,ILX 
         IF (NOCONV(I) .EQ. 1) THEN
           KCBLMX = AMAX0(KLPBL(I),KCBLMX)
           RZ     = (ZF(i,KCBL(I)) - ZF(i,1)) * DZHI(i,1)
           DTLIM(I)  = AMIN1(XX / (MBARKS(1,I) * RZ),DTLIM(I))
         ENDIF
       ENDDO
       
       DO K = kts,KL
         DO I = its,ILX
           VCI(1,I,K) = THETA(I,K)
           DO L = 2, NSPX
             VCI(L,I,K) = QXS(I,K,L-1)
           END DO
         ENDDO
       ENDDO
       
       DO I = its,ILX
         FS(1,I) = -UST(I) * TST(I)
         DO L = 2, NSPX
           FS(L,I) = 0.0      ! SURFACE FLUXES OF VERTICALLY-DIFFUSED TRACERS = 0, except for water vapor, below
         END DO
         FS(NTQVX+1,I) = -UST(I) * QST(I)
       ENDDO
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO I = its,ILX      
 
         NLP   = INT(DTPBL / DTLIM(I) + 1.0)
         DTS   = (DTPBL / NLP)
         DTRAT = DTS / DTPBL
         DO NL = 1,NLP           ! LOOP OVER SUB TIME LOOP              
 
 !-- COMPUTE ARRAY ELEMENTS THAT ARE INDEPENDANT OF SPECIES
 
           DO K = kts,kte
             AI(K) = 0.0
             BI(K) = 0.0
             CI(K) = 0.0
             EI(K) = 0.0
           ENDDO
 
           DO K = 2, KCBL(I)
             EI(K-1) = -CRANKP * MDWN(K,I) * DTS * DZH(i,K) * DZHI(i,K-1)
             BI(K)   = 1.0 + CRANKP * MDWN(K,I) * DTS
             AI(K)   = -CRANKP * MBARKS(K,I) * DTS
           ENDDO
 
           EI(1) = EI(1) -EDDYZ(I,1) * CRANKP * DZHI(i,1) * DZFI(i,1) * DTS
           AI(2) = AI(2) -EDDYZ(I,1) * CRANKP * DZHI(i,2) * DZFI(i,1) * DTS
 
           DO K =  KCBL(I)+1, KL
             BI(K) = 1.0
           ENDDO
 
           DO K = 2,KL
             XPLUS(K)  = EDDYZ(I,K) * DZHI(i,K) * DZFI(i,K) * DTS
             XMINUS(K) = EDDYZ(I,K-1) * DZHI(i,K) * DZFI(i,K-1) * DTS
             CI(K)     = - XMINUS(K) * CRANKP
             EI(K)     = EI(K) - XPLUS(K) * CRANKP
             BI(K)     = BI(K) + XPLUS(K) * CRANKP + XMINUS(K) * CRANKP
           ENDDO
 
           IF (NOCONV(I) .EQ. 1) THEN
             BI(1) = 1.0 + CRANKP * MBARKS(1,I) * (PBL(I) - ZF(i,1)) * DTS   &
                   * DZHI(i,1) + EDDYZ(I,1) * CRANKP * DZHI(i,1) * DZFI(i,1) * DTS
           ELSE
             BI(1) = 1.0  + EDDYZ(I,1) * CRANKP * DZHI(i,1) * DZFI(i,1) * DTS
           ENDIF
 
 
           DO K = 1,KL
             DO L = 1,NSPX                    
               DI(L,K) = 0.0
             ENDDO
           ENDDO
 !
 !**   COMPUTE TENDENCY OF CBL CONCENTRATIONS - SEMI-IMPLICIT SOLUTION
           DO K = 2,KCBL(I)
             DO L = 1,NSPX                    
               DELC = DTS * (MBARKS(K,I) * VCI(L,I,1) - MDWN(K,I) *          &
                  VCI(L,I,K) + DZH(i,K+1) * DZHI(i,K) *                      &
                         MDWN(K+1,I) * VCI(L,I,K+1))
               DI(L,K)   = VCI(L,I,K) + (1.0 - CRANKP) * DELC
             ENDDO
           ENDDO
 
           DO K = KCBL(I)+1, KL
             DO L = 1,NSPX                    
               DI(L,K) = VCI(L,I,K)
             ENDDO
           ENDDO
 
           DO K = 2,KL
             IF (K .EQ. KL) THEN
               DO L = 1,NSPX                    
                 DI(L,K) = DI(L,K)  - (1.0 - CRANKP) * XMINUS(K) *           &
                           (VCI(L,I,K) - VCI(L,I,K-1))
               ENDDO
             ELSE
               DO L = 1,NSPX                    
                 DI(L,K) = DI(L,K) + (1.0 - CRANKP) * XPLUS(K) *             &
                           (VCI(L,I,K+1) - VCI(L,I,K))  -                    &
                           (1.0 - CRANKP) * XMINUS(K) *                      &
                           (VCI(L,I,K) - VCI(L,I,K-1))
               ENDDO
             ENDIF
           ENDDO
 
           IF (NOCONV(I) .EQ. 1) THEN
             DO L = 1,NSPX                    
               DI(L,1) = VCI(L,I,1) + (FS(L,I) - (1.0 - CRANKP)              &
                         * (MBARKS(1,I) *                                    &
                         (PBL(I) - ZF(i,1)) * VCI(L,I,1) -                   &
                         MDWN(2,I) * VCI(L,I,2) * DZH(i,2))) * DZHI(i,1) * DTS
             ENDDO
           ELSE
             DO L = 1,NSPX                    
               DI(L,1) = VCI(L,I,1) + FS(L,I) * DZHI(i,1) * DTS
             ENDDO
           ENDIF
           DO L = 1,NSPX                    
             DI(L,1) = DI(L,1) + (1.0 - CRANKP) * EDDYZ(I,1) * DZHI(i,1)     &
                     * DZFI(i,1) * DTS * (VCI(L,I,2) - VCI(L,I,1))
           ENDDO
           IF ( NOCONV(I) .EQ. 1 ) THEN
             CALL MATRIX (AI, BI, CI, DI, EI, UI, KL, NSPX)
           ELSE
             CALL TRI (CI, BI, EI, DI, UI, KL, NSPX)
           END IF
 !
 !-- COMPUTE NEW THETAV AND Q
           DO K = 1,KL
             DO L = 1,NSPX                    
               VCI(L,I,K) = UI(L,K)
             ENDDO
           ENDDO
 
         ENDDO                   ! END I LOOP
       ENDDO                     ! END SUB TIME LOOP
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !
       DO K = kts,KL
         DO I = its,ILX
           THETAX(I,K) = VCI(1,I,K)
           DO L = 2, NSPX
             !if (L == 5) then
               !write(*,*) i,k,l,VCI(L,I,K)
             !end if
             QXX(I,K,L-1) = VCI(L,I,K)
           END DO
         END DO
       END DO
       
 
       DEALLOCATE (DI)       
       DEALLOCATE (UI)  
       DEALLOCATE (FS)
       DEALLOCATE (VCI)
 
    END SUBROUTINE ACM
    
    SUBROUTINE ACMM (DTPBL,  NOCONV, ZF, DZH, DZHI,      &
                    KLPBL, PBL, DZFI, MOL,                     &
                    USTM, EDDYZM,                 &
                    US,    VS,                                      &
                    UX,    VX,                                      &
                    ids,ide, kds,kde,                      &
                    ims,ime, kms,kme,                      &
                    its,ite, kts,kte)
 !**********************************************************************
 !   PBL model called the Asymmetric Convective Model, Version 2 (ACM2) 
 !   -- See top of module for summary and references
 !
 !---- REVISION HISTORY:
 !   AX     3/2005 - developed WRF version based on ACM2 in the MM5 PX LSM
 !   JP and RG 8/2006 - updates
 !
 !**********************************************************************
 !  ARGUMENTS:
 !-- DTPBL           PBL time step
 !-- NOCONV          If free convection =0, no; =1, yes
 !-- ZF              Height for full layer
 !-- DZH             Layer thickness --> ZF(K) - ZF(K-1)
 !-- DZHI            Inverse of layer thickness
 !-- KLPBL           PBL level at K index
 !-- PBL             PBL height in m
 !-- DZFI            Inverse layer thickness --> 1/(Z(K+1)-Z(K))
 !-- MOL             Monin-Obukhov length in 1D form
 !-- USTM            U* for computation of momemtum flux 
 !-- EDDYZM          eddy diffusivity for momentum KM
 !-- US              U wind 
 !-- VS              V wind
 !-- THETA           potential temperature
 !-- QVS             water vapor mixing ratio (Kg/Kg)
 !-- QCS             cloud mixing ratio (Kg/Kg)
 !-- QIS             ice mixing ratio (Kg/Kg)
 !-- UX              new U wind 
 !-- VX              new V wind
 !-- THETAX          new potential temperature
 !-- QVX             new water vapor mixing ratio (Kg/Kg)
 !-- QCX             new cloud mixing ratio (Kg/Kg)
 !-- QIX             new ice mixing ratio (Kg/Kg)
 !-----------------------------------------------------------------------
 !-----------------------------------------------------------------------
 
       IMPLICIT NONE
 
 !.......Arguments
 
 !... Integer
       INTEGER,  INTENT(IN)      ::      ids,ide,kds,kde, &
                                         ims,ime,kms,kme, &
                                         its,ite,kts,kte
       INTEGER,  DIMENSION( its:ite ), INTENT(IN)  :: NOCONV
       INTEGER,  DIMENSION( ims:ime ), INTENT(IN)  :: KLPBL
 
 !... Real
       REAL , DIMENSION( ims:ime ),          INTENT(IN)  :: PBL
       REAL ,                                INTENT(IN)  :: DTPBL
       REAL , DIMENSION( its:ite ),          INTENT(IN)  :: MOL, USTM
       REAL , DIMENSION( its:ite, kts:kte ), INTENT(IN)  :: DZHI, DZH, DZFI
       REAL , DIMENSION( its:ite, 0:kte ),   INTENT(IN)  :: ZF
       REAL , DIMENSION( its:ite, kts:kte ), INTENT(INOUT)  :: EDDYZM
       REAL , DIMENSION( ims:ime, kms:kme ), INTENT(IN)  :: US, VS
       REAL , DIMENSION( its:ite, kts:kte ), INTENT(OUT) :: UX, VX
 !.......Local variables
 
 !... Parameters
       INTEGER, PARAMETER :: NSP   = 2
 !
       REAL,    PARAMETER :: XX    = 0.5          ! FACTOR APPLIED TO CONV MIXING TIME STEP
 
 !... Integer
       INTEGER :: ILX, KL, KLM, I, K, NSPX, NLP, NL, JJ, L
       INTEGER :: KCBLMX
       INTEGER, DIMENSION( its:ite ) :: KCBL
 
 !... Real
       REAL                               :: MBMAX, HOVL, MEDDY, MBAR
       REAL                               :: EKZ, RZ, FM, WSPD, DTS, DTRAT, F1
       REAL, DIMENSION( its:ite )         :: FSACM, DTLIM
       REAL, DIMENSION( kts:kte, its:ite) :: MBARKS, MDWN
       REAL, DIMENSION( 1:NSP, its:ite )  :: FS
       REAL, DIMENSION( kts:kte )         :: XPLUS, XMINUS
       REAL  DELC
       REAL, DIMENSION( 1:NSP,its:ite,kts:kte  ) :: VCI
 
       REAL, DIMENSION( kts:kte )               :: AI, BI, CI, EI !, Y
       REAL, DIMENSION( 1:NSP,kts:kte )         :: DI, UI    
 !
 !--Start Exicutable ----
 
       ILX = ite
       KL  = kte
       KLM = kte - 1
 
       KCBLMX = 0
       MBMAX  = 0.0
 
 !---COMPUTE ACM MIXING RATE
       DO I = its, ILX
         DTLIM(I)  = DTPBL
         KCBL(I)   = 1
         FSACM(I)  = 0.0
 
         IF (NOCONV(I) .EQ. 1) THEN
           KCBL(I) = KLPBL(I)
 
 !-------MBARKS IS UPWARD MIXING RATE; MDWN IS DOWNWARD MIXING RATE
 !--New couple ACM & EDDY-------------------------------------------------------------
           HOVL     = -PBL(I) / MOL(I)
           FSACM(I) = 1./(1.+((KARMAN/(HOVL))**0.3333)/(0.72*KARMAN))
           MEDDY    = EDDYZM(I,1) * DZFI(i,1) / (PBL(I) - ZF(i,1))
           MBAR     = MEDDY * FSACM(I)
           DO K = kts,KCBL(I)-1
             EDDYZM(I,K) = EDDYZM(I,K) * (1.0 - FSACM(I))
           ENDDO
 
           MBMAX = AMAX1(MBMAX,MBAR)
           DO K = kts+1,KCBL(I)
             MBARKS(K,I) = MBAR
             MDWN(K,I)   = MBAR * (PBL(I) - ZF(i,K-1)) * DZHI(i,K)
           ENDDO
           MBARKS(1,I) = MBAR
           MBARKS(KCBL(I),I) = MDWN(KCBL(I),I)
           MDWN(KCBL(I)+1,I) = 0.0
         ENDIF
       ENDDO                              ! end of I loop
 
       DO K = kts,KLM
         DO I = its,ILX
           EKZ   = EDDYZM(I,K) * DZFI(i,K) * DZHI(i,K)
           DTLIM(I) = AMIN1(0.75 / EKZ,DTLIM(I))
         ENDDO
       ENDDO
        
       DO I = its,ILX 
         IF (NOCONV(I) .EQ. 1) THEN
           KCBLMX = AMAX0(KLPBL(I),KCBLMX)
           RZ     = (ZF(i,KCBL(I)) - ZF(i,1)) * DZHI(i,1)
           DTLIM(I)  = AMIN1(XX / (MBARKS(1,I) * RZ),DTLIM(I))
         ENDIF
       ENDDO
 
       DO K = kts,KL
         DO I = its,ILX
           VCI(1,I,K) = US(I,K)
           VCI(2,I,K) = VS(I,K)
         ENDDO
       ENDDO
 
       NSPX=2
 
       DO I = its,ILX
         FM      = -USTM(I) * USTM(I)
         WSPD    = SQRT(US(I,1) * US(I,1) + VS(I,1) * VS(I,1)) + 1.E-9
         FS(1,I) = FM * US(I,1) / WSPD
         FS(2,I) = FM * VS(I,1) / WSPD
       ENDDO
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO I = its,ILX      
 
         NLP   = INT(DTPBL / DTLIM(I) + 1.0)
         DTS   = (DTPBL / NLP)
         DTRAT = DTS / DTPBL
         DO NL = 1,NLP           ! LOOP OVER SUB TIME LOOP              
 
 !-- COMPUTE ARRAY ELEMENTS THAT ARE INDEPENDANT OF SPECIES
 
           DO K = kts,KL
             AI(K) = 0.0
             BI(K) = 0.0
             CI(K) = 0.0
             EI(K) = 0.0
           ENDDO
 
           DO K = 2, KCBL(I)
             EI(K-1) = -CRANKP * MDWN(K,I) * DTS * DZH(i,K) * DZHI(i,K-1)
             BI(K)   = 1.0 + CRANKP * MDWN(K,I) * DTS
             AI(K)   = -CRANKP * MBARKS(K,I) * DTS
           ENDDO
 
           EI(1) = EI(1) -EDDYZM(I,1) * CRANKP * DZHI(i,1) * DZFI(i,1) * DTS
           AI(2) = AI(2) -EDDYZM(I,1) * CRANKP * DZHI(i,2) * DZFI(i,1) * DTS
 
           DO K =  KCBL(I)+1, KL
             BI(K) = 1.0
           ENDDO
 
           DO K = 2,KL
             XPLUS(K)  = EDDYZM(I,K) * DZHI(i,K) * DZFI(i,K) * DTS
             XMINUS(K) = EDDYZM(I,K-1) * DZHI(i,K) * DZFI(i,K-1) * DTS
             CI(K)     = - XMINUS(K) * CRANKP
             EI(K)     = EI(K) - XPLUS(K) * CRANKP
             BI(K)     = BI(K) + XPLUS(K) * CRANKP + XMINUS(K) * CRANKP
           ENDDO
 
           IF (NOCONV(I) .EQ. 1) THEN
             BI(1) = 1.0 + CRANKP * MBARKS(1,I) * (PBL(I) - ZF(i,1)) * DTS   &
                   * DZHI(i,1) + EDDYZM(I,1) * DZHI(i,1) * DZFI(i,1) * CRANKP * DTS
           ELSE
             BI(1) = 1.0  + EDDYZM(I,1) * DZHI(i,1) * DZFI(i,1) * CRANKP * DTS
           ENDIF
 
 
           DO K = 1,KL
             DO L = 1,NSPX                    
               DI(L,K) = 0.0
             ENDDO
           ENDDO
 !
 !**   COMPUTE TENDENCY OF CBL CONCENTRATIONS - SEMI-IMPLICIT SOLUTION
           DO K = 2,KCBL(I)
             DO L = 1,NSPX                    
               DELC = DTS * (MBARKS(K,I) * VCI(L,I,1) - MDWN(K,I) *          &
                  VCI(L,I,K) + DZH(i,K+1) * DZHI(i,K) *                      &
                         MDWN(K+1,I) * VCI(L,I,K+1))
               DI(L,K)   = VCI(L,I,K) + (1.0 - CRANKP) * DELC
             ENDDO
           ENDDO
 
           DO K = KCBL(I)+1, KL
             DO L = 1,NSPX                    
               DI(L,K) = VCI(L,I,K)
             ENDDO
           ENDDO
 
           DO K = 2,KL
             IF (K .EQ. KL) THEN
               DO L = 1,NSPX                    
                 DI(L,K) = DI(L,K)  - (1.0 - CRANKP) * XMINUS(K) *           &
                           (VCI(L,I,K) - VCI(L,I,K-1))
               ENDDO
             ELSE
               DO L = 1,NSPX                    
                 DI(L,K) = DI(L,K) + (1.0 - CRANKP) * XPLUS(K) *             &
                           (VCI(L,I,K+1) - VCI(L,I,K))  -                    &
                           (1.0 - CRANKP) * XMINUS(K) *                      &
                           (VCI(L,I,K) - VCI(L,I,K-1))
               ENDDO
             ENDIF
           ENDDO
 
           IF (NOCONV(I) .EQ. 1) THEN
             DO L = 1,NSPX                    
               DI(L,1) = VCI(L,I,1) + (FS(L,I) - (1.0 - CRANKP)              &
                         * (MBARKS(1,I) *                                    &
                         (PBL(I) - ZF(i,1)) * VCI(L,I,1) -                   &
                         MDWN(2,I) * VCI(L,I,2) * DZH(i,2))) * DZHI(i,1) * DTS
             ENDDO
           ELSE
             DO L = 1,NSPX                    
               DI(L,1) = VCI(L,I,1) + FS(L,I) * DZHI(i,1) * DTS
             ENDDO
           ENDIF
           DO L = 1,NSPX                    
             DI(L,1) = DI(L,1) + (1.0 - CRANKP) * EDDYZM(I,1) * DZHI(i,1)    &
                     * DZFI(i,1) * DTS * (VCI(L,I,2) - VCI(L,I,1))
           ENDDO
           IF ( NOCONV(I) .EQ. 1 ) THEN
             CALL MATRIX (AI, BI, CI, DI, EI, UI, KL, NSPX)
           ELSE
             CALL TRI (CI, BI, EI, DI, UI, KL, NSPX)
           END IF
 !
 !-- COMPUTE NEW THETAV AND Q
           DO K = 1,KL
             DO L = 1,NSPX                    
               VCI(L,I,K) = UI(L,K)
             ENDDO
           ENDDO
 
         ENDDO                   ! END I LOOP
       ENDDO                     ! END SUB TIME LOOP
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !
       DO K = kts,KL
         DO I = its,ILX
           UX(I,K)     = VCI(1,I,K)
           VX(I,K)     = VCI(2,I,K)
         ENDDO
       ENDDO
 
    END SUBROUTINE ACMM
    
    SUBROUTINE MATRIX(A,B,C,D,E,X,KL,NSP)
       
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
       IMPLICIT NONE
    !
    !-- Bordered band diagonal matrix solver for ACM2
    
    !-- ACM2 Matrix is in this form:
    !   B1 E1
    !   A2 B2 E2
    !   A3 C3 B3 E3
    !   A4    C4 B4 E4
    !   A5       C5 B5 E5
    !   A6          C6 B6
    
    !--Upper Matrix is
    !  U11 U12
    !      U22 U23
    !          U33 U34
    !              U44 U45
    !                  U55 U56
    !                      U66
    
    !--Lower Matrix is:
    !  1
    ! L21  1
    ! L31 L32  1
    ! L41 L42 L43  1
    ! L51 L52 L53 L54  1
    ! L61 L62 L63 L64 L65 1
    !---------------------------------------------------------
    !...Arguments
          INTEGER, INTENT(IN)   :: KL
          INTEGER, INTENT(IN)   :: NSP
          REAL A(KL),B(KL),E(KL)
          REAL C(KL),D(NSP,KL),X(NSP,KL)
    
    !...Locals
          REAL Y(NSP,KL),AIJ,SUM
          REAL L(KL,KL),UII(KL),UIIP1(KL),RUII(KL)
          INTEGER I,J,V
    
    !-- Define Upper and Lower matrices
          L(1,1) = 1.
          UII(1) = B(1)
          RUII(1) = 1./UII(1)
          DO I = 2, KL
    	      L(I,I) = 1.
    	      L(I,1) = A(I)/B(1)
            UIIP1(I-1)=E(I-1)
    	      IF(I.GE.3) THEN
    	        DO J = 2,I-1
    	          IF(I.EQ.J+1) THEN
    	            AIJ = C(I)
    	          ELSE
    	            AIJ = 0.
    	          ENDIF
    	          L(I,J) = (AIJ-L(I,J-1)*E(J-1))/      &
                          (B(J)-L(J,J-1)*E(J-1))
    	        ENDDO
    	      ENDIF
          ENDDO
    	  
          DO I = 2,KL
            UII(I) = B(I)-L(I,I-1)*E(I-1)
            RUII(I) = 1./UII(I)
          ENDDO
      
    !-- Forward sub for Ly=d
          DO V= 1, NSP
            Y(V,1) = D(V,1)
            DO I=2,KL
    	        SUM = D(V,I)
    	        DO J=1,I-1
    	          SUM = SUM - L(I,J)*Y(V,J)
    	        ENDDO
    	        Y(V,I) = SUM
            ENDDO
          ENDDO
    
    !-- Back sub for Ux=y
    
          DO V= 1, NSP
            X(V,KL) = Y(V,KL)*RUII(KL)
          ENDDO
          DO I = KL-1,1,-1
            DO V= 1, NSP
             X(V,I) = (Y(V,I)-UIIP1(I)*X(V,I+1))*RUII(I)
            ENDDO
          ENDDO
    
    END SUBROUTINE MATRIX
    
    SUBROUTINE TRI ( L, D, U, B, X,KL,NSP)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!  FUNCTION:
!    Solves tridiagonal system by Thomas algorithm. 
!   The associated tri-diagonal system is stored in 3 arrays
!   D : diagonal
!   L : sub-diagonal
!   U : super-diagonal
!   B : right hand side function
!   X : return solution from tridiagonal solver

!     [ D(1) U(1) 0    0    0 ...       0     ]
!     [ L(2) D(2) U(2) 0    0 ...       .     ]
!     [ 0    L(3) D(3) U(3) 0 ...       .     ]
!     [ .       .     .     .           .     ] X(i) = B(i)
!     [ .             .     .     .     0     ]
!     [ .                   .     .     .     ]
!     [ 0                           L(n) D(n) ]

!-----------------------------------------------------------------------

      IMPLICIT NONE

! Arguments:

      INTEGER, INTENT(IN)   :: KL
      INTEGER, INTENT(IN)   :: NSP

      REAL        L( KL )               ! subdiagonal
      REAL        D(KL)   ! diagonal
      REAL        U( KL )               ! superdiagonal
      REAL        B(NSP,KL )   ! R.H. side
      REAL        X( NSP,KL )   ! solution

! Local Variables:

      REAL        GAM( KL )
      REAL        BET
      INTEGER     V, K

! Decomposition and forward substitution:
      BET = 1.0 / D( 1 )
      DO V = 1, NSP
         X( V,1 ) = BET * B(V,1 )
      ENDDO

      DO K = 2, KL
        GAM(K ) = BET * U( K-1 )
        BET = 1.0 / ( D( K ) - L( K ) * GAM( K ) )
        DO V = 1, NSP
           X( V, K ) = BET * ( B( V,K ) - L( K ) * X( V,K-1 ) )
        ENDDO
      ENDDO

! Back-substitution:

      DO K = KL - 1, 1, -1
        DO V = 1, NSP
          X( V,K ) = X( V,K ) - GAM( K+1 ) * X( V,K+1 )
        ENDDO
      ENDDO
   
    END SUBROUTINE TRI
!-------------------------------------------------------------------------------
end module bl_acm
!-------------------------------------------------------------------------------
