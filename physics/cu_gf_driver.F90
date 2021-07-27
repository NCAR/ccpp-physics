!>\file cu_gf_driver.F90
!! This file is scale-aware Grell-Freitas cumulus scheme driver.


module cu_gf_driver

   ! DH* TODO: replace constants with arguments to cu_gf_driver_run
   use physcons  , g => con_g, cp => con_cp, xlv => con_hvap, r_v => con_rv
   use machine   , only: kind_phys
   use cu_gf_deep, only: cu_gf_deep_run,neg_check,autoconv,aeroevap,fct1d3
   use cu_gf_sh  , only: cu_gf_sh_run

   implicit none

   private

   public :: cu_gf_driver_init, cu_gf_driver_run, cu_gf_driver_finalize

contains

!> \brief Brief description of the subroutine
!!
!! \section arg_table_cu_gf_driver_init Argument Table
!! \htmlinclude cu_gf_driver_init.html
!!
      subroutine cu_gf_driver_init(imfshalcnv, imfshalcnv_gf, imfdeepcnv, &
                          imfdeepcnv_gf,mpirank, mpiroot, errmsg, errflg)

         implicit none
         
         integer,                   intent(in) :: imfshalcnv, imfshalcnv_gf
         integer,                   intent(in) :: imfdeepcnv, imfdeepcnv_gf       
         integer,                   intent(in)    :: mpirank
         integer,                   intent(in)    :: mpiroot
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg
         
         ! initialize ccpp error handling variables
         errmsg = ''
         errflg = 0
         
         ! DH* temporary
         if (mpirank==mpiroot) then
            write(0,*) ' -----------------------------------------------------------------------------------------------------------------------------'
            write(0,*) ' --- WARNING --- the CCPP Grell Freitas convection scheme is currently under development, use at your own risk --- WARNING ---'
            write(0,*) ' -----------------------------------------------------------------------------------------------------------------------------'
         end if
         ! *DH temporary

         ! Consistency checks
         if (.not. (imfshalcnv == imfshalcnv_gf .or.                       &
        &        imfdeepcnv == imfdeepcnv_gf)) then
           write(errmsg,'(*(a))') 'Logic error: namelist choice of',       &
        &    ' convection is different from Grell-Freitas scheme'
           errflg = 1
           return
         end if

      end subroutine cu_gf_driver_init

      subroutine cu_gf_driver_finalize()
      end subroutine cu_gf_driver_finalize
!
! t2di is temp after advection, but before physics
! t = current temp (t2di + physics up to now)
!===================

!> \defgroup cu_gf_group Grell-Freitas Convection Scheme Module
!! This is the Grell-Freitas scale and aerosol aware scheme.
!>\defgroup cu_gf_driver  Grell-Freitas Convection Scheme Driver Module
!> \ingroup cu_gf_group
!! This is the Grell-Freitas convection scheme driver module.
!! \section arg_table_cu_gf_driver_run Argument Table
!! \htmlinclude cu_gf_driver_run.html
!!
!>\section gen_gf_driver GSD GF Cumulus Scheme General Algorithm
!> @{
      subroutine cu_gf_driver_run(ntracer,garea,im,km,dt,cactiv,                &
               forcet,forceqv_spechum,phil,raincv,qv_spechum,t,cld1d,           &
               us,vs,t2di,w,qv2di_spechum,p2di,psuri,                           &
               hbot,htop,kcnv,xland,hfx2,qfx2,cliw,clcw,                        &
               pbl,ud_mf,dd_mf,dt_mf,cnvw_moist,cnvc,imfshalcnv,                &
               flag_for_scnv_generic_tend,flag_for_dcnv_generic_tend,           &
               dtend,dtidx,ntqv,ntiw,ntcw,index_of_temperature,index_of_x_wind, &
               index_of_y_wind,index_of_process_scnv,index_of_process_dcnv,     &
               ldiag3d,qci_conv,errmsg,errflg)
!-------------------------------------------------------------
      implicit none
      integer, parameter :: maxiens=1
      integer, parameter :: maxens=1
      integer, parameter :: maxens2=1
      integer, parameter :: maxens3=16
      integer, parameter :: ensdim=16
      integer, parameter :: imid_gf=1    ! testgf2 turn on middle gf conv.
      integer, parameter :: ideep=1
      integer, parameter :: ichoice=0	! 0 2 5 13 8
     !integer, parameter :: ichoicem=5	! 0 2 5 13
      integer, parameter :: ichoicem=13	! 0 2 5 13
      integer, parameter :: ichoice_s=3	! 0 1 2 3
      real(kind=kind_phys), parameter :: aodccn=0.1
      real(kind=kind_phys) :: dts,fpi,fp
      integer, parameter :: dicycle=0 ! diurnal cycle flag
      integer, parameter :: dicycle_m=0 !- diurnal cycle flag
      integer            :: ishallow_g3 ! depend on imfshalcnv
!-------------------------------------------------------------
   integer      :: its,ite, jts,jte, kts,kte 
   integer, intent(in   ) :: im,km,ntracer
   logical, intent(in   ) :: flag_for_scnv_generic_tend,flag_for_dcnv_generic_tend
   logical, intent(in   ) :: ldiag3d

   real(kind=kind_phys), optional, intent(inout)            :: dtend(:,:,:)
   integer, intent(in)                                      :: dtidx(:,:), &
        index_of_x_wind, index_of_y_wind, index_of_temperature,            &
        index_of_process_scnv, index_of_process_dcnv, ntqv, ntcw, ntiw
   
   real(kind=kind_phys),  dimension( : , : ), intent(in    ) :: forcet,forceqv_spechum,w,phil
   real(kind=kind_phys),  dimension( : , : ), intent(inout ) :: t,us,vs
   real(kind=kind_phys),  dimension( : , : ), intent(inout ) :: qci_conv
   real(kind=kind_phys),  dimension( : , : ), intent(out   ) :: cnvw_moist,cnvc
   real(kind=kind_phys),  dimension( : , : ), intent(inout ) :: cliw, clcw

   real(kind=kind_phys), allocatable :: clcw_save(:,:), cliw_save(:,:)

   integer, dimension (:), intent(out) :: hbot,htop,kcnv
   integer, dimension (:), intent(in)  :: xland
   real(kind=kind_phys),    dimension (:), intent(in) :: pbl
   integer, dimension (im) :: tropics
!  ruc variable
   real(kind=kind_phys), dimension (:),   intent(in)  :: hfx2,qfx2,psuri
   real(kind=kind_phys), dimension (:,:), intent(out) :: ud_mf,dd_mf,dt_mf
   real(kind=kind_phys), dimension (:),   intent(out) :: raincv,cld1d
   real(kind=kind_phys), dimension (:,:), intent(in)  :: t2di,p2di
   ! Specific humidity from FV3
   real(kind=kind_phys), dimension (:,:), intent(in) :: qv2di_spechum
   real(kind=kind_phys), dimension (:,:), intent(inout) :: qv_spechum
   ! Local water vapor mixing ratios and cloud water mixing ratios
   real(kind=kind_phys), dimension (im,km) :: qv2di, qv, forceqv, cnvw
   !
   real(kind=kind_phys), dimension(:),intent(in) :: garea
   real(kind=kind_phys), intent(in   ) :: dt 

   integer, intent(in   ) :: imfshalcnv
   integer, dimension(:), intent(inout) :: cactiv

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!  local variables
   integer, dimension(im) :: k22_shallow,kbcon_shallow,ktop_shallow
   real(kind=kind_phys), dimension (im)    :: rand_mom,rand_vmas
   real(kind=kind_phys), dimension (im,4)  :: rand_clos
   real(kind=kind_phys), dimension (im,km,11) :: gdc,gdc2
   real(kind=kind_phys), dimension (im)    :: ht
   real(kind=kind_phys), dimension (im)    :: dx
   real(kind=kind_phys), dimension (im,km) :: outt,outq,outqc,phh,subm,cupclw,cupclws
   real(kind=kind_phys), dimension (im,km) :: dhdt,zu,zus,zd,phf,zum,zdm,outum,outvm
   real(kind=kind_phys), dimension (im,km) :: outts,outqs,outqcs,outu,outv,outus,outvs
   real(kind=kind_phys), dimension (im,km) :: outtm,outqm,outqcm,submm,cupclwm
   real(kind=kind_phys), dimension (im,km) :: cnvwt,cnvwts,cnvwtm
   real(kind=kind_phys), dimension (im,km) :: hco,hcdo,zdo,zdd,hcom,hcdom,zdom
   real(kind=kind_phys), dimension    (km) :: zh
   real(kind=kind_phys), dimension (im)    :: tau_ecmwf,edt,edtm,edtd,ter11,aa0,xlandi
   real(kind=kind_phys), dimension (im)    :: pret,prets,pretm,hexec
   real(kind=kind_phys), dimension (im,10) :: forcing,forcing2

   integer, dimension (im) :: kbcon, ktop,ierr,ierrs,ierrm,kpbli
   integer, dimension (im) :: k22s,kbcons,ktops,k22,jmin,jminm
   integer, dimension (im) :: kbconm,ktopm,k22m

   integer :: iens,ibeg,iend,jbeg,jend,n
   integer :: ibegh,iendh,jbegh,jendh
   integer :: ibegc,iendc,jbegc,jendc,kstop
   real(kind=kind_phys), dimension(im,km) :: rho_dryar
   real(kind=kind_phys) :: pten,pqen,paph,zrho,pahfs,pqhfl,zkhvfl,pgeoh
   integer, parameter :: ipn = 0

!
! basic environmental input includes moisture convergence (mconv)
! omega (omeg), windspeed (us,vs), and a flag (ierr) to turn off
! convection for this call only and at that particular gridpoint
!
   real(kind=kind_phys), dimension (im,km) :: qcheck,zo,t2d,q2d,po,p2d,rhoi
   real(kind=kind_phys), dimension (im,km) :: tn,qo,tshall,qshall,dz8w,omeg
   real(kind=kind_phys), dimension (im)    :: ccn,z1,psur,cuten,cutens,cutenm
   real(kind=kind_phys), dimension (im)    :: umean,vmean,pmean
   real(kind=kind_phys), dimension (im)    :: xmbs,xmbs2,xmb,xmbm,xmb_dumm,mconv

   integer :: i,j,k,icldck,ipr,jpr,jpr_deep,ipr_deep,uidx,vidx,tidx,qidx
   integer :: itf,jtf,ktf,iss,jss,nbegin,nend,cliw_idx,clcw_idx
   integer :: high_resolution
   real(kind=kind_phys)    :: clwtot,clwtot1,excess,tcrit,tscl_kf,dp,dq,sub_spread,subcenter
   real(kind=kind_phys)    :: dsubclw,dsubclws,dsubclwm,dtime_max,ztm,ztq,hfm,qfm,rkbcon,rktop
   real(kind=kind_phys), dimension(km)   :: massflx,trcflx_in1,clw_in1,clw_ten1,po_cup
!  real(kind=kind_phys), dimension(km)   :: trcflx_in2,clw_in2,clw_ten2
   real(kind=kind_phys), dimension (im)  :: flux_tun,tun_rad_mid,tun_rad_shall,tun_rad_deep
   character*50 :: ierrc(im),ierrcm(im)
   character*50 :: ierrcs(im)
!  ruc variable
!  hfx2 -- sensible heat flux (k m/s), positive upward from sfc
!  qfx2 -- latent heat flux (kg/kg m/s), positive upward from sfc 
!  gf needs them in w/m2. define hfx and qfx after simple unit conversion
   real(kind=kind_phys), dimension (im)  :: hfx,qfx
   real(kind=kind_phys) tem,tem1,tf,tcr,tcrf
   real(kind=kind_phys) :: cliw_shal,clcw_shal,tem_shal, cliw_both, weight_sum
   real(kind=kind_phys) :: cliw_deep,clcw_deep,tem_deep, clcw_both
   integer :: cliw_deep_idx, clcw_deep_idx, cliw_shal_idx, clcw_shal_idx

  !parameter (tf=243.16, tcr=270.16, tcrf=1.0/(tcr-tf)) ! FV3 original
  !parameter (tf=263.16, tcr=273.16, tcrf=1.0/(tcr-tf))
  !parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
  parameter (tf=258.16, tcr=273.16, tcrf=1.0/(tcr-tf)) ! as fim, HCB tuning
  ! initialize ccpp error handling variables
     errmsg = ''
     errflg = 0

     if(ldiag3d) then
       if(flag_for_dcnv_generic_tend) then
         cliw_deep_idx=0
         clcw_deep_idx=0
       else
         cliw_deep_idx=dtidx(100+ntiw,index_of_process_dcnv)
         clcw_deep_idx=dtidx(100+ntcw,index_of_process_dcnv)
       endif
       if(flag_for_scnv_generic_tend) then
         cliw_shal_idx=0
         clcw_shal_idx=0
       else
         cliw_shal_idx=dtidx(100+ntiw,index_of_process_scnv)
         clcw_shal_idx=dtidx(100+ntcw,index_of_process_scnv)
       endif
       if(cliw_deep_idx>=1 .or. clcw_deep_idx>=1 .or. &
            cliw_shal_idx>=1 .or.  clcw_shal_idx>=1) then
         allocate(clcw_save(im,km), cliw_save(im,km))
         clcw_save=clcw
         cliw_save=cliw
       endif
     endif

!
! Scale specific humidity to dry mixing ratio
!
     ! state in before physics
     qv2di = qv2di_spechum/(1.0_kind_phys-qv2di_spechum)
     ! forcing by dynamics, based on state in
     forceqv = forceqv_spechum/(1.0_kind_phys-qv2di_spechum)
     ! current state (updated by preceeding physics)
     qv = qv_spechum/(1.0_kind_phys-qv_spechum)
!
!
! these should be coming in from outside
!
!    cactiv(:)      = 0
     rand_mom(:)    = 0.
     rand_vmas(:)   = 0.
     rand_clos(:,:) = 0.
!
     its=1
     ite=im
     itf=ite
     jts=1
     jte=1
     jtf=jte
     kts=1
     kte=km
     ktf=kte-1
! 
     tropics(:)=0
!
!> - Set tuning constants for radiation coupling
!
     tun_rad_shall(:)=.02
     tun_rad_mid(:)=.15
     tun_rad_deep(:)=.13
     edt(:)=0.
     edtm(:)=0.
     edtd(:)=0.
     zdd(:,:)=0.
     flux_tun(:)=5.
! 10/11/2016 dx and tscl_kf are replaced with input dx(i), is dlength. 
! dx for scale awareness
!    dx=40075000./float(lonf)
!    tscl_kf=dx/25000.
     ccn(its:ite)=150.
  
     if (imfshalcnv == 3) then
      ishallow_g3 = 1
     else
      ishallow_g3 = 0
     end if
     high_resolution=0
     subcenter=0.
     iens=1
!
! these can be set for debugging
!
     ipr=0
     jpr=0
     ipr_deep=0
     jpr_deep= 0 !53322 ! 528196 !0 ! 1136 !0 !421755 !3536
!
!
     ibeg=its
     iend=ite
     tcrit=258.

     ztm=0.
     ztq=0.
     hfm=0.
     qfm=0.
     ud_mf =0.
     dd_mf =0.
     dt_mf =0.
     tau_ecmwf(:)=0.
!
     j=1
     ht(:)=phil(:,1)/g
     do i=its,ite
      cld1d(i)=0.
      zo(i,:)=phil(i,:)/g
      dz8w(i,1)=zo(i,2)-zo(i,1)
      zh(1)=0.
      kpbli(i)=2
      do k=kts+1,ktf
       dz8w(i,k)=zo(i,k+1)-zo(i,k)
      enddo
      do k=kts+1,ktf
       zh(k)=zh(k-1)+dz8w(i,k-1)
       if(zh(k).gt.pbl(i))then
        kpbli(i)=max(2,k)
        exit
       endif
      enddo
     enddo

     do i= its,itf
      forcing(i,:)=0.
      forcing2(i,:)=0.
      ccn(i)=100.
      hbot(i)  =kte
      htop(i)  =kts
      raincv(i)=0.
      xlandi(i)=real(xland(i))
!     if(abs(xlandi(i)-1.).le.1.e-3) tun_rad_shall(i)=.15     
!     if(abs(xlandi(i)-1.).le.1.e-3) flux_tun(i)=1.5     
     enddo
     do i= its,itf
      mconv(i)=0.
     enddo
     do k=kts,kte
      do i= its,itf
       omeg(i,k)=0.
       zu(i,k)=0.
       zum(i,k)=0.
       zus(i,k)=0.
       zd(i,k)=0.
       zdm(i,k)=0.
      enddo
     enddo

     psur(:)=0.01*psuri(:)
     do i=its,itf
      ter11(i)=max(0.,ht(i))
     enddo
     do k=kts,kte
      do i=its,ite
       cnvw(i,k)=0.
       cnvc(i,k)=0.
       gdc(i,k,1)=0.
       gdc(i,k,2)=0.
       gdc(i,k,3)=0.
       gdc(i,k,4)=0.
       gdc(i,k,7)=0.
       gdc(i,k,8)=0.
       gdc(i,k,9)=0.
       gdc(i,k,10)=0.
       gdc2(i,k,1)=0.
      enddo
     enddo
     ierr(:)=0
     ierrm(:)=0
     ierrs(:)=0
     cuten(:)=0.
     cutenm(:)=0.
     cutens(:)=0.
     ierrc(:)=" "

     kbcon(:)=0
     kbcons(:)=0
     kbconm(:)=0

     ktop(:)=0
     ktops(:)=0
     ktopm(:)=0

     xmb(:)=0.
     xmb_dumm(:)=0.
     xmbm(:)=0.
     xmbs(:)=0.
     xmbs2(:)=0.

     k22s(:)=0
     k22m(:)=0
     k22(:)=0

     jmin(:)=0
     jminm(:)=0

     pret(:)=0.
     prets(:)=0.
     pretm(:)=0.

     umean(:)=0.
     vmean(:)=0.
     pmean(:)=0.

     cupclw(:,:)=0.
     cupclwm(:,:)=0.
     cupclws(:,:)=0.

     cnvwt(:,:)=0.
     cnvwts(:,:)=0.

     hco(:,:)=0.
     hcom(:,:)=0.
     hcdo(:,:)=0.
     hcdom(:,:)=0.

     outt(:,:)=0.
     outts(:,:)=0.
     outtm(:,:)=0.

     outu(:,:)=0.
     outus(:,:)=0.
     outum(:,:)=0.

     outv(:,:)=0.
     outvs(:,:)=0.
     outvm(:,:)=0.

     outq(:,:)=0.
     outqs(:,:)=0.
     outqm(:,:)=0.

     outqc(:,:)=0.
     outqcs(:,:)=0.
     outqcm(:,:)=0.

     subm(:,:)=0.
     dhdt(:,:)=0.
     
     do k=kts,ktf
      do i=its,itf
        p2d(i,k)=0.01*p2di(i,k)
        po(i,k)=p2d(i,k) !*.01
        rhoi(i,k) = 100.*p2d(i,k)/(287.04*(t2di(i,k)*(1.+0.608*qv2di(i,k))))
        qcheck(i,k)=qv(i,k)
        tn(i,k)=t(i,k)!+forcet(i,k)*dt
        qo(i,k)=max(1.e-16,qv(i,k))!+forceqv(i,k)*dt
        t2d(i,k)=t2di(i,k)-forcet(i,k)*dt
        q2d(i,k)=max(1.e-16,qv2di(i,k)-forceqv(i,k)*dt)
        if(qo(i,k).lt.1.e-16)qo(i,k)=1.e-16
        tshall(i,k)=t2d(i,k)
        qshall(i,k)=q2d(i,k)
      enddo
     enddo
123  format(1x,i2,1x,2(1x,f8.0),1x,2(1x,f8.3),3(1x,e13.5))
     do i=its,itf
      do k=kts,kpbli(i)
         tshall(i,k)=t(i,k)
         qshall(i,k)=max(1.e-16,qv(i,k))
      enddo
     enddo
!
! converting hfx2 and qfx2 to w/m2
!    hfx=cp*rho*hfx2
!    qfx=xlv*qfx2
     do i=its,itf
      hfx(i)=hfx2(i)*cp*rhoi(i,1)
      qfx(i)=qfx2(i)*xlv*rhoi(i,1)
      dx(i) = sqrt(garea(i))
     enddo

     do i=its,itf
      do k=kts,kpbli(i)
       tn(i,k)=t(i,k) 
       qo(i,k)=max(1.e-16,qv(i,k))
      enddo
     enddo
     nbegin=0
     nend=0
     do i=its,itf
      do k=kts,kpbli(i)
       dhdt(i,k)=cp*(forcet(i,k)+(t(i,k)-t2di(i,k))/dt) +  & 
                 xlv*(forceqv(i,k)+(qv(i,k)-qv2di(i,k))/dt) 
!      tshall(i,k)=t(i,k) 
!      qshall(i,k)=qv(i,k) 
      enddo
     enddo
     do k=  kts+1,ktf-1
      do i = its,itf
       if((p2d(i,1)-p2d(i,k)).gt.150.and.p2d(i,k).gt.300)then
         dp=-.5*(p2d(i,k+1)-p2d(i,k-1))
         umean(i)=umean(i)+us(i,k)*dp
         vmean(i)=vmean(i)+vs(i,k)*dp
         pmean(i)=pmean(i)+dp
       endif
      enddo
     enddo
     do k=kts,ktf-1
      do i = its,itf
        omeg(i,k)= w(i,k) !-g*rhoi(i,k)*w(i,k)
!       dq=(q2d(i,k+1)-q2d(i,k))
!       mconv(i)=mconv(i)+omeg(i,k)*dq/g
      enddo
     enddo
     do i = its,itf
      if(mconv(i).lt.0.)mconv(i)=0.
     enddo
!
!---- call cumulus parameterization
!
       if(ishallow_g3.eq.1)then

          do i=its,ite
           ierrs(i)=0
           ierrm(i)=0
          enddo
!
!> - Call shallow: cu_gf_sh_run()
!
          call cu_gf_sh_run (us,vs,                                              &
! input variables, must be supplied
                         zo,t2d,q2d,ter11,tshall,qshall,p2d,psur,dhdt,kpbli,     &
                         rhoi,hfx,qfx,xlandi,ichoice_s,tcrit,dt,                 &
! input variables. ierr should be initialized to zero or larger than zero for
! turning off shallow convection for grid points
                         zus,xmbs,kbcons,ktops,k22s,ierrs,ierrcs,                &
! output tendencies
                         outts,outqs,outqcs,outus,outvs,cnvwt,prets,cupclws,     &
! dimesnional variables
                         itf,ktf,its,ite, kts,kte,ipr,tropics)


          do i=its,itf
           if(xmbs(i).gt.0.)cutens(i)=1.
          enddo
!> - Call neg_check() for GF shallow convection
          call neg_check('shallow',ipn,dt,qcheck,outqs,outts,outus,outvs,   &
                                 outqcs,prets,its,ite,kts,kte,itf,ktf,ktops)
       endif

       ipr=0
       jpr_deep=0 !340765
!> - Call cu_gf_deep_run() for middle GF convection
      if(imid_gf == 1)then
       call cu_gf_deep_run(        &
               itf,ktf,its,ite, kts,kte  &
              ,dicycle_m       &
              ,ichoicem       &
              ,ipr           &
              ,ccn           &
              ,dt            &
              ,imid_gf       &
              ,kpbli         &
              ,dhdt          &
              ,xlandi        &

              ,zo            &
              ,forcing2      &
              ,t2d           &
              ,q2d           &
              ,ter11         &
              ,tshall        &
              ,qshall        &
              ,p2d          &
              ,psur          &
              ,us            &
              ,vs            &
              ,rhoi          &
              ,hfx           &
              ,qfx           &
              ,dx            & !hj dx(im)
              ,mconv         &
              ,omeg          &

              ,cactiv        &
              ,cnvwtm        &
              ,zum           &
              ,zdm           & ! hli
              ,zdd           &
              ,edtm          &
              ,edtd          & ! hli
              ,xmbm          &
              ,xmb_dumm      &
              ,xmbs          &
              ,pretm         &
              ,outum         &
              ,outvm         &
              ,outtm         &
              ,outqm         &
              ,outqcm        &
              ,kbconm        &
              ,ktopm         &
              ,cupclwm       &
              ,ierrm         &
              ,ierrcm        &
!    the following should be set to zero if not available
              ,rand_mom      & ! for stochastics mom, if temporal and spatial patterns exist
              ,rand_vmas     & ! for stochastics vertmass, if temporal and spatial patterns exist
              ,rand_clos     & ! for stochastics closures, if temporal and spatial patterns exist
              ,0             & ! flag to what you want perturbed
                               ! 1 = momentum transport 
                               ! 2 = normalized vertical mass flux profile
                               ! 3 = closures
                               ! more is possible, talk to developer or
                               ! implement yourself. pattern is expected to be
                               ! betwee -1 and +1
#if ( wrf_dfi_radar == 1 )
              ,do_capsuppress,cap_suppress_j &
#endif
              ,k22m          &
              ,jminm,tropics)

            do i=its,itf
             do k=kts,ktf
              qcheck(i,k)=qv(i,k) +outqs(i,k)*dt
             enddo
            enddo
!> - Call neg_check() for middle GF convection
      call neg_check('mid',ipn,dt,qcheck,outqm,outtm,outum,outvm,   &
                     outqcm,pretm,its,ite,kts,kte,itf,ktf,ktopm)
     endif
!> - Call cu_gf_deep_run() for deep GF convection
     if(ideep.eq.1)then
      call cu_gf_deep_run(        &
               itf,ktf,its,ite, kts,kte  &

              ,dicycle       &
              ,ichoice       &
              ,ipr           &
              ,ccn           &
              ,dt            &
              ,0             &

              ,kpbli         &
              ,dhdt          &
              ,xlandi        &

              ,zo            &
              ,forcing       &
              ,t2d           &
              ,q2d           &
              ,ter11         &
              ,tn            &
              ,qo            &
              ,p2d           &
              ,psur          &
              ,us            &
              ,vs            &
              ,rhoi          &
              ,hfx           &
              ,qfx           &
              ,dx            & !hj replace dx(im)
              ,mconv         &
              ,omeg          &

              ,cactiv       &
              ,cnvwt        &
              ,zu           &
              ,zd           &
              ,zdm          & ! hli
              ,edt          &
              ,edtm         & ! hli
              ,xmb          &
              ,xmbm         &
              ,xmbs         &
              ,pret         &
              ,outu         &
              ,outv         &
              ,outt         &
              ,outq         &
              ,outqc        &
              ,kbcon        &
              ,ktop         &
              ,cupclw       &
              ,ierr         &
              ,ierrc        &
!    the following should be set to zero if not available
              ,rand_mom      & ! for stochastics mom, if temporal and spatial patterns exist
              ,rand_vmas     & ! for stochastics vertmass, if temporal and spatial patterns exist
              ,rand_clos     & ! for stochastics closures, if temporal and spatial patterns exist
              ,0             & ! flag to what you want perturbed
                               ! 1 = momentum transport 
                               ! 2 = normalized vertical mass flux profile
                               ! 3 = closures
                               ! more is possible, talk to developer or
                               ! implement yourself. pattern is expected to be
                               ! betwee -1 and +1
#if ( wrf_dfi_radar == 1 )
              ,do_capsuppress,cap_suppress_j &
#endif
              ,k22          &
              ,jmin,tropics)
          jpr=0
          ipr=0
          do i=its,itf
           do k=kts,ktf
            qcheck(i,k)=qv(i,k) +(outqs(i,k)+outqm(i,k))*dt
           enddo
          enddo
!> - Call neg_check() for deep GF convection
       call neg_check('deep',ipn,dt,qcheck,outq,outt,outu,outv,   &
                      outqc,pret,its,ite,kts,kte,itf,ktf,ktop)
!
      endif
!            do i=its,itf
!              kcnv(i)=0  
!              if(pret(i).gt.0.)then
!                 cuten(i)=1.
!                 kcnv(i)= 1 !jmin(i) 
!              else 
!                 kbcon(i)=0
!                 ktop(i)=0
!                 cuten(i)=0.
!              endif   ! pret > 0
!              if(pretm(i).gt.0.)then
!                 kcnv(i)= 1 !jmin(i)  
!                 cutenm(i)=1.
!              else 
!                 kbconm(i)=0
!                 ktopm(i)=0
!                 cutenm(i)=0.
!              endif   ! pret > 0
!            enddo
            do i=its,itf
              kcnv(i)=0
              if(pretm(i).gt.0.)then
                 kcnv(i)= 1 !jmin(i)  
                 cutenm(i)=1.
              else
                 kbconm(i)=0
                 ktopm(i)=0
                 cutenm(i)=0.
              endif   ! pret > 0

              if(pret(i).gt.0.)then
                 cuten(i)=1.
                 cutenm(i)=0.
                 pretm(i)=0.
                 kcnv(i)= 1 !jmin(i) 
                 ktopm(i)=0
                 kbconm(i)=0
              else
                 kbcon(i)=0
                 ktop(i)=0
                 cuten(i)=0.
              endif   ! pret > 0
            enddo
!
            do i=its,itf
            massflx(:)=0.
            trcflx_in1(:)=0.
            clw_in1(:)=0.
            clw_ten1(:)=0.
            po_cup(:)=0.
            kstop=kts
            if(ktopm(i).gt.kts .or. ktop(i).gt.kts)kstop=max(ktopm(i),ktop(i))
            if(ktops(i).gt.kts)kstop=max(kstop,ktops(i))
            if(kstop.gt.2)then
            htop(i)=kstop
            if(kbcon(i).gt.2 .or. kbconm(i).gt.2)then
               hbot(i)=max(kbconm(i),kbcon(i)) !jmin(i)
            endif

            dtime_max=dt
            do k=kts,kstop
               cnvc(i,k) = 0.04 * log(1. + 675. * zu(i,k) * xmb(i)) +   &
                           0.04 * log(1. + 675. * zum(i,k) * xmbm(i)) + &
                           0.04 * log(1. + 675. * zus(i,k) * xmbs(i))
               cnvc(i,k) = min(cnvc(i,k), 0.6)
               cnvc(i,k) = max(cnvc(i,k), 0.0)
               cnvw(i,k)=cnvwt(i,k)*xmb(i)*dt+cnvwts(i,k)*xmbs(i)*dt+cnvwtm(i,k)*xmbm(i)*dt
               ud_mf(i,k)=cuten(i)*zu(i,k)*xmb(i)*dt
               dd_mf(i,k)=cuten(i)*zd(i,k)*edt(i)*xmb(i)*dt
               t(i,k)=t(i,k)+dt*(cutens(i)*outts(i,k)+cutenm(i)*outtm(i,k)+outt(i,k)*cuten(i))
               qv(i,k)=max(1.e-16,qv(i,k)+dt*(cutens(i)*outqs(i,k)+cutenm(i)*outqm(i,k)+outq(i,k)*cuten(i)))
               gdc(i,k,7)=sqrt(us(i,k)**2 +vs(i,k)**2)
               us(i,k)=us(i,k)+outu(i,k)*cuten(i)*dt +outum(i,k)*cutenm(i)*dt +outus(i,k)*cutens(i)*dt
               vs(i,k)=vs(i,k)+outv(i,k)*cuten(i)*dt +outvm(i,k)*cutenm(i)*dt +outvs(i,k)*cutens(i)*dt

               gdc(i,k,1)= max(0.,tun_rad_shall(i)*cupclws(i,k)*cutens(i))      ! my mod
               gdc2(i,k,1)=max(0.,tun_rad_deep(i)*(cupclwm(i,k)*cutenm(i)+cupclw(i,k)*cuten(i)))
               qci_conv(i,k)=gdc2(i,k,1)
               gdc(i,k,2)=(outt(i,k))*86400.
               gdc(i,k,3)=(outtm(i,k))*86400. 
               gdc(i,k,4)=(outts(i,k))*86400.
               gdc(i,k,7)=-(gdc(i,k,7)-sqrt(us(i,k)**2 +vs(i,k)**2))/dt
              !gdc(i,k,8)=(outq(i,k))*86400.*xlv/cp
               gdc(i,k,8)=(outqm(i,k)+outqs(i,k)+outq(i,k))*86400.*xlv/cp 
               gdc(i,k,9)=gdc(i,k,2)+gdc(i,k,3)+gdc(i,k,4)
!
!> - Calculate subsidence effect on clw
!
!              dsubclw=0.
!              dsubclwm=0.
!              dsubclws=0.
!              dp=100.*(p2d(i,k)-p2d(i,k+1))
!              if (clcw(i,k) .gt. -999.0 .and. clcw(i,k+1) .gt. -999.0 )then
!                 clwtot = cliw(i,k) + clcw(i,k)
!                 clwtot1= cliw(i,k+1) + clcw(i,k+1)
!                 dsubclw=((-edt(i)*zd(i,k+1)+zu(i,k+1))*clwtot1   &
!                      -(-edt(i)*zd(i,k)  +zu(i,k))  *clwtot  )*g/dp
!                 dsubclwm=((-edtm(i)*zdm(i,k+1)+zum(i,k+1))*clwtot1   &
!                      -(-edtm(i)*zdm(i,k)  +zum(i,k))  *clwtot  )*g/dp
!                 dsubclws=(zus(i,k+1)*clwtot1-zus(i,k)*clwtot)*g/dp
!                 dsubclw=dsubclw+(zu(i,k+1)*clwtot1-zu(i,k)*clwtot)*g/dp 
!                 dsubclwm=dsubclwm+(zum(i,k+1)*clwtot1-zum(i,k)*clwtot)*g/dp 
!                 dsubclws=dsubclws+(zus(i,k+1)*clwtot1-zus(i,k)*clwtot)*g/dp 
!              endif
!              tem  = dt*(outqcs(i,k)*cutens(i)+outqc(i,k)*cuten(i)       &
!                    +outqcm(i,k)*cutenm(i)                           &
!                     +dsubclw*xmb(i)+dsubclws*xmbs(i)+dsubclwm*xmbm(i) &
!                    )
!              tem1 = max(0.0, min(1.0, (tcr-t(i,k))*tcrf))
!              if (clcw(i,k) .gt. -999.0) then
!               cliw(i,k) = max(0.,cliw(i,k) + tem * tem1)            ! ice
!               clcw(i,k) = max(0.,clcw(i,k) + tem *(1.0-tem1))       ! water
!              else
!                cliw(i,k) = max(0.,cliw(i,k) + tem)
!              endif
!
!            enddo

!> - FCT treats subsidence effect to cloud ice/water (begin)
               dp=100.*(p2d(i,k)-p2d(i,k+1))
               dtime_max=min(dtime_max,.5*dp)
               po_cup(k)=.5*(p2d(i,k)+p2d(i,k+1))
               if (clcw(i,k) .gt. -999.0 .and. clcw(i,k+1) .gt. -999.0 )then
                  clwtot = cliw(i,k) + clcw(i,k)
                  if(clwtot.lt.1.e-32)clwtot=0.
                  clwtot1= cliw(i,k+1) + clcw(i,k+1)
                  if(clwtot1.lt.1.e-32)clwtot1=0.
                  clw_in1(k)=clwtot
                  massflx(k)=-(xmb(i) *( zu(i,k)- edt(i)* zd(i,k)))   &
                             -(xmbm(i)*(zdm(i,k)-edtm(i)*zdm(i,k)))   &
                             -(xmbs(i)*zus(i,k))
                  trcflx_in1(k)=massflx(k)*.5*(clwtot+clwtot1)
               endif
             enddo

             massflx   (1)=0.
             trcflx_in1(1)=0.
             call fct1d3 (kstop,kte,dtime_max,po_cup,                  &
                            clw_in1,massflx,trcflx_in1,clw_ten1,g)

             do k=1,kstop
               tem  = dt*(outqcs(i,k)*cutens(i)+outqc(i,k)*cuten(i)    &
                      +outqcm(i,k)*cutenm(i)                           &
                      +clw_ten1(k)                                     &
                         )
               tem1 = max(0.0, min(1.0, (tcr-t(i,k))*tcrf))
               if (clcw(i,k) .gt. -999.0) then
                cliw(i,k) = max(0.,cliw(i,k) + tem * tem1)            ! ice
                clcw(i,k) = max(0.,clcw(i,k) + tem *(1.0-tem1))       ! water
               else
                cliw(i,k) = max(0.,cliw(i,k) + tem)
               endif

             enddo

            gdc(i,1,10)=forcing(i,1)
            gdc(i,2,10)=forcing(i,2)
            gdc(i,3,10)=forcing(i,3)
            gdc(i,4,10)=forcing(i,4)
            gdc(i,5,10)=forcing(i,5)
            gdc(i,6,10)=forcing(i,6)
            gdc(i,7,10)=forcing(i,7)
            gdc(i,8,10)=forcing(i,8)
            gdc(i,10,10)=xmb(i)
            gdc(i,11,10)=xmbm(i)
            gdc(i,12,10)=xmbs(i)
            gdc(i,13,10)=hfx(i)
            gdc(i,15,10)=qfx(i)
            gdc(i,16,10)=pret(i)*3600.
            if(ktop(i).gt.2 .and.pret(i).gt.0.)dt_mf(i,ktop(i)-1)=ud_mf(i,ktop(i))
            endif
            enddo
            do i=its,itf
              if(pret(i).gt.0.)then
                 cactiv(i)=1
                 raincv(i)=.001*(cutenm(i)*pretm(i)+cutens(i)*prets(i)+cuten(i)*pret(i))*dt
              else
                 cactiv(i)=0
                 if(pretm(i).gt.0)raincv(i)=.001*cutenm(i)*pretm(i)*dt
              endif   ! pret > 0
            enddo
 100    continue
!
! Scale dry mixing ratios for water wapor and cloud water to specific humidy / moist mixing ratios
!
        qv_spechum = qv/(1.0_kind_phys+qv)
        cnvw_moist = cnvw/(1.0_kind_phys+qv)
!
! Diagnostic tendency updates
!
        if(ldiag3d) then
          if(ishallow_g3.eq.1 .and. .not.flag_for_scnv_generic_tend) then
            uidx=dtidx(index_of_x_wind,index_of_process_scnv)
            vidx=dtidx(index_of_y_wind,index_of_process_scnv)
            tidx=dtidx(index_of_temperature,index_of_process_scnv)
            qidx=dtidx(100+ntqv,index_of_process_scnv)
            if(uidx>=1) then
              do k=kts,ktf
                dtend(:,k,uidx) = dtend(:,k,uidx) + cutens(:)*outus(:,k) * dt
              enddo
            endif
            if(vidx>=1) then
              do k=kts,ktf
                dtend(:,k,vidx) = dtend(:,k,vidx) + cutens(:)*outvs(:,k) * dt
              enddo
            endif
            if(tidx>=1) then
              do k=kts,ktf
                dtend(:,k,tidx) = dtend(:,k,tidx) + cutens(:)*outts(:,k) * dt
              enddo
            endif
            if(qidx>=1) then
              do k=kts,ktf
                do i=its,itf
                  tem = cutens(i)*outqs(i,k)* dt
                  tem = tem/(1.0_kind_phys+tem)
                  dtend(i,k,qidx) = dtend(i,k,qidx) + tem
                enddo
              enddo
            endif
          endif
          if((ideep.eq.1. .or. imid_gf.eq.1) .and. .not.flag_for_dcnv_generic_tend) then
            uidx=dtidx(index_of_x_wind,index_of_process_dcnv)
            vidx=dtidx(index_of_y_wind,index_of_process_dcnv)
            tidx=dtidx(index_of_temperature,index_of_process_dcnv)
            if(uidx>=1) then
              do k=kts,ktf
                dtend(:,k,uidx) = dtend(:,k,uidx) + (cuten*outu(:,k)+cutenm*outum(:,k)) * dt
              enddo
            endif
            if(vidx>=1) then
              do k=kts,ktf
                dtend(:,k,vidx) = dtend(:,k,vidx) + (cuten*outv(:,k)+cutenm*outvm(:,k)) * dt
              enddo
            endif
            if(tidx>=1) then
              do k=kts,ktf
                dtend(:,k,tidx) = dtend(:,k,tidx) + (cuten*outt(:,k)+cutenm*outtm(:,k)) * dt
              enddo
            endif

            qidx=dtidx(100+ntqv,index_of_process_dcnv)
            if(qidx>=1) then
              do k=kts,ktf
                do i=its,itf
                  tem = (cuten(i)*outq(i,k) + cutenm(i)*outqm(i,k))* dt
                  tem = tem/(1.0_kind_phys+tem)
                  dtend(i,k,qidx) = dtend(i,k,qidx) + tem
                enddo
              enddo
            endif
          endif
          if(allocated(clcw_save)) then
            do k=kts,ktf
              do i=its,itf
                tem_shal = dt*(outqcs(i,k)*cutens(i)+outqcm(i,k)*cutenm(i))
                tem_deep = dt*(outqc(i,k)*cuten(i)+clw_ten1(k))
                tem  = tem_shal+tem_deep
                tem1 = max(0.0, min(1.0, (tcr-t(i,k))*tcrf))
                weight_sum = abs(tem_shal)+abs(tem_deep)
                if(weight_sum<1e-12) then
                  cycle
                endif
                
                if (clcw_save(i,k) .gt. -999.0) then
                  cliw_both = max(0.,cliw_save(i,k) + tem * tem1) - cliw_save(i,k)
                  clcw_both = max(0.,clcw_save(i,k) + tem) - clcw_save(i,k)
                else if(cliw_idx>=1) then
                  cliw_both = max(0.,cliw_save(i,k) + tem) - cliw_save(i,k)
                  clcw_both = 0
                endif
                if(cliw_deep_idx>=1) then
                  dtend(i,k,cliw_deep_idx) = dtend(i,k,cliw_deep_idx) + abs(tem_deep)/weight_sum*cliw_both
                endif
                if(clcw_deep_idx>=1) then
                  dtend(i,k,clcw_deep_idx) = dtend(i,k,clcw_deep_idx) + abs(tem_deep)/weight_sum*clcw_both
                endif
                if(cliw_shal_idx>=1) then
                  dtend(i,k,cliw_shal_idx) = dtend(i,k,cliw_shal_idx) + abs(tem_shal)/weight_sum*cliw_both
                endif
                if(clcw_shal_idx>=1) then
                  dtend(i,k,clcw_shal_idx) = dtend(i,k,clcw_shal_idx) + abs(tem_shal)/weight_sum*clcw_both
                endif
              enddo
            enddo
          endif
        endif
   end subroutine cu_gf_driver_run
!> @}
end module cu_gf_driver
