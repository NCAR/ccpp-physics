!>\file sgscloud_radpre.F90
!!  Contains the preliminary (interstitial) work to the call to the radiation schemes:
!!    1) Backs up the original qc & qi
!!    2) Adds the partioning of convective condensate into liqice/ice for effective radii
!!    3) Adds the subgrid clouds mixing ratio and cloud fraction to the original (resolved-
!!       scale) qc, qi and cloud fraction coming from the microphysics scheme.
!!    4) Recompute the diagnostic high, mid, low, total and bl clouds to be consistent with radiation

      module sgscloud_radpre

      contains

!> \defgroup sgsradpre_group sgscloud_radpre_run Module
!! This interstitial code adds the subgrid clouds to the resolved-scale clouds
!! if there is no resolved-scale clouds in that particular grid box. It can also
!! specify a cloud fraction for resolved-scale clouds as is done currently when
!! using MYNN-EDMF. For clouds coming from the convection schemes (in this case
!! only used by GF scheme), two cloud fraction options are available:
!! Xu-Randall (XR1996) or Chaboureau and Bechtold (CB2005), chosen by the
!! switch "conv_cf_opt" = 0: CB2005, 1: XR1996.
!!
!> \section arg_table_sgscloud_radpre_run Argument Table
!! \htmlinclude sgscloud_radpre_run.html
!!
!!    cloud array description:                                          !
!!          clouds(:,:,1)  -  layer total cloud fraction                !
!!          clouds(:,:,2)  -  layer cloud liq water path                !
!!          clouds(:,:,3)  -  mean effective radius for liquid cloud    !
!!          clouds(:,:,4)  -  layer cloud ice water path                !
!!          clouds(:,:,5)  -  mean effective radius for ice cloud       !
!!          clouds(:,:,6)  -  layer rain drop water path                !
!!          clouds(:,:,7)  -  mean effective radius for rain drop       !
!!          clouds(:,:,8)  -  layer snow flake water path               !
!!          clouds(:,:,9)  -  mean effective radius for snow flake
!!
!>\section sgscloud_radpre_mod  SGS Cloud Scheme Pre General Algorithm
      subroutine sgscloud_radpre_run(    &
           im,dt,fhswr,levs,             &
           flag_init,flag_restart,       &
           con_g, con_pi, eps, epsm1,    &
           r_v, cpv, rcp,                &
           xlv, xlf, cp,                 &
           do_mynnedmf,                  &
           qc, qi, qv, T3D, P3D, exner,  &
           qr, qs, qg,                   &
           qci_conv,qlc,qli,ud_mf,       &
           imfdeepcnv, imfdeepcnv_gf,    &
           imfdeepcnv_unified,           &
           imfdeepcnv_sas,               &
           qc_save, qi_save, qs_save,    &
           qc_bl,qi_bl,cldfra_bl,        &
           delp,clouds1,clouds2,clouds3, &
           clouds4,clouds5,              &
           clouds8,clouds9,slmsk,        &
           nlay, plyr, xlat, dz,de_lgth, &
           cldsa,mtopa,mbota,            &
           imp_physics, imp_physics_gfdl,&
           imp_physics_fa,               &
           iovr,                         &
           errmsg, errflg                )

      use machine , only : kind_phys
      use module_radiation_clouds, only : gethml
      use radcons,  only: qmin          ! Minimum values for various calculations
      use funcphys, only: fpvs          ! Function to compute sat. vapor pressure over liq.
!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 
      ! Interface variables
      real(kind=kind_phys), intent(in) :: con_g, con_pi, eps, epsm1
      real(kind=kind_phys), intent(in) :: r_v, cpv, rcp
      real(kind=kind_phys), intent(in) :: xlv, xlf, cp
      real(kind=kind_phys), intent(in) :: dt,fhswr
      real :: xls, xlvcp, xlscp !derived below
      real(kind=kind_phys)             :: gfac
      integer,             intent(in)  :: im, levs, imfdeepcnv, imfdeepcnv_gf, &
           &  nlay, imfdeepcnv_sas, imfdeepcnv_unified, imp_physics, & 
           &  imp_physics_gfdl, imp_physics_fa
      logical,             intent(in)  :: flag_init, flag_restart, do_mynnedmf

      real(kind=kind_phys), dimension(:,:), intent(inout) :: qc, qi
      real(kind=kind_phys), dimension(:,:), intent(inout) :: qr, qs, qg
      ! note: qci_conv only allocated if GF is used
      real(kind=kind_phys), dimension(:,:), intent(inout) :: qci_conv
      real(kind=kind_phys), dimension(:,:), intent(inout) :: qlc, qli !for SAS
      real(kind=kind_phys), dimension(:,:), intent(in)    :: ud_mf
      real(kind=kind_phys), dimension(:,:), intent(in)    :: T3D,delp
      real(kind=kind_phys), dimension(:,:), intent(in)    :: qv,P3D,exner
      real(kind=kind_phys), dimension(:,:), intent(inout) ::  &
           &         clouds1,clouds2,clouds3,clouds4,clouds5, &
           &         clouds8,clouds9
      real(kind=kind_phys), dimension(:,:), intent(inout) :: qc_save, qi_save, qs_save
      real(kind=kind_phys), dimension(:,:), intent(in)    :: qc_bl, qi_bl, cldfra_bl
      real(kind=kind_phys), dimension(:),   intent(in)    :: slmsk, xlat, de_lgth
      real(kind=kind_phys), dimension(:,:), intent(in)    :: plyr, dz      
      real(kind=kind_phys), dimension(:,:), intent(inout) :: cldsa
      integer,              dimension(:,:), intent(inout) :: mbota, mtopa
      integer,                              intent(in)    :: iovr
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      ! pressure limits of cloud domain interfaces (low,mid,high) in mb (0.1kPa)
      real(kind=kind_phys) :: ptop1(im,3+1)  !< pressure limits of cloud domain interfaces
      real(kind=kind_phys) :: ptopc(3+1,2 )  !< pressure limits of cloud domain interfaces
                                             !! (low, mid, high) in mb (0.1kPa)
      data ptopc / 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /
      real(kind=kind_phys), dimension(im,nlay) :: cldcnv
      real(kind=kind_phys), dimension(im)      :: rxlat
      real(kind=kind_phys) :: Tc, Tk, liqfrac, iwc, ice_frac, snow_frac
      integer              :: i, k, id
      ! DH* 20200723 - see comment at the end of this routine around 'gethml'
      real(kind=kind_phys), dimension(im,nlay) :: alpha_dummy
      ! *DH

      ! PARAMETERS FOR RANDALL AND XU (1996) CLOUD FRACTION
      real, parameter  ::  coef_p = 0.25, coef_gamm = 0.49, coef_alph = 100.
      real :: rhgrid,h2oliq,qsat,tem1,tem2,clwt,es,onemrh,value

      !Chaboureau and Bechtold (2002 and 2005)
      real :: a, f, sigq, qmq, qt, xl, th, thl, rsl, cpm, cb_cf
      real(kind=kind_phys) :: tlk

      !Option to convective cloud fraction
      integer, parameter :: conv_cf_opt = 0  !0: C-B, 1: X-R

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! some derived variables from incoming constants:
      xls=xlv+xlf
      xlvcp=xlv/cp
      xlscp=(xlv+xlf)/cp

      !write(0,*)"=============================================="
      !write(0,*)"in SGSCLoud_RadPre"
      gfac=1.0e5/con_g
      if (flag_init .and. (.not. flag_restart)) then
        !write (0,*) 'Skip this flag_init = ', flag_init
        ! return
        ! Need default cloud fraction when MYNN is not used: Resort to
        ! Xu-Randall (1996).
        ! cloud fraction =
        ! {1-exp[-100.0*qc/((1-RH)*qsat)**0.49]}*RH**0.25
        do k = 1, levs
          do i = 1, im
            if ( qi(i,k) > 1E-7 .OR. qc(i,k) > 1E-7 ) then
              es     = min( p3d(i,k),  fpvs( t3d(i,k) ) )  ! fpvs and prsl in pa
              qsat   = max( QMIN, eps * es / (p3d(i,k) + epsm1*es) )
              rhgrid = max( 0., min( 1., qv(i,k)/qsat ) )
              h2oliq = qc(i,k) + qi(i,k) + qr(i,k) + qs(i,k) + qg(i,k)   ! g/kg
              clwt   = 1.0e-6 * (p3d(i,k)*0.00001)

              if (h2oliq > clwt) then
                onemrh= max( 1.e-10, 1.0-rhgrid )
                tem1  = min(max((onemrh*qsat)**0.49,0.0001),1.0)  !jhan
                tem1  = 100.0 / tem1
                value = max( min( tem1*(h2oliq-clwt), 50.0 ), 0.0 )
                tem2  = sqrt( sqrt(rhgrid) )

                clouds1(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
              endif

            endif
          enddo
        enddo

      else ! timestep > 1 or restart

        ! Back-up microphysics cloud information:
        do k = 1, levs
          do i = 1, im
            qc_save(i,k) = qc(i,k)
            qi_save(i,k) = qi(i,k)
            qs_save(i,k) = qs(i,k)
          end do
        end do

        if ( do_mynnedmf ) then

          ! add boundary layer clouds - Note: now the temperature-dependent sorting of
          ! ice and water subgrid-scale clouds is done inside the MYNN-EDMF

          do k = 1, levs
            do i = 1, im

              !if (imp_physics == imp_physics_gfdl) then
              !  ! only complement the GFDL cloud fractions
              !  if (clouds1(i,k) < 0.01 .and. cldfra_bl(i,k) > 0.01) then
              !    clouds1(i,k) = cldfra_bl(i,k)
              !  endif
              !else
                 clouds1(i,k) = cldfra_bl(i,k)
              !endif

              if (qc(i,k) < 1.e-6 .and. cldfra_bl(i,k)>0.001) then
                qc(i,k) = qc_bl(i,k)

                !eff radius cloud water (microns) from Miles et al. (2007)
                if (nint(slmsk(i)) == 1) then !land
                  if(qc(i,k)>1.E-8)clouds3(i,k)=5.4
                else
                  if(qc(i,k)>1.E-8)clouds3(i,k)=9.6
                endif

                !calculate the liquid water path using additional BL clouds
                clouds2(i,k) = max(0.0, qc(i,k) * gfac * delp(i,k))
              endif

              Tc = T3D(i,k) - 273.15
              !crudely split frozen species into 50% ice and 50% snow below
              !~700 mb and decrease snow to zero by ~300 mb 
              snow_frac = min(0.5, max((p3d(i,k)-30000.0),0.0)/140000.0)
              ice_frac  = 1.0 - snow_frac
              if (qi(i,k) < 1.e-9 .and. cldfra_bl(i,k)>0.001) then
                qi(i,k) = ice_frac*qi_bl(i,k)

                !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 6b)
                if(qi(i,k)>1.E-8)clouds5(i,k)=max(173.45 + 2.14*Tc, 20.)
                !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 8b)
                !iwc = qi(i,k)*1.0e6*rho(i,k)
                !IF(qi(i,k)>1.E-8)clouds5(i,k)=MAX(139.7 + 1.76*Tc + 13.49*LOG(iwc), 20.)

                !calculate the ice water path using additional BL clouds
                clouds4(i,k) = max(0.0, qi(i,k) * gfac * delp(i,k))
              endif

              if (qs(i,k) < 1.e-9 .and. cldfra_bl(i,k)>0.001) then
                qs(i,k) = snow_frac*qi_bl(i,k)

                !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 6b)
                if(qs(i,k)>1.E-8)clouds9(i,k)=max(2.*(173.45 + 2.14*Tc), 50.)

                !calculate the snow water path using additional BL clouds
                clouds8(i,k) = max(0.0, qs(i,k) * gfac * delp(i,k))
              endif
            enddo
          enddo

        elseif (imp_physics /= imp_physics_gfdl) then 

          ! Non-MYNN cloud fraction AND non-GFDL microphysics, since both
          ! have their own cloud fractions. In this case, we resort to
          ! Xu-Randall (1996).
          ! cloud fraction =
          ! {1-exp[-100.0*qc/((1-RH)*qsat)**0.49]}*RH**0.25
          do k = 1, levs
            do i = 1, im
              if ( qi(i,k) > 1E-7 .OR. qc(i,k) > 1E-7 ) then

                es     = min( p3d(i,k),  fpvs( t3d(i,k) ) )  ! fpvs and prsl in pa
                qsat   = max( QMIN, eps * es / (p3d(i,k) + epsm1*es) )
                rhgrid = max( 0., min( 1., qv(i,k)/qsat ) )
                h2oliq = qc(i,k) + qi(i,k) + qr(i,k) + qs(i,k) + qg(i,k) ! g/kg
                clwt   = 1.0e-6 * (p3d(i,k)*0.00001)

                if (h2oliq > clwt) then
                  onemrh= max( 1.e-10, 1.0-rhgrid )
                  tem1  = min(max((onemrh*qsat)**0.49,0.0001),1.0)  !jhan
                  tem1  = 100.0 / tem1
                  value = max( min( tem1*(h2oliq-clwt), 50.0 ), 0.0 )
                  tem2  = sqrt( sqrt(rhgrid) )

                  clouds1(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
                endif

              endif
            enddo
          enddo

        endif ! end MYNN or OTHER choice for background clouds fractions

        ! At this point, we have cloud properties for all non-deep convective clouds.
        ! So now we add the convective clouds:

        if (imfdeepcnv == imfdeepcnv_gf .or. imfdeepcnv == imfdeepcnv_unified) then
          do k = 1, levs
            do i = 1, im
              if ( qci_conv(i,k) > 0. ) then
                Tk = T3D(i,k)
                Tc = Tk - 273.15

                !Partition the convective clouds into water & frozen species
                liqfrac = min(1., max(0., (Tk-244.)/29.))
                qc(i,k) = qc(i,k)+qci_conv(i,k)*liqfrac
                !split ice & snow 50-50%
                qi(i,k) = qi(i,k)+0.5*qci_conv(i,k)*(1. - liqfrac)
                qs(i,k) = qs(i,k)+0.5*qci_conv(i,k)*(1. - liqfrac)

                !eff radius cloud water (microns)
                if (nint(slmsk(i)) == 1) then !land
                  if(qc(i,k)>1.E-8)clouds3(i,k)=5.4
                else
                  !from Miles et al.
                  if(qc(i,k)>1.E-8)clouds3(i,k)=9.6
                endif
                !from Mishra et al. (2014, JGR Atmos), assume R_sno = 2*R_ice
                if(qi(i,k)>1.e-8)clouds5(i,k)=max(     173.45 + 2.14*Tc , 20.)
                if(qs(i,k)>1.e-8)clouds9(i,k)=max(2.0*(173.45 + 2.14*Tc), 50.)

                if ( conv_cf_opt .eq. 0 ) then
                   !print *,'Chab-Bechtold cloud fraction used'
                   !  clouds1(i,k) = cldfra_bl(i,k)

                   !Alternatively, use Chaboureau-Bechtold (CB) convective component
                   !Based on both CB2002 and CB2005.
                   xl  = xlv*liqfrac + xls*(1.-liqfrac)  ! blended heat capacity
                   tlk = t3d(i,k) - xlvcp/exner(i,k)*qc(i,k) &
                       &          - xlscp/exner(i,k)*qi(i,k)! liquid temp
                   ! get saturation water vapor mixing ratio at tl and p
                   es  = min( p3d(i,k), fpvs( tlk ) )   ! fpvs and prsl in pa
                   qsat= max( QMIN, eps*es / (p3d(i,k) + epsm1*es) )
                   rsl = xl*qsat / (r_v*tlk**2)   ! slope of C-C curve at t = tl 
                                                  ! CB02, Eqn. 4
                   qt  = qc(i,k) + qi(i,k) + qv(i,k) !total water
                   cpm = cp + qt*cpv              ! CB02, sec. 2, para. 1
                   a   = 1./(1. + xl*rsl/cpm)     ! CB02 variable "a"
                   !Now calculate convective component of the cloud fraction:
                   if (a > 0.0) then
                      f = min(1.0/a, 4.0)         ! f is the vertical profile
                   else                           ! scaling function (CB2005)
                      f = 1.0
                   endif
                   sigq = 1.5E-3 * ud_mf(i,k)/dt * f
                   !sigq = 3.E-3 * ud_mf(i,k)/dt * f
                   sigq = SQRT(sigq**2 + 1e-10)   ! combined conv + background components
                   qmq  = a * (qt - qsat)         ! saturation deficit/excess;
                                                  !   the numerator of Q1
                   cb_cf= min(max(0.5 + 0.36 * atan(1.55*(qmq/sigq)),0.0),0.99)
                   if (qci_conv(i,k) .lt. 1e-9) cb_cf = 0.0
                   if (do_mynnedmf .and. qmq .ge. 0.0) then
                      ! leverage C-B stratus clouds from MYNN in saturated conditions
                      if (cb_cf .gt. 0.0) then
                         clouds1(i,k) = 0.5*(clouds1(i,k) + cb_cf)
                      else
                         !default to MYNN clouds - already specified
                      endif
                   else                           ! unsaturated
                      clouds1(i,k) = cb_cf
                   endif
                else
                   !print *,'GF with Xu-Randall cloud fraction'
                   ! Xu-Randall (1996) cloud fraction
                   es     = min( p3d(i,k),  fpvs( t3d(i,k) ) )  ! fpvs and prsl in pa
                   qsat   = max( QMIN, eps*es / (p3d(i,k) + epsm1*es) )
                   rhgrid = max( 0., min( 1.00, qv(i,k)/qsat ) )
                   h2oliq = qc(i,k) + qi(i,k) + qr(i,k) + qs(i,k) + qg(i,k) ! g/kg
                   clwt   = 1.0e-6 * (p3d(i,k)*0.00001)

                   if (h2oliq > clwt) then
                      onemrh= max( 1.e-10, 1.0-rhgrid )
                      tem1  = min(max((onemrh*qsat)**0.49,0.0001),1.0)  !jhan
                      tem1  = 100.0 / tem1
                      value = max( min( tem1*(h2oliq-clwt), 50.0 ), 0.0 )
                      tem2  = sqrt( sqrt(rhgrid) )

                      clouds1(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
                   else
                      clouds1(i,k) = 0.0
                   endif
                   !print*,"XuRandla- cf:",clouds1(i,k)," rh:",rhgrid," qt:",h2oliq
                   !print*,"XuRandlb- clwt:",clwt," qsat:",qsat," p:",p3d(i,k)
                endif ! end convective cf choice
              endif ! qci_conv
            enddo
          enddo

        elseif (imfdeepcnv == imfdeepcnv_sas) then

          do k = 1, levs
            do i = 1, im
              h2oliq  = qlc(i,k)+qli(i,k)
              if ( h2oliq > 0. ) then
                Tk = T3D(i,k)
                Tc = Tk - 273.15

                !Partition the convective clouds into water & frozen species
                liqfrac = min(1., max(0., (Tk-244.)/29.))

                qc(i,k) = qc(i,k)+qlc(i,k)
                !split ice & snow 50-50%
                qi(i,k) = qi(i,k)+0.5*qli(i,k)
                qs(i,k) = qs(i,k)+0.5*qli(i,k)

                !eff radius cloud water (microns)
                if (nint(slmsk(i)) == 1) then !land
                  if(qc(i,k)>1.E-8)clouds3(i,k)=5.4
                else
                  !from Miles et al.
                  if(qc(i,k)>1.E-8)clouds3(i,k)=9.6
                endif
                !from Mishra et al. (2014, JGR Atmos), assume R_sno = 2*R_ice
                if(qi(i,k)>1.e-8)clouds5(i,k)=max(     173.45 + 2.14*Tc , 20.)
                if(qs(i,k)>1.e-8)clouds9(i,k)=max(2.0*(173.45 + 2.14*Tc), 50.)

                if ( conv_cf_opt .eq. 0 ) then
                   !print *,'Chab-Bechtold cloud fraction used'
                   !Alternatively, use Chaboureau-Bechtold (CB) convective component
                   !Based on both CB2002 and CB2005.
                   xl  = xlv*liqfrac + xls*(1.-liqfrac)  ! blended heat capacity
                   tlk = t3d(i,k) - xlvcp/exner(i,k)*qc(i,k) &
                       &          - xlscp/exner(i,k)*qi(i,k)! liquid temp
                   ! get saturation water vapor mixing ratio at tl and p
                   es  = min( p3d(i,k), fpvs( tlk ) )   ! fpvs and prsl in pa
                   qsat= max( QMIN, eps*es / (p3d(i,k) + epsm1*es) )
                   rsl = xl*qsat / (r_v*tlk**2)   ! slope of C-C curve at t = tl
                                                  ! CB02, Eqn. 4
                   qt  = qc(i,k) + qi(i,k) + qv(i,k) !total water
                   cpm = cp + qt*cpv              ! CB02, sec. 2, para. 1
                   a   = 1./(1. + xl*rsl/cpm)     ! CB02 variable "a"
                   !Now calculate convective component of the cloud fraction:
                   if (a > 0.0) then
                      f = min(1.0/a, 4.0)         ! f is the vertical profile
                   else                           ! scaling function (CB2005)
                      f = 1.0
                   endif
                   sigq = 1.5E-3 * ud_mf(i,k)/dt * f
                   !sigq = 3.E-3 * ud_mf(i,k)/dt * f
                   sigq = SQRT(sigq**2 + 1e-10)   ! combined conv + background components
                   qmq  = a * (qt - qsat)         ! saturation deficit/excess;
                                                  !   the numerator of Q1
                   cb_cf= min(max(0.5 + 0.36 * atan(1.55*(qmq/sigq)),0.0),0.99)
                   if (h2oliq .lt. 1e-9) cb_cf = 0.0
                   if (do_mynnedmf .and. qmq .ge. 0.0) then
                      ! leverage C-B stratus clouds from MYNN in saturated conditions
                      if (cb_cf .gt. 0.0) then
                         clouds1(i,k) = 0.5*(clouds1(i,k) + cb_cf)
                      else
                         !default to MYNN clouds - already specified
                      endif
                   else                           ! unsaturated
                      clouds1(i,k) = cb_cf
                   endif
                else
                   !print *,'SAS with Xu-Randall cloud fraction'
                   ! Xu-Randall (1996) cloud fraction
                   es     = min( p3d(i,k),  fpvs( t3d(i,k) ) )  ! fpvs and prsl in pa
                   qsat   = max( QMIN, eps*es / (p3d(i,k) + epsm1*es) )
                   rhgrid = max( 0., min( 1.00, qv(i,k)/qsat ) )
                   h2oliq = qc(i,k) + qi(i,k) + qr(i,k) + qs(i,k) + qg(i,k) ! g/kg
                   clwt   = 1.0e-6 * (p3d(i,k)*0.00001)

                   if (h2oliq > clwt) then
                      onemrh= max( 1.e-10, 1.0-rhgrid )
                      tem1  = min(max((onemrh*qsat)**0.49,0.0001),1.0)  !jhan
                      tem1  = 100.0 / tem1
                      value = max( min( tem1*(h2oliq-clwt), 50.0 ), 0.0 )
                      tem2  = sqrt( sqrt(rhgrid) )

                      clouds1(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
                   else
                      clouds1(i,k) = 0.0
                   endif
                   !print*,"XuRandla- cf:",clouds1(i,k)," rh:",rhgrid," qt:",h2oliq
                   !print*,"XuRandlb- clwt:",clwt," qsat:",qsat," p:",p3d(i,k)
                endif ! end convective cf choice
              endif ! qlc/qli check
            enddo
          enddo

        endif ! convection scheme check

      endif ! timestep > 1

      end subroutine sgscloud_radpre_run
      end module sgscloud_radpre
