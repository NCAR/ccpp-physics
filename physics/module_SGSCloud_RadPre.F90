!>\file module_SGSCloud_RadPre.F90
!!  Contains the preliminary (interstitial) work to the call to the radiation schemes:
!!    1) Backs up the original qc & qi
!!    2) Adds the partioning of convective condensate into liqice/ice for effective radii
!!    3) Adds the subgrid clouds mixing ratio and cloud fraction to the original (resolved-
!!       scale) qc, qi and cloud fraction coming from the microphysics scheme.
!!    4) Recompute the diagnostic high, mid, low, total and bl clouds to be consistent with radiation

      module sgscloud_radpre

      contains

      subroutine sgscloud_radpre_init ()
      end subroutine sgscloud_radpre_init

      subroutine sgscloud_radpre_finalize ()
      end subroutine sgscloud_radpre_finalize

!> \defgroup sgsrad_group GSD sgscloud_radpre_run Module
!> \ingroup sgscloud_radpre
!! This interstitial code adds the subgrid clouds to the resolved-scale clouds 
!! if there is no resolved-scale clouds in that particular grid box. It can also 
!! specify a cloud fraction for resolved-scale clouds, using Xu-Randall (1996),
!! if desired.
!> \section arg_table_sgscloud_radpre_run Argument Table
!! \htmlinclude sgscloud_radpre_run.html
!!
!!    cloud array description:                                          !
!!          clouds(:,:,1)  -  layer total cloud fraction                !
!!          clouds(:,:,2)  -  layer cloud liq water path                !
!!          clouds(:,:,3)  -  mean effective radius for liquid cloud    !
!!          clouds(:,:,4)  -  layer cloud ice water path                !
!!          clouds(:,:,5)  -  mean effective radius for ice cloud       !
!!
!>\section sgscloud_radpre GSD SGS Scheme General Algorithm
!> @{
      subroutine sgscloud_radpre_run(    &
           im,levs,                      &
           flag_init,flag_restart,       &
           do_mynnedmf,                  &
           qc, qi, qv, T3D, P3D,         &
           qr, qs, qg,                   &
           qci_conv,                     &
           imfdeepcnv, imfdeepcnv_gf,    &
           qc_save, qi_save,             &
           qc_bl,qi_bl,cldfra_bl,        &
           delp,clouds1,clouds2,clouds3, &
           clouds4,clouds5,slmsk,        &
           nlay, plyr, xlat, dz,de_lgth, &
           cldsa,mtopa,mbota,            &
           imp_physics, imp_physics_gfdl,&
           iovr,                         &
           errmsg, errflg                )

! should be moved to inside the mynn:
      use machine , only : kind_phys
      use physcons, only : con_g, con_pi, &
                        eps   => con_eps, & ! Rd/Rv
                      epsm1 => con_epsm1    ! Rd/Rv-1
      use module_radiation_clouds, only : gethml
      use radcons, only: qmin               ! Minimum vlaues for varius calculations
      use funcphys, only: fpvs              ! Function ot compute sat. vapor pressure over liq.
!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 
      ! Interface variables
      real (kind=kind_phys), parameter :: gfac=1.0e5/con_g
      integer,          intent(in)  :: im, levs, imfdeepcnv, imfdeepcnv_gf, &
           &               nlay, imp_physics, imp_physics_gfdl
      logical,          intent(in)  :: flag_init, flag_restart, do_mynnedmf
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: qc, qi
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: qr, qs, qg
      ! qci_conv only allocated if GF is used
      real(kind=kind_phys), dimension(:,:),     intent(inout) :: qci_conv
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: T3D,delp, &
           &                                                     qv,P3D
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: &
           &         clouds1,clouds2,clouds3,clouds4,clouds5
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: qc_save, qi_save
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: qc_bl, qi_bl, cldfra_bl
      real(kind=kind_phys), dimension(im),      intent(in)    :: slmsk, xlat, de_lgth
      real(kind=kind_phys), dimension(im,nlay), intent(in)    :: plyr, dz      
      real(kind=kind_phys), dimension(im,5),    intent(inout) :: cldsa
      integer,              dimension(im,3),    intent(inout) :: mbota, mtopa
      integer,                                  intent(in)    :: iovr
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      ! pressure limits of cloud domain interfaces (low,mid,high) in mb (0.1kPa)
      real (kind=kind_phys) :: ptop1(im,3+1)  !< pressure limits of cloud domain interfaces
      real (kind=kind_phys) :: ptopc(3+1,2 )  !< pressure limits of cloud domain interfaces
                                              !! (low, mid, high) in mb (0.1kPa)
      data ptopc / 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /
      real(kind=kind_phys), dimension(im,nlay) :: cldcnv
      real(kind=kind_phys), dimension(im)      :: rxlat
      real (kind=kind_phys):: Tc, iwc
      integer              :: i, k, id
      ! DH* 20200723 - see comment at the end of this routine around 'gethml'
      real(kind=kind_phys), dimension(im,nlay) :: alpha_dummy
      ! *DH

      ! PARAMETERS FOR RANDALL AND XU (1996) CLOUD FRACTION
      REAL, PARAMETER  ::  coef_p = 0.25, coef_gamm = 0.49, coef_alph = 100.
      REAL :: rhgrid,h2oliq,qsat,tem1,tem2,clwt,es,onemrh,value

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !write(0,*)"=============================================="
      !write(0,*)"in SGSCLoud_RadPre"

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

              !if( qr(i,k) > 1.0e-7 .OR. qs(i,k) > 1.0e-7.or.qci_conv(i,k)>1.0e-7)THEN
                 !Keep Xu-RandalL clouds fraction - do not overwrite
              !else
              !    clouds1(i,k) = cldfra_bl(i,k)
              !endif

              if (qc(i,k) < 1.e-6 .and. cldfra_bl(i,k)>0.001) then
                qc(i,k) = qc_bl(i,k)*cldfra_bl(i,k)
                if (nint(slmsk(i)) == 1) then !land
                  if(qc(i,k)>1.E-8)clouds3(i,k)=5.4                !eff radius cloud water (microns)
                else
                  !eff radius cloud water (microns), from Miles et al.
                  if(qc(i,k)>1.E-8)clouds3(i,k)=9.6
                endif
                !calculate the liquid water path using additional BL clouds
                clouds2(i,k) = max(0.0, qc(i,k) * gfac * delp(i,k))
              endif
              if (qi(i,k) < 1.e-8 .and. cldfra_bl(i,k)>0.001) then
                qi(i,k) = qi_bl(i,k)*cldfra_bl(i,k)
                Tc = T3D(i,k) - 273.15
                !iwc = qi(i,k)*1.0e6*rho(i,k)
                if (nint(slmsk(i)) == 1) then !land
                  !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 6b)
                  if(qi(i,k)>1.E-8)clouds5(i,k)=max(173.45 + 2.14*Tc, 20.)
                else
                  if(qi(i,k)>1.E-8)clouds5(i,k)=max(173.45 + 2.14*Tc, 20.)
                  !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 8b)
                  !IF(qi(i,k)>1.E-8)clouds5(i,k)=MAX(139.7 + 1.76*Tc + 13.49*LOG(iwc), 20.)
                endif
                !calculate the ice water path using additional BL clouds
                clouds4(i,k) = max(0.0, qi(i,k) * gfac * delp(i,k))
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
        ! So now we add the convective clouds,

        if (imfdeepcnv == imfdeepcnv_gf) then
          do k = 1, levs
            do i = 1, im
              !if ( qci_conv(i,k) > 0. .AND. (qi(i,k) < 1E-7 .AND. qc(i,k) < 1E-7 ) ) then
              if ( qci_conv(i,k) > 0. ) then
                !Partition the convective clouds into water & ice according to a linear
                qc(i,k) = qc(i,k)+qci_conv(i,k)*(min(1., max(0., (T3D(i,k)-244.)/25.)))
                qi(i,k) = qi(i,k)+qci_conv(i,k)*(1. - min(1., max(0., (T3D(i,k)-244.)/25.)))

                Tc = T3D(i,k) - 273.15

                if (nint(slmsk(i)) == 1) then !land
                  if(qc(i,k)>1.E-8)clouds3(i,k)=5.4                !eff radius cloud water (microns)
                  !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos)
                  if(qi(i,k)>1.e-8)clouds5(i,k)=max(173.45 + 2.14*Tc, 20.)
                else
                  !eff radius cloud water (microns), from Miles et al. 
                  if(qc(i,k)>1.E-8)clouds3(i,k)=9.6
                  !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 6b)
                  if(qi(i,k)>1.E-8)clouds5(i,k)=max(173.45 + 2.14*Tc, 20.)
                endif

                if ( do_mynnedmf .or. (imp_physics == imp_physics_gfdl) ) then
                  !print *,'MYNN PBL or GFDL MP cldcov used'
                else
                  !print *,'GF with Xu-Randall cloud fraction'
                  ! Xu-Randall (1996) cloud fraction
                  es     = min( p3d(i,k),  fpvs( t3d(i,k) ) )  ! fpvs and prsl in pa
                  qsat   = max( QMIN, eps * es / (p3d(i,k) + epsm1*es) )
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
                endif ! not MYNN PBL or GFDL MP
              endif ! qci_conv
            enddo
          enddo
        endif ! imfdeepcnv_gf

      endif ! timestep > 1

!> - Compute SFC/low/middle/high cloud top pressure for each cloud domain for given latitude.

      do i =1, im
        rxlat(i) = abs( xlat(i) / con_pi )      ! if xlat in pi/2 -> -pi/2 range
!       rxlat(i) = abs(0.5 - xlat(i)/con_pi)    ! if xlat in 0 -> pi range
      enddo

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)
        do i =1, im
          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*rxlat(i)-1.0 )
        enddo
      enddo

      cldcnv = 0.

! DH* 20200723
! iovr == 4 or 5 requires alpha, which is computed in GFS_rrmtg_pre,
! which comes after SGSCloud_RadPre. Computing alpha here requires
! a lot more input variables and computations (dzlay etc.), and
! recomputing it in GFS_rrmtg_pre is a waste of time. Workaround:
! pass a dummy array initialized to zero to gethml for other values of iovr.
      if ( iovr == 4 .or. iovr == 5 ) then
        errmsg = 'Logic error in sgscloud_radpre: iovr==4 or 5 not implemented'
        errflg = 1
        return
      end if
!! Call subroutine get_alpha_exp to define alpha parameter for EXP and ER cloud overlap options
!      if ( iovr == 4 .or. iovr == 5 ) then 
!        call get_alpha_exp                                              &
!!  ---  inputs:
!             (im, nlay, dzlay, iovr, latdeg, julian, yearlen, clouds1,  &
!!  ---  outputs:
!              alpha                                                     &
!            )
!      endif
      alpha_dummy = 0.0
! *DH 2020723

!> - Recompute the diagnostic high, mid, low, total and bl cloud fraction
      call gethml                                                       &
!  ---  inputs:
           ( plyr, ptop1, clouds1, cldcnv, dz, de_lgth, alpha_dummy,    &
!  ---  outputs:
             im, nlay, cldsa, mtopa, mbota)

       !print*,"===Finished adding subgrid clouds to the resolved-scale clouds"
       !print*,"qc_save:",qc_save(1,1)," qi_save:",qi_save(1,1)

      end subroutine sgscloud_radpre_run

      end module sgscloud_radpre
