!>\file module_SGSCloud_RadPre.F90
!!  Contains the preliminary (interstitial) work to the call to the radiation schemes:
!!    1) Backs up the original qc & qi
!!    2) Adds the partioning of convective condensate into liqice/ice for effective radii
!!    3) Adds the subgrid clouds mixing ratio and cloud fraction to the original qc, qi and cloud fraction coming from the microphysics scheme.
!!    4) Recompute the diagnostic high, mid, low, total and bl clouds to be consistent with radiation

      module sgscloud_radpre

      contains

      subroutine sgscloud_radpre_init ()
      end subroutine sgscloud_radpre_init

      subroutine sgscloud_radpre_finalize ()
      end subroutine sgscloud_radpre_finalize

!> \defgroup sgsrad_group GSD sgscloud_radpre_run Module
!> \ingroup sgscloud_radpre
!! This interstitial code adds the subgrid clouds to the resolved-scale clouds if there is no resolved-scale clouds in that particular grid box.
!> \section arg_table_sgscloud_radpre_run Argument Table
!! \htmlinclude sgscloud_radpre_run.html
!!
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
           ix,im,levs,                   &
           flag_init,flag_restart,       &
           do_mynnedmf,                  &
           qc, qi, T3D,                  &
           qr, qs,                       &
           qci_conv,                     &
           imfdeepcnv, imfdeepcnv_gf,    &
           qc_save, qi_save,             &
           qc_bl,cldfra_bl,              &
           delp,clouds1,clouds2,clouds3, &
           clouds4,clouds5,slmsk,        &
           nlay, plyr, xlat, dz,de_lgth, &
           cldsa,mtopa,mbota,            &
           errmsg, errflg                )

! should be moved to inside the mynn:
      use machine , only : kind_phys
      use physcons,                only : con_g, con_pi
      use module_radiation_clouds, only : gethml

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 
      ! Interface variables
      real (kind=kind_phys), parameter :: gfac=1.0e5/con_g
      integer, intent(in)  :: ix, im, levs, imfdeepcnv, imfdeepcnv_gf, nlay
      logical,          intent(in)  :: flag_init, flag_restart, do_mynnedmf
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: qc, qi
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: qr, qs
      ! qci_conv only allocated if GF is used
      real(kind=kind_phys), dimension(:,:),     intent(inout) :: qci_conv
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: T3D,delp
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: &
           &         clouds1,clouds2,clouds3,clouds4,clouds5
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: qc_save, qi_save
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: qc_bl, cldfra_bl
      real(kind=kind_phys), dimension(im),      intent(in)    :: slmsk, xlat, de_lgth
      real(kind=kind_phys), dimension(im,nlay), intent(in)    :: plyr, dz      
      real(kind=kind_phys), dimension(im,5),    intent(inout) :: cldsa
      integer,              dimension(im,3),    intent(inout) :: mbota, mtopa
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
      real (kind=kind_phys):: Tc, iwc, tem1
      integer              :: i, k, id

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !write(0,*)"=============================================="
      !write(0,*)"in mynn rad pre"

      if (flag_init .and. (.not. flag_restart)) then
       !write (0,*) 'Skip MYNNrad_pre flag_init = ', flag_init
        return
      endif
      ! Back-up microphysics cloud information:
      do k = 1, levs
        do i = 1, im
          qc_save(i,k) = qc(i,k)
          qi_save(i,k) = qi(i,k)
        end do
      end do
 
      ! add boundary layer clouds
      if (do_mynnedmf) then
        do k = 1, levs
          do i = 1, im
            clouds1(i,k) = cldfra_bl(i,k)

            !if( qr(i,k) > 1.0e-7 .OR. qs(i,k) > 1.0e-7.or.qci_conv(i,k)>1.0e-7)THEN
                !Keep Xu-RandalL clouds fraction - do not overwrite
            !else
            !    clouds1(i,k) = cldfra_bl(i,k)
            !endif

            if (qc(i,k) < 1.e-6 .and. qi(i,k) < 1.e-8 .and. cldfra_bl(i,k)>0.001) then
              !Partition the BL clouds into water & ice according to a linear
              !approximation of Hobbs et al. (1974). This allows us to only use
              !one 3D array for both cloud water & ice.
              !Wice = 1. - MIN(1., MAX(0., (t(i,k)-254.)/15.))
              !Wh2o = 1. - Wice
              !clouds1(i,k)=MAX(clouds1(i,k),CLDFRA_BL(i,k))
              !clouds1(i,k)=MAX(0.0,MIN(1.0,clouds1(i,k)))
              qc(i,k) = qc_bl(i,k)*(min(1., max(0., (T3D(i,k)-244.)/25.)))*cldfra_bl(i,k)
              qi(i,k) = qc_bl(i,k)*(1. - min(1., max(0., (T3D(i,k)-244.)/25.)))*cldfra_bl(i,k)
              Tc = T3D(i,k) - 273.15
              !iwc = qi(i,k)*1.0e6*rho(i,k)

              if (nint(slmsk(i)) == 1) then !land
                if(qc(i,k)>1.e-8)clouds3(i,k)=5.4                !eff radius cloud water (microns)
                !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos)
                if(qi(i,k)>1.E-8)clouds5(i,k)=max(173.45 + 2.14*Tc, 20.)
              else
                !eff radius cloud water (microns), from Miles et al. 
                if(qc(i,k)>1.E-8)clouds3(i,k)=9.6
                !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 6b)
                if(qi(i,k)>1.E-8)clouds5(i,k)=max(173.45 + 2.14*Tc, 20.)
                !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 8b)
                !IF(qi(i,k)>1.E-8)clouds5(i,k)=MAX(139.7 + 1.76*Tc + 13.49*LOG(iwc), 20.)
              endif
              !calculate water and ice paths for additional BL clouds
              clouds2(i,k) = max(0.0, qc(i,k) * gfac * delp(i,k))
              clouds4(i,k) = max(0.0, qi(i,k) * gfac * delp(i,k))
            endif

          enddo
        enddo
      endif ! do_mynnedmf

      ! add convective clouds
      if (imfdeepcnv == imfdeepcnv_gf) then
        do k = 1, levs
          do i = 1, im
            if ( qci_conv(i,k) > 0.) then
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
            endif
          enddo
        enddo
      endif
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

!> - Recompute the diagnostic high, mid, low, total and bl cloud fraction
      call gethml                                                       &
!  ---  inputs:
           ( plyr, ptop1, clouds1, cldcnv, dz, de_lgth, im, nlay,       &
!  ---  outputs:
             cldsa, mtopa, mbota)

       !print*,"===Finished adding subgrid clouds to the resolved-scale clouds"
       !print*,"qc_save:",qc_save(1,1)," qi_save:",qi_save(1,1)

      end subroutine sgscloud_radpre_run

      end module sgscloud_radpre
