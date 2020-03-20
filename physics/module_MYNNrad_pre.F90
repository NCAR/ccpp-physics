!> \file module_MYNNrad_pre.F90
!!  Contains the preliminary (interstitial) work to the call to the radiation schemes:
!!    1) Backs up the original qc & qi
!!    2) Adds the subgrid clouds mixing ratio and cloud fraction to the original qc, qi and cloud fraction coming from the microphysics scheme.

      MODULE mynnrad_pre

      contains

      subroutine mynnrad_pre_init ()
      end subroutine mynnrad_pre_init

      subroutine mynnrad_pre_finalize ()
      end subroutine mynnrad_pre_finalize

!>\defgroup gsd_mynnrad_pre GSD mynnrad_pre_run Module
!>\ingroup gsd_mynn_edmf
!! This interstitial code adds the subgrid clouds to the resolved-scale clouds if there is no resolved-scale clouds in that particular grid box.
#if 0
!> \section arg_table_mynnrad_pre_run Argument Table
!! \htmlinclude mynnrad_pre_run.html
!!
#endif
!
!    cloud array description:                                          !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud liq water path                !
!          clouds(:,:,3)  -  mean effective radius for liquid cloud    !
!          clouds(:,:,4)  -  layer cloud ice water path                !
!          clouds(:,:,5)  -  mean effective radius for ice cloud       !
!
!###===================================================================
SUBROUTINE mynnrad_pre_run(                &
     &     ix,im,levs,                     &
     &     flag_init,flag_restart,         &
     &     qc, qi, T3D,                    &
     &     qc_save, qi_save,               &
     &     qc_bl,cldfra_bl,                &
     &     delp,clouds1,clouds2,clouds3,   &
     &     clouds4,clouds5,slmsk,          &
     &     errmsg, errflg                  )

! should be moved to inside the mynn:
      use machine , only : kind_phys
      ! DH* TODO - input argument, not constant
      use physcons,            only : con_g

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 
      ! Interface variables
      real (kind=kind_phys), parameter :: gfac=1.0e5/con_g
      integer, intent(in)  :: ix, im, levs
      logical,          intent(in)  :: flag_init, flag_restart
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: qc, qi
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: T3D,delp
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: &
           &         clouds1,clouds2,clouds3,clouds4,clouds5
      real(kind=kind_phys), dimension(im,levs), intent(out)   :: qc_save, qi_save
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: qc_bl, cldfra_bl
      ! DH* TODO add intent() information for delp,clouds1,clouds2,clouds3,clouds4,clouds5
      real(kind=kind_phys), dimension(im), intent(in)         :: slmsk
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer              :: i, k
      real                 :: Tc, iwc

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !write(0,*)"=============================================="
      !write(0,*)"in mynn rad pre"

      if (flag_init .and. (.not. flag_restart)) then
       !write (0,*) 'Skip MYNNrad_pre flag_init = ', flag_init
        return
      endif
     ! Add subgrid cloud information:
        do k = 1, levs
           do i = 1, im

              qc_save(i,k) = qc(i,k)
              qi_save(i,k) = qi(i,k)
              clouds1(i,k) = CLDFRA_BL(i,k)

              IF (qc(i,k) < 1.E-6 .AND. qi(i,k) < 1.E-8 .AND. CLDFRA_BL(i,k)>0.001) THEN
                !Partition the BL clouds into water & ice according to a linear
                !approximation of Hobbs et al. (1974). This allows us to only use
                !one 3D array for both cloud water & ice.
!               Wice = 1. - MIN(1., MAX(0., (t(i,k)-254.)/15.))
!               Wh2o = 1. - Wice
                !clouds1(i,k)=MAX(clouds1(i,k),CLDFRA_BL(i,k))
                !clouds1(i,k)=MAX(0.0,MIN(1.0,clouds1(i,k)))
                qc(i,k) = QC_BL(i,k)*(MIN(1., MAX(0., (T3D(i,k)-254.)/15.)))*CLDFRA_BL(i,k)
                qi(i,k) = QC_BL(i,k)*(1. - MIN(1., MAX(0., (T3D(i,k)-254.)/15.)))*CLDFRA_BL(i,k)

                Tc = T3D(i,k) - 273.15
                !iwc = qi(i,k)*1.0e6*rho(i,k)

                IF (nint(slmsk(i)) == 1) then !land
                  IF(qc(i,k)>1.E-8)clouds3(i,k)=5.4                !eff radius cloud water (microns)
                  !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos)
                  IF(qi(i,k)>1.E-8)clouds5(i,k)=MAX(173.45 + 2.14*Tc, 20.)
                ELSE
                  !eff radius cloud water (microns), from Miles et al. 
                  IF(qc(i,k)>1.E-8)clouds3(i,k)=9.6
                  !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 6b)
                  IF(qi(i,k)>1.E-8)clouds5(i,k)=MAX(173.45 + 2.14*Tc, 20.)
                  !eff radius cloud ice (microns), from Mishra et al. (2014, JGR Atmos, fig 8b)
                  !IF(qi(i,k)>1.E-8)clouds5(i,k)=MAX(139.7 + 1.76*Tc + 13.49*LOG(iwc), 20.)
                ENDIF

                !water and ice paths
                clouds2(i,k) = max(0.0, qc(i,k) * gfac * delp(i,k))
                clouds4(i,k) = max(0.0, qi(i,k) * gfac * delp(i,k))

              ENDIF

           enddo
        enddo

        !print*,"===Finished adding subgrid clouds to the resolved-scale clouds"
        !print*,"qc_save:",qc_save(1,1)," qi_save:",qi_save(1,1)

  END SUBROUTINE mynnrad_pre_run

!###=================================================================

END MODULE mynnrad_pre
