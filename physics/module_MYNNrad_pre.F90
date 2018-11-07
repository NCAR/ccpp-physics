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

!!
!> \brief This interstitial code adds the subgrid clouds to the resolved-scale clouds if there is no resolved-scale clouds in that particular grid box.
#if 0
!! \section arg_table_mynnrad_pre_run Argument Table
!! | local_name          | standard_name                                                               | long_name                                                                  | units   | rank | type      |    kind   | intent | optional |
!! |---------------------|-----------------------------------------------------------------------------|----------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | ix                  | horizontal_dimension                                                        | horizontal dimension                                                       | count   |    0 | integer   |           | in     | F        |
!! | im                  | horizontal_loop_extent                                                      | horizontal loop extent                                                     | count   |    0 | integer   |           | in     | F        |
!! | levs                | vertical_dimension                                                          | vertical layer dimension                                                   | count   |    0 | integer   |           | in     | F        |
!! | qc                  | cloud_condensed_water_mixing_ratio                                          | moist (dry+vapor, no condensates) mixing ratio of cloud water (condensate) | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | qi                  | ice_water_mixing_ratio                                                      | moist (dry+vapor, no condensates) mixing ratio of ice water                | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | T3D                 | air_temperature                                                             | layer mean air temperature                                                 | K       |    2 | real      | kind_phys | in     | F        |
!! | CLDFRA              | total_cloud_fraction                                                        | layer total cloud fraction                                                 | frac    |    2 | real      | kind_phys | inout  | F        |
!! | qc_save             | cloud_liquid_water_mixing_ratio_save                                        | cloud liquid water mixing ratio before entering a physics scheme           | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | qi_save             | cloud_ice_water_mixing_ratio_save                                           | cloud ice water mixing ratio before entering a physics scheme              | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | QC_BL               | subgrid_cloud_mixing ratio_pbl                                              | subgrid cloud cloud mixing ratio from PBL scheme                           | kg kg-1 |    2 | real      | kind_phys | in     | F        |
!! | CLDFRA_BL           | subgrid_cloud_fraction_pbl                                                  | subgrid cloud fraction from PBL scheme                                     | frac    |    2 | real      | kind_phys | in     | F        |
!! | errmsg              | ccpp_error_message                                                          | error message for error handling in CCPP                                   | none    |    0 | character | len=*     | out    | F        |
!! | errflg              | ccpp_error_flag                                                             | error flag for error handling in CCPP                                      | flag    |    0 | integer   |           | out    | F        |
!!
#endif
!###===================================================================
SUBROUTINE mynnrad_pre_run(                &
     &     ix,im,levs,                     &
     &     qc,qi,T3D,cldfra,               &
     &     qc_save, qi_save,               &
     &     qc_bl,cldfra_bl,                &
     &     errmsg, errflg                  )

! should be moved to inside the mynn:
      use machine , only : kind_phys

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 
      ! Interface variables
      integer, intent(in)  :: ix, im, levs
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: qc, qi
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: T3D
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: cldfra
      real(kind=kind_phys), dimension(im,levs), intent(out)   :: qc_save, qi_save
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: qc_bl, cldfra_bl
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer              :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      write(0,*)"=============================================="
      write(0,*)"in mynn rad pre"

     ! Add subgrid cloud information:
        do k = 1, levs
           do i = 1, im

              qc_save(i,k) = qc(i,k)
              qi_save(i,k) = qi(i,k)

              IF (qc(i,k) < 1.E-6 .AND. qi(i,k) < 1.E-6 .AND. CLDFRA_BL(i,k)>0.001) THEN
                !Partition the BL clouds into water & ice according to a linear
                !approximation of Hobbs et al. (1974). This allows us to only use
                !one 3D array for both cloud water & ice.
!               Wice = 1. - MIN(1., MAX(0., (t(i,k)-254.)/15.))
!               Wh2o = 1. - Wice
                CLDFRA(i,k)=MAX(CLDFRA(i,k),CLDFRA_BL(i,k))
                CLDFRA(i,k)=MAX(0.0,MIN(1.0,CLDFRA(i,k)))
                qc(i,k)=qc(i,k) + QC_BL(i,k)*(MIN(1., MAX(0., (T3D(i,k)-254.)/15.)))*CLDFRA_BL(i,k)
                qi(i,k)=qi(i,k) + QC_BL(i,k)*(1. - MIN(1., MAX(0., (T3D(i,k)-254.)/15.)))*CLDFRA_BL(i,k)
              ENDIF

           enddo
        enddo


        print*,"===Finished adding subgrid clouds to the resolved-scale clouds"
        !print*,"qc_save:",qc_save(1,1)," qi_save:",qi_save(1,1)

  END SUBROUTINE mynnrad_pre_run

!###=================================================================

END MODULE mynnrad_pre
