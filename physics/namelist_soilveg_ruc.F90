!>\file namelist_soilveg_ruc.F90
!>\ingroup RUC_lsm

      module namelist_soilveg_ruc

      use machine ,   only : kind_phys

      implicit none
      save

      INTEGER MAX_SLOPETYP
      INTEGER MAX_SOILTYP
      INTEGER MAX_VEGTYP

      PARAMETER(MAX_SLOPETYP = 30)
      PARAMETER(MAX_SOILTYP = 30)
      PARAMETER(MAX_VEGTYP = 30)

      real(kind_phys) SLOPE_DATA(MAX_SLOPETYP)
!> vegetation
      real(kind_phys) ALBTBL(MAX_VEGTYP)
      real(kind_phys) Z0TBL(MAX_VEGTYP)
      real(kind_phys) LEMITBL(MAX_VEGTYP)
      real(kind_phys) PCTBL(MAX_VEGTYP)
      real(kind_phys) SHDTBL(MAX_VEGTYP)
      INTEGER IFORTBL(MAX_VEGTYP)
      real(kind_phys) RSTBL(MAX_VEGTYP)
      real(kind_phys) RGLTBL(MAX_VEGTYP)
      real(kind_phys) HSTBL(MAX_VEGTYP)
      real(kind_phys) SNUPTBL(MAX_VEGTYP)
      real(kind_phys) LAITBL(MAX_VEGTYP)
      real(kind_phys) MAXALB(MAX_VEGTYP)
      real(kind_phys) MFSNO(MAX_VEGTYP)
      real(kind_phys) SNCOVFAC(MAX_VEGTYP)
      LOGICAL LPARAM
      real(kind_phys) TOPT_DATA
      real(kind_phys) CMCMAX_DATA
      real(kind_phys) CFACTR_DATA
      real(kind_phys) RSMAX_DATA
      INTEGER BARE
      INTEGER GLACIER
      INTEGER NATURAL
      INTEGER CROP
      INTEGER URBAN
      INTEGER DEFINED_VEG
      INTEGER DEFINED_SOIL
      INTEGER DEFINED_SLOPE
!> -- soils
      real(kind_phys) BB(MAX_SOILTYP)
      real(kind_phys) DRYSMC(MAX_SOILTYP)
      real(kind_phys) HC(MAX_SOILTYP)
      real(kind_phys) MAXSMC(MAX_SOILTYP)
      real(kind_phys) REFSMC(MAX_SOILTYP)
      real(kind_phys) SATPSI(MAX_SOILTYP)
      real(kind_phys) SATDK(MAX_SOILTYP)
      real(kind_phys) SATDW(MAX_SOILTYP)
      real(kind_phys) WLTSMC(MAX_SOILTYP)
      real(kind_phys) QTZ(MAX_SOILTYP)
      real(kind_phys) REFSMCnoah(MAX_SOILTYP)
      real(kind_phys) WLTSMCnoah(MAX_SOILTYP)
      real(kind_phys) BBnoah(MAX_SOILTYP)
      real(kind_phys) SATDKnoah(MAX_SOILTYP)
      real(kind_phys) SATPSInoah(MAX_SOILTYP)
      real(kind_phys) MAXSMCnoah(MAX_SOILTYP)
      end module namelist_soilveg_ruc
