!>\file module_get_aerosols_for_mp.F90
!! This contains routines used for to convert MERRA2 aerosol to water-friendly
!! and ice-friendly concentrations for TEMPO and Thompson microphysics schemes. 

module module_get_aerosols_for_mp

    use machine, only : kind_phys

    implicit none

    private
    public :: get_niwfa

    contains
    subroutine get_niwfa(aerfld, nifa, nwfa, ncol, nlev)
         ! To calculate nifa and nwfa from bins of aerosols.
         ! In GOCART and MERRA2, aerosols are given as mixing ratio (kg/kg). To
         ! convert from kg/kg to #/kg, the "unit mass" (mass of one particle)
         ! within the mass bins is calculated. A lognormal size distribution
         ! within aerosol bins is used to find the size based upon the median
         ! mass. NIFA is mainly summarized over five dust bins and NWFA over the
         ! other 10 bins. The parameters besides each bins are carefully tuned
         ! for a good performance of the scheme.
         !
         ! The fields for the last index of the aerfld array
         ! are specified as below.
         ! 1: dust bin 1,                     0.1 to 1.0  micrometers
         ! 2: dust bin 2,                     1.0 to 1.8  micrometers
         ! 3: dust bin 3,                     1.8 to 3.0  micrometers
         ! 4: dust bin 4,                     3.0 to 6.0  micrometers
         ! 5: dust bin 5,                     6.0 to 10.0 micrometers
         ! 6: sea salt bin 1,                 0.03 to 0.1 micrometers
         ! 7: sea salt bin 2,                 0.1 to 0.5  micrometers
         ! 8: sea salt bin 3,                 0.5 to 1.5  micrometers 
         ! 9: sea salt bin 4,                 1.5 to 5.0  micrometers
         ! 10: sea salt bin 5,                5.0 to 10.0 micrometers
         ! 11: Sulfate,                       0.35 (mean) micrometers
         ! 15: water-friendly organic carbon, 0.35 (mean) micrometers
         !
         ! Bin densities are as follows:
         ! 1:    dust bin 1:         2500 kg/m2
         ! 2-5:  dust bin 2-5:       2650 kg/m2
         ! 6-10: sea salt bins 6-10: 2200 kg/m2
         ! 11:   sulfate:            1700 kg/m2
         ! 15:   organic carbon:     1800 kg/m2
         
         integer, intent(in)::ncol, nlev
         real (kind=kind_phys), dimension(:,:,:), intent(in)  :: aerfld
         real (kind=kind_phys), dimension(:,:),   intent(out ):: nifa, nwfa

         nifa=(aerfld(:,:,1)/4.0737762+aerfld(:,:,2)/30.459203+aerfld(:,:,3)/153.45048+ &
              aerfld(:,:,4)/1011.5142+ aerfld(:,:,5)/5683.3501)*1.e15

         nwfa=((aerfld(:,:,6)/0.0045435214+aerfld(:,:,7)/0.2907854+aerfld(:,:,8)/12.91224+ &
              aerfld(:,:,9)/206.2216+ aerfld(:,:,10)/4326.23)*9.+aerfld(:,:,11)/0.3053104*5+ &
              aerfld(:,:,15)/0.3232698*8)*1.e15
    end subroutine get_niwfa

end module module_get_aerosols_for_mp