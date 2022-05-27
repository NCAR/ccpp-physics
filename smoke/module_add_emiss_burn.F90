!>\file  module_add_emiss_burn.F90
!! This file adds the biomass burning emissions to the smoke field.

module module_add_emiss_burn
!RAR: significantly modified for the new BB emissions
  use machine ,        only : kind_phys
  use rrfs_smoke_data
  use rrfs_smoke_config
CONTAINS
  subroutine add_emis_burn(data,dtstep,ktau,dz8w,rho_phy,rel_hum,    &
                           chem,julday,gmt,xlat,xlong,               &
                           !luf_igbp,lu_fire1,                       &
                           vegtype,vfrac,peak_hr,                    &
                           time_int,ebu,                             &   ! RAR
                           r_q,fhist,aod3d_smoke,aod3d_dust,         &
 !                         nwfa,nifa,                                &
                           rainc,rainnc, swdown,smoke_forecast,      &
                           ids,ide, jds,jde, kds,kde,                &
                           ims,ime, jms,jme, kms,kme,                &
                           its,ite, jts,jte, kts,kte                 )

!   USE module_configure, only: grid_config_rec_type
!   USE module_state_description
   IMPLICIT NONE
   type(smoke_data), intent(inout) :: data

!   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: ktau, julday,        &
                                  ids,ide, jds,jde, kds,kde,      &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   real(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem

   real(kind_phys), DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         INTENT(IN) ::                                   ebu

   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  xlat,xlong, rainc,rainnc,swdown, peak_hr, vfrac
   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(OUT)    ::  r_q    ! RAR:
   real(kind_phys), DIMENSION(ims:ime,jms:jme), INTENT(INOUT)  ::  fhist  ! RAR:
   real(kind_phys), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT) ::  aod3d_smoke, aod3d_dust    ! RAR:
   integer, DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  vegtype

   real(kind_phys), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: dz8w,rho_phy,rel_hum
!   real(kind_phys), DIMENSION(ims:ime,1:nlcat,jms:jme), INTENT(IN) :: luf_igbp

!   real(kind_phys),  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
!          OPTIONAL, INTENT(INOUT   ) ::                 nwfa,nifa                   ! RAR:

    real(kind_phys), INTENT(IN) ::  dtstep, gmt
    real(kind_phys), INTENT(IN) ::  time_int       ! RAR: time in seconds since start of simulation
    logical,         INTENT(IN) ::  smoke_forecast

    integer :: i,j,k,n,m
    real(kind_phys) :: conv_rho, conv, ext2, dm_smoke, daero_num_wfa, daero_num_ifa !, lu_sum1_5, lu_sum12_14
    !real(kind_phys) :: ebumax
!    CHARACTER (LEN=80) :: message

    INTEGER, PARAMETER :: kfire_max=35    ! max vertical level for BB plume rise
    ! Diameters and standard deviations for emissions
    ! the diameters are the volume (mass) geometric mean diameters, following MADE_SORGAM
    real(kind_phys), PARAMETER :: dgvem_i= 0.08E-6 !0.03E-6 ! [ m ]
    real(kind_phys), PARAMETER :: sgem_i = 1.8     !1.7

    ! *** Accumulation mode:
    real(kind_phys), PARAMETER :: dgvem_j= 0.3E-6 ! [ m ]
    real(kind_phys), PARAMETER :: sgem_j = 2.0

    ! *** Coarse mode
    real(kind_phys), PARAMETER :: dgvem_c= 6.0E-6 ! [ m ] 
    real(kind_phys), PARAMETER :: sgem_c=  2.2
    real(kind_phys), PARAMETER :: pic= 3.14159

    ! RAR: factors for getting number emissions rate from mass emissions rate following made_sorgam
    real(kind_phys), PARAMETER :: fact_numn= 1.e-9*6.0/pic*exp(4.5*log(sgem_i)**2)/dgvem_i**3       ! Aitken mode
    real(kind_phys), PARAMETER :: fact_numa= 1.e-9*6.0/pic*exp(4.5*log(sgem_j)**2)/dgvem_j**3       ! accumulation mode
    real(kind_phys), PARAMETER :: fact_numc= 1.e-9*6.0/pic*exp(4.5*log(sgem_c)**2)/dgvem_c**3       ! coarse mode

    real(kind_phys), PARAMETER :: dens_oc_aer=1.4e3, dens_ec_aer=1.7e3  ! kg/m3
!    real(kind_phys), PARAMETER :: rinti=2.1813936e-8, cx=2.184936* 3600, timeq_max=3600.*24.   ! constants for the diurnal cycle calculations
    real(kind_phys), PARAMETER :: ax1=531., cx1=7800.  ! For cropland, urban and small fires
!    real(kind_phys), PARAMETER :: rinti=2.1813936e-8, ax2=3200., const2=100., coef2=10.6712963e-4, cx=2.184936* 3600, timeq_max=3600.*24.
    real(kind_phys), PARAMETER :: rinti=2.1813936e-8, ax2=3400., const2=130., coef2=10.6712963e-4, cx2=7200., timeq_max=3600.*24.    ! New parameters
    real(kind_phys), PARAMETER :: sc_me= 4.0, ab_me=0.5     ! m2/g, scattering and absorption efficiency for smoke

!   Parameters used for the wfa and ifa in mp physics per Trude E. (NCAR)
!   Water friendly: radius: 0.04 micron, standard deviation: 1.8, kappa (for hygroscopic growth): 0.2, real index of refraction: 1.53, imaginary index of refraction: 1e-7
!   Ice friendly: radius: 0.4 micron, standard deviation: 1.8, kappa : 0.04, real index of refraction: 1.56, imaginary index of refraction: 3e-3

    !    real, parameter :: cx        =  2.184936 * 3600., rinti     =  2.1813936e-8    , ax        =  2000.6038
    ! bx_bburn  =  20.041288 * 3600.,  RAR: this depends on the vegetation class, location (local time) etc.
    real(kind_phys) :: timeq, dt1,dt2,dtm         ! For BB emis. diurnal cycle calculation

    timeq= gmt*3600. + real(time_int,4)
    timeq= mod(timeq,timeq_max)

! Main loops to add BB emissions
    do j=jts,jte
       do i=its,ite
          !if( luf_igbp(i,17,j)>0.99 .OR. ebu(i,1,j,p_ebu_smoke) < 1.e-6) cycle       ! no BB emissions or water pixels
          if( (1.-vfrac (i,j))>0.99 .OR. ebu(i,1,j) < 1.e-6) cycle       ! no BB emissions or water pixels

          ! RAR: the decrease in the BB emissions after >18 hrs of forecast, the decrease occurs at night. The decrease occurs at night.
          IF (time_int>64800. .AND. swdown(i,j)<.1 .AND. fhist(i,j)>.75 ) THEN
              fhist(i,j)= 0.75
          ENDIF

          IF (time_int>129600. .AND. swdown(i,j)<.1 .AND. fhist(i,j)>.5 ) THEN      ! After 36 hr forecast 
              fhist(i,j)= 0.5
          ENDIF

          IF ( (rainc(i,j) + rainnc(i,j))>=10. .AND. fhist(i,j)>.3 ) THEN    ! If it rains more than 1cm, then the BB emissions are reduced
              fhist(i,j)= 0.3
          ENDIF

! RAR: Grasslands (29% of ther western HRRR CONUS domain) probably also need to be added below, check this later
! RAR: In the HRRR CONUS domain (western part) crop 11%, 2% cropland/natural vegetation and 0.4% urban of pixels
!.OR. lu_index(i,j)==14)  then ! Croplands(12), Urban and Built-Up(13), cropland/natural vegetation (14) mosaic in MODI-RUC vegetation classes
! Peak hours for the fire activity depending on the latitude
!                   if (xlong(i,j)<-130.) then  max_ti= 24.041288* 3600.    ! peak at 24 UTC, fires in  Alaska
!                   elseif (xlong(i,j)<-100.) then   max_ti= 22.041288* 3600.    ! peak at 22 UTC, fires in the western US
!                   elseif (xlong(i,j)<-70.) then   ! peak at 20 UTC, fires in the eastern US,   max_ti= 20.041288* 3600.
!                   else   max_ti= 18.041288* 3600.
!                   endif

           !IF ( lu_fire1(i,j)>0.9 ) then    !Ag, urban fires, bare land etc.
           IF ( vegtype(i,j)==12 .or. vegtype(i,j)==13 ) then    !Ag, urban fires, bare land etc.
           !   these fires will have exponentially decreasing diurnal cycle, these fires decrease 55% in 2 hours, end in 5 hours
               r_q(i,j) = rinti* ax1 * exp(- (time_int**2)/(cx1**2) )
           ELSE
               ! RAR: Gaussian profile for wildfires
               dt1= abs(timeq - peak_hr(i,j))
               dt2= timeq_max - peak_hr(i,j) + timeq   ! peak hour is always <86400.
               dtm= MIN(dt1,dt2)
               r_q(i,j) = rinti*( ax2 * exp(- dtm**2/(2.*cx2**2) ) + const2 - coef2*timeq )
           ENDIF

           r_q(i,j) = fhist(i,j)* max(0.,r_q(i,j)*timeq_max)

           !IF (swdown(i,j)<.1) THEN
           !    r_q(i,j)= MIN(0.5,r_q(i,j))   ! lower BB emissions at night
           !ENDIF

          !IF (.NOT. config_flags%bb_dcycle) THEN
           !IF (.NOT. bb_dcycle) THEN
           !    r_q(i,j)= fhist(i,j)    ! no diurnal cycle
           !END IF

          !IF (.NOT. smoke_forecast) THEN
               r_q(i,j)= 1.
          !END IF

           do k=kts,kfire_max
              conv= r_q(i,j)*dtstep/(rho_phy(i,k,j)* dz8w(i,k,j))

              ! RAR: in this case tracer_1 is fire emitted CO
              ! conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
              ! chem(i,k,j,p_tracer_1) = chem(i,k,j,p_tracer_1) + ebu(i,k,j,p_ebu_co)*conv_rho

!              dm_oc_bb = conv* ebu(i,k,j,p_ebu_oc)     ! Assume that BB primary PM25 is mostly OC, 1.25 is OM/OC ratio
!              dm_p25_bb= conv* ebu(i,k,j,p_ebu_pm25)
!              dm_ec_bb = conv* ebu(i,k,j,p_ebu_bc)
!              dm_smk = conv* ebu(i,k,j,p_ebu_smoke)
              !IF (k==kts) THEN                         ! Partition takes place here to avoid double counting of smold. and flam. BB emiss.
              !   C11= (1.-flam_frac(i,j))*r_q(i,j)
              !ELSE
              !   C11= flam_frac(i,j)*r_q(i,j)
              !ENDIF
              dm_smoke= conv*ebu(i,k,j)
!              print*,'hli dm_smoke',dm_smoke,conv,ebu(i,k,j,p_ebu_smoke)

              chem(i,k,j,p_smoke) = chem(i,k,j,p_smoke) + dm_smoke
              chem(i,k,j,p_smoke) = MIN(chem(i,k,j,p_smoke),5.e+3)

             if (ktau<1000 .and. dbg_opt) then
             !  if ( k==kts ) then
             !    WRITE(6,*) 'add_emiss_burn: ktau,gmt,dtstep,time_int ',ktau,gmt,dtstep,time_int
             !    WRITE(*,*) 'add_emiss_burn: i,j,xlat(i,j),xlong(i,j) ',i,j,xlat(i,j),xlong(i,j)
                 !WRITE(*,*) 'add_emiss_burn: luf_igbp(i,:,j) ',luf_igbp(i,:,j)
                 !WRITE(*,*) 'add_emiss_burn: lu_fire1(i,j) ',lu_fire1(i,j)
             !    WRITE(6,*) 'add_emiss_burn: timeq,peak_hr(i,j),fhist(i,j),r_q(i,j) ',timeq,peak_hr(i,j),fhist(i,j),r_q(i,j)
             !    WRITE(*,*) 'add_emiss_burn: rainc(i,j),rainnc(i,j) ', rainc(i,j),rainnc(i,j)
             !  endif
               if ( k==kts .OR. k==kfire_max ) then
                 WRITE(6,*) 'add_emiss_burn: i,j,k ',i,j,k
                 WRITE(6,*) 'add_emiss_burn: rho_phy(i,k,j),dz8w(i,k,j),conv ',rho_phy(i,k,j),dz8w(i,k,j),conv
                 WRITE(6,*) 'add_emiss_burn: ebu(i,k,j),dm_smoke ', ebu(i,k,j),dm_smoke
               endif
             endif

              enddo
            enddo
          enddo

          ext2= sc_me + ab_me   
          do j=jts,jte
           do k=kts,kte
            do i=its,ite
 
               ! Check for NaNs, negative and too large numbers
               IF (.NOT. (chem(i,k,j,p_smoke)>=0. .AND. chem(i,k,j,p_smoke)<1.1e+4)) THEN
                   chem(i,k,j,p_smoke)=1.e-16
               END IF

               aod3d_smoke(i,k,j)= 1.e-6* ext2* chem(i,k,j,p_smoke )*rho_phy(i,k,j)*dz8w(i,k,j)
               aod3d_dust (i,k,j)= 1.e-6* ext2* chem(i,k,j,p_dust_1)*rho_phy(i,k,j)*dz8w(i,k,j)
            enddo
           enddo
          enddo

     IF ( ktau<2000 .and. dbg_opt ) then
         WRITE(*,*) 'add_emis_burn: i,j,k,ext2 ',i,j,k,ext2
         WRITE(*,*) 'add_emis_burn: rel_hum(its,kts,jts),rel_hum(ite,kfire_max,jte) ',rel_hum(its,kts,jts),rel_hum(ite,kfire_max,jte)
         WRITE(*,*) 'add_emis_burn: aod3d_smoke(its,kts,jts),aod3d_smoke(ite,kfire_max,jte) ',aod3d_smoke(its,kts,jts),aod3d_smoke(ite,kfire_max,jte)
         WRITE(*,*) 'add_emis_burn: aod3d_dust(its,kts,jts),aod3d_dust(ite,kfire_max,jte) ',aod3d_dust(its,kts,jts),aod3d_dust(ite,kfire_max,jte)
     END IF

!     CASE DEFAULT
!       call wrf_debug(15,'nothing done with burn emissions for chem array')
!    END SELECT emiss_select

    END subroutine add_emis_burn

END module module_add_emiss_burn
