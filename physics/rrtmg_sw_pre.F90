!>\file rrtmg_sw_pre.f90
!! This file contains a subroutine to module_radiation_surface::setalb() to
!! setup surface albedo for SW radiation.
      module rrtmg_sw_pre
      contains

!>\defgroup rrtmg_sw_pre GFS RRTMG scheme Pre
!! @{
      subroutine rrtmg_sw_pre_init ()
      end subroutine rrtmg_sw_pre_init

!> \section arg_table_rrtmg_sw_pre_run Argument Table
!! \htmlinclude rrtmg_sw_pre_run.html
!!
      subroutine rrtmg_sw_pre_run (im, lndp_type, n_var_lndp, lsswr, lndp_var_list, lndp_prt_list, tsfg, tsfa, coszen,     &
        alb1d, slmsk, snowd, sncovr, snoalb, zorl, hprime, alvsf, alnsf, alvwf,&
        alnwf, facsf, facwf, fice, tisfc, sfalb, nday, idxday, sfcalb1,        &
        sfcalb2, sfcalb3, sfcalb4, errmsg, errflg)

      use machine,                   only: kind_phys

      use module_radiation_surface,  only: NF_ALBD, setalb

      implicit none

      integer,                              intent(in)    :: im, lndp_type, n_var_lndp
      character(len=3)    , dimension(:),   intent(in)    :: lndp_var_list
      logical,                              intent(in)    :: lsswr
      real(kind=kind_phys), dimension(:),   intent(in)    :: lndp_prt_list
      real(kind=kind_phys), dimension(im),  intent(in)    :: tsfg, tsfa, coszen
      real(kind=kind_phys), dimension(im),  intent(in)    :: alb1d
      real(kind=kind_phys), dimension(im),  intent(in)    :: slmsk, snowd,     &
                                                             sncovr, snoalb,   &
                                                             zorl, hprime,     &
                                                             alvsf, alnsf,     &
                                                             alvwf, alnwf,     &
                                                             facsf, facwf,     &
                                                             fice, tisfc
      real(kind=kind_phys), dimension(im),  intent(inout) :: sfalb
      integer,                              intent(out)   :: nday
      integer, dimension(im),               intent(out)   :: idxday
      real(kind=kind_phys), dimension(im),  intent(out)   :: sfcalb1, sfcalb2, &
                                                             sfcalb3, sfcalb4
      character(len=*),                     intent(out)   :: errmsg
      integer,                              intent(out)   :: errflg
      ! Local variables
      integer :: i
      real(kind=kind_phys), dimension(im,NF_ALBD) :: sfcalb

      real(kind=kind_phys) :: lndp_alb

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  --- ...  start radiation calculations
!           remember to set heating rate unit to k/sec!
!> -# Start SW radiation calculations
      if (lsswr) then

!>  - Check for daytime points for SW radiation.
        nday = 0
        idxday = 0
        do i = 1, IM
          if (coszen(i) >= 0.0001) then
            nday = nday + 1
            idxday(nday) = i
          endif
        enddo

! set albedo pert, if requested.
        lndp_alb = -999.
        if (lndp_type==1) then
          do i =1,n_var_lndp
            if (lndp_var_list(i) == 'alb') then
                lndp_alb = lndp_prt_list(i)
            endif
          enddo
        endif

!>  - Call module_radiation_surface::setalb() to setup surface albedo.
!!  for SW radiation.

        call setalb (slmsk, snowd, sncovr, snoalb, zorl,  coszen, tsfg, tsfa,  &  !  ---  inputs
                     hprime, alvsf, alnsf, alvwf, alnwf, facsf, facwf, fice,   &
                     tisfc, IM, alb1d, lndp_alb,                               &  !  mg, sfc-perts
                     sfcalb)                                           !  ---  outputs

!> -# Approximate mean surface albedo from vis- and nir-  diffuse values.
        sfalb(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
      else
        nday   = 0
        idxday = 0
        sfcalb = 0.0
      endif

      do i = 1, im
        sfcalb1(i) = sfcalb(i,1)
        sfcalb2(i) = sfcalb(i,2)
        sfcalb3(i) = sfcalb(i,3)
        sfcalb4(i) = sfcalb(i,4)
      enddo

      end subroutine rrtmg_sw_pre_run

      subroutine rrtmg_sw_pre_finalize ()
      end subroutine rrtmg_sw_pre_finalize

!! @}
      end module rrtmg_sw_pre
