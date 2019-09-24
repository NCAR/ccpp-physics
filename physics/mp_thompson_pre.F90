!>\file mp_thompson_pre.F90
!!

! CCPP license goes here, as well as further documentation
!>\ingroup aathompson
module mp_thompson_pre

      use machine, only : kind_phys

      use module_mp_thompson, only : naIN0, naIN1, naCCN0, naCCN1, eps

      use module_mp_thompson_make_number_concentrations, only: make_IceNumber, make_DropletNumber, make_RainNumber

      implicit none

      public :: mp_thompson_pre_init, mp_thompson_pre_run, mp_thompson_pre_finalize

      private

   contains

      subroutine mp_thompson_pre_init()
      end subroutine mp_thompson_pre_init

#if 0
!! \section arg_table_mp_thompson_pre_run Argument Table
!! \htmlinclude mp_thompson_pre_run.html
!!
#endif
      subroutine mp_thompson_pre_run(ncol, nlev, kdt, con_g, con_rd,    &
                                  spechum, qc, qr, qi, qs, qg, ni, nr,       &
                                  is_aerosol_aware, nc, nwfa, nifa, nwfa2d,  &
                                  nifa2d, tgrs, tgrs_save, prsl, phil, area, &
                                  mpirank, mpiroot, blkno, errmsg, errflg)

         implicit none

         ! Interface variables
         ! Dimensions and constants
         integer,                   intent(in   ) :: ncol
         integer,                   intent(in   ) :: nlev
         integer,                   intent(in   ) :: kdt
         real(kind_phys),           intent(in   ) :: con_g
         real(kind_phys),           intent(in   ) :: con_rd
         ! Hydrometeors
         real(kind_phys),           intent(inout) :: spechum(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qc(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qr(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qi(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qs(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: qg(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: ni(1:ncol,1:nlev)
         real(kind_phys),           intent(inout) :: nr(1:ncol,1:nlev)
         ! Aerosols
         logical,                   intent(in   ) :: is_aerosol_aware
         real(kind_phys), optional, intent(inout) :: nc(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nwfa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nifa(1:ncol,1:nlev)
         real(kind_phys), optional, intent(inout) :: nwfa2d(1:ncol)
         real(kind_phys), optional, intent(inout) :: nifa2d(1:ncol)
         ! State variables and timestep information
         real(kind_phys),           intent(in   ) :: tgrs(1:ncol,1:nlev)
         real(kind_phys),           intent(  out) :: tgrs_save(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: prsl(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: phil(1:ncol,1:nlev)
         real(kind_phys),           intent(in   ) :: area(1:ncol)
         ! MPI information
         integer,                   intent(in   ) :: mpirank
         integer,                   intent(in   ) :: mpiroot
         ! Blocking information
         integer,                   intent(in   ) :: blkno
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Local variables
         real (kind=kind_phys) :: hgt(1:ncol,1:nlev)  ! m
         real (kind=kind_phys) :: rho(1:ncol,1:nlev)  ! kg m-3
         real (kind=kind_phys) :: orho(1:ncol,1:nlev) ! m3 kg-1
         real (kind=kind_phys) :: h_01, airmass, niIN3, niCCN3
         integer :: i, k

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Save current air temperature for tendency limiters in mp_thompson_post
         tgrs_save = tgrs

         ! Return if not first timestep
         if (kdt > 1) return

         ! Consistency check
         if (is_aerosol_aware     .and. &
             (.not.present(nc)     .or. &
              .not.present(nwfa2d) .or. &
              .not.present(nifa2d) .or. &
              .not.present(nwfa)   .or. &
              .not.present(nifa)        )) then
             write(errmsg,fmt='(*(a))') 'Logic error in mp_thompson_pre_run:',                 &
                                        ' aerosol-aware microphysics require all of the following', &
                                        ' optional arguments: nc, nwfa2d, nifa2d, nwfa, nifa'
             errflg = 1
             return
          end if

         ! Fix initial values of hydrometeors
         where(spechum<0) spechum = 0.0
         where(qc<0)      qc = 0.0
         where(qr<0)      qr = 0.0
         where(qi<0)      qi = 0.0
         where(qs<0)      qs = 0.0
         where(qg<0)      qg = 0.0
         where(ni<0)      ni = 0.0
         where(nr<0)      nr = 0.0

         if (is_aerosol_aware) then
            ! Fix initial values of aerosols
            where(nc<0)     nc = 0.0
            where(nwfa<0)   nwfa = 0.0
            where(nifa<0)   nifa = 0.0
            where(nwfa2d<0) nwfa2d = 0.0
            where(nifa2d<0) nifa2d = 0.0
         end if

         ! Geopotential height in m2 s-2 to height in m
         hgt = phil/con_g

         ! Density of air in kg m-3 and inverse density of air
         rho = prsl/(con_rd*tgrs)
         orho = 1.0/rho

         !  Prior to calling the functions: make_DropletNumber, make_IceNumber, make_RainNumber,
         !  the incoming mixing ratios should be converted to units of mass/num per cubic meter
         !  rather than per kg of air.  So, to pass back to the model state variables,
         !  they also need to be switched back to mass/number per kg of air, because
         !  what is returned by the functions is in units of number per cubic meter.

         ! If qi is in boundary conditions but ni is not, calculate ni from qi, rho and tgrs
         if (maxval(qi)>0.0 .and. maxval(ni)==0.0) then
             ni = make_IceNumber(qi*rho, tgrs) * orho
         end if

         ! If ni is in boundary conditions but qi is not, reset ni to zero
         if (maxval(ni)>0.0 .and. maxval(qi)==0.0) ni = 0.0

         ! If qr is in boundary conditions but nr is not, calculate nr from qr, rho and tgrs
         if (maxval(qr)>0.0 .and. maxval(nr)==0.0) then
             nr = make_RainNumber(qr*rho, tgrs) * orho
         end if

         ! If nr is in boundary conditions but qr is not, reset nr to zero
         if (maxval(nr)>0.0 .and. maxval(qr)==0.0) nr = 0.0

         ! Return if aerosol-aware option is not used
         if (.not. is_aerosol_aware) return

!..Check for existing aerosol data, both CCN and IN aerosols.  If missing
!.. fill in just a basic vertical profile, somewhat boundary-layer following.

!.. CCN
         if (MAXVAL(nwfa) .lt. eps) then
            if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial CCN aerosols.'
            do i = 1, ncol
               if (hgt(i,1).le.1000.0) then
                  h_01 = 0.8
               elseif (hgt(i,1).ge.2500.0) then
                  h_01 = 0.01
               else
                  h_01 = 0.8*cos(hgt(i,1)*0.001 - 1.0)
               endif
               niCCN3 = -1.0*ALOG(naCCN1/naCCN0)/h_01
               nwfa(i,1) = naCCN1+naCCN0*exp(-((hgt(i,2)-hgt(i,1))/1000.)*niCCN3)
               airmass = 1./orho(i,1) * (hgt(i,2)-hgt(i,1))*area(i) ! kg
               nwfa2d(i) = nwfa(i,1) * 0.000196 * (airmass*2.E-10)
               do k = 2, nlev
                  nwfa(i,k) = naCCN1+naCCN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niCCN3)
               enddo
            enddo
         else
            if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently initial CCN aerosols are present.'
            if (MAXVAL(nwfa2d) .lt. eps) then
! Hard-coded switch between new (from WRFv4.0, top) and old (until WRFv3.9.1.1, bottom) surface emission rate calculations
#if 0
               !+---+-----------------------------------------------------------------+
               !..Scale the lowest level aerosol data into an emissions rate.  This is
               !.. very far from ideal, but need higher emissions where larger amount
               !.. of (climo) existing and lesser emissions where there exists fewer to
               !.. begin as a first-order simplistic approach.  Later, proper connection to
               !.. emission inventory would be better, but, for now, scale like this:
               !.. where: Nwfa=50 per cc, emit 0.875E4 aerosols per second per grid box unit
               !..        that was tested as ~(20kmx20kmx50m = 2.E10 m**-3)
               !+---+-----------------------------------------------------------------+
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial CCN aerosol surface emission rates.'
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Use new (WRFv4+) formula to calculate CCN surface emission rates.'
               do i = 1, ncol
                  airmass = 1./orho(i,1) * (hgt(i,2)-hgt(i,1))*area(i) ! kg
                  nwfa2d(i) = nwfa(i,1) * 0.000196 * (airmass*2.E-10)
               enddo
#else
               !+---+-----------------------------------------------------------------+
               !..Scale the lowest level aerosol data into an emissions rate.  This is
               !.. very far from ideal, but need higher emissions where larger amount
               !.. of existing and lesser emissions where not already lots of aerosols
               !.. for first-order simplistic approach.  Later, proper connection to
               !.. emission inventory would be better, but, for now, scale like this:
               !.. where: Nwfa=50 per cc, emit 0.875E4 aerosols per kg per second
               !..        Nwfa=500 per cc, emit 0.875E5 aerosols per kg per second
               !..        Nwfa=5000 per cc, emit 0.875E6 aerosols per kg per second
               !.. for a grid with 20km spacing and scale accordingly for other spacings.
               !+---+-----------------------------------------------------------------+
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial CCN aerosol surface emission rates.'
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Use old (pre WRFv4) formula to calculate CCN surface emission rates.'
               do i = 1, ncol
                  if (SQRT(area(i))/20000.0 .ge. 1.0) then
                     h_01 = 0.875
                  else
                     h_01 = (0.875 + 0.125*((20000.-SQRT(area(i)))/16000.)) * SQRT(area(i))/20000.
                  endif
                  nwfa2d(i) = 10.0**(LOG10(nwfa(i,1)*1.E-6)-3.69897)
                  nwfa2d(i) = nwfa2d(i)*h_01 * 1.E6
               enddo
#endif
            else
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently initial CCN aerosol surface emission rates are present.'
            endif
         endif

!.. IN
         if (MAXVAL(nifa) .lt. eps) then
            if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial IN aerosols.'
            do i = 1, ncol
               if (hgt(i,1).le.1000.0) then
                  h_01 = 0.8
               elseif (hgt(i,1).ge.2500.0) then
                  h_01 = 0.01
               else
                  h_01 = 0.8*cos(hgt(i,1)*0.001 - 1.0)
               endif
               niIN3 = -1.0*ALOG(naIN1/naIN0)/h_01
               nifa(i,1) = naIN1+naIN0*exp(-((hgt(i,2)-hgt(i,1))/1000.)*niIN3)
               nifa2d(i) = 0.
               do k = 2, nlev
                  nifa(i,k) = naIN1+naIN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niIN3)
               enddo
            enddo
         else
            if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently initial IN aerosols are present.' 
            if (MAXVAL(nifa2d) .lt. eps) then
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently there are no initial IN aerosol surface emission rates, set to zero.'
               ! calculate IN surface flux here, right now just set to zero
               nifa2d = 0.
            else
               if (mpirank==mpiroot .and. blkno==1) write(*,*) ' Apparently initial IN aerosol surface emission rates are present.'
            endif
         endif

         ! If qc is in boundary conditions but nc is not, calculate nc from qc, rho and nwfa
         if (maxval(qc)>0.0 .and. maxval(nc)==0.0) then
            nc = make_DropletNumber(qc*rho, nwfa) * orho
         end if

         ! If nc is in boundary conditions but qc is not, reset nc to zero
         if (maxval(nc)>0.0 .and. maxval(qc)==0.0) nc = 0.0

      end subroutine mp_thompson_pre_run

      subroutine mp_thompson_pre_finalize()
      end subroutine mp_thompson_pre_finalize

end module mp_thompson_pre
