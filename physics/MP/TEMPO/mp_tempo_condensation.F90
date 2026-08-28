!>\file mp_tempo_condensation.F90
!! This file contains tempo condensation.


!>\defgroup tempo condensation
!! This module contains tempo condensation
module mp_tempo_condensation

      use mpi_f08
      use machine, only : kind_phys, kind_dyn, kind_dbl_prec

      use module_mp_tempo_params, only : roverrv, rv, rdry, eps, r1, t0, cp, initialize_parameters
      use module_mp_tempo_cfgs, only : ty_tempo_cfgs
      use module_mp_tempo_main, only : cloud_check_and_update
      use module_mp_tempo_utils, only : calc_rslf, get_constant_cloud_number

      implicit none

      public :: mp_tempo_condensation_init, mp_tempo_condensation_run, mp_tempo_condensation_final

      private

   contains

!! \section arg_table_mp_tempo_condensation_init Argument Table
!! \htmlinclude mp_tempo_condensation_init.html
!!
      subroutine mp_tempo_condensation_init(errmsg, errflg)

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         errmsg = ''
         errflg = 0
      end subroutine mp_tempo_condensation_init


!> \section arg_table_mp_tempo_condensation_run Argument Table
!! \htmlinclude mp_tempo_condensation_run.html
!!
      subroutine mp_tempo_condensation_run(mpicomm, mpirank, mpiroot, &
           tgrs, prsl, spechum, qc, mdt, errmsg, errflg, &
           is, ie, js, je, isc1, iec1, jsc1, jec1, &
           isc2, iec2, jsc2, jec2, isd, ied, jsd, jed, km)

         integer,                   intent(in   ) :: isc1, iec1, jsc1, jec1
         integer,                   intent(in   ) :: isd, ied, jsd, jed, is, ie, js, je, km, isc2, iec2, jsc2, jec2
         type(MPI_Comm),            intent(in   ) :: mpicomm
         integer,                   intent(in   ) :: mpirank
         integer,                   intent(in   ) :: mpiroot
         real(kind_dyn),            intent(in   ) :: mdt
         real(kind_dyn),            intent(inout) :: spechum(isd:ied, jsd:jed, 1:km)
         real(kind_dyn),            intent(inout) :: qc(isd:ied, jsd:jed, 1:km)
         real(kind_dyn),            intent(inout) :: tgrs(isd:ied, jsd:jed, 1:km)
         real(kind_dyn),            intent(in   ) :: prsl(is:ie, 1:km+1, js:je)

         real(kind_phys) :: qv(is:ie, js:je, 1:km)
         real(kind_phys) :: qc_mixing_ratio(is:ie, js:je, 1:km)
         real(kind_phys) :: nc3d(is:ie, js:je, 1:km)                  
         real(kind_phys) :: qvs(is:ie, js:je, 1:km)         
         real(kind_phys) :: rho(is:ie, js:je, 1:km)
         real(kind_phys) :: temp(is:ie, js:je, 1:km)
         real(kind_phys) :: satw(1:km)
         real(kind_phys) :: ssatw(1:km)
         real(kind_phys) :: lvap(1:km)
         real(kind_phys) :: tcond(1:km)
         real(kind_phys) :: diffu(1:km)
         real(kind_phys) :: ocp(1:km)
         real(kind_phys) :: lvt2(1:km)                                             
         real(kind_phys) :: rc(1:km)
         real(kind_phys) :: nc(1:km)
         real(kind_dbl_prec) :: ilamc(1:km)
         real(kind_phys) :: mvd_c(1:km)
         real(kind_phys) :: qcten(1:km)
         real(kind_phys) :: ncten(1:km)
         logical         :: l_qc(1:km)                                                                        
         real(kind_phys) :: condensation(is:ie, js:je, 1:km)         
         real(kind_phys) :: orho, clap, fcd, dfcd, xrc, odt
         
         ! CCPP error handling
         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg
         integer :: i, j, k, n
         logical, save :: need_tempo_params = .true.
        
         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0
         odt = 1._kind_phys/real(mdt,kind=kind_phys)

         if (need_tempo_params) then
            call initialize_parameters()
            need_tempo_params = .false.
         endif

         if (mpirank == mpiroot) write(*,*) 'tempo condensation on re-mapping timestep'
         do k = 1, km
            do j = js, je
               do i = is, ie
                  spechum(i,j,k) = max(1.e-10, spechum(i,j,k))
                  qv(i,j,k) = spechum(i,j,k)/(1.0_kind_phys-spechum(i,j,k))
                  qc_mixing_ratio(i,j,k) = real(qc(i,j,k),kind=kind_phys) / (1._kind_phys-real(spechum(i,j,k),kind=kind_phys)) ! convert from moist to dry rho
                  temp(i,j,k) = tgrs(i,j,k) / ((0.608 * spechum(i,j,k)) + 1.) ! tgrs is virtual temperature
                  rho(i,j,k) = roverrv * exp(prsl(i,k,j)) / (rdry*temp(i,j,k) * (qv(i,j,k)+roverrv))
                  qvs(i,j,k) = calc_rslf(real(exp(prsl(i,k,j)),kind=kind_phys), real(temp(i,j,k),kind=kind_phys)) ! prsl is log pressure
               enddo
            enddo
         enddo

         condensation = 0.         
         do j = js, je
            do i = is, ie
               l_qc = .false.
               rc = 0.
               nc = 0.
               qcten = 0.
               ncten = 0.

               call get_constant_cloud_number(nc=nc3d(i,j,:))

               ! returns cloud mass concentration
               call cloud_check_and_update(dt=real(mdt,kind=kind_phys), odt=odt, rho=rho(i,j,:), l_qc=l_qc, &
                    qc1d=qc_mixing_ratio(i,j,:), nc1d=nc3d(i,j,:), rc=rc, nc=nc, qcten=qcten, ncten=ncten, ilamc=ilamc, mvd_c=mvd_c)
                
               do k = 1, km
                  satw(k) = qv(i,j,k)/qvs(i,j,k)
                  ssatw(k) = satw(k) - 1.
                  if (abs(ssatw(k)) < eps) ssatw(k) = 0.
                  lvap(k) = 2.5e6 + (2106.0 - 4218.0)*(temp(i,j,k)-t0)
                  tcond(k) = (5.69 + 0.0168*(temp(i,j,k)-t0))*1.0e-5 * 418.936
                  diffu(k) = 2.11e-5*(temp(i,j,k)/t0)**1.94 * (101325./exp(prsl(i,k,j)))
                  ocp(k) = 1./(cp*(1.+0.887*qv(i,j,k)))         
                  lvt2(k) = lvap(k)*lvap(k)*ocp(k)*(1./rv)*(1./temp(i,j,k))*(1./temp(i,j,k))
                  
                  if (abs(ssatw(k)) >= eps) then
                     orho = 1./rho(i,j,k)
                     clap = (qv(i,j,k)-qvs(i,j,k))/(1. + lvt2(k)*qvs(i,j,k))
                     do n = 1, 3
                        fcd = qvs(i,j,k)*exp(lvt2(k)*clap) - qv(i,j,k) + clap
                        dfcd = qvs(i,j,k)*lvt2(k)*exp(lvt2(k)*clap) + 1.
                        clap = clap - fcd/dfcd
                     enddo
                     xrc = rc(k) + clap*rho(i,j,k)
                     
                     if (xrc > r1) then
                        condensation(i,j,k) = clap / mdt
                        
                        if (l_qc(k) .and. ssatw(k) < -1.e-6 .and. clap < -eps) then ! evaporation
                           condensation(i,j,k) = max(-rc(k)*0.99*orho/mdt, condensation(i,j,k))
                        endif
                     else
                        condensation(i,j,k) = -rc(k)*orho/mdt
                     endif
                  endif
                  
                  qv(i,j,k) = qv(i,j,k) - condensation(i,j,k)*mdt
                  qc_mixing_ratio(i,j,k) = qc_mixing_ratio(i,j,k) + condensation(i,j,k)*mdt
                  temp(i,j,k) = temp(i,j,k) + lvap(k)*ocp(k)*condensation(i,j,k)*mdt
               enddo

               call cloud_check_and_update(dt=real(mdt,kind=kind_phys), odt=odt, rho=rho(i,j,:), l_qc=l_qc, &
                    qc1d=qc_mixing_ratio(i,j,:), nc1d=nc3d(i,j,:), rc=rc, nc=nc, qcten=qcten, ncten=ncten, ilamc=ilamc, mvd_c=mvd_c)
            enddo
         enddo

         do k = 1, km
            do j = js, je
               do i = is, ie
                  spechum(i,j,k) = max(1.e-10, real(qv(i,j,k),kind=kind_dyn)/(1. + real(qv(i,j,k), kind=kind_dyn)))
                  qc(i,j,k) = real(qc_mixing_ratio(i,j,k), kind=kind_dyn)/(1. + real(qv(i,j,k), kind=kind_dyn))
                  tgrs(i,j,k) = temp(i,j,k)*((0.608*spechum(i,j,k)) + 1.)
               enddo
            enddo
         enddo
         
      end subroutine mp_tempo_condensation_run

!> \section arg_table_mp_tempo_condensation_final Argument Table
!! \htmlinclude mp_tempo_condensation_final.html
!!
      subroutine mp_tempo_condensation_final(errmsg, errflg)

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         errmsg = ''
         errflg = 0
      end subroutine mp_tempo_condensation_final

end module mp_tempo_condensation
