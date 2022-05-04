!>\file get_phi_fv3.F90
!! This file contains a subroutine to calculate geopotential from within physics.

module get_phi_fv3

   use machine,  only: kind_phys
   use physcons, only: con_fvirt

!--- public declarations
   public get_phi_fv3_run

!--- local variables
   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
   real(kind=kind_phys), parameter :: half = 0.5_kind_phys
   real(kind=kind_phys), parameter :: one  = 1.0_kind_phys

contains

!! \section arg_table_get_phi_fv3_run Argument Table
!! \htmlinclude get_phi_fv3_run.html
!!
   subroutine get_phi_fv3_run(ix, levs, con_fvirt, gt0, gq01, del_gz, phii, phil, errmsg, errflg)

     implicit none

     ! Interface variables
     integer, intent(in) :: ix, levs
     real(kind=kind_phys), intent(in) :: con_fvirt
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: gt0
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: gq01
     real(kind=kind_phys), dimension(:,:),     intent(inout) :: del_gz
     real(kind=kind_phys), dimension(:,:),     intent(out)   :: phii
     real(kind=kind_phys), dimension(:,:),     intent(out)   :: phil
     character(len=*),                         intent(out)   :: errmsg
     integer,                                  intent(out)   :: errflg

     ! Local variables
     integer :: i, k

     ! Initialize CCPP error handling variables
     errmsg = ''
     errflg = 0

! SJL: Adjust the height hydrostatically in a way consistent with FV3 discretization
     do i=1,ix
        phii(i,1) = zero
     enddo
     do k=1,levs
       do i=1,ix
         del_gz(i,k) = del_gz(i,k)*gt0(i,k) *                          &
     &                 (one + con_fvirt*max(zero,gq01(i,k)))
         phii(i,k+1) = phii(i,k) + del_gz(i,k)
         phil(i,k)   = half*(phii(i,k) + phii(i,k+1))
       enddo
     enddo

   end subroutine get_phi_fv3_run

end module get_phi_fv3