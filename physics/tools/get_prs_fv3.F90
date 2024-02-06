!>\file get_prs_fv3.F90
!! This file contains a subroutine to "adjust the geopotential height hydrostatically in a way consistent with FV3 discretization," 
!! according to SJ Lin.

module get_prs_fv3

   use machine,  only: kind_phys

!--- public declarations
   public get_prs_fv3_run

!--- local variables
   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
   real(kind=kind_phys), parameter :: one  = 1.0_kind_phys

contains

!> \section arg_table_get_prs_fv3_run Argument Table
!! \htmlinclude get_prs_fv3_run.html
!!
   subroutine get_prs_fv3_run(ix, levs, con_fvirt, phii, prsi, tgrs, qgrs1, del, del_gz, errmsg, errflg)

     implicit none

     ! Interface variables
     integer, intent(in) :: ix, levs
     real(kind=kind_phys), intent(in) :: con_fvirt
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: phii
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: prsi
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: tgrs
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: qgrs1
     real(kind=kind_phys), dimension(:,:),     intent(out)   :: del
     real(kind=kind_phys), dimension(:,:),     intent(out)   :: del_gz
     character(len=*),                         intent(out)   :: errmsg
     integer,                                  intent(out)   :: errflg

     ! Local variables
     integer :: i, k

     ! Initialize CCPP error handling variables
     errmsg = ''
     errflg = 0

! SJL: Adjust the geopotential height hydrostatically in a way consistent with FV3 discretization
! del_gz is a temp array recording the old info before (t,q) are adjusted
     do k=1,levs
       do i=1,ix
            del(i,k) = prsi(i,k) - prsi(i,k+1)
         del_gz(i,k) = (phii(i,k+1) - phii(i,k)) /                    &
                        (tgrs(i,k)*(one + con_fvirt*max(zero,qgrs1(i,k))))
       enddo
     enddo

   end subroutine get_prs_fv3_run

end module get_prs_fv3
