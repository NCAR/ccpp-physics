
!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Auto-generated cap module for the CCPP suite
!!
!
module ccpp_suite_cap

   use ccpp_group_time_vary_cap, only: time_vary_init_cap,time_vary_run_cap,time_vary_finalize_cap
   use ccpp_group_radiation_cap, only: radiation_init_cap,radiation_run_cap,radiation_finalize_cap
   use ccpp_group_physics_cap, only: physics_init_cap,physics_run_cap,physics_finalize_cap
   use ccpp_group_stochastics_cap, only: stochastics_init_cap,stochastics_run_cap,stochastics_finalize_cap


   implicit none

   private
   public :: suite_init_cap,suite_run_cap,suite_finalize_cap

   contains

   function suite_init_cap(cdata,GFS_Interstitial,GFS_Control,GFS_Data,CCPP_shared) result(ierr)

      use ccpp_types, only: ccpp_t
      use GFS_typedefs, only: GFS_interstitial_type
      use GFS_typedefs, only: GFS_control_type
      use GFS_typedefs, only: GFS_data_type
      use CCPP_typedefs, only: CCPP_shared_type

      implicit none

      integer :: ierr
      type(ccpp_t), intent(inout) :: cdata
      type(GFS_interstitial_type), intent(inout) :: GFS_Interstitial(:)
      type(GFS_control_type), intent(inout) :: GFS_Control
      type(GFS_data_type), intent(inout) :: GFS_Data(:)
      type(CCPP_shared_type), intent(in) :: CCPP_shared(:)

      ierr = 0


      ierr = time_vary_init_cap(GFS_Interstitial=GFS_Interstitial,cdata=cdata,GFS_Data=GFS_Data,GFS_Control=GFS_Control, &
                  CCPP_shared=CCPP_shared)
      if (ierr/=0) return

      ierr = radiation_init_cap()
      if (ierr/=0) return

      ierr = physics_init_cap(cdata=cdata,GFS_Control=GFS_Control)
      if (ierr/=0) return

      ierr = stochastics_init_cap()
      if (ierr/=0) return


   end function suite_init_cap

   function suite_run_cap(con_t0c,GFS_Data,GFS_Control,con_epsm1,con_hvap,GFS_Interstitial,con_rd, &
                  con_rv,con_g,con_cp,con_cvap,CCPP_shared,LTP,con_eps,con_pi,cdata,con_fvirt, &
                  con_cliq) result(ierr)

      use machine, only: kind_phys
      use GFS_typedefs, only: GFS_data_type
      use GFS_typedefs, only: GFS_control_type
      use GFS_typedefs, only: GFS_interstitial_type
      use CCPP_typedefs, only: CCPP_shared_type
      use ccpp_types, only: ccpp_t

      implicit none

      integer :: ierr
      real(kind_phys), intent(in) :: con_t0c
      type(GFS_data_type), intent(inout) :: GFS_Data(:)
      type(GFS_control_type), intent(inout) :: GFS_Control
      real(kind_phys), intent(in) :: con_epsm1
      real(kind_phys), intent(in) :: con_hvap
      type(GFS_interstitial_type), intent(inout) :: GFS_Interstitial(:)
      real(kind_phys), intent(in) :: con_rd
      real(kind_phys), intent(in) :: con_rv
      real(kind_phys), intent(in) :: con_g
      real(kind_phys), intent(in) :: con_cp
      real(kind_phys), intent(in) :: con_cvap
      type(CCPP_shared_type), intent(in) :: CCPP_shared(:)
      integer, intent(in) :: LTP
      real(kind_phys), intent(in) :: con_eps
      real(kind_phys), intent(in) :: con_pi
      type(ccpp_t), intent(inout) :: cdata
      real(kind_phys), intent(in) :: con_fvirt
      real(kind_phys), intent(in) :: con_cliq

      ierr = 0


      ierr = time_vary_run_cap(cdata=cdata,GFS_Data=GFS_Data,GFS_Control=GFS_Control,CCPP_shared=CCPP_shared)
      if (ierr/=0) return

      ierr = radiation_run_cap(cdata=cdata,GFS_Interstitial=GFS_Interstitial,GFS_Data=GFS_Data,GFS_Control=GFS_Control, &
                  LTP=LTP)
      if (ierr/=0) return

      ierr = physics_run_cap(GFS_Control=GFS_Control,GFS_Data=GFS_Data,con_epsm1=con_epsm1,con_hvap=con_hvap, &
                  GFS_Interstitial=GFS_Interstitial,con_rd=con_rd,con_rv=con_rv,con_g=con_g, &
                  con_cp=con_cp,con_cvap=con_cvap,con_t0c=con_t0c,con_eps=con_eps,con_pi=con_pi, &
                  cdata=cdata,con_fvirt=con_fvirt,con_cliq=con_cliq)
      if (ierr/=0) return

      ierr = stochastics_run_cap(cdata=cdata,GFS_Data=GFS_Data,GFS_Control=GFS_Control)
      if (ierr/=0) return


   end function suite_run_cap

   function suite_finalize_cap(cdata) result(ierr)

      use ccpp_types, only: ccpp_t

      implicit none

      integer :: ierr
      type(ccpp_t), intent(inout) :: cdata

      ierr = 0


      ierr = time_vary_finalize_cap(cdata=cdata)
      if (ierr/=0) return

      ierr = radiation_finalize_cap()
      if (ierr/=0) return

      ierr = physics_finalize_cap(cdata=cdata)
      if (ierr/=0) return

      ierr = stochastics_finalize_cap()
      if (ierr/=0) return


   end function suite_finalize_cap

end module ccpp_suite_cap
