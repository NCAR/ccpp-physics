
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
!! @brief Auto-generated cap module for the CCPP time_vary group
!!
!
module ccpp_group_time_vary_cap

   use GFS_time_vary_pre, only: GFS_time_vary_pre_init
   use GFS_rrtmg_setup, only: GFS_rrtmg_setup_init
   use GFS_phys_time_vary, only: GFS_phys_time_vary_init
   use stochastic_physics, only: stochastic_physics_init
   use stochastic_physics_sfc, only: stochastic_physics_sfc_init
   use GFS_time_vary_pre, only: GFS_time_vary_pre_run
   use GFS_rrtmg_setup, only: GFS_rrtmg_setup_run
   use GFS_rad_time_vary, only: GFS_rad_time_vary_run
   use GFS_phys_time_vary, only: GFS_phys_time_vary_run
   use stochastic_physics, only: stochastic_physics_run
   use GFS_time_vary_pre, only: GFS_time_vary_pre_finalize
   use GFS_rrtmg_setup, only: GFS_rrtmg_setup_finalize
   use GFS_phys_time_vary, only: GFS_phys_time_vary_finalize
   use stochastic_physics, only: stochastic_physics_finalize


   implicit none

   private
   public :: time_vary_init_cap,time_vary_run_cap,time_vary_finalize_cap

   logical, save :: initialized = .false.

   contains

   function time_vary_init_cap(GFS_Interstitial,cdata,GFS_Data,GFS_Control,CCPP_shared) result(ierr)

      use GFS_typedefs, only: GFS_interstitial_type
      use ccpp_types, only: ccpp_t
      use GFS_typedefs, only: GFS_data_type
      use GFS_typedefs, only: GFS_control_type
      use CCPP_typedefs, only: CCPP_shared_type

      implicit none

      integer                     :: ierr

      type(GFS_interstitial_type), intent(inout) :: GFS_Interstitial(:)
      type(ccpp_t), intent(inout) :: cdata
      type(GFS_data_type), intent(inout) :: GFS_Data(:)
      type(GFS_control_type), intent(inout) :: GFS_Control
      type(CCPP_shared_type), intent(in) :: CCPP_shared(:)

      ierr = 0


      if (initialized) return



      call GFS_time_vary_pre_init(errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_time_vary_pre_init"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_rrtmg_setup_init(si=GFS_Control%si,levr=GFS_Control%levr,ictm=GFS_Control%ictm,isol=GFS_Control%isol, &
                  ico2=GFS_Control%ico2,iaer=GFS_Control%iaer,ialb=GFS_Control%ialb,iems=GFS_Control%iems, &
                  ntcw=GFS_Control%ntcw,num_p2d=GFS_Control%num_p2d,num_p3d=GFS_Control%num_p3d, &
                  npdf3d=GFS_Control%npdf3d,ntoz=GFS_Control%ntoz,iovr_sw=GFS_Control%iovr_sw, &
                  iovr_lw=GFS_Control%iovr_lw,isubc_sw=GFS_Control%isubc_sw,isubc_lw=GFS_Control%isubc_lw, &
                  icliq_sw=GFS_Control%icliq_sw,crick_proof=GFS_Control%crick_proof,ccnorm=GFS_Control%ccnorm, &
                  imp_physics=GFS_Control%imp_physics,norad_precip=GFS_Control%norad_precip, &
                  idate=GFS_Control%idate,iflip=GFS_Control%iflip,im=GFS_Control%blksz(cdata%blk_no), &
                  faerlw=GFS_Interstitial(cdata%thrd_no)%faerlw,faersw=GFS_Interstitial(cdata%thrd_no)%faersw, &
                  aerodp=GFS_Interstitial(cdata%thrd_no)%aerodp,me=GFS_Control%me,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_rrtmg_setup_init"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_phys_time_vary_init(Data=GFS_Data,Model=GFS_Control,Interstitial=GFS_Interstitial,nthrds=CCPP_shared(cdata%thrd_no)%nthreads, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_phys_time_vary_init"
        ierr=cdata%errflg
        return
      end if

    
      call stochastic_physics_init(Model=GFS_Control,nthreads=CCPP_shared(cdata%thrd_no)%nthreads,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in stochastic_physics_init"
        ierr=cdata%errflg
        return
      end if

    
      call stochastic_physics_sfc_init(Model=GFS_Control,Data=GFS_Data,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in stochastic_physics_sfc_init"
        ierr=cdata%errflg
        return
      end if

    


      initialized = .true.


   end function time_vary_init_cap

   function time_vary_run_cap(cdata,GFS_Data,GFS_Control,CCPP_shared) result(ierr)

      use ccpp_types, only: ccpp_t
      use GFS_typedefs, only: GFS_data_type
      use GFS_typedefs, only: GFS_control_type
      use CCPP_typedefs, only: CCPP_shared_type

      implicit none

      integer                     :: ierr

      type(ccpp_t), intent(inout) :: cdata
      type(GFS_data_type), intent(inout) :: GFS_Data(:)
      type(GFS_control_type), intent(inout) :: GFS_Control
      type(CCPP_shared_type), intent(in) :: CCPP_shared(:)

      ierr = 0


      if (.not.initialized) then
        write(cdata%errmsg,'(*(a))') 'time_vary_run called before time_vary_init'
        cdata%errflg = 1
        return
      end if



      call GFS_time_vary_pre_run(Model=GFS_Control,errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_time_vary_pre_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_rrtmg_setup_run(idate=GFS_Control%idat,jdate=GFS_Control%jdat,deltsw=GFS_Control%fhswr, &
                  deltim=GFS_Control%dtf,lsswr=GFS_Control%lsswr,me=GFS_Control%me,slag=GFS_Control%slag, &
                  sdec=GFS_Control%sdec,cdec=GFS_Control%cdec,solcon=GFS_Control%solcon,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_rrtmg_setup_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_rad_time_vary_run(Model=GFS_Control,Data=GFS_Data,nthrds=CCPP_shared(cdata%thrd_no)%nthreads, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_rad_time_vary_run"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_phys_time_vary_run(Data=GFS_Data,Model=GFS_Control,nthrds=CCPP_shared(cdata%thrd_no)%nthreads, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_phys_time_vary_run"
        ierr=cdata%errflg
        return
      end if

    
      call stochastic_physics_run(Model=GFS_Control,Data=GFS_Data,nthreads=CCPP_shared(cdata%thrd_no)%nthreads, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in stochastic_physics_run"
        ierr=cdata%errflg
        return
      end if

    



   end function time_vary_run_cap

   function time_vary_finalize_cap(cdata) result(ierr)

      use ccpp_types, only: ccpp_t

      implicit none

      integer                     :: ierr

      type(ccpp_t), intent(inout) :: cdata

      ierr = 0


      if (.not.initialized) return



      call GFS_time_vary_pre_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_time_vary_pre_finalize"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_rrtmg_setup_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_rrtmg_setup_finalize"
        ierr=cdata%errflg
        return
      end if

    
      call GFS_phys_time_vary_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in GFS_phys_time_vary_finalize"
        ierr=cdata%errflg
        return
      end if

    
      call stochastic_physics_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
      if (cdata%errflg/=0) then
        write(cdata%errmsg,'(a)') "An error occured in stochastic_physics_finalize"
        ierr=cdata%errflg
        return
      end if

    


      initialized = .false.


   end function time_vary_finalize_cap

end module ccpp_group_time_vary_cap
