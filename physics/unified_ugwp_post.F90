!>  \file unified_ugwp_post.F90
!! This file contains
module unified_ugwp_post

contains

!>\defgroup unified_ugwp_post unified_UGWP Scheme Post
!! @{

    subroutine unified_ugwp_post_init ()
    end subroutine unified_ugwp_post_init

!>@brief The subroutine initializes the unified UGWP

!> \section arg_table_unified_ugwp_post_run Argument Table
!! \htmlinclude unified_ugwp_post_run.html
!!
     subroutine unified_ugwp_post_run (ldiag_ugwp, dtf, im, levs,     &
         gw_dtdt, gw_dudt, gw_dvdt, tau_tofd, tau_mtb, tau_ogw,     &
         tau_ngw, zmtb, zlwb, zogw, dudt_mtb, dudt_ogw, dudt_tms,   &
         tot_zmtb, tot_zlwb, tot_zogw,                              &
         tot_tofd, tot_mtb, tot_ogw, tot_ngw,                       &
         du3dt_mtb,du3dt_ogw, du3dt_tms, du3dt_ngw, dv3dt_ngw,      &
         dtdt, dudt, dvdt, errmsg, errflg)

        use machine,                only: kind_phys

        implicit none

        ! Interface variables
        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtf
        logical,              intent(in) :: ldiag_ugwp      !< flag for CIRES UGWP Diagnostics

        real(kind=kind_phys), intent(in),    dimension(:)   :: zmtb, zlwb, zogw
        real(kind=kind_phys), intent(in),    dimension(:)   :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_mtb, tot_ogw, tot_tofd, tot_ngw
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_zmtb, tot_zlwb, tot_zogw
        real(kind=kind_phys), intent(in),    dimension(:,:) :: gw_dtdt, gw_dudt, gw_dvdt, dudt_mtb, dudt_ogw, dudt_tms
        real(kind=kind_phys), intent(inout), dimension(:,:) :: du3dt_mtb, du3dt_ogw, du3dt_tms, du3dt_ngw, dv3dt_ngw
        real(kind=kind_phys), intent(inout), dimension(:,:) :: dtdt, dudt, dvdt

        character(len=*),        intent(out) :: errmsg
        integer,                 intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (ldiag_ugwp) then
          tot_zmtb =  tot_zmtb + dtf *zmtb
          tot_zlwb =  tot_zlwb + dtf *zlwb
          tot_zogw =  tot_zogw + dtf *zogw
    
          tot_tofd  = tot_tofd + dtf *tau_tofd
          tot_mtb   = tot_mtb +  dtf *tau_mtb
          tot_ogw   = tot_ogw +  dtf *tau_ogw
          tot_ngw   = tot_ngw +  dtf *tau_ngw
    
          du3dt_mtb = du3dt_mtb + dtf *dudt_mtb
          du3dt_tms = du3dt_tms + dtf *dudt_tms
          du3dt_ogw = du3dt_ogw + dtf *dudt_ogw
          du3dt_ngw = du3dt_ngw + dtf *gw_dudt
          dv3dt_ngw = dv3dt_ngw + dtf *gw_dvdt
        endif

        dtdt = dtdt + gw_dtdt
        dudt = dudt + gw_dudt
        dvdt = dvdt + gw_dvdt

      end subroutine unified_ugwp_post_run

      subroutine unified_ugwp_post_finalize ()
      end subroutine unified_ugwp_post_finalize

!! @}
end module unified_ugwp_post
