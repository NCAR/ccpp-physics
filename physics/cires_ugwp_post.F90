!>  \file cires_ugwp_post.F90
!! This file contains
module cires_ugwp_post

contains

!>\defgroup cires_ugwp_post CIRES UGWP Scheme Post
!! @{
!> \section arg_table_cires_ugwp_post_init Argument Table
!!
    subroutine cires_ugwp_post_init ()
    end subroutine cires_ugwp_post_init

!>@brief The subroutine initializes the CIRES UGWP
#if 0
!> \section arg_table_cires_ugwp_post_run Argument Table
!! \htmlinclude cires_ugwp_post_run.html
!!
#endif


     subroutine cires_ugwp_post_run (ldiag_ugwp, dtf, im, levs,     &
         gw_dudt, tau_tofd, tau_mtb, tau_ogw, tau_ngw,              &
         zmtb, zlwb, zogw, dudt_mtb, dudt_ogw, dudt_tms,            &
         tot_zmtb, tot_zlwb, tot_zogw,                              &
         tot_tofd, tot_mtb, tot_ogw, tot_ngw,                       &
         du3dt_mtb,du3dt_ogw, du3dt_tms, du3dt_ngw,                 &
         cnvgwd, errmsg, errflg)

        use machine,                only: kind_phys

        implicit none

        ! Interface variables
        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtf
        logical,              intent(in) :: ldiag_ugwp      !< flag for CIRES UGWP Diagnostics
        logical,              intent(inout) :: cnvgwd       !< flag to turn on/off convective gwd

        real(kind=kind_phys), intent(in),    dimension(im)       :: zmtb, zlwb, zogw
        real(kind=kind_phys), intent(in),    dimension(im)       :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
        real(kind=kind_phys), intent(inout), dimension(im)       :: tot_mtb, tot_ogw, tot_tofd, tot_ngw
        real(kind=kind_phys), intent(inout), dimension(im)       :: tot_zmtb, tot_zlwb, tot_zogw
        real(kind=kind_phys), intent(in),    dimension(im, levs) :: gw_dudt, dudt_mtb, dudt_ogw, dudt_tms 
        real(kind=kind_phys), intent(inout), dimension(im, levs) :: du3dt_mtb, du3dt_ogw, du3dt_tms, du3dt_ngw

        character(len=*),        intent(out) :: errmsg
        integer,                 intent(out) :: errflg


        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (.not. (ldiag_ugwp)) return


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
         endif          


         cnvgwd = .false.

      end subroutine cires_ugwp_post_run

!> \section arg_table_cires_ugwp_post_finalize Argument Table
!!
      subroutine cires_ugwp_post_finalize ()
      end subroutine cires_ugwp_post_finalize

!! @}
end module cires_ugwp_post
