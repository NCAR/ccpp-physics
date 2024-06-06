!>  \file unified_ugwp_post.F90
!! This file saves CIRES UGWP diagnostics.
module unified_ugwp_post

contains

!>\defgroup unified_ugwp_post unified_UGWP Scheme Post
!> @{
!! The subroutine saves CIRES UGWP diagnostics.
!! \section arg_table_unified_ugwp_post_run Argument Table
!! \htmlinclude unified_ugwp_post_run.html
!!
     subroutine unified_ugwp_post_run (ldiag3d, ldiag_ugwp,         &
         dtf, im, levs,                                             &
         gw_dtdt, gw_dudt, gw_dvdt, tau_tofd, tau_mtb, tau_ogw,     &
         tau_ngw, zmtb, zlwb, zogw, dudt_mtb, dudt_ogw, dvdt_ogw,   &
         dudt_tms, tot_zmtb, tot_zlwb, tot_zogw,                    &
         tot_tofd, tot_mtb, tot_ogw, tot_ngw,                       &
         du3dt_mtb,du3dt_ogw, du3dt_tms, du3dt_ngw, dv3dt_ngw,      &
         ldu3dt_ogw, ldu3dt_obl, ldu3dt_oss, ldu3dt_ofd,            &
         dudt_ngw, dvdt_ngw, dtdt_ngw,                              &
         ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw, dudt_obl, dvdt_obl,    &
         dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd, dws3dt_ogw,        &
         dws3dt_obl, dws3dt_oss, dws3dt_ofd, du_ogwcol, dv_ogwcol,  &
         du_oblcol, dv_oblcol, du_osscol, dv_osscol, du_ofdcol,     &
         dv_ofdcol, du3_ogwcol, dv3_ogwcol, du3_oblcol, dv3_oblcol, &
         du3_osscol, dv3_osscol, du3_ofdcol, dv3_ofdcol,            & 
         dtdt, dudt, dvdt, errmsg, errflg)

        use machine,                only: kind_phys

        implicit none

        ! Interface variables
        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtf
        logical,              intent(in) :: ldiag_ugwp      !< flag for CIRES UGWP Diagnostics
        logical,              intent(in) :: ldiag3d

        real(kind=kind_phys), intent(in),    dimension(:)   :: zmtb, zlwb, zogw
        real(kind=kind_phys), intent(in),    dimension(:)   :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_mtb, tot_ogw, tot_tofd, tot_ngw
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_zmtb, tot_zlwb, tot_zogw
        real(kind=kind_phys), intent(in),    dimension(:,:) :: gw_dtdt, gw_dudt, gw_dvdt, dudt_mtb
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_ogw, dvdt_ogw
        real(kind=kind_phys), intent(in),    dimension(:,:) :: dudt_tms
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: du3dt_mtb, du3dt_ogw, du3dt_tms, du3dt_ngw, dv3dt_ngw
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: ldu3dt_ogw, ldu3dt_obl, ldu3dt_oss, ldu3dt_ofd
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_ngw, dvdt_ngw, dtdt_ngw
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_obl, dvdt_obl
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: dws3dt_obl, dws3dt_ogw
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: dws3dt_oss, dws3dt_ofd
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: du_ogwcol, dv_ogwcol
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: du_oblcol, dv_oblcol
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: du_osscol, dv_osscol
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: du_ofdcol, dv_ofdcol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_ogwcol, dv3_ogwcol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_oblcol, dv3_oblcol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_osscol, dv3_osscol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_ofdcol, dv3_ofdcol
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

          dws3dt_ogw = dws3dt_ogw + dtf *sqrt(dudt_ogw**2+dvdt_ogw**2)
          dws3dt_obl = dws3dt_obl + dtf *sqrt(dudt_obl**2+dvdt_obl**2)
          dws3dt_oss = dws3dt_oss + dtf *sqrt(dudt_oss**2+dvdt_oss**2)
          dws3dt_ofd = dws3dt_ofd + dtf *sqrt(dudt_ofd**2+dvdt_ofd**2)
          ldu3dt_ogw  = ldu3dt_ogw  + dtf *dudt_ogw
          ldu3dt_obl  = ldu3dt_obl  + dtf *dudt_obl
          ldu3dt_oss  = ldu3dt_oss  + dtf *dudt_oss
          ldu3dt_ofd  = ldu3dt_ofd  + dtf *dudt_ofd
          du3_ogwcol = du3_ogwcol + dtf *du_ogwcol
          dv3_ogwcol = dv3_ogwcol + dtf *dv_ogwcol
          du3_oblcol = du3_oblcol + dtf *du_oblcol
          dv3_oblcol = dv3_oblcol + dtf *dv_oblcol
          du3_osscol = du3_osscol + dtf *du_osscol
          dv3_osscol = dv3_osscol + dtf *dv_osscol
          du3_ofdcol = du3_ofdcol + dtf *du_ofdcol
          dv3_ofdcol = dv3_ofdcol + dtf *dv_ofdcol
          ! Special treatment for non-stationary GWD diagnostics
          ldu3dt_ngw = ldu3dt_ngw + dtf *dudt_ngw
          ldv3dt_ngw = ldv3dt_ngw + dtf *dvdt_ngw
          ldt3dt_ngw = ldt3dt_ngw + dtf *dtdt_ngw
        end if

        dtdt = dtdt + gw_dtdt
        dudt = dudt + gw_dudt
        dvdt = dvdt + gw_dvdt

      end subroutine unified_ugwp_post_run

!> @}
end module unified_ugwp_post
