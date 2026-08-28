!> \file GFS_photochemistry.F90
!!

module GFS_photochemistry
  use machine,        only: kind_phys
  use module_ozphys,  only: ty_ozphys
  use module_h2ophys, only: ty_h2ophys
  implicit none
contains

!> \section arg_table_GFS_photochemistry_init Argument Table
!! \htmlinclude GFS_photochemistry_init.html
!!
  subroutine GFS_photochemistry_init(oz_phys_2006, oz_phys_2015, h2o_phys, errmsg, errflg)
    logical, intent(in) :: &
         oz_phys_2015, & !< Do ozone photochemistry? (2015)
         oz_phys_2006, & !< Do ozone photochemistry? (2006)
         h2o_phys        !< Do stratospheric h2o photochemistry?
    character(len=*), intent(out) :: &
         errmsg          !< CCPP Error message.
    integer,  intent(out) :: &
         errflg          !< CCPP Error flag.

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! If no photchemical scheme is on, but SDF has this module, report an error?
    if ((.not. oz_phys_2006) .and. (.not. oz_phys_2015) .and. (.not. h2o_phys)) then
       write (errmsg,'(*(a))') 'Logic error: One of [oz_phys_2006, oz_phys_2015, or h2o_phys] must == .true. '
       errflg = 1
       return
    endif
    
    ! Only one ozone scheme can be on. Otherwise, return and report error.
    if (oz_phys_2006 .and. oz_phys_2015) then
       write (errmsg,'(*(a))') 'Logic error: Only one ozone scheme can be enabled at a time'
       errflg = 1
       return
    endif
    
  end subroutine GFS_photochemistry_init

!> \section arg_table_GFS_photochemistry_run Argument Table
!! \htmlinclude GFS_photochemistry_run.html
!!
! #########################################################################################
  subroutine GFS_photochemistry_run (dtp, ntqv, ntoz, im, levs, ozphys, oz_phys_2015, oz_phys_2006, con_1ovg,   &
       prsl, dp, ozpl, h2o_phys, h2ophys, h2opl, gq0, gt0, ten_q, ten_u, ten_v, ten_t, do3_dt_prd, do3_dt_ozmx, &
       do3_dt_temp, do3_dt_ohoz, dqv_dt_prd, dqv_dt_qvmx, errmsg, errflg)
    
    ! Inputs
    real(kind=kind_phys), intent(in) :: &
         dtp,          & ! Model timestep
         con_1ovg        ! Physical constant (1./gravity)
    integer, intent(in) :: &
         ntqv,           &! index for specific humidity in the tracer array
         ntoz,           &! index for ozone in the the tracer array
         im,             &! horizontal loop extent
         levs            ! vertical dimension
    real(kind=kind_phys), intent(in), dimension(:,:) :: &
         prsl,         & ! Air pressure (Pa)
         dp,           & ! Pressure thickness (Pa)
         gt0             ! Air temperature (K)
    real(kind=kind_phys), intent(in), dimension(:,:,:) :: &
         gq0             ! tracer concentration
    real(kind=kind_phys), intent(in), dimension(:,:,:) :: &
         ozpl,         & ! Ozone data for current model timestep.
         h2opl           ! h2o data for curent model timestep.
    logical,              intent(in) :: &
         oz_phys_2015, & ! Do ozone photochemistry? (2015)
         oz_phys_2006, & ! Do ozone photochemistry? (2006)
         h2o_phys        ! Do stratospheric h2o photochemistry?
    type(ty_ozphys), intent(in) :: &
         ozphys          ! DDT with ozone photochemistry scheme/data.
    type(ty_h2ophys), intent(in) :: &
         h2ophys         ! DDT with h2o photochemistry scheme/data.

    ! Outputs (optional)
    real(kind=kind_phys), intent(inout), dimension(:,:), optional :: &
         do3_dt_prd,   & ! Physics tendency: production and loss effect
         do3_dt_ozmx,  & ! Physics tendency: ozone mixing ratio effect
         do3_dt_temp,  & ! Physics tendency: temperature effect
         do3_dt_ohoz,  & ! Physics tendency: overhead ozone effect
         dqv_dt_prd,   & ! Physics tendency: Climatological net production effect
         dqv_dt_qvmx     ! Physics tendency: specific humidity effect

    ! Outputs
    real(kind=kind_phys), intent(out), dimension(:,:) :: &
         ten_u, ten_v, ten_t
    real(kind=kind_phys), intent(out), dimension(:,:,:) :: &
         ten_q             ! tendency of tracer concentration
    character(len=*), intent(out) :: &
         errmsg          ! CCPP Error message.
    integer,  intent(out) :: &
         errflg          ! CCPP Error flag.
    
    ! Locals
    integer :: i,k
    real(kind=kind_phys), dimension(im,levs) :: &
         init_oz0,  oz0,          & ! initial and updated local ozone concentration
         init_h2o0, h2o0            ! initial and updated local h2o concentration
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    ten_u(:,:) = 0.0_kind_phys
    ten_v(:,:) = 0.0_kind_phys
    ten_t(:,:) = 0.0_kind_phys
    ten_q(:,:,:) = 0.0_kind_phys
    if (oz_phys_2015) then
       init_oz0  = gq0(:,:,ntoz)
       oz0 = init_oz0
       call ozphys%run_o3prog_2015(con_1ovg, dtp, prsl, gt0, dp, ozpl, oz0, do3_dt_prd,    &
            do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz)
       do i=1, im
         do k=1, levs
           ten_q(i,k,ntoz) = (oz0(i,k) - init_oz0(i,k))/dtp
         end do
       end do
    endif
    if (oz_phys_2006) then
       init_oz0  = gq0(:,:,ntoz)
       oz0 = init_oz0
       call ozphys%run_o3prog_2006(con_1ovg, dtp, prsl, gt0, dp, ozpl, oz0, do3_dt_prd,    &
            do3_dt_ozmx, do3_dt_temp, do3_dt_ohoz)
       do i=1, im
         do k=1, levs
           ten_q(i,k,ntoz) = (oz0(i,k) - init_oz0(i,k))/dtp
         end do
       end do
    endif
    if (h2o_phys) then
       init_h2o0 = gq0(:,:,ntqv)
       h2o0 = init_h2o0
       call h2ophys%run(dtp, prsl, h2opl, h2o0, dqv_dt_prd, dqv_dt_qvmx)
       do i=1, im
         do k=1, levs
           ten_q(i,k,ntqv) = (h2o0(i,k) - init_h2o0(i,k))/dtp
         end do
       end do
    endif

  end subroutine GFS_photochemistry_run

end module GFS_photochemistry
