#define CCPP
!>  \file noahmpdrv_time_vary.F90
!!  This file contains the IAU Updates for the NoahMP land surface scheme driver.

!>\defgroup NoahMP_LSM NoahMP LSM Model
!! \brief This is the NoahMP LSM the IAU Updates module

!> This module contains the CCPP-compliant IAU Update module for NoahMP land surface model driver.
!> The noahmpdrv_time_vary module is an alternative to calling the IAU updates directly from within the noahmpdrv module 
!> The current "CCPP_driver" module's ccpp_step(step="timestep_init") function call only handles group="time_vary" and not "physics"
!
module noahmpdrv_time_vary

    !   use module_sf_noahmplsm
    ! 3.5.24 for use in IAU
    use lnd_iau_mod,  only: lnd_iau_control_type, lnd_iau_external_data_type,&
                            lnd_iau_mod_set_control, lnd_iau_mod_init,lnd_iau_mod_getiauforcing 

    implicit none

    private

    public :: noahmpdrv_time_vary_init, noahmpdrv_time_vary_timestep_init  !, noahmpdrv_time_vary_run
!   public :: noahmpdrv_time_vary_timestep_finalize, noahmpdrv_time_vary_finalize

    ! IAU data and control
    type (lnd_iau_control_type)                  :: LND_IAU_Control
    type (lnd_iau_external_data_type)            :: LND_IAU_Data      !(number of blocks):each proc holds nblks

    contains

!> \ingroup NoahMP_LSM
!! \brief This subroutine is called during the CCPP initialization phase to  
!! initialize Land IAU Control and Land IAU Data structures.
!! \section arg_table_noahmpdrv_time_vary_init Argument Table
!! \htmlinclude noahmpdrv_time_vary_init.html
!!
  subroutine noahmpdrv_time_vary_init(lsm, lsm_noahmp, me, mpi_root,                        &
                                fn_nml, input_nml_file, isc, jsc, ncols, nx, ny, nblks,     &
                                blksz, xlon, xlat,                                          &     
                                lsoil, lsnow_lsm, dtp, fhour, errmsg, errflg)

    use machine,          only: kind_phys
    !use GFS_typedefs, only: GFS_control_type
    ! use GFS_typedefs, only: GFS_data_type

    implicit none

    integer,              intent(in) :: lsm
    integer,              intent(in) :: lsm_noahmp    
    integer,              intent(in) :: me         !  mpi_rank       
    integer,                       intent(in) :: mpi_root   ! = GFS_Control%master        
    character(*),                  intent(in) :: fn_nml
    character(len=:), intent(in), dimension(:), pointer :: input_nml_file 
    integer,                       intent(in) :: isc, jsc, ncols, nx, ny, nblks      !=GFS_Control%ncols, %nx, %ny, nblks
    integer, dimension(:),         intent(in) :: blksz   !(one:) !GFS_Control%blksz
    real(kind_phys), dimension(:), intent(in) :: xlon    ! longitude !GFS_Data(cdata%blk_no)%Grid%xlon
    real(kind_phys), dimension(:), intent(in) :: xlat    ! latitude
    integer,                       intent(in) :: lsoil, lsnow_lsm
    real(kind=kind_phys),          intent(in) :: dtp, fhour
    character(len=*),     intent(out) :: errmsg
    integer,              intent(out) :: errflg

    ! type(gfs_data_type), dimension(:), intent(inout)          :: GFS_Data ! !(one:)
    ! type(gfs_control_type), intent(in)          :: GFS_Control

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! 3.7.24 init iau for land
    call lnd_iau_mod_set_control(LND_IAU_Control, fn_nml, input_nml_file, me, mpi_root, isc,jsc, nx, ny, nblks, blksz,  &
                                    lsoil, lsnow_lsm, dtp, fhour, errmsg, errflg)
!        print*, 'proc errmsg, errflg after set control', me, errmsg, errflg
!        print*, 'proc iau_control isc, nx, dtp fhour', me, LND_IAU_Control%isc, LND_IAU_Control%nx, &
!                LND_IAU_Control%dtp, LND_IAU_Control%fhour
!        print*, 'proc iau_control incfiles(1)', me, LND_IAU_Control%iau_inc_files(1)
!        stop
    call lnd_iau_mod_init (LND_IAU_Control, LND_IAU_Data, xlon, xlat, errmsg, errflg)
    !print*, 'proc errmsg, errflg interval after lnd_iau_init ', me,trim(errmsg), errflg, LND_IAU_Data%in_interval
    print*, 'proc nblks blksize(1) after lnd_iau_mod_init ', me,LND_IAU_Control%nblks, LND_IAU_Control%blksz(1)  

  end subroutine noahmpdrv_time_vary_init

!> \ingroup NoahMP_LSM
!! \brief This subroutine is called before noahmpdrv_run timestep to update
!!  states with iau increments
!! \section arg_table_noahmpdrv_time_vary_timestep_init Argument Table
!! \htmlinclude noahmpdrv_time_vary_timestep_init.html
!!
  subroutine noahmpdrv_time_vary_timestep_init (itime, fhour, delt, km,  &      !me, mpi_root,
                                      stc, slc, errmsg, errflg)       ! smc, t2mmp, q2mp,

    use machine,          only: kind_phys

    implicit none

    ! integer,                                  intent(in) :: me         !mpi_rank
    ! integer,                                  intent(in) :: mpi_root   ! = GFS_Control%master
    integer                                   , intent(in) :: itime      !current forecast iteration
    real(kind=kind_phys)                      , intent(in) :: fhour      !current forecast time (hr)
    real(kind=kind_phys)                      , intent(in) :: delt       ! time interval [s]
    integer                                   , intent(in) :: km         !vertical soil layer dimension
    real(kind=kind_phys), dimension(:,:)   , intent(inout) :: stc        ! soiltemp [K]
    real(kind=kind_phys), dimension(:,:)   , intent(inout) :: slc        !liquid soil moisture [m3/m3]'
    character(len=*),                          intent(out) :: errmsg
    integer,                                   intent(out) :: errflg

    !  ---  local variable
    ! integer :: nb, im        ! vertical soil layer dimension

    ! IAU update
    real,allocatable :: stc_inc_flat(:,:)
    real,allocatable :: slc_inc_flat(:,:)
    ! real,allocatable :: tmp2m_inc_flat(:)
    ! real,allocatable :: spfh2m_inc_flat(:)
    integer :: j, k, ib
  !  --- end declaration

  !  --- Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    !> update current forecast hour
    ! GFS_control%jdat(:) = jdat(:)
    LND_IAU_Control%fhour=fhour

    if(LND_IAU_Control%me == LND_IAU_Control%mpi_root) then
      print*,"itime ",itime," GFScont%fhour ",fhour," IauCon%fhour",LND_IAU_Control%fhour,   &
             " delt ",delt," IauCont%dtp",LND_IAU_Control%dtp
    endif

    !> 3.7.24 read iau increments
    call lnd_iau_mod_getiauforcing(LND_IAU_Control, LND_IAU_Data, errmsg, errflg)   !call getiauforcing(GFS_control,IAU_data)
    if (errflg .ne. 0) then
      if(LND_IAU_Control%me == LND_IAU_Control%mpi_root) then
        print*, "noahmpdrv_timestep_init: lnd_iau_mod_getiauforcing returned nonzero value"
        print*, errmsg
      endif
      return
    endif

    !> update with iau increments
    if (LND_IAU_Data%in_interval) then
        if(LND_IAU_Control%me == LND_IAU_Control%mpi_root) then
          print*, "adding land iau increments "
        endif

      if (LND_IAU_Control%lsoil .ne. km) then
        write(errmsg,*) 'noahmpdrv_timestep_init: LND_IAU_Data%lsoil ',LND_IAU_Control%lsoil,' not equal to km ',km
        errflg = 1
        return
      endif

      ! local variable to copy blocked data LND_IAU_Data%stc_inc
      allocate(stc_inc_flat(LND_IAU_Control%nx * LND_IAU_Control%ny, km))  !GFS_Control%ncols
      allocate(slc_inc_flat(LND_IAU_Control%nx * LND_IAU_Control%ny, km))  !GFS_Control%ncols
      ! allocate(tmp2m_inc_flat(LND_IAU_Control%nx * LND_IAU_Control%ny))  !GFS_Control%ncols
      ! allocate(spfh2m_inc_flat(LND_IAU_Control%nx * LND_IAU_Control%ny))  !GFS_Control%ncols
      ib = 1
      do j = 1, LND_IAU_Control%ny  !ny
        do k = 1, km
          stc_inc_flat(ib:ib+LND_IAU_Control%nx-1, k) =LND_IAU_Data%stc_inc(:,j, k)
          slc_inc_flat(ib:ib+LND_IAU_Control%nx-1, k) = LND_IAU_Data%slc_inc(:,j, k)
        enddo
      ! ib = 1
      ! do j = 1, LND_IAU_Control%ny  !ny
        ! tmp2m_inc_flat(ib:ib+LND_IAU_Control%nx-1) =LND_IAU_Data%tmp2m_inc(:,j, 1)
        ! spfh2m_inc_flat(ib:ib+LND_IAU_Control%nx-1)=LND_IAU_Data%spfh2m_inc(:,j, 1)

        ib = ib + LND_IAU_Control%nx  !nlon
      enddo

      ! delt=GFS_Control%dtf
      if ((LND_IAU_Control%dtp - delt) > 0.0001) then
        if(LND_IAU_Control%me == LND_IAU_Control%mpi_root) then
          print*, "Warning noahmpdrv_timevary_tstep delt ",delt,"different from LND_IAU_Control%dtp ",LND_IAU_Control%dtp
        endif
      endif
      !IAU increments are in units of 1/sec     !LND_IAU_Control%dtp
      do k = 1, km
        stc(:,k) = stc(:,k) + stc_inc_flat(:,k)*delt !LND_IAU_Control%dtp
        slc(:,k) = slc(:,k) + slc_inc_flat(:,k)*delt !LND_IAU_Control%dtp
      enddo
      ! t2mmp = t2mmp +  &
      ! tmp2m_inc_flat(LND_IAU_Control%blk_strt_indx(nb):LND_IAU_Control%blk_strt_indx(nb) + im-1)*delt !LND_IAU_Control%dtp
      ! q2mp = q2mp +   &
      ! spfh2m_inc_flat(LND_IAU_Control%blk_strt_indx(nb):LND_IAU_Control%blk_strt_indx(nb)+ im-1)*delt !LND_IAU_Control%dtp

      deallocate(stc_inc_flat, slc_inc_flat)   !, tmp2m_inc_flat,spfh2m_inc_flat)

    endif

  end subroutine noahmpdrv_time_vary_timestep_init


! !> \ingroup NoahMP_LSM
! !! \brief 
! !! \section arg_table_noahmpdrv_time_vary_run Argument Table
! !! \htmlinclude noahmpdrv_time_vary_run.html
! !!
! !! \section general_noahmpdrv_time_vary_run 
! !!  @{
! !!    - Initialize CCPP error handling variables.

!   subroutine noahmpdrv_time_vary_run(nb, im, km, lsnowl, itime, fhour, errmsg, errflg)                               
! ! !  ---  inputs:
! ! !  ---  in/outs:
! !       weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
! ! ! --- Noah MP specific
! ! !  ---  outputs:
! !       )

!   use machine ,   only : kind_phys

!   implicit none
      
! !
! !  ---  CCPP interface fields (in call order)
! !
!   integer                                , intent(in)    :: nb         !=cdata%blk_no, 
!   integer                                , intent(in)    :: im         ! horiz dimension and num of used pts
!   integer                                , intent(in)    :: km         ! vertical soil layer dimension
!   integer                                , intent(in)    :: lsnowl     ! lower bound for snow level arrays
!   integer                                , intent(in)    :: itime      ! NOT USED current forecast iteration
!   real(kind=kind_phys)                   , intent(in)    :: fhour      ! currentforecast time (hr)

! !   real(kind=kind_phys), dimension(:)     , intent(inout) :: weasd      ! water equivalent accumulated snow depth [mm]
! !   real(kind=kind_phys), dimension(:)     , intent(inout) :: snwdph     ! snow depth [mm]
! !   real(kind=kind_phys), dimension(:)     , intent(inout) :: tskin      ! ground surface skin temperature [K]
! !   real(kind=kind_phys), dimension(:)     , intent(inout) :: tprcp      ! total precipitation [m]
! !   real(kind=kind_phys), dimension(:)     , intent(inout) :: srflag     ! snow/rain flag for precipitation
! !   real(kind=kind_phys), dimension(:,:)   , intent(inout) :: smc        ! total soil moisture content [m3/m3]
! !   real(kind=kind_phys), dimension(:,:)   , intent(inout) :: stc        ! soil temp [K]
! !   real(kind=kind_phys), dimension(:,:)   , intent(inout) :: slc        ! liquid soil moisture [m3/m3]
! !   real(kind=kind_phys), dimension(:)     , intent(inout) :: canopy     ! canopy moisture content [mm]
! !   real(kind=kind_phys), dimension(:)     , intent(inout) :: trans      ! total plant transpiration [m/s]
! !   real(kind=kind_phys), dimension(:)     , intent(inout) :: tsurf      ! surface skin temperature [K]
! !   real(kind=kind_phys), dimension(:)     , intent(inout) :: zorl       ! surface roughness [cm]

!   character(len=*)    ,                    intent(out)   :: errmsg
!   integer             ,                    intent(out)   :: errflg
! !
! !  --- end declaration
! !     

! !
! !  --- Initialize CCPP error handling variables
! !
!   errmsg = ''
!   errflg = 0

! !   if(LND_IAU_Control%me == LND_IAU_Control%mpi_root) then 
! !     print*,"nb ",nb," itime ",itime," GFScont%fhour ",fhour," iauCon%fhour",LND_IAU_Control%fhour," delt ",delt," iauCont%dtp",LND_IAU_Control%dtp
! !   endif
! !   ! 3.7.24 read iau increments 
! !   call lnd_iau_mod_getiauforcing(LND_IAU_Control, LND_IAU_Data, errmsg, errflg) 
! !   if (errflg .ne. 0) return
! !   ! update with iau increments
! !   if (LND_IAU_Data%in_interval) then
! !     if (LND_IAU_Control%lsoil .ne. km) then
! !       write(errmsg, *)'in noahmpdrv_run, lnd_iau_mod update increments:LND_IAU_Control%lsoil ',LND_IAU_Control%lsoil,' not equal to km ',km
! !       errflg = 1
! !       return
! !     endif
! !     ! LND_IAU_Data%stc_inc(is:ie, js:je, km))  size of (nx, ny)
! !     ! xlatin(im) stc(im,km) slc() t2mmp(:) q2mp(im) km=n_soill, im =
! !     ! GFS_Control%blksz(cdata%blk_no)
! !     ! >> need to get (cdata%blk_no from function call 

! !     ! local variable to copy blocked data LND_IAU_Data%stc_inc
! !     allocate(stc_inc_flat(LND_IAU_Control%nx * LND_IAU_Control%ny, km))  !GFS_Control%ncols
! !     allocate(slc_inc_flat(LND_IAU_Control%nx * LND_IAU_Control%ny, km))  !GFS_Control%ncols
! !     allocate(tmp2m_inc_flat(LND_IAU_Control%nx * LND_IAU_Control%ny))  !GFS_Control%ncols
! !     allocate(spfh2m_inc_flat(LND_IAU_Control%nx * LND_IAU_Control%ny))  !GFS_Control%ncols
! !     ib = 1
! !     do j = 1, LND_IAU_Control%ny  !ny 
! !       do k = 1, km    
! !         stc_inc_flat(ib:ib+LND_IAU_Control%nx-1, k) = LND_IAU_Data%stc_inc(:,j,k)  
! !         slc_inc_flat(ib:ib+LND_IAU_Control%nx-1, k) = LND_IAU_Data%slc_inc(:,j,k) 
! !       enddo
! !     ! ib = 1
! !     ! do j = 1, LND_IAU_Control%ny  !ny
! !       tmp2m_inc_flat(ib:ib+LND_IAU_Control%nx-1) = LND_IAU_Data%tmp2m_inc(:,j,1)  
! !       spfh2m_inc_flat(ib:ib+LND_IAU_Control%nx-1) = LND_IAU_Data%spfh2m_inc(:,j,1) 
! !       ib = ib + LND_IAU_Control%nx  !nlon    
! !     enddo

! !     !IAU increments are in units of 1/sec     !LND_IAU_Control%dtp
! !     ! delt=GFS_Control%dtf
! !     if ((LND_IAU_Control%dtp - delt) > 0.0001) then 
! !       if(LND_IAU_Control%me == LND_IAU_Control%mpi_root) then 
! !         print*, "Warning time step used in noahmpdrv_run delt ",delt," different from LND_IAU_Control%dtp ",LND_IAU_Control%dtp
! !       endif
! !     endif
! !     do k = 1, km 
! !     stc(:,k)=stc(:,k)+stc_inc_flat(LND_IAU_Control%blk_strt_indx(nb):LND_IAU_Control%blk_strt_indx(nb)+im-1, k)*delt !LND_IAU_Control%dtp
! !     slc(:,k)=slc(:,k)+slc_inc_flat(LND_IAU_Control%blk_strt_indx(nb):LND_IAU_Control%blk_strt_indx(nb)+im-1, k)*delt !LND_IAU_Control%dtp    
! !     enddo
! !     t2mmp = t2mmp+tmp2m_inc_flat(LND_IAU_Control%blk_strt_indx(nb):LND_IAU_Control%blk_strt_indx(nb)+im-1)*delt !LND_IAU_Control%dtp   
! !     q2mp = q2mp +spfh2m_inc_flat(LND_IAU_Control%blk_strt_indx(nb):LND_IAU_Control%blk_strt_indx(nb)+im-1)*delt !LND_IAU_Control%dtp 

! !     deallocate(stc_inc_flat, slc_inc_flat, tmp2m_inc_flat, spfh2m_inc_flat)
  
! !   end if  
!   end subroutine noahmpdrv_time_vary_run

!   subroutine noahmpdrv_time_vary_timestep_finalize (errmsg, errflg)       ! smc, t2mmp, q2mp,
                                  

!     use machine,          only: kind_phys

!     implicit none

!     character(len=*),                          intent(out) :: errmsg
!     integer,                                   intent(out) :: errflg

!   !  --- Initialize CCPP error handling variables
!     errmsg = ''
!     errflg = 0

!   end subroutine noahmpdrv_time_vary_timestep_finalize

!   subroutine noahmpdrv_time_vary_finalize (errmsg, errflg)       ! smc, t2mmp, q2mp,
                                  

!     use machine,          only: kind_phys

!     implicit none

!     character(len=*),                          intent(out) :: errmsg
!     integer,                                   intent(out) :: errflg

!   !  --- Initialize CCPP error handling variables
!     errmsg = ''
!     errflg = 0

!   end subroutine noahmpdrv_time_vary_finalize

end module noahmpdrv_time_vary
