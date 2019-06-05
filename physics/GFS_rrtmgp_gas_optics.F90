module GFS_rrtmgp_gas_optics
  use machine,      only: kind_phys
  use GFS_typedefs, only: GFS_control_type,GFS_radtend_type
contains
  
!! \section arg_table_GFS_rrtmgp_gas_optics_init Argument Table
!! | local_name   | standard_name                     | long_name                                               | units    | rank |  type                 |   kind    | intent | optional |
!! |--------------|-----------------------------------|---------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model        | GFS_control_type_instance         | Fortran DDT containing FV3-GFS model control parameters | DDT      |    0 | GFS_control_type      |           | in     | F        |
!! | Radtend      | GFS_radtend_type_instance         | Fortran DDT containing FV3-GFS radiation tendencies     | DDT      |    0 | GFS_radtend_type      |           | inout  | F        |
!! | errmsg       | ccpp_error_message                | error message for error handling in CCPP                | none     |    0 | character             | len=*     | out    | F        |
!! | errflg       | ccpp_error_flag                   | error flag for error handling in CCPP                   | flag     |    0 | integer               |           | out    | F        |
!!
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_gas_optics_init()
  ! #########################################################################################
  subroutine GFS_rrtmgp_gas_optics_init(Model, Radtend, errmsg, errflg)
    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model      ! DDT containing model control parameters
   type(GFS_radtend_type), intent(inout) :: &
        Radtend     ! Fortran DDT containing FV3-GFS radiation tendencies 
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg     ! Error message
    integer, intent(out) :: &  
         errflg     ! Error flag

    ! Local variables
    character(len=1) :: tempstr
    integer :: ij, count
    integer,dimension(Model%ngases,2) :: gasIndices

    ! Initialize
    errmsg = ''
    errflg = 0


    ! Which gases are active? Provided via physics namelist.
    if (len(Model%active_gases) .gt. 0) then

       ! Pull out gas names from list...
       ! First grab indices in character array corresponding to start:end of gas name.
       gasIndices(1,1)=1
       count=1
       do ij=1,len(Model%active_gases)
          tempstr=trim(Model%active_gases(ij:ij))
          if (tempstr .eq. '_') then
             gasIndices(count,2)=ij-1
             gasIndices(count+1,1)=ij+1
             count=count+1
          endif
       enddo
       gasIndices(Model%ngases,2)=len(trim(Model%active_gases))
       ! Now extract the gas names
       do ij=1,Model%ngases
          Radtend%active_gases(ij) = Model%active_gases(gasIndices(ij,1):gasIndices(ij,2))
       enddo
    endif
  end subroutine GFS_rrtmgp_gas_optics_init
  !
  subroutine GFS_rrtmgp_gas_optics_run()
  end subroutine GFS_rrtmgp_gas_optics_run
  !
  subroutine GFS_rrtmgp_gas_optics_finalize()
  end subroutine GFS_rrtmgp_gas_optics_finalize
  !
end module GFS_rrtmgp_gas_optics
