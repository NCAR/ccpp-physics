module tracer_sanitizer

  use machine, only : kind_phys

  implicit none

  private

  public :: tracer_sanitizer_init, tracer_sanitizer_run, tracer_sanitizer_finalize

  real(kind=kind_phys), parameter :: zero  = 0.0_kind_phys
  real(kind=kind_phys), parameter :: qvmin = 1.0E-6_kind_phys

contains

  subroutine tracer_sanitizer_init()
  end subroutine tracer_sanitizer_init

!> \section arg_table_tracer_sanitizer_run Argument Table
!! \htmlinclude tracer_sanitizer_run.html
!!
  subroutine tracer_sanitizer_run(tracers, ntqv, ntcw, ntiw, ntrw, ntsw, ntgl, &
                              ntlnc, ntinc, ntrnc, ntsnc, ntgnc, errmsg, errflg)

    ! Interface variables
    integer,              intent(in   ) :: ntqv, ntcw, ntiw, ntrw, ntsw, ntgl, &
                                           ntlnc, ntinc, ntrnc, ntsnc, ntgnc
    real(kind=kind_phys), intent(inout) :: tracers(:,:,:)
    character(len=*),     intent(  out) :: errmsg
    integer,              intent(  out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Water vapor specific humidity
    if (ntqv>0) then
      where (tracers(:,:,ntqv)<qvmin)
        tracers(:,:,ntqv)=qvmin
      end where
    end if

    ! Cloud water
    if (ntcw>0) then
      where (tracers(:,:,ntcw)<zero)
        tracers(:,:,ntcw)=zero
      end where
      ! Adjust second moments
      if (ntlnc>0) then
        where (tracers(:,:,ntlnc)==zero)
          tracers(:,:,ntlnc)=zero
        end where
      end if
    end if

    ! Ice water
    if (ntiw>0) then
      where (tracers(:,:,ntiw)<zero)
        tracers(:,:,ntiw)=zero
      end where
      ! Adjust second moments
      if (ntinc>0) then
        where (tracers(:,:,ntinc)==zero)
          tracers(:,:,ntinc)=zero
        end where
      end if
    end if

    ! Rain water
    if (ntrw>0) then
      where (tracers(:,:,ntrw)<zero)
        tracers(:,:,ntrw)=zero
      end where
      ! Adjust second moments
      if (ntrnc>0) then
        where (tracers(:,:,ntrnc)==zero)
          tracers(:,:,ntrnc)=zero
        end where
      end if
    end if

    ! Snow
    if (ntsw>0) then
      where (tracers(:,:,ntsw)<zero)
        tracers(:,:,ntsw)=zero
      end where
      ! Adjust second moments
      if (ntsnc>0) then
        where (tracers(:,:,ntsnc)==zero)
          tracers(:,:,ntsnc)=zero
        end where
      end if
    end if

    ! Graupel
    if (ntgl>0) then
      where (tracers(:,:,ntgl)<zero)
        tracers(:,:,ntgl)=zero
      end where
      ! Adjust second moments
      if (ntgnc>0) then
        where (tracers(:,:,ntgnc)==zero)
          tracers(:,:,ntgnc)=zero
        end where
      end if
    end if

  end subroutine tracer_sanitizer_run

  subroutine tracer_sanitizer_finalize()
  end subroutine tracer_sanitizer_finalize

end module tracer_sanitizer