module w3emc
  use machine, only: kind_sngl_prec, kind_dbl_prec
  implicit none
  public :: w3fs26

  interface w3difdat
     module procedure :: w3difdat32
     module procedure :: w3difdat64
  end interface w3difdat

  interface w3movdat
     module procedure :: w3movdat32
     module procedure :: w3movdat64
  end interface w3movdat

  interface w3reddat
     module procedure :: w3reddat32
     module procedure :: w3reddat64
  end interface w3reddat
  
contains  

  !  @brief Computes julian day number from year (4 digits), month, and day.
  !  @author Ralph Jones @date 1987-03-29
  !  Computes julian day number from year (4 digits), month,
  !  and day. iw3jdn is valid for years 1583 a.d. to 3300 a.d.
  !  Julian day number can be used to compute day of week, day of
  !  year, record numbers in an archive, replace day of century,
  !  find the number of days between two dates.
  !  @param[in] IYEAR Integer year (4 Digits)
  !  @param[in] MONTH Integer month of year (1 - 12)
  !  @param[in] IDAY Integer day of month (1 - 31)
  !  @return IW3JDN Integer Julian day number
  !  - Jan 1, 1960 is Julian day number 2436935
  !  - Jan 1, 1987 is Julian day number 2446797

  !  @note Julian period was devised by joseph scaliger in 1582.
  !  Julian day number #1 started on Jan. 1,4713 B.C. Three major
  !  chronological cycles begin on the same day. A 28-year solar
  !  cycle, a 19-year luner cycle, a 15-year indiction cycle, used
  !  in ancient rome to regulate taxes. It will take 7980 years
  !  to complete the period, the product of 28, 19, and 15.
  !  scaliger named the period, date, and number after his father
  !  Julius (not after the julian calendar). This seems to have
  !  caused a lot of confusion in text books. Scaliger name is
  !  spelled three different ways. Julian date and Julian day
  !  number are interchanged. A Julian date is used by astronomers
  !  to compute accurate time, it has a fraction. When truncated to
  !  an integer it is called an Julian day number. This function
  !  was in a letter to the editor of the communications of the acm
  !  volume 11 / number 10 / october 1968. The Julian day number
  !  can be converted to a year, month, day, day of week, day of
  !  year by calling subroutine w3fs26.
  !  @author Ralph Jones @date 1987-03-29
  function iw3jdn(iyear, month, iday)
    integer, intent(in) :: iyear, month, iday
    integer :: iw3jdn
    iw3jdn  =    iday - 32075 &
                 + 1461 * (iyear + 4800 + (month - 14) / 12) / 4 &
                 + 367 * (month - 2 - (month -14) / 12 * 12) / 12 &
                 - 3 * ((iyear + 4900 + (month - 14) / 12) / 100) / 4
  end function iw3jdn

  
  
  !> @brief Return the real kind and integer kind used in w3 lib.
  !> @author Jun Wang @date 2011-06-24
  !> This subprogram returns the real kind and the integer kind that the w3 lib
  !> is compiled with.
  !> @param[out] KINDREAL Kind of real number in w3 lib
  !> @param[out] KINDINT Kind of integer number in w3 lib
  !>
  !> @author Jun Wang @date 2011-06-24
  subroutine w3kind(kindreal, kindint)
    implicit none
    integer,intent(out) :: kindreal,kindint
    !  get real kind from a real number
    kindreal=kind(1.0)
    kindint=kind(1)
  end subroutine w3kind


  !> @brief Returns the integer day of week, the day
  !> of year, and julian day given an NCEP absolute date and time.
  !> @author Mark Iredell @date 1998-01-05
  
  !> @param[in] IDAT Integer (8) NCEP absolute date and time
  !> (year, month, day, time zone, hour, minute, second, millisecond)
  !> @param[out] JDOW Integer day of week (1-7, where 1 is sunday)
  !> @param[out] JDOY Integer day of year (1-366, where 1 is january 1)
  !> @param[out] JDAY Integer julian day (day number from jan. 1,4713 b.c.)
  !>
  !> @author Mark Iredell @date 1998-01-05
  subroutine w3doxdat(idat, jdow, jdoy, jday)
    integer :: idat(8), jdow, jdoy, jday
    integer :: jy, jm, jd

    !  get julian day and then get day of week and day of year
    jday=iw3jdn(idat(1), idat(2), idat(3))
    call w3fs26(jday, jy, jm, jd, jdow, jdoy)
  end subroutine w3doxdat

  
  !> @brief Return a time interval between two dates.
  !> @author Mark Iredell @date 1998-01-05
  !> Returns the elapsed time interval from
  !> an NCEP absolute date and time given in the second argument until
  !> an NCEP absolute date and time given in the first argument.
  !> The output time interval is in one of seven canonical forms
  !> of the ncep relative time interval data structure.
  !> @param[in] JDAT Integer (8) ncep absolute date and time
  !> (year, month, day, time zone, hour, minute, second, millisecond)
  !> @param[in] IDAT Integer (8) ncep absolute date and time
  !> (year, month, day, time zone, hour, minute, second, millisecond)
  !> @param[in] IT Integer relative time interval format type
  !> (-1 for first reduced type (hours always positive),
  !> 0 for second reduced type (hours can be negative),
  !> 1 for days only, 2 for hours only, 3 for minutes only,
  !> 4 for seconds only, 5 for milliseconds only)
  !> @param[out] RINC Real (5) ncep relative time interval
  !> (days, hours, minutes, seconds, milliseconds)
  !> (time interval is positive if jdat is later than idat.)
  !>
  !> @author Mark Iredell @date 1998-01-05
  subroutine w3difdat32(jdat, idat, it, rinc)
    integer :: it
    integer :: jdat(8),idat(8)
    real(kind_sngl_prec) :: rinc(5)
    real(kind_sngl_prec) :: rinc1(5)

    !  difference the days and time and put into canonical form
    rinc1(1)=iw3jdn(jdat(1),jdat(2),jdat(3)) - &
             iw3jdn(idat(1),idat(2),idat(3))
    rinc1(2:5)=jdat(5:8)-idat(5:8)
    call w3reddat(it,rinc1,rinc)
  end subroutine w3difdat32

  subroutine w3difdat64(jdat, idat, it, rinc)
    integer :: it
    integer :: jdat(8),idat(8)
    real(kind_dbl_prec) :: rinc(5)
    real(kind_dbl_prec) :: rinc1(5)

    !  difference the days and time and put into canonical form
    rinc1(1)=iw3jdn(jdat(1),jdat(2),jdat(3)) - &
             iw3jdn(idat(1),idat(2),idat(3))
    rinc1(2:5)=jdat(5:8)-idat(5:8)
    call w3reddat(it,rinc1,rinc)
  end subroutine w3difdat64

  

  ! @brief Year, month, day from julian day number
  ! @author Ralph Jones @date 1987-03-29
  subroutine w3fs26(jldayn, iyear, month, iday, idaywk, idayyr)
    integer :: jldayn, iyear, month, iday, idaywk, idayyr
    integer :: i, j, l
    real :: n
    l      = jldayn + 68569
    n      = 4 * l / 146097
    l      = l - (146097 * n + 3) / 4
    i      = 4000 * (l + 1) / 1461001
    l      = l - 1461 * i / 4 + 31
    j      = 80 * l / 2447
    iday   = l - 2447 * j / 80
    l      = j / 11
    month  = j + 2 - 12 * l
    iyear  = 100 * (n - 49) + i + l
    idaywk = mod((jldayn + 1),7) + 1
    idayyr = jldayn - &
         (-31739 +1461 * (iyear+4799) / 4 - 3 * ((iyear+4899)/100)/4)
  end subroutine w3fs26

  !> @file
  !> @brief Reduce a time interval to a canonical form.
  !> @author Mark Iredell @date 1998-01-05
  !> ### Program History Log:
  !> Date | Programmer | Comment
  !> -----|------------|--------
  !> 1998-01-05 | Mark Iredell | Initial.  
  !>
  !> @param[in] IT Relative time interval format type
  !> - (-1 for first reduced type (hours always positive),
  !> - 0 for second reduced type (hours can be negative),
  !> - 1 for days only, 2 for hours only, 3 for minutes only,
  !> - 4 for seconds only, 5 for milliseconds only)
  !> @param[in] RINC NCEP relative time interval (days, hours, minutes, seconds,
  !> milliseconds)
  !> @param[out] DINC NCEP relative time interval (days, hours, minutes,
  !> seconds, milliseconds)
  !>
  !> @author Mark Iredell @date 1998-01-05
  subroutine w3reddat32(it, rinc, dinc)
    integer :: it
    real(kind_sngl_prec) :: rinc(5), dinc(5)
    !  parameters for number of units in a day
    !  and number of milliseconds in a unit
    !  and number of next smaller units in a unit, respectively
    integer, dimension(5), parameter :: itd=(/1,24,1440,86400,86400000/), &
         itm=itd(5)/itd
    integer, dimension(4), parameter :: itn=itd(2:5)/itd(1:4)
    integer, parameter :: np=16
    integer :: iinc(4), jinc(5), kinc(5)
    integer :: ms

    !  first reduce to the first reduced form
    iinc=floor(rinc(1:4))
    !  convert all positive fractional parts to milliseconds
    !  and determine canonical milliseconds
    jinc(5)=nint(dot_product(rinc(1:4)-iinc,real(itm(1:4)))+rinc(5))
    kinc(5)=modulo(jinc(5),itn(4))
    !  convert remainder to seconds and determine canonical seconds
    jinc(4)=iinc(4)+(jinc(5)-kinc(5))/itn(4)
    kinc(4)=modulo(jinc(4),itn(3))
    !  convert remainder to minutes and determine canonical minutes
    jinc(3)=iinc(3)+(jinc(4)-kinc(4))/itn(3)
    kinc(3)=modulo(jinc(3),itn(2))
    !  convert remainder to hours and determine canonical hours
    jinc(2)=iinc(2)+(jinc(3)-kinc(3))/itn(2)
    kinc(2)=modulo(jinc(2),itn(1))
    !  convert remainder to days and compute milliseconds of the day
    kinc(1)=iinc(1)+(jinc(2)-kinc(2))/itn(1)
    ms=dot_product(kinc(2:5),itm(2:5))

    !  next reduce to either single value canonical form
    !  or to one of the two reduced forms
    if(it.ge.1.and.it.le.5) then
       !  ensure that exact multiples of 1./np       dinc(it)=real(kinc(1))*itd(it)+rp/np
    else
       !  the reduced form is done except the second reduced form is modified
       !  for negative time intervals with fractional days
       dinc=kinc
       if(it.eq.0.and.kinc(1).lt.0.and.ms.gt.0) then
          dinc(1)=dinc(1)+1
          dinc(2:5)=mod(ms-itm(1),itm(1:4))/itm(2:5)
       endif
    endif
  end subroutine w3reddat32
  
  subroutine w3reddat64(it, rinc, dinc)
    integer :: it
    real(kind_dbl_prec) :: rinc(5), dinc(5)
    !  parameters for number of units in a day
    !  and number of milliseconds in a unit
    !  and number of next smaller units in a unit, respectively
    integer, dimension(5), parameter :: itd=(/1,24,1440,86400,86400000/), &
         itm=itd(5)/itd
    integer, dimension(4), parameter :: itn=itd(2:5)/itd(1:4)
    integer, parameter :: np=16
    integer :: iinc(4), jinc(5), kinc(5)
    integer :: ms

    !  first reduce to the first reduced form
    iinc=floor(rinc(1:4))
    !  convert all positive fractional parts to milliseconds
    !  and determine canonical milliseconds
    jinc(5)=nint(dot_product(rinc(1:4)-iinc,real(itm(1:4)))+rinc(5))
    kinc(5)=modulo(jinc(5),itn(4))
    !  convert remainder to seconds and determine canonical seconds
    jinc(4)=iinc(4)+(jinc(5)-kinc(5))/itn(4)
    kinc(4)=modulo(jinc(4),itn(3))
    !  convert remainder to minutes and determine canonical minutes
    jinc(3)=iinc(3)+(jinc(4)-kinc(4))/itn(3)
    kinc(3)=modulo(jinc(3),itn(2))
    !  convert remainder to hours and determine canonical hours
    jinc(2)=iinc(2)+(jinc(3)-kinc(3))/itn(2)
    kinc(2)=modulo(jinc(2),itn(1))
    !  convert remainder to days and compute milliseconds of the day
    kinc(1)=iinc(1)+(jinc(2)-kinc(2))/itn(1)
    ms=dot_product(kinc(2:5),itm(2:5))

    !  next reduce to either single value canonical form
    !  or to one of the two reduced forms
    if(it.ge.1.and.it.le.5) then
       !  ensure that exact multiples of 1./np       dinc(it)=real(kinc(1))*itd(it)+rp/np
    else
       !  the reduced form is done except the second reduced form is modified
       !  for negative time intervals with fractional days
       dinc=kinc
       if(it.eq.0.and.kinc(1).lt.0.and.ms.gt.0) then
          dinc(1)=dinc(1)+1
          dinc(2:5)=mod(ms-itm(1),itm(1:4))/itm(2:5)
       endif
    endif
  end subroutine w3reddat64

  
  !> @file
  !> @brief Return a date from a time interval and date
  !> @author Mark Iredell @date 1998-08-01
  !> This subprogram returns the date and time that is a given
  !> NCEP relative time interval from an NCEP absolute date and time.
  !> The output is in the NCEP absolute date and time data structure.
  !>
  !> ### Program History Log:
  !> Date | Programmer | Comment
  !> -----|------------|--------
  !> 1998-01-05 | Mark Iredell | Initial.
  !>
  !> @param[in] RINC NCEP relative time interval (days, hours, minutes, seconds
  !> milliseconds)
  !> @param[in] IDAT NCEP absolute date and time (year, month, day, time zone,
  !> hour, minute, second, millisecond)
  !> @param[in] JDAT NCEP absolute date and time (year, month, day, time zone,
  !> hour, minute, second, millisecond) (jdat is later than idat if time
  !> interval is positive.)
  !>
  !> @author Mark Iredell @date 1998-08-01
  subroutine w3movdat32(rinc, idat, jdat)
    real(kind_sngl_prec) :: rinc(5)
    integer :: idat(8), jdat(8), jdow, jdoy, jldayn
    real(kind_sngl_prec) :: rinc1(5), rinc2(5)

    !  add the interval to the input time of day and put into reduced form
    !  and then compute new date using julian day arithmetic.
    rinc1(1)=rinc(1)
    rinc1(2:5)=rinc(2:5)+idat(5:8)
    call w3reddat(-1,rinc1,rinc2)
    jldayn=iw3jdn(idat(1),idat(2),idat(3))+nint(rinc2(1))
    call w3fs26(jldayn,jdat(1),jdat(2),jdat(3),jdow,jdoy)
    jdat(4)=idat(4)
    jdat(5:8)=nint(rinc2(2:5))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine w3movdat32

  subroutine w3movdat64(rinc, idat, jdat)
    real(kind_dbl_prec) :: rinc(5)
    integer :: idat(8), jdat(8), jdow, jdoy, jldayn
    real(kind_dbl_prec) :: rinc1(5), rinc2(5)

    !  add the interval to the input time of day and put into reduced form
    !  and then compute new date using julian day arithmetic.
    rinc1(1)=rinc(1)
    rinc1(2:5)=rinc(2:5)+idat(5:8)
    call w3reddat(-1,rinc1,rinc2)
    jldayn=iw3jdn(idat(1),idat(2),idat(3))+nint(rinc2(1))
    call w3fs26(jldayn,jdat(1),jdat(2),jdat(3),jdow,jdoy)
    jdat(4)=idat(4)
    jdat(5:8)=nint(rinc2(2:5))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine w3movdat64
  
end module w3emc
