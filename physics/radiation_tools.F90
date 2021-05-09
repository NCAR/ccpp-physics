module radiation_tools
  use machine, only: &
       kind_phys                   ! Working type
  implicit none
contains

  ! ######################################################################################### 
  ! ######################################################################################### 
  subroutine cmp_tlev(nCol,nLev,minP,minT,maxT,p_lay,t_lay,p_lev,tsfc,t_lev)
    ! Inputs
    integer, intent(in) :: &
         nCol,nLev
    real(kind_phys),intent(in) :: &
         minP,minT,maxT
    real(kind_phys),dimension(nCol),intent(in) :: &
         tsfc
    real(kind_phys),dimension(nCol,nLev),intent(in) :: &
         p_lay,t_lay
    real(kind_phys),dimension(nCol,nLev+1),intent(in) :: &
         p_lev

    ! Outputs
    real(kind_phys),dimension(nCol,nLev+1),intent(out) :: &
         t_lev

    ! Local
    integer :: iCol,iLay, iSFC, iTOA
    logical :: top_at_1
    real(kind_phys), dimension(nCol,nLev) :: tem2da, tem2db

    top_at_1 = (p_lev(1,1) .lt.  p_lev(1, nLev))
    if (top_at_1) then
       iSFC = nLev
       iTOA = 1
    else
       iSFC = 1
       iTOA = nLev
    endif

    if (iTOA .eq. 1) then
       tem2da(1:nCol,2:iSFC) = log(p_lay(1:nCol,2:iSFC))
       tem2db(1:nCol,2:iSFC) = log(p_lev(1:nCol,2:iSFC))
       do iCol = 1, nCol
           tem2da(iCol,1)    = log(p_lay(iCol,1) )
           tem2db(iCol,1)    = log(max(minP, p_lev(iCol,1)) )
           tem2db(iCol,iSFC) = log(p_lev(iCol,iSFC) )
       enddo
       t_lev(1:NCOL,1)      = t_lay(1:NCOL,iTOA)
       do iLay = 2, iSFC
          do iCol = 1, nCol
            t_lev(iCol,iLay) = t_lay(iCol,iLay) + (t_lay(iCol,iLay-1) - t_lay(iCol,iLay))&
                     * (tem2db(iCol,iLay)   - tem2da(iCol,iLay))                   &
                     / (tem2da(iCol,iLay-1) - tem2da(iCol,iLay))
           enddo
        enddo
       t_lev(1:NCOL,iSFC+1) = tsfc(1:NCOL)
    else
       tem2da(1:nCol,2:iTOA) = log(p_lay(1:nCol,2:iTOA))
       tem2db(1:nCol,2:iTOA) = log(p_lev(1:nCol,2:iTOA))
       do iCol = 1, nCol
           tem2da(iCol,1)    = log(p_lay(iCol,1))
           tem2db(iCol,1)    = log(p_lev(iCol,1))
           tem2db(iCol,iTOA) = log(max(minP, p_lev(iCol,iTOA)) )
       enddo
 
       t_lev(1:NCOL,1)      = tsfc(1:NCOL)
       do iLay = 1, iTOA-1
          do iCol = 1, nCol
            t_lev(iCol,iLay+1) = t_lay(iCol,iLay) + (t_lay(iCol,iLay+1) - t_lay(iCol,iLay))&
                     * (tem2db(iCol,iLay+1) - tem2da(iCol,iLay))                   &
                     / (tem2da(iCol,iLay+1) - tem2da(iCol,iLay))
           enddo
        enddo
       t_lev(1:NCOL,iTOA+1) = t_lay(1:NCOL,iTOA)
    endif

    ! Bound temperature at layer interfaces
    do iCol=1,NCOL
       do iLay=1,nLev+1
          if (t_lev(iCol,iLay) .le. minT) then
             t_lev(iCol,iLay) = minT + epsilon(minT)
          endif
          if (t_lev(iCol,iLay) .ge. maxT) then
             t_lev(iCol,iLay) = maxT - epsilon(maxT)
          endif
       enddo
    enddo

  end subroutine cmp_tlev

  ! #########################################################################################
  ! SUBROUTINE check_error_msg
  ! #########################################################################################
  subroutine check_error_msg(routine_name, error_msg)
    character(len=*), intent(in) :: &
         error_msg, routine_name
    
    if(error_msg /= "") then
       print*,"ERROR("//trim(routine_name)//"): "
       print*,trim(error_msg)
       return
    end if
  end subroutine check_error_msg  

end module radiation_tools
