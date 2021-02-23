module cires_tauamf_data

  use machine, only: kind_phys
!...........................................................................................
! tabulated GW-sources: GRACILE/Ern et al., 2018 and/or Resolved GWs from C384-Annual run
!...........................................................................................
implicit none

   integer           :: ntau_d1y, ntau_d2t  
   real(kind=kind_phys), allocatable :: ugwp_taulat(:)
   real(kind=kind_phys), allocatable :: tau_limb(:,:), days_limb(:)
   logical           :: flag_alloctau = .false.          
   character(len=255):: ugwp_taufile =  'ugwp_limb_tau.nc' 

   public :: read_tau_amf, cires_indx_ugwp, tau_amf_interp 

contains
  
   subroutine read_tau_amf(me, master, errmsg, errflg)
  
    use  netcdf
    integer, intent(in) ::  me, master   
    integer :: ncid,  iernc, vid, dimid, status         
    integer :: k
    
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg    
!
      
      iernc=NF90_OPEN(trim(ugwp_taufile), nf90_nowrite, ncid)
     
       if(iernc.ne.0) then         
          write(errmsg,'(*(a))') "read_tau_amf: cannot open file_limb_tab data-file ",  &
                                    trim(ugwp_taufile)
	    print *, 'cannot open ugwp-v1 tau-file=',trim(ugwp_taufile)			    
          errflg = 1
          return
        else


       status = nf90_inq_dimid(ncid, "lat", DimID)
!      if (status /= nf90_noerr) call handle_err(status)
!
       status = nf90_inquire_dimension(ncid, DimID,  len =ntau_d1y )
       
       status = nf90_inq_dimid(ncid, "days", DimID)
       status = nf90_inquire_dimension(ncid, DimID,  len =ntau_d2t )
       
           if (me == master)  print *, ntau_d1y, ntau_d2t, ' dimd of tau_ngw ugwp-v1 '
	   if (ntau_d2t .le. 0 .or. ntau_d1y .le. 0) then 
	       print *, 'ugwp-v1 tau-file=',    trim(ugwp_taufile)	   
	       print *, '  ugwp-v1: ', 'ntau_d2t=',ntau_d2t, 'ntau_d2t=',ntau_d1y
	       stop
	   endif
	   	   
        if (.not.allocated(ugwp_taulat))  allocate (ugwp_taulat(ntau_d1y ))
        if (.not.allocated(days_limb))    allocate (days_limb(ntau_d2t))
        if (.not.allocated(tau_limb))     allocate (tau_limb(ntau_d1y, ntau_d2t ))   	   
              
	iernc=nf90_inq_varid( ncid, 'DAYS', vid )
        iernc= nf90_get_var( ncid, vid, days_limb)
	iernc=nf90_inq_varid( ncid, 'LATS', vid )
        iernc= nf90_get_var( ncid, vid, ugwp_taulat)
	iernc=nf90_inq_varid( ncid, 'ABSMF', vid )
        iernc= nf90_get_var( ncid, vid, tau_limb)
			
	iernc=nf90_close(ncid)
	
	endif    
	
  end  subroutine read_tau_amf  
  
    subroutine cires_indx_ugwp (npts, me, master, dlat,j1_tau,j2_tau, w1_j1tau, w2_j2tau)
     
    use machine, only: kind_phys
    		 
    implicit none
    
      integer, intent(in)                                      ::   npts, me, master
      real(kind=kind_phys) ,   dimension(npts), intent(in)     ::   dlat 
           
      integer, dimension(npts), intent(inout)                  ::  j1_tau,   j2_tau
      real(kind=kind_phys) ,   dimension(npts), intent(inout)  ::  w1_j1tau, w2_j2tau
      
!locals

      integer :: i,j, j1, j2     
!     
      do j=1,npts
        j2_tau(j) = ntau_d1y
        do i=1,ntau_d1y
          if (dlat(j) < ugwp_taulat(i)) then
            j2_tau(j) = i
            exit
          endif
        enddo
	
      
        j2_tau(j) = min(j2_tau(j),ntau_d1y)
        j1_tau(j) = max(j2_tau(j)-1,1)	
	
        if (j1_tau(j) /= j2_tau(j) ) then
          w2_j2tau(j) = (dlat(j)  - ugwp_taulat(j1_tau(j))) &
                 / (ugwp_taulat(j2_tau(j))-ugwp_taulat(j1_tau(j)))       	 
        else
          w2_j2tau(j) = 1.0
        endif
          w1_j1tau(j) = 1.0 -	w2_j2tau(j)	
      enddo
      return
    end subroutine cires_indx_ugwp   
    
    subroutine tau_amf_interp(me, master, im, idate, fhour, j1_tau,j2_tau, ddy_j1, ddy_j2, tau_ddd)    
    use machine, only: kind_phys	           
    implicit none
    
!input    
    integer, intent(in)               :: me, master
    integer, intent(in)               :: im, idate(4)
    real(kind=kind_phys), intent(in)  :: fhour
      
    real(kind=kind_phys), intent(in), dimension(im) ::  ddy_j1, ddy_j2
    integer             , intent(in), dimension(im) ::  j1_tau,j2_tau        
!ouput    
    real(kind=kind_phys),  dimension(im)   ::  tau_ddd
!locals

    integer :: i, j1, j2, it1, it2 , iday
    integer :: ddd    
    real(kind=kind_phys)  :: tx1, tx2, w1, w2, fddd 
!
! define day of year ddd ..... from the old-fashioned "GFS-style"
! 
         call gfs_idate_calendar(idate, fhour, ddd, fddd)  
    
            it1 = 2
         do iday=1, ntau_d2t
	    if (fddd .lt. days_limb(iday) ) then
	    it2 = iday
	    exit
	    endif
	 enddo
	 
	 it2 = min(it2,ntau_d2t)	 
	 it1 = max(it2-1,1)
	 if (it2 > ntau_d2t ) then
	  print *, ' Error in time-interpolation for tau_amf_interp '	 
	  print *, ' it1, it2, ntau_d2t ', it1, it2, ntau_d2t
	  print *, ' Error in time-interpolation see cires_tauamf_data.F90 '	  
	  stop
	 endif
	 
	 w2 = (fddd-days_limb(it1))/(days_limb(it2)-days_limb(it1))
	 w1 = 1.0-w2     
       
      do i=1, im	 
	 j1 = j1_tau(i)
	 j2 = j2_tau(i)
	 tx1 = tau_limb(j1, it1)*ddy_j1(i)+tau_limb(j2, it1)*ddy_j2(i)
	 tx2 = tau_limb(j1, it2)*ddy_j1(i)+tau_limb(j2, it2)*ddy_j2(i)	 
	 tau_ddd(i) =  tx1*w1 + w2*tx2
      enddo
             
    end subroutine tau_amf_interp  
    
    subroutine gfs_idate_calendar(idate, fhour, ddd, fddd) 
    
    use machine, only: kind_phys    		 
    implicit none  
! input     
    integer, intent(in)                 :: idate(4)
    real(kind=kind_phys), intent(in)   :: fhour
!out    
    integer, intent(out)                :: ddd    
    real(kind=kind_phys), intent(out)  :: fddd  
!
!locals
!
      real(kind=kind_phys) :: rinc(5), rjday
      integer              :: jdow, jdoy, jday
      real(4)              :: rinc4(5)
      integer              :: w3kindreal, w3kindint
      
      integer ::  iw3jdn
      integer :: jd1, jddd
      
      integer  idat(8),jdat(8)  
          
       
      idat(1:8)    = 0
      idat(1) = idate(4)
      idat(2) = idate(2)
      idat(3) = idate(3)
      idat(5) = idate(1)
      rinc(1:5)    = 0.
      rinc(2) = fhour
!    
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal==4) then
        rinc4 = rinc
        call w3movdat(rinc4, idat,jdat)
      else
        call w3movdat(rinc,  idat,jdat)
      endif           
!     jdate(8)- date and time (yr, mo, day, [tz], hr, min, sec)
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow, ddd, jday)
      fddd = float(ddd) + jdat(5) / 24.        
    end  subroutine gfs_idate_calendar    
    
end  module cires_tauamf_data
