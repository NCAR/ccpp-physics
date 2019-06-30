!===============================
!  Part -3  init  wave solvers
!===============================

  module ugwp_lsatdis_init 
     implicit none
     
      integer  :: nwav, nazd
      integer  :: nst
      real     :: eff 
      integer, parameter  :: incdim = 4
      integer, parameter  :: iazdim = 4            
!    
     contains
           
     subroutine initsolv_lsatdis(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw) 
       
     implicit none
! 
     integer  :: me, master
     integer  :: nwaves, nazdir
     integer  :: nstoch
     real     :: effac  
     logical  :: do_physb
     real     :: kxw     
!     
!locals: define azimuths and Ch(nwaves) - domain when physics-based soureces
!                                          are not actibve
!
     integer :: inc, jk, jl, iazi, i, j, k    
            
     if( nwaves == 0 .or. nstoch == 1 ) then
!                                redefine from the default       
       nwav = incdim
       nazd = iazdim
       nst  = 0
       eff  = 1.0
     else
!                                from input_nml multi-wave spectra   
       nwav = nwaves
       nazd = nazdir
       nst  = nstoch
       eff  = effac  
     endif        
!
!       
       
       
     end subroutine initsolv_lsatdis
!     
  end module ugwp_lsatdis_init     
!
