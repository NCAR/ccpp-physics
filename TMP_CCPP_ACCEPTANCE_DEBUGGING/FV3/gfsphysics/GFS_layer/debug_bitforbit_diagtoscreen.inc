#ifdef DHDEBUG
   write(0,*) "Calling GFS_diagtoscreen_run"
   call GFS_diagtoscreen_run(Model, Statein, Stateout, Sfcprop, Coupling, &
                             Grid, Tbd, Cldprop, Radtend, Diag, Tbd%blkno)
#endif