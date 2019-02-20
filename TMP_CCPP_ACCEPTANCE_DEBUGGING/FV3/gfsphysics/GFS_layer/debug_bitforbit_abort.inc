#if !defined(CCPP) && defined(TRANSITION) && defined(DBGPRNT)
   write(0,*) "Calling GFS_abort_run"
   call GFS_abort_run(Model, Tbd%blkno)
#endif
