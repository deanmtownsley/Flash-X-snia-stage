subroutine Driver_evolveAll()
   use Driver_interface, only: Driver_getMype
   use Timers_interface, only: Timers_getSummary, Timers_start, Timers_stop
   use Logfile_interface, only: Logfile_stamp

   use Spacetime_interface, only: Spacetime_unitTest

   use iso_fortran_env, only: stdout => output_unit

#include "constants.h"

   implicit none

   logical ::  perfect = .true.

   integer :: procID

   call Driver_getMype(GLOBAL_COMM, procID)

   write (stdout, "(a,i0.4,a)") "Proc #", procID, " - Spacetime unit test started"

   perfect = .true.

   call Timers_start("testing")
   call Spacetime_unitTest(stdout, perfect)
   call Timers_stop("testing")

   if (perfect) then
      write (stdout, "(a,i0.4,a)") "Proc #", procID, " - Spacetime unit test PASSED"
   else
      write (stdout, "(a,i0.4,a)") "Proc #", procID, " - Spacetime unit test FAILED"
   end if

   call Timers_getSummary(0)

   call Logfile_stamp("Flash-X Spacetime unit test run complete.", "LOGFILE_END")
end subroutine Driver_evolveAll
