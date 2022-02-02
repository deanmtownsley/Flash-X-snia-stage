!!****if* source/Driver/DriverMain/dr_wallClockLimitExceeded
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!
!!  dr_wallClockLimitExceeded
!!
!! SYNOPSIS
!!
!!  call dr_wallClockLimitExceeded(logical(OUT)  :: endRunWallClock)
!!
!! DESCRIPTION
!!
!!  Determine whether the wall clock time limit has been exceeded.
!!
!! ARGUMENTS
!!
!!   endRunWallClock : the answer is returned here.
!!
!! NOTES
!!
!!  This routine must be called collectively by all MPI tasks.
!!  The elapsed time is determiend on the master PE and then
!!  broadcast to all tasks.
!!
!!  This routines is usually called from Driver_evolveAll.
!!
!!***

subroutine dr_wallClockLimitExceeded(endRunWallClock)
  use Driver_interface, ONLY: Driver_getElapsedWCTime
  use Driver_data, ONLY: dr_globalMe, dr_globalComm, &
                         dr_elapsedWCTime, dr_wallClockTimeLimit
#include "Flashx_mpi_implicitNone.fh"
  logical, intent(OUT) :: endRunWallClock

#include "constants.h"
  integer :: ierr

  endRunWallClock = .false.
  if (dr_wallClockTimeLimit .GE. 0.0) then
     call Driver_getElapsedWCTime(dr_elapsedWCTime)
     if (dr_globalMe == MASTER_PE) then
        if (dr_elapsedWCTime >  dr_wallClockTimeLimit) endRunWallClock = .true.
     end if
     call MPI_Bcast(endRunWallClock, 1, FLASH_LOGICAL, MASTER_PE, dr_globalComm, ierr)
  end if

end subroutine dr_wallClockLimitExceeded
