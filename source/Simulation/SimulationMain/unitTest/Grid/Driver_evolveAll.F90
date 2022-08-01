!!****if* source/Simulation/SimulationMain/unitTest/Grid/Driver_evolveAll
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
!!  Driver_evolveAll
!!
!! SYNOPSIS
!!
!!  Driver_evolveAll()
!!
!! DESCRIPTION
!!
!!  This is the global driver specifically for the Grid unit test
!!  It opens a file for each of the processor that it is running
!!  on and call the Grid_unitTest with the file unit. The Grid_unitTest
!!  functions carries out the testing of the unit, and reports success or
!!  failure through a logical argument "perfect". Upon return from 
!!  Grid_unitTest, the Driver_evolveAll function makes appropriate notification 
!!  in the file which the test suite can parse to determine if the test was
!!  successful.
!!
!!
!!***

subroutine Driver_evolveAll()

#include "constants.h"

  use Driver_data, ONLY: dr_globalMe 
  use Grid_interface, ONLY : Grid_unitTest
  use Timers_interface, ONLY : Timers_getSummary,Timers_start,Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

  logical ::  perfect = .true.

  character(len=20) :: fileName
  integer, parameter        :: fileUnit = 2
  integer,dimension(4) :: prNum
  integer :: temp,i

  ! stays true if no errors are found
  perfect = .true.
  
  temp = dr_globalMe

  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))

  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe

  print *, "Preparing to call Grid_unitTest"
  call Timers_start("testing")
  call Grid_unitTest(fileUnit,perfect)
  call Timers_stop("testing")
  print *, "    Returned from Grid_unitTest"

  call Timers_start("writing")
  if (perfect) then
    write(fileUnit,'("all results conformed with expected values.")')
  endif

  close(fileUnit)
  call Timers_stop("writing")

  call Timers_getSummary(0)
  call Logfile_stamp("Flash-X unitTest run complete.", "LOGFILE_END")

  return
  
end subroutine Driver_evolveAll
