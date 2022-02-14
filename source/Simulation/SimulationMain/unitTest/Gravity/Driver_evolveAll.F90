!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Driver_evolveAll
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
!!  call Driver_evolveAll()
!!
!! DESCRIPTION
!!
!!  This is a modification of the standard Driver_evolveAll for the unitTest for 
!!  the Gravity unit.
!!
!! NOTES
!!
!!  Some Gravity tests include a more specific implementation of Driver_evolveAll.
!!  If an implementation Driver_evolveAll.F90 exists in a subdirectory of
!!  SimulationMain/unitTest/Gravity, that one will normally be used in preference
!!  over this less specifc implementation.
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveAll()
       
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use IO_interface, ONLY :IO_writeCheckpoint, IO_writePlotfile
  use Gravity_interface, ONLY: Gravity_potential, Gravity_unitTest

  use Driver_data, ONLY: dr_globalMe, dr_nstep

  implicit none

#include "constants.h"
#include "Simulation.h"

  integer   :: i


  
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr


  logical ::  perfect = .true.

  character(len=20) :: fileName
  integer, parameter        :: fileUnit = 2
  integer,dimension(4) :: prNum
  integer :: temp


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



  call Logfile_stamp( 'Starting Calculation' , '[Driver_evolveAll]')
  print*,'started calculation'
  call Timers_start("calculation")


  call Gravity_potential()      ! no op for potential implementations that are constant in time
  print*,'Gravity_potential has been called.'

  call Gravity_unitTest(fileUnit,perfect)

  call Timers_stop("calculation")

  call Logfile_stamp( 'Ending Calculation' , '[Driver_evolveAll]')

  call IO_writeCheckpoint()

  call IO_writePlotfile( )

  call Timers_getSummary( dr_nstep)

  
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")

  call Logfile_close()


  !finish unit test write out file

  if (perfect) then
    write(fileUnit,'("all results conformed with expected values.")')
  endif

  close(fileUnit)



  return
  
end subroutine Driver_evolveAll



