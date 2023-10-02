!!****if* source/physics/Eos/unitTest/Hybrid/eos_initTest
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
!!  eos_initTest
!!
!!
!! SYNOPSIS
!!
!!  call eos_initTest()
!!
!!
!! DESCRIPTION
!!
!!  Initialize Eos unitTest scope variables which are typically from runtime parameters.
!!  The part of the unitTest located in the Simulation unit tree should arrange for this
!!  initialization routine to be called.
!!
!! ARGUMENTS
!!
!!
!!  none
!!
!!
!! PARAMETERS
!!
!!
!!    eos_testPresMode [STRING]
!!    eos_testEintMode [STRING]
!!    eos_testTempMode [STRING]
!!***

subroutine eos_initTest()

   use RuntimeParameters_interface, only: RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt

   use eos_testData, only: eos_testPresMode, eos_testPresModeStr, &
                           eos_testEintMode, eos_testEintModeStr, &
                           eos_testTempMode, eos_testTempModeStr, &
                           eos_test1allB, eos_test2allB, &
                           eos_test3allB, eos_test4allB, &
                           eos_testTolerance

#include "constants.h"

   implicit none

   character(len=MAX_STRING_LENGTH) :: testPresModeStr, &
                                       testEintModeStr, &
                                       testTempModeStr

   call RuntimeParameters_get("eos_testPresMode", testPresModeStr)
   call RuntimeParameters_get("eos_testEintMode", testEintModeStr)
   call RuntimeParameters_get("eos_testTempMode", testTempModeStr)

   eos_testPresModeStr = trim(testPresModeStr)
   eos_testEintModeStr = trim(testEintModeStr)
   eos_testTempModeStr = trim(testTempModeStr)

   call RuntimeParameters_mapStrToInt(eos_testPresModeStr, eos_testPresMode)
   call RuntimeParameters_mapStrToInt(eos_testEintModeStr, eos_testEintMode)
   call RuntimeParameters_mapStrToInt(eos_testTempModeStr, eos_testTempMode)

   call RuntimeParameters_get("eos_testTolerance", eos_testTolerance)

   print *, "DENS_PRES:", eos_testPresModeStr, eos_testPresMode
   print *, "DENS_EINT:", eos_testEintModeStr, eos_testEintMode
   print *, "DENS_TEMP:", eos_testTempModeStr, eos_testTempMode

   eos_test1allB = .true.
   eos_test2allB = .true.
   eos_test3allB = .true.
   eos_test4allB = .true.

end subroutine eos_initTest
