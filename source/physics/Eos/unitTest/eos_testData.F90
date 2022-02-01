!!****ih* source/physics/Eos/unitTest/eos_testData
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
!!  eos_testData
!!
!! 
!! SYNOPSIS
!!
!! use eos_testData
!!
!! DESCRIPTION
!!
!!***

module eos_testData

  implicit none

#include "constants.h"

  character(len=MAX_STRING_LENGTH), save :: eos_testPresModeStr, eos_testEintModeStr, eos_testTempModeStr
  integer, save :: eos_testPresMode, eos_testEintMode, eos_testTempMode

  real, save    :: eos_testTolerance

  logical, save :: eos_test1allB,eos_test2allB,eos_test3allB,eos_test4allB !for all blocks

end module eos_testData
