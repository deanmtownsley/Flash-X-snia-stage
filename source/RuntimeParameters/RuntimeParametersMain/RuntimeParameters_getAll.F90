!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_getAll
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
!!
!! NAME
!!  RuntimeParameters_getAll
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_getAll(          integer(inout) :: num, 
!!            character(len=MAX_STRING_LENGTH)(inout) :: names(num), 
!!                            real/int/str/log(inout) :: values(num), 
!!                                 logical(inout)     :: changed(num))
!!
!! DESCRIPTION
!!
!!  This function gets all the parameters of a given type from the
!!  RuntimeParameters database, and indicates whether they have changed
!!  during the current restart. 
!!  
!!  The routines below are implemented; they are not overloaded,
!!    so call them directly. 
!!  RuntimeParameter_getAllReal
!!  RuntimeParameter_getAllInt
!!  RuntimeParameter_getAllStr
!!  RuntimeParameter_getAllLog
!!
!! ARGUMENTS
!!
!!
!! num:        number of parameter to get
!! names:      names of parameters
!! values:     values of parameters
!! changed:    logicals indicating if parameter changed or not 
!!             from the initial run
!!
!!
!!
!!***


   
subroutine RuntimeParameters_getAllReal (num, names, values, changed)

  use RuntimeParameters_data, ONLY : parameter

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)          :: names(num)
  real, intent(inout)                                   :: values(num)
  logical, intent(inout)                                :: changed(num)


  call nameValueLL_getAllReal(parameter, num, names, values, changed)

  return
  
end subroutine RuntimeParameters_getAllReal


   
subroutine RuntimeParameters_getAllInt (num, names, values, changed)

  use RuntimeParameters_data, ONLY : parameter

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)       :: names(num)
  integer, intent(inout)                                :: values(num)
  logical, intent(inout)                                :: changed(num)

  call nameValueLL_getAllInt(parameter, num, names, values, changed)

  return
  
end subroutine RuntimeParameters_getAllInt


   
subroutine RuntimeParameters_getAllStr (num, names, values, changed)

  use RuntimeParameters_data, ONLY : parameter

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)       :: names(num), values(num)
  logical, intent(inout)                                :: changed(num)


  call nameValueLL_getAllStr(parameter, num, names, values, changed)

  return
  
end subroutine RuntimeParameters_getAllStr


   
subroutine RuntimeParameters_getAllLog (num, names, values, changed)

  use RuntimeParameters_data, ONLY : parameter

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)          :: names(num)
  logical, intent(inout)                                   :: values(num)
  logical, intent(inout)                                :: changed(num)


  call nameValueLL_getAllLog(parameter, num, names, values, changed)

  return
  
end subroutine RuntimeParameters_getAllLog



