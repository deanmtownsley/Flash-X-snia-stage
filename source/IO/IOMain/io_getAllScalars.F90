!!****if* source/IO/IOMain/io_getAllScalars
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
!!  io_getAllScalars
!!
!! SYNOPSIS
!!
!!  io_getAllScalars(integer(inout)        :: num, 
!!                   char(len=*)(in)       :: names(num),
!!                   real/int/str/log(out) :: values(num),
!!                   logical               :: changed(num))
!!                      
!!
!! DESCRIPTION
!!
!! This function gets all the scalars of a given type and returns
!! an array 'changed' indicating if 
!! they have changed since most recent restart.
!! Uses linked list implemented under the hood.  
!!
!! ARGUMENTS
!!
!! num:         number of scalars to get
!! names:       names of scalar
!! values:      scalar values
!! changed:     returned array indicating if values have changed since 
!!              restart
!!
!!
!!
!!***


   
subroutine io_getAllScalarsReal (num, names, values, changed)

  use IO_data, ONLY : io_scalar

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)          :: names(num)
  real, intent(inout)                                   :: values(num)
  logical, intent(inout)                                :: changed(num)


  call NameValueLL_getAllReal(io_scalar, num, names, values, changed)

  return
  
end subroutine io_getAllScalarsReal


   
subroutine io_getAllScalarsInt (num, names, values, changed)

  use IO_data, ONLY : io_scalar

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)       :: names(num)
  integer, intent(inout)                                :: values(num)
  logical, intent(inout)                                :: changed(num)



  call NameValueLL_getAllInt(io_scalar, num, names, values, changed)

  return
  
end subroutine io_getAllScalarsInt


   
subroutine io_getAllScalarsStr (num, names, values, changed)

  use IO_data, ONLY : io_scalar

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)       :: names(num), values(num)
  logical, intent(inout)                                :: changed(num)


  call NameValueLL_getAllStr(io_scalar, num, names, values, changed)

  return
  
end subroutine io_getAllScalarsStr


   
subroutine io_getAllScalarsLog (num, names, values, changed)

  use IO_data, ONLY : io_scalar

implicit none
#include "constants.h"

  integer, intent(inout)                                :: num
  character(len=MAX_STRING_LENGTH), intent(inout)          :: names(num)
  logical, intent(inout)                                   :: values(num)
  logical, intent(inout)                                :: changed(num)


  call NameValueLL_getAllLog(io_scalar, num, names, values, changed)

  return
  
end subroutine io_getAllScalarsLog



