!!****if* source/monitors/Logfile/LogfileMain/log_getUnitsArr
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
!!  log_getUnitsArr
!!
!! SYNOPSIS
!!
!!  log_getUnitsArr(integer(in) :: arrLen,
!!                  character(len=MAX_STRING_LENGTH)(inout) :: strArr)
!!
!! DESCRIPTION
!!
!!  This is a helper routine which formats the list of all the units
!!  included in a given simulation so that they can be written to the
!!  Logfile.
!!
!! ARGUMENTS
!!
!!   arrLen : number of units included in the simulation
!!
!!   strArr : returned value holding the formatted output 
!!
!!
!!
!!***


subroutine log_getUnitsArr(arrLen, strArr)

  use Logfile_data, ONLY : log_endOfLine, log_strArr

  implicit none


#include "constants.h"
      
  integer, intent(in)                                :: arrLen
  character(len=MAX_STRING_LENGTH), intent(inout)    :: strArr(arrLen)
  integer                                            :: i, j
  character(len=MAX_STRING_LENGTH)                   :: flash_mods(arrLen)
    


  call log_allocateStrArr(arrLen, 1)
  call setup_getFlashUnits(log_strArr)
  
  do i = 1, arrLen
     do j = 1, MIN(len_trim(log_strArr(i,1)), (MAX_STRING_LENGTH-1))
        strArr(i)(j:) = log_strArr(i,1)(j:j)
     end do
     strArr(i)(j:) = log_endOfLine
  end do
  deallocate(log_strArr)

end subroutine log_getUnitsArr
