!!****ih* source/flashUtilities/nameValueLL/nameValueLL_logRules
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
!!  NAME
!!    nameValueLL_logRules
!!
!!  SYNOPSIS
!!
!!    nameValueLL_logRules(  character(len=*), intent(in)      :: name,
!!                           integer(in)                       :: numValues,
!!                           integer/Real, dimension(numValues):: minValues,
!!                           integer/Real, dimension(numValues):: maxValues)
!! 
!!  DESCRIPTION
!!    Prints the rules for given variable in log file
!!
!!  ARGUMENTS
!!
!!    name:         name of parameter
!!    numValues:    number of valid values
!!    minValues:    array of given size, provides minimum of valid values
!!    maxValues:    array of given size, provides maximum of valid values
!!
!!***

subroutine nameValueLL_logRulesInt ( name, numValues, minValues, maxValues)

  implicit none

#include "constants.h"

   character(len=*), intent(in) :: name
   integer, intent(in) :: numValues
   integer,dimension(numValues),intent(in):: minValues,maxValues
   integer :: ctr

   ! print the message to screen directly (no Logfile setup yet)
   print *,"[nameValueLL] Rules for integer Parameter ",trim(name),": Union of following closed intervals"
   do ctr =1, numValues
      print *, "[nameValueLL] ",minValues(ctr)," ... ", maxValues(ctr)
   enddo

end subroutine nameValueLL_logRulesInt
   
subroutine nameValueLL_logRulesReal ( name, numValues, minValues, maxValues)

  implicit none

#include "constants.h"

   character(len=*), intent(in) :: name
   integer, intent(in) :: numValues
   real,dimension(numValues),intent(in):: minValues,maxValues
   integer :: ctr

   print *, "[nameValueLL] Rules for real Parameter ", trim(name),": Union of following intervals"
   do ctr =1, numValues
      print *, "[nameValueLL] ",minValues(ctr)," ... ", maxValues(ctr)
   enddo

end subroutine nameValueLL_logRulesReal

subroutine nameValueLL_logRulesStr ( name, numValues, validValues )

  implicit none
#include "constants.h"

   character(len=*), intent(in) :: name
   integer, intent(in) :: numValues
   character(len=*),dimension(numValues),intent(in):: validValues
   character (len=MAX_STRING_LENGTH) :: buf
   integer :: ctr

   print *,"[nameValueLL] Rules for string Parameter ",trim(name),": One of the following values"
   do ctr =1, numValues
      print *,"[nameValueLL] '",trim(validValues(ctr)),"'"
   enddo

end subroutine nameValueLL_logRulesStr
