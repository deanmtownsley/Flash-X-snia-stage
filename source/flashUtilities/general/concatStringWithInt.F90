!!****if* source/flashUtilities/general/concatStringWithInt
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
!!    concatStringWithInt
!!
!! SYNOPSIS
!!
!!    concatStringWithInt()
!!
!! DESCRIPTION
!!
!!  This routine is used for generating a new string name where
!!  the base name is appended with some integer value at the end
!!  For example if stringBaseName is "refine_var_", and the value
!!  of the integer is "1", then the routine will return "refine_var_1"
!!
!!
!!***

subroutine concatStringWithInt(inString, ind, outString )
  
#include "constants.h"

  implicit none
  character(len=*),intent(IN) :: inString
  integer,intent(IN) :: ind
  character(len=MAX_STRING_LENGTH),intent(OUT) :: outString
  
  
  character(len=MAX_STRING_LENGTH) :: parname
  integer :: pos


  pos=len_trim(inString)+1
  parname= ""
  write(parname, *) ind
  parname = ADJUSTL(parname)
  outString=inString(:pos-1)//parname
  

end subroutine concatStringWithInt
