!!****if* source/flashUtilities/nameValueLL/makeLowercase
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
!!    makeLowercase
!! 
!!  SYNOPSIS
!!    makeLowercase(character(len=*)(INOUT) :: str)
!!
!!  DESCRIPTION 
!!    Convert a string to lowercase.
!!
!!  ARGUMENTS
!!    str:    string to be converted 
!!
!!***

subroutine makeLowercase (str)
  
implicit none
  character(len=*),intent(INOUT) :: str
  integer         :: i
  
  do i = 1, len_trim(str)
     if (lge(str(i:i), 'A') .and. lle(str(i:i), 'Z')) & 
          &      str(i:i) = achar( iachar(str(i:i)) + 32 )
  enddo
  
  return
end subroutine makeLowercase
