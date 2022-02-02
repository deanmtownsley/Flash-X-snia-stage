!!****f* source/Grid/Grid_parseNonRep
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
!!  Grid_parseNonRep
!!
!!
!! SYNOPSIS
!!
!!  Grid_parseNonRep(character(IN) :: strlwr(*),
!!                     integer(OUT):: nonrep,
!!                     integer(OUT):: idx)
!!
!!
!! DESCRIPTION
!!
!!  Given the string name of a nonreplicated variable array element, returns the array id and index.
!!
!!
!! ARGUMENTS
!!  
!!   strlwr - the all lowercase string name of the array element, MUST BE ALL lowercase
!!   nonrep - receives the array id
!!   idx - receives the array index
!!
!!***
subroutine Grid_parseNonRep(strlwr, nonrep, idx)
    implicit none
    character(*), intent(in) :: strlwr
    integer, intent(out) :: nonrep, idx
    ! signal no match
    nonrep = 0
    idx = 0
end subroutine
