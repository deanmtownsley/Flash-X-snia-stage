!!****ih* source/flashUtilities/general/ut_generalInterface
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
!!  ut_generalInterface
!!
!! SYNOPSIS
!!
!!  use ut_generalInterface
!!
!! DESCRIPTION
!!
!!  Interface module for some general utilities.
!!
!!***

! Modification history:
!     Created   March 2016  KW

module ut_generalInterface

  interface
     integer function ut_getFreeFileUnit ()
     end function ut_getFreeFileUnit
  end interface

  interface
     subroutine ut_printMatrix (fileUnit,                   &
          title,                      &
          rowMinMatrix, rowMaxMatrix, &
          colMinMatrix, colMaxMatrix, &
          rowMinPrint , rowMaxPrint,  &
          colMinPrint , colMaxPrint,  &
          matrix                      )
       implicit none
       integer,           intent (in) :: fileUnit
       character (len=*), intent (in) :: title
       integer,           intent (in) :: rowMinMatrix, rowMaxMatrix
       integer,           intent (in) :: colMinMatrix, colMaxMatrix
       integer,           intent (in) :: rowMinPrint , rowMaxPrint
       integer,           intent (in) :: colMinPrint , colMaxPrint
       real,              intent (in) :: matrix (rowMinMatrix:rowMaxMatrix , colMinMatrix:colMaxMatrix)
     end subroutine ut_printMatrix
  end interface


end module ut_generalInterface
