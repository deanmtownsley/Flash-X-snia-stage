!!****if* source/Grid/GridSolvers/Grid_setSolverDbgContextInfo
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
!!  Grid_setSolverDbgContextInfo
!!
!! SYNOPSIS
!!
!!  call Grid_setSolverDbgContextInfo(integer(in),OPTIONAL  :: component,
!!                                    integer(in),OPTIONAL  :: group)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   component : 
!!
!!   group : 
!!
!!
!!
!!***

subroutine Grid_setSolverDbgContextInfo(component,group)
  use gr_solversData, ONLY: gr_solversDbgContext
  implicit none
  integer,intent(in),OPTIONAL :: component, group

  if (.NOT. present(group) .AND. .NOT. present(component)) then
     ! clear the info
     gr_solversDbgContext%component = 0
     gr_solversDbgContext%group     = 0
  else if (present(group)) then
     gr_solversDbgContext%group     = group
     if (present(component)) then
        gr_solversDbgContext%component = component
     else
        gr_solversDbgContext%component = 3
     end if
  else
     gr_solversDbgContext%component = component
  end if

  gr_solversDbgContext%libErrCode = 0
  gr_solversDbgContext%flashErrCode = 0

end subroutine Grid_setSolverDbgContextInfo
