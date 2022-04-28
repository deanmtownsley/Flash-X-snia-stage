!!****if* source/Grid/GridMain/Grid_getNumVars
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
!!  Grid_getNumVars
!! 
!! SYNOPSIS
!!  call Grid_getNumVars(integer(in)  :: gridStruct,
!!                       integer(out) :: nVar)
!!
!! DESCRIPTION
!!  Returns the number of variables in each mesh data
!!  structure.  Allows us to write more flexible code that
!!  is not dependent on NUNK_VARS and NFACE_VARS, etc.
!!
!! ARGUMENTS
!!  gridStruct -- Integer value representing mesh data structure;
!!                one of CENTER,FACEX,FACEY,FACEZ,SCRATCH_CTR,etc.
!!  nVar -- Number of mesh variables in gridStruct
!!
!! NOTES
!!
!!  The symbols CENTER,FACEX,FACEY,FACEZ,SCRATCH_CTR, etc., for use in
!!  in the gridStruct argument, are defined in header file constants.h .
!!
!!***
#include "constants.h"
#include "Simulation.h"

subroutine Grid_getNumVars(gridStruct, nVar)  
  use Driver_interface, ONLY: Driver_abort
  implicit none
  integer, intent(in) :: gridStruct
  integer, intent(out) :: nVar

  select case (gridStruct)
  case(CENTER)
     nVar = NUNK_VARS
  case(FACEX,FACEY,FACEZ)
     nVar = NFACE_VARS
  case(SCRATCH)
     nVar = NSCRATCH_GRID_VARS
  case(SCRATCH_CTR)
     nVar = NSCRATCH_CENTER_VARS
  case(SCRATCH_FACEX)
     nVar = NSCRATCH_FACEX_VARS
  case(SCRATCH_FACEY)
     nVar = NSCRATCH_FACEY_VARS
  case(SCRATCH_FACEZ)
     nVar = NSCRATCH_FACEZ_VARS
  case DEFAULT
     call Driver_abort("[Grid_getNumVars]: Invalid data structure")
  end select
end subroutine Grid_getNumVars
