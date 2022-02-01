!!****f* source/Driver/Driver_getComm
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
!!  Driver_getComm
!!
!! SYNOPSIS
!!
!!  Driver_getComm(integer(IN) :: communicatorType,
!!                 integer(OUT):: communicator,
!!       optional, integer(IN) :: axis)
!!               
!!  
!! DESCRIPTION 
!!
!!  All units can use this interface to query the 
!!  communicators in use. At initialization all units call this interface.
!!  One can ask for either global communicator, or mesh communicators
!!  or axis communicators used in some implementations of Uniform mesh
!!  When there are no duplicate copies of the mesh, the mesh communicator
!!  defaults to MPI_COMM_WORLD
!! 
!!
!! ARGUMENTS
!!
!!  communicatorType - Input argument indicating whether to return the global 
!!                     communicator (GLOBAL_COMM), or mesh communicator
!!                     (MESH_COMM) that allows duplication of mesh or
!!                     it is communicator that handles rows and colums of
!!                     processors (AXIS_COMM)
!! communicator      - output argument in which to return the communicator
!! axis              - this optional argument is included only if 
!!                     communicatorType is AXIS_COMM, the valid 
!!                     values are IAXIS, JAXIS or KAXIS
!!
!!***

#include "constants.h"

subroutine Driver_getComm(communicatorType, communicator, axis)

  implicit none
  integer, INTENT(IN) :: communicatorType
  integer, INTENT(OUT) :: communicator
  integer, optional, INTENT(IN) :: axis
  communicator = NO_COMM
end subroutine Driver_getComm
