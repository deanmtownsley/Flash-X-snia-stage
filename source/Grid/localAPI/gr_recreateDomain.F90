!!****if* source/Grid/localAPI/gr_recreateDomain
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
!!  gr_recreateDomain
!!
!! 
!! SYNOPSIS
!!
!!  call gr_recreateDomain()
!!
!!
!! DESCRIPTION
!!
!!  Construct the top-level block structure.  Normally, when the Grid
!!  is created from scratch,   
!!  gr_createDomain() sets up an Nblockx * Nblocky * Nblockz array
!!  of top-level blocks.
!!  This is a pared-down variant of gr_createDomain appropriate when
!!  restarting from a checkpoint.  It is only really needed in simulations
!!  where Simulation_defineDomain is used to add obstacle blocks.
!!
!!
!! ARGUMENTS
!!
!!  NONE
!!
!! NOTES
!! 
!!  The user-defined Simulation_defineDomain should define the obstacles
!!  in the same ways as in previous runs up to the checkpoint from which
!!  the simulation is restarting. If Simulation_defineDomain is changed
!!  in mid-run, the resulting set of obstacle blocks may not be consistent
!!  any more with the properties, especially neighbors, of blocks read
!!  from the checkpoint file, and all bets are off.
!!
!!***

subroutine gr_recreateDomain

  implicit none
end subroutine gr_recreateDomain
