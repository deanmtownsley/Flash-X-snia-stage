!!****f* source/Simulation/Simulation_initBlock
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
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            integer(IN)  :: blockDesc  )
!!
!!
!!
!! DESCRIPTION
!!  This routine applies initial conditions of a specific simulation
!!  to the specified block.
!!
!! 
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  blockDesc -        describes the block to initialize
!!
!! PARAMETERS
!!
!!  eosModeInit -     after this routine sets up initial conditions,
!!                    the grid package calls Eos to insure the
!!                    values are thermodynamically consistent.  This
!!                    parameter controls the mode of application of
!!                    the Eos.  Its default is "dens_ie", and it can
!!                    be one of "dens_ie", "dens_pres", "dens_temp".
!!                    Setting this value to dens_temp, for instance,
!!                    would make it possible to leave this routine
!!                    having just initialized density and temperature,
!!                    leaving Eos to calculate the rest.
!!
!! SEE ALSO
!!
!!  Eos_multiDim
!!***

subroutine Simulation_initBlock(solnData, tileDesc)
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

  real,dimension(:,:,:,:),pointer :: solnData
  type(Grid_tile_t), intent(in) :: tileDesc

  return
end subroutine Simulation_initBlock
