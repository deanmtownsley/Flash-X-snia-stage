!!****if* source/IO/IOMain/IO_updateScalars
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
!!  IO_updateScalars
!!
!!
!! SYNOPSIS
!!
!!  IO_updateScalars() 
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  It calls all of the Unit output routines that allow a user to write out
!!  scalar data.  The purpose is to send each unit's scalars to a scalar list
!!  to be output together in a checkpoint file.
!!
!! 
!! ARGUMENTS
!!
!!
!! NOTES
!!
!!
!!
!! SIDE EFFECTS
!!
!!
!!***

subroutine IO_updateScalars()
  use Hydro_interface, ONLY : Hydro_sendOutputData
  use Simulation_interface, ONLY : Simulation_sendOutputData
  use Driver_interface, ONLY : Driver_sendOutputData
  use Particles_interface, ONLY : Particles_sendOutputData
  use Grid_interface, ONLY : Grid_sendOutputData
  use IO_interface, ONLY : IO_sendOutputData
implicit none

  call Grid_sendOutputData()

  call Driver_sendOutputData()

  call IO_sendOutputData()

  call Hydro_sendOutputData()

  call Particles_sendOutputData()

  !!call IO_outputScalars()
  

  !! add in other calls to output
  call Simulation_sendOutputData()

  
end subroutine IO_updateScalars
