!!****if* source/Simulation/SimulationMain/EosGraph/Flash
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!! NAME
!!
!!  Flash
!!
!!
!! SYNOPSIS
!!
!!  N/A
!!
!!
!! DESCRIPTION
!!
!!  The source file Flash.F90 in the Simulation unit contains the Fortran
!!  PROGRAM. As such it can be considered the top-level "driver" of an application.
!!  By default it is set up to drive the simulation of a time-dependent
!!  problem by calling:
!!  - Driver_initFlash  for initializations,
!!  - Driver_evolveFlash  for managing the computation, and
!!  - Driver_finalizeFlash  for cleaning up.
!!
!! SEE ALSO
!!
!!  Driver_initFlash
!!  Driver_evolveFlash
!!  Driver_finalizeFlash
!!
!!***

program Flash

  use Driver_interface, ONLY : Driver_initParallel, Driver_initFlash,&
       Driver_evolveFlash

  implicit none

  call Driver_initParallel()

  call Driver_initFlash()

  call Driver_evolveFlash( )


end program Flash
