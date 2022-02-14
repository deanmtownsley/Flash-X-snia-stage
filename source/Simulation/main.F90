!!****f* source/Simulation/main
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
!!  main
!!
!!
!! SYNOPSIS
!!
!!  N/A
!!
!!
!! DESCRIPTION
!!
!!  The source file main.F90 in the Simulation unit contains the Fortran
!!  PROGRAM. As such it can be considered the top-level "driver" of an application.
!!  By default it is set up to drive the simulation of a time-dependent
!!  problem by calling:
!!  - Driver_initAll  for initializations,
!!  - Driver_evolveAll  for managing the computation, and
!!  - Driver_finalizeAll  for cleaning up.
!!
!! SEE ALSO
!!
!!  Driver_initAll
!!  Driver_evolveAll
!!  Driver_finalizeAll
!!
!!
!!***

program Flashx

  use Driver_interface, ONLY : Driver_initParallel, Driver_initAll,&
       Driver_evolveAll, Driver_finalizeAll

  implicit none

  call Driver_initParallel()

  call Driver_initAll()

  call Driver_evolveAll( )

  call Driver_finalizeAll ( )
  

end program Flashx
