!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/Drive_evolveAll
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
!!  Drive_evolveAll
!!
!! SYNOPSIS
!!
!!  Drive_evolveAll ()
!!
!! DESCRIPTION
!!
!!  Simple stripped down version for testing single units.
!!
!! NOTES
!!
!!***

subroutine Drive_evolveAll ()

  implicit none

  real :: told, tnew

  call cpu_time (told)

  call sim_RungeKuttaTest ()

  call cpu_time (tnew)

  write (*,*) ' Cpu time = ',tnew - told

  return
end subroutine Drive_evolveAll
