!!****if* source/Driver/DriverMain/Driver_initNumericalTools
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
!!  Driver_initNumericalTools
!!
!! SYNOPSIS
!!
!!  Driver_initNumericalTools ()
!!
!! DESCRIPTION
!!
!!  Initializes all numerical tool Units by calling their respective
!!  initialization routines viz. Roots_init, RungeKutta_init, etc.
!!  
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine Driver_initNumericalTools ()

  use Roots_interface,      ONLY:  Roots_init
  use RungeKutta_interface, ONLY:  RungeKutta_init
  use Stencils_interface,   ONLY:  Stencils_init

  implicit none

  call Roots_init      ()
  call RungeKutta_init ()
  call Stencils_init   ()

end subroutine Driver_initNumericalTools
