!!****if* source/Driver/DriverMain/MoL/dr_molPreEvolve
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  NAME
!!
!!      dr_molPreEvolve
!!
!!  SYNOPSIS
!!
!!      call dr_molPreEvolve(real, intent(in) :: t)
!!
!!  DESCRIPTION
!!
!!      Perform any pre-evolution work that must occur after all *_init
!!      and intBlock calls (e.g. setting evolved variables from primitives)
!!
!!
!!  ARGUMENTS
!!
!!      t  : Current time
!!
!!***
subroutine dr_molPreEvolve(t)
   use Hydro_interface, only: Hydro_molPreEvolve
   use RadTrans_interface, only: RadTrans_molPreEvolve
   use Simulation_interface, only: Simulation_molPreEvolve

   implicit none

   real, intent(in) :: t

   call Hydro_molPreEvolve(t)
   call RadTrans_molPreEvolve(t)
   call Simulation_molPreEvolve(t)
end subroutine dr_molPreEvolve
