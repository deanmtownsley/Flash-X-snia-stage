!!****if* source/Driver/DriverMain/MoL/dr_molFastRHS
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
!!  NAME
!!
!!      dr_molFastRHS
!!
!!  SYNOPSIS
!!
!!      call dr_molFastRHS(real, intent(in) :: t)
!!
!!  DESCRIPTION
!!
!!      Calculate fast RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      t : Current time
!!
!!***
subroutine dr_molFastRHS(t)
   ! use Spacetime_interface,  only: Spacetime_molFastRHS
   use Hydro_interface, only: Hydro_molFastRHS
   use RadTrans_interface, only: RadTrans_molFastRHS
   use Simulation_interface, only: Simulation_molFastRHS

   implicit none

   real, intent(in) :: t

   ! call Spacetime_molFastRHS(t)
   call Hydro_molFastRHS(t)
   call RadTrans_molFastRHS(t)
   call Simulation_molFastRHS(t)
end subroutine dr_molFastRHS
