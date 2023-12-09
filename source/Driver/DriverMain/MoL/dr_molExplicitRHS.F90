!!****if* source/Driver/DriverMain/MoL/dr_molExplicitRHS
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
!!      dr_molExplicitRHS
!!
!!  SYNOPSIS
!!
!!      call dr_molExplicitRHS(real,    intent(in) :: t,
!!                             integer, intent(in) :: activeRHS,
!!                             real,    intent(in) :: dtWeight)
!!
!!  DESCRIPTION
!!
!!      Calculate explicit RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      t          : Current time
!!      activeRHS : The RHS data struct to fill
!!      dtWeight  : Weighted timestep for the current stage (e.g. for flux corrections)
!!
!!***
subroutine dr_molExplicitRHS(t, activeRHS, dtWeight)
   use Hydro_interface, only: Hydro_molExplicitRHS
   use RadTrans_interface, only: RadTrans_molExplicitRHS
   use Simulation_interface, only: Simulation_molExplicitRHS

   implicit none

   real, intent(in) :: t
   integer, intent(in) :: activeRHS
   real, intent(in) :: dtWeight

   call Hydro_molExplicitRHS(t, activeRHS, dtWeight)
   call RadTrans_molExplicitRHS(t, activeRHS, dtWeight)
   call Simulation_molExplicitRHS(t, activeRHS, dtWeight)
end subroutine dr_molExplicitRHS
