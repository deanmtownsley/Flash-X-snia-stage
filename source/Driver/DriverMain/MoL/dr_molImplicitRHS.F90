!!****if* source/Driver/DriverMain/MoL/dr_molImplicitRHS
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
!!      dr_molImplicitRHS
!!
!!  SYNOPSIS
!!
!!      call dr_molImplicitRHS(real,    intent(in) :: t,
!!                             integer, intent(in) :: activeRHS,
!!                             real,    intent(in) :: dtWeight)
!!
!!  DESCRIPTION
!!
!!      Calculate implicit RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      t          : Current time
!!      activeRHS : The RHS data struct to fill
!!      dtWeight  : Weighted timestep for the current stage (e.g. for flux corrections)
!!
!!***
subroutine dr_molImplicitRHS(t, activeRHS, dtWeight)
   use RadTrans_interface, only: RadTrans_molImplicitRHS
   use Simulation_interface, only: Simulation_molImplicitRHS

   implicit none

   real, intent(in) :: t
   integer, intent(in) :: activeRHS
   real, intent(in) :: dtWeight

   call RadTrans_molImplicitRHS(t, activeRHS, dtWeight)
   call Simulation_molImplicitRHS(t, activeRHS, dtWeight)
end subroutine dr_molImplicitRHS
