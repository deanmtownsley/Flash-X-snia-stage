!!****f* source/Simulation/Simulation_molFastRHS
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
!!      Simulation_molFastRHS
!!
!!  SYNOPSIS
!!
!!      call Simulation_molFastRHS(real,    intent(in) :: t,
!                                  integer, intent(in) :: activeRHS
!!                                 real,    intent(in) :: dtWeight)
!!
!!  DESCRIPTION
!!
!!      Calculate fast RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      t         : Current time
!!      activeRHS : RHS data struct to fill
!!      dtWeight  : Weighted timestep (e.g. for flux corrections)
!!
!!***
subroutine Simulation_molFastRHS(t, activeRHS, dtWeight)

   implicit none

   real, intent(in) :: t
   integer, intent(in) :: activeRHS
   real, intent(in) :: dtWeight

   return
end subroutine Simulation_molFastRHS
