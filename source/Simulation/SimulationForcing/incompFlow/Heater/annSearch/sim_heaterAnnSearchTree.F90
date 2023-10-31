!!****if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterAnnSearchTree
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
!!***

#include "constants.h"
#include "Simulation.h"

subroutine sim_heaterAnnSearchTree(heater, queryPt, annElems, annIdx)

   use ANN_mod
   use ANN_types_mod
   use sim_heaterData, ONLY: sim_heaterType
   implicit none

   type(sim_heaterType), intent(IN)  :: heater
   integer, intent(IN) :: annElems
   real, dimension(:), target, intent(IN) :: queryPt ! query point
   integer, dimension(:), target, intent(OUT):: annIdx ! indices of nearest neighbors

   ! local variables
   real :: eps = 0.
   real, dimension(:), allocatable, target :: annDists

   allocate (annDists(annElems))
   call ann_kSearch(c_loc(queryPt), heater%dims, annElems, c_loc(annIdx), c_loc(annDists), eps, heater%kdTree)
   deallocate (annDists)

end subroutine sim_heaterAnnSearchTree
