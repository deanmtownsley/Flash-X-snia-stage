!!****if* source/physics/sourceTerms/Heater/localAPI/htr_annSearchTree
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

subroutine htr_annSearchTree(heater, queryPt, annElems, annIdx)
   use Heater_type, ONLY: Heater_type_t
   implicit none
   type(Heater_type_t), intent(IN)  :: heater
   integer, intent(IN) :: annElems
   real, dimension(:), target, intent(IN) :: queryPt ! query point
   integer, dimension(:), target, intent(OUT):: annIdx ! indices of nearest neighbors
end subroutine htr_annSearchTree
