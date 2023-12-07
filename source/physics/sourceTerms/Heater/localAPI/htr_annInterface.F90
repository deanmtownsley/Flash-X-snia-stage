!!****if* source/physics/sourceTerms/Heater/localAPI/htr_annInterface
!!
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
!!
!! SYNOPSIS
!!  htr_annInterface()
!!
!! DESCRIPTION
!!  This is an interface specific for heater geometry and specifications
!!
!!***

#include "constants.h"
#include "Simulation.h"

Module htr_annInterface

   implicit none

   interface
      subroutine htr_annBuildTree(heater)
         use Heater_type, ONLY: Heater_type_t
         type(Heater_type_t), intent(INOUT)  :: heater
      end subroutine htr_annBuildTree
   end interface

   interface
      subroutine htr_annSearchTree(heater, queryPt, annElems, annIdx)
         use Heater_type, ONLY: Heater_type_t
         type(Heater_type_t), intent(IN)  :: heater
         integer, intent(IN) :: annElems
         ! query point
         real, dimension(:), target, intent(IN) :: queryPt
         ! indices of nearest neighbors
         integer, dimension(:), target, intent(OUT):: annIdx
      end subroutine htr_annSearchTree
   end interface

End module htr_annInterface
