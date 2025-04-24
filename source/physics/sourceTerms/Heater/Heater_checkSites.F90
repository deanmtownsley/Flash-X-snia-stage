!!****f* source/physics/sourceTerms/Heater/Heater_checkSites
!!
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

#include "Simulation.h"
#include "constants.h"

subroutine Heater_checkSites(solnData,del,lo,hi, lblock)

!----------------------------------------------------------------------------------------
   implicit none
   integer, intent(in) :: lblock
   real, pointer, dimension(:, :, :, :) :: solnData
   integer, dimension(MDIM),intent(IN)          :: lo, hi
   real,intent(IN)    :: del(MDIM)
end subroutine Heater_checkSites
