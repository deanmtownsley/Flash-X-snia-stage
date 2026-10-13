!!***if* source/physics/sourceTerms/Outlet/Outlet_setForcing
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
!!REORDER(4): solnData, face[xyz]Data

#include "Simulation.h"
#include "constants.h"

subroutine Outlet_setForcing(solnData, facexData, faceyData, facezData,&
                             xC, yC, zC, boundBox, del, lo, hi, dt)
  
  real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
  real, dimension(:),intent(IN) :: xC,yC,zC
  integer, dimension(MDIM),intent(IN)   :: lo, hi
  real,dimension(MDIM),intent(IN)    :: del
  real, intent(in) :: dt
  real, dimension(LOW:HIGH, 1:MDIM),intent(IN)    :: boundBox
end subroutine Outlet_setForcing
