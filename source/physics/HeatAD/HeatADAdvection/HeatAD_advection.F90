!!***if* source/physics/HeatAD/HeatADAdvection/HeatAD_advection
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
!!
!!
!!***
!!REORDER(4): solnData,face[xyz]Data

#include "constants.h"
#include "HeatAD.h"
#include "Simulation.h"

subroutine HeatAD_advection(solnData, facexData, faceyData, facezData, del, lo, hi)

   use HeatAD_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Stencils_interface, ONLY: Stencils_advectWeno

!--------------------------------------------------------------------------------------------
   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real,dimension(MDIM),intent(IN) :: del
   integer, dimension(MDIM+1) :: face
   integer, dimension(MDIM),intent(IN) :: lo, hi
   integer, parameter :: center=1, facex=2,facey=3,facez=4

!---------------------------------------------------------------------------------------------


   face=0

   call Timers_start("HeatAD_advection")

   face(center)=1
   call Stencils_advectWeno(solnData(HTN0_VAR, :, :, :), &
                              solnData(TEMP_VAR, :, :, :), &
                              facexData(ht_iVelFVar, :, :, :), &
                              faceyData(ht_iVelFVar, :, :, :), &
                              facezData(ht_iVelFvar, :, :, :), &
                              del, lo, hi, face)


   call Timers_stop("HeatAD_advection")

end subroutine HeatAD_advection
