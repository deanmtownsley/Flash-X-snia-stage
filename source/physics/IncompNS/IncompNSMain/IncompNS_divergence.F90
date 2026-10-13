!!****if* source/physics/IncompNS/IncompNSMain/IncompNS_divergence
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
!!
!!***
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_divergence(solnData, facexData, faceyData, facezData, del, lo, hi)

   use ins_interface, ONLY: ins_divergence
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use IncompNS_data

!-----------------------------------------------------------------------------------------
   implicit none

   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   integer, dimension(MDIM) :: lo, hi
   real,dimension(MDIM),intent(IN) :: del

!----------------------------------------------------------------------------------------

   ! compute divergence of intermediate velocities
   call ins_divergence(facexData(VELC_FACE_VAR, :, :, :), &
                       faceyData(VELC_FACE_VAR, :, :, :), &
                       facezData(VELC_FACE_VAR, :, :, :), &
                       lo( IAXIS), hi( IAXIS), &
                       lo( JAXIS), hi( JAXIS), &
                       lo( KAXIS), hi( KAXIS), &
                       del(DIR_X), del(DIR_Y), del(DIR_Z), &
                       solnData(DUST_VAR, :, :, :))

   ! Release pointers:

   call Timers_stop("IncompNS_divergence")

   return
end subroutine IncompNS_divergence
