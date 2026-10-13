!!****if* source/physics/IncompNS/IncompNSMain/varDens/IncompNS_corrector
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

subroutine IncompNS_corrector(solnData, facexData, faceyData, facezData, del, lo, hi, dt)

   use ins_interface, ONLY: ins_corrector_vardens
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use IncompNS_data

   implicit none
   !-----Argument List-----
   real, dimension(:,:,:,:), pointer :: solnData, facexData, faceyData, facezData
   real,dimension(MDIM),intent(IN) :: del
   integer, dimension(MDIM),intent(IN) :: lo, hi
   real, INTENT(IN) :: dt
   

   call Timers_start("IncompNS_corrector")

   ! update divergence-free velocities (not on block boundary)
   call ins_corrector_vardens(facexData(VELC_FACE_VAR, :, :, :), &
                              faceyData(VELC_FACE_VAR, :, :, :), &
                              facezData(VELC_FACE_VAR, :, :, :), &
                              !------------------------------!
                              facexData(SIGM_FACE_VAR, :, :, :), &
                              faceyData(SIGM_FACE_VAR, :, :, :), &
                              facezData(SIGM_FACE_VAR, :, :, :), &
                              !------------------------------!
                              facexData(PGN1_FACE_VAR, :, :, :), &
                              faceyData(PGN1_FACE_VAR, :, :, :), &
                              facezData(PGN1_FACE_VAR, :, :, :), &
                              !------------------------------!
                              facexData(PGN2_FACE_VAR, :, :, :), &
                              faceyData(PGN2_FACE_VAR, :, :, :), &
                              facezData(PGN2_FACE_VAR, :, :, :), &
                              !------------------------------!
                              facexData(RHOF_FACE_VAR, :, :, :), &
                              faceyData(RHOF_FACE_VAR, :, :, :), &
                              facezData(RHOF_FACE_VAR, :, :, :), &
                              !------------------------------!
                              solnData(PRES_VAR, :, :, :), &
                              ins_rhoGas, dt, del(DIR_X), del(DIR_Y), del(DIR_Z), &
                              lo(IAXIS), hi(IAXIS), &
                              lo(JAXIS), hi(JAXIS), &
                              lo(KAXIS), hi(KAXIS))


   call Timers_stop("IncompNS_corrector")

   return
end subroutine IncompNS_corrector
