!!****if* source/physics/IncompNS/IncompNSMain/constDens/IncompNS_predictor
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
!!***
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_predictor(solnData, facexData, faceyData, facezData, del, lo, hi, dt)

   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_abort
   use Stencils_interface, ONLY: Stencils_integrateEuler, Stencils_integrateAB2
   use IncompNS_data

   implicit none

   !-----Argument-List-----!
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real, dimension(MDIM), intent(in) :: del
   integer,dimension(MDIM), intent(in) :: lo,hi
   real, INTENT(IN) :: dt

   !------------------------------------------------------------------------------------------
   integer, dimension(MDIM) :: hi1

   !
   call Timers_start("IncompNS_predictor")

   if (ins_intSchm /= 1 .and. ins_intSchm /= 2) then
      call Driver_abort("[IncompNS_predictor] ins_intSchm should be 1 or 2 for constant density configuration")
   end if

   !
   !------------------------------------------------------------------------------------------------------
   ! COMPUTE RIGHT HAND SIDE AND PREDICTOR STEP:
   ! ------- ----- ---- ---- --- --------- ----
   !------------------------------------------------------------------------------------------------------
   hi1(:) = hi(:)
   
   if (ins_intSchm == 1) then
      hi1(IAXIS)=hi1(IAXIS)+1
      call Stencils_integrateEuler(facexData(VELC_FACE_VAR, :, :, :), &
                                   facexData(HVN0_FACE_VAR, :, :, :), &
                                   dt,lo, hi1, &
                                   iSource=ins_prescoeff*facexData(PGN1_FACE_VAR, :, :, :) &
                                   + facexData(VFRC_FACE_VAR, :, :, :) &
                                   - ins_dpdx + ins_gravX)

      hi1(IAXIS)=hi1(IAXIS)-1; hi1(JAXIS)=hi1(JAXIS)+1
      call Stencils_integrateEuler(faceyData(VELC_FACE_VAR, :, :, :), &
                                   faceyData(HVN0_FACE_VAR, :, :, :), &
                                   dt, lo, hi1,&
                                   iSource=ins_prescoeff*faceyData(PGN1_FACE_VAR, :, :, :) &
                                   + faceyData(VFRC_FACE_VAR, :, :, :) &
                                   - ins_dpdy + ins_gravY)

#if NDIM == 3
      hi1(JAXIS)=hi1(JAXIS)-1; hi1(KAXIS)=hi1(KAXIS)+1
      call Stencils_integrateEuler(facezData(VELC_FACE_VAR, :, :, :), &
                                   facezData(HVN0_FACE_VAR, :, :, :), &
                                   dt, lo, hi1, &
                                   iSource=ins_prescoeff*facezData(PGN1_FACE_VAR, :, :, :) &
                                   + facezData(VFRC_FACE_VAR, :, :, :) &
                                   - ins_dpdz + ins_gravZ)
#endif

   else if (ins_intSchm == 2) then
      hi1(IAXIS)=hi1(IAXIS)+1
      call Stencils_integrateAB2(facexData(VELC_FACE_VAR, :, :, :), &
                                 facexData(HVN0_FACE_VAR, :, :, :), &
                                 facexData(HVN1_FACE_VAR, :, :, :), &
                                 dt, lo, hi1, &
                                 iSource=ins_prescoeff*facexData(PGN1_FACE_VAR, :, :, :) &
                                 + facexData(VFRC_FACE_VAR, :, :, :) &
                                 - ins_dpdx + ins_gravX)

      facexData(HVN1_FACE_VAR, :, :, :) = facexData(HVN0_FACE_VAR, :, :, :)
      hi1(IAXIS)=hi1(IAXIS)-1; hi1(JAXIS)=hi1(JAXIS)+1
      call Stencils_integrateAB2(faceyData(VELC_FACE_VAR, :, :, :), &
                                 faceyData(HVN0_FACE_VAR, :, :, :), &
                                 faceyData(HVN1_FACE_VAR, :, :, :), &
                                 dt, lo, hi1, &
                                 iSource=ins_prescoeff*faceyData(PGN1_FACE_VAR, :, :, :) &
                                 + faceyData(VFRC_FACE_VAR, :, :, :) &
                                 - ins_dpdy + ins_gravY)

      faceyData(HVN1_FACE_VAR, :, :, :) = faceyData(HVN0_FACE_VAR, :, :, :)

#if NDIM == 3
      hi1(JAXIS)=hi1(JAXIS)-1; hi1(KAXIS)=hi1(KAXIS)+1
      call Stencils_integrateAB2(facezData(VELC_FACE_VAR, :, :, :), &
                                 facezData(HVN0_FACE_VAR, :, :, :), &
                                 facezData(HVN1_FACE_VAR, :, :, :), &
                                 dt,lo, hi1, &
                                 iSource=ins_prescoeff*facezData(PGN1_FACE_VAR, :, :, :) &
                                 + facezData(VFRC_FACE_VAR, :, :, :) &
                                 - ins_dpdz + ins_gravZ)

      facezData(HVN1_FACE_VAR, :, :, :) = facezData(HVN0_FACE_VAR, :, :, :)
#endif

   end if

   call Timers_stop("IncompNS_predictor")

   return
end subroutine IncompNS_predictor
