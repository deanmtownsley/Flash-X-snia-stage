!!****if* source/physics/IncompNS/IncompNSMain/constDens/IncompNS_advection
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

subroutine IncompNS_advection(solnData, facexData, faceyData, facezData, del, lo, hi)

!!$   use Grid_tile, ONLY: Grid_tile_t
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_abort
   use Stencils_interface, ONLY: Stencils_advectWeno, Stencils_advectCentral
   use IncompNS_data

!------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
!!$   type(Grid_tile_t), intent(in) :: tileDesc

   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real, dimension(MDIM), intent(in) :: del
   integer,dimension(MDIM), intent(in) :: lo,hi
!------------------------------------------------------------------------------------------
   integer,dimension(MDIM) :: hi1
   integer, dimension(MDIM+1) :: face
   integer, parameter :: center=1, facex=2, facey=3, facez=4
!------------------------------------------------------------------------------------------

   call Timers_start("IncompNS_advection")

   if (ins_advSchm /= 2) then
      call Driver_abort("[IncompNS_advection] ins_intSchm should be 2 for constant density configuration")
   end if
   hi1(:) = hi(:)
   !
   ! compute RHS of momentum equation
   face=0
   hi1(IAXIS)=hi1(IAXIS)+1
   face(facex)=1
   call Stencils_advectCentral(facexData(HVN0_FACE_VAR, :, :, :), &
                                 facexData(VELC_FACE_VAR, :, :, :), &
                                 facexData(VELC_FACE_VAR, :, :, :), &
                                 faceyData(VELC_FACE_VAR, :, :, :), &
                                 facezData(VELC_FACE_VAR, :, :, :), &
                                 del, lo, hi1, face)

   hi1(IAXIS)=hi1(IAXIS)-1; hi1(JAXIS)=hi1(JAXIS)+1
   face(facex)=0; face(facey)=1
   call Stencils_advectCentral(faceyData(HVN0_FACE_VAR, :, :, :), &
                                 faceyData(VELC_FACE_VAR, :, :, :), &
                                 facexData(VELC_FACE_VAR, :, :, :), &
                                 faceyData(VELC_FACE_VAR, :, :, :), &
                                 facezData(VELC_FACE_VAR, :, :, :), &
                                 del, lo, hi1, face)

#if NDIM==3
   hi1(JAXIS)=hi1(JAXIS)-1; hi1(KAXIS)=hi1(KAXIS)+1
   face(facey)=0; face(facez)=1
   call Stencils_advectCentral(facezData(HVN0_FACE_VAR, :, :, :), &
                                 facezData(VELC_FACE_VAR, :, :, :), &
                                 facexData(VELC_FACE_VAR, :, :, :), &
                                 faceyData(VELC_FACE_VAR, :, :, :), &
                                 facezData(VELC_FACE_VAR, :, :, :), &
                                 del, lo, hi1, face)

#endif

   call Timers_stop("IncompNS_advection")

   return
end subroutine IncompNS_advection
