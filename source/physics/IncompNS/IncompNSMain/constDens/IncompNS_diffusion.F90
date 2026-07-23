!!****if* source/physics/IncompNS/IncompNSMain/constDens/IncompNS_diffusion
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

subroutine IncompNS_diffusion(solnData, facexData, faceyData, facezData, del, lo, hi)

   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Stencils_interface, ONLY: Stencils_diffusion
   use IncompNS_data

!------------------------------------------------------------------------------------------
   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real, dimension(MDIM), intent(in) :: del
   integer,dimension(MDIM), intent(in) :: lo,hi
!------------------------------------------------------------------------------------------
   integer,dimension(MDIM) :: hi1   
   !
   call Timers_start("IncompNS_diffusion")
   !
   hi1(:)=hi(:)
   ! compute RHS of momentum equation
   hi1(IAXIS)=hi1(IAXIS)+1   
   call Stencils_diffusion(facexData(HVN0_FACE_VAR, :, :, :), &
                             facexData(VELC_FACE_VAR, :, :, :), &
                             del,&
                             ins_invReynolds, &
                             lo, hi1)

   hi1(IAXIS)=hi1(IAXIS)-1; hi1(JAXIS)=hi1(JAXIS)+1
   call Stencils_diffusion(faceyData(HVN0_FACE_VAR, :, :, :), &
                             faceyData(VELC_FACE_VAR, :, :, :), &
                             del,&
                             ins_invReynolds, &
                             lo, hi1)
#if NDIM == 3
   hi1(JAXIS)=hi1(JAXIS)-1; hi1(KAXIS)=hi1(KAXIS)+1
   call Stencils_diffusion(facezData(HVN0_FACE_VAR, :, :, :), &
                             facezData(VELC_FACE_VAR, :, :, :), &
                             del,&
                             ins_invReynolds, &
                             lo, hi1)

#endif

   call Timers_stop("IncompNS_diffusion")

   return
end subroutine IncompNS_diffusion
