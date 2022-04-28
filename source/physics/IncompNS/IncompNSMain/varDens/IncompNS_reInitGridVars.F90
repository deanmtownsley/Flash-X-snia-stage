!!****if* source/physics/IncompNS/IncompNSMain/varDens/IncompNS_reInitGridVars
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
!! IncompNS_reInitGridVars
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "constants.h"
#include "IncompNS.h"
#include "Simulation.h"

subroutine IncompNS_reInitGridVars(tileDesc)

   use Grid_tile, ONLY: Grid_tile_t
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use IncompNS_data

!------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData

!------------------------------------------------------------------------------------------
   nullify (solnData, facexData, faceyData, facezData)

   call Timers_start("IncompNS_reInitGridVars")

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)
   solnData(VISC_VAR, :, :, :) = 1.
   facexData(RHOF_FACE_VAR, :, :, :) = 1.
   faceyData(RHOF_FACE_VAR, :, :, :) = 1.
   facexData(SIGM_FACE_VAR, :, :, :) = 0.
   faceyData(SIGM_FACE_VAR, :, :, :) = 0.
   facexData(HVN0_FACE_VAR, :, :, :) = 0.
   faceyData(HVN0_FACE_VAR, :, :, :) = 0.
   facexData(VFRC_FACE_VAR, :, :, :) = 0.
   faceyData(VFRC_FACE_VAR, :, :, :) = 0.
   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)
#if NDIM ==3
   call tileDesc%getDataPtr(facezData, FACEZ)
   facezData(RHOF_FACE_VAR, :, :, :) = 1.
   facezData(SIGM_FACE_VAR, :, :, :) = 0.
   facezData(HVN0_FACE_VAR, :, :, :) = 0.
   facezData(VFRC_FACE_VAR, :, :, :) = 0.
   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

   call Timers_stop("IncompNS_reInitGridVars")

   return
end subroutine IncompNS_reInitGridVars
