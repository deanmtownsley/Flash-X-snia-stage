!!***if* source/physics/sourceTerms/Outlet/OutletMain/phaseForcing/Outlet_setForcing
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
!!***
!!REORDER(4): solnData, face[xyz]Data

#include "Simulation.h"
#include "constants.h"

subroutine Outlet_setForcing(solnData, facexData, faceyData, facezData,&
                             xC, yC,zC, boundBox, del, lo, hi, dt)
  
  use Outlet_data, ONLY: out_sink, out_flag, &
       out_buffer, out_growthRate, &
       out_QAuxLiq, out_QAuxGas, out_volAuxLiq, &
       out_volAuxGas, out_QOutLiq, out_QOutGas, out_velRefScale, &
       out_meshMe, out_xMin, out_xMax, out_yMin, out_yMax
#if NDIM == MDIM
  use Outlet_data, ONLY: out_zMin, out_zMax
#endif
  
  use out_interface, ONLY: out_lsDamping, out_velFrcPhased
  
  use IncompNS_data, ONLY: ins_gravX, ins_gravY, ins_gravZ
  use IncompNS_interface, ONLY: IncompNS_setVectorProp
  use Timers_interface, ONLY: Timers_start, Timers_stop
  
  implicit none
  !----------------------------------------------------------------------------------------
  real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
  real, dimension(:),intent(IN) :: xC,yC,zC
  integer, dimension(MDIM),intent(IN)   :: lo, hi
  real,dimension(MDIM),intent(IN)    :: del
  real, intent(in) :: dt
  real, dimension(LOW:HIGH, 1:MDIM),intent(IN)    :: boundBox
  
  integer :: ierr

!----------------------------------------------------------------------------------------

   call Timers_start("Outlet_setForcing")


#if NDIM < MDIM
   call out_lsDamping(solnData(DFRC_VAR, :, :, :), &
                      solnData(DFUN_VAR, :, :, :), &
                      xC, yC, zC, boundBox, &
                      dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                      lo(IAXIS), hi(IAXIS), &
                      lo(JAXIS), hi(JAXIS), &
                      lo(KAXIS), hi(KAXIS), &
                      out_flag, out_sink, out_buffer, &
                      out_growthRate, &
                     out_xMin, out_xMax, out_yMin, out_yMax, 0., 0.)

   call out_velFrcPhased(facexData(VELC_FACE_VAR, :, :, :), &
                         facexData(VFRC_FACE_VAR, :, :, :), &
                         facexData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xC-del(IAXIS)/2, yC, zC, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         lo(IAXIS), hi(IAXIS)+1, &
                         lo(JAXIS), hi(JAXIS), &
                         lo(KAXIS), hi(KAXIS), &
                         out_xMin, out_xMax, out_yMin, out_yMax, 0., 0., &
                         out_flag, out_buffer, out_growthRate, &
                         IAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

   call out_velFrcPhased(faceyData(VELC_FACE_VAR, :, :, :), &
                         faceyData(VFRC_FACE_VAR, :, :, :), &
                         faceyData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xC, yC-del(JAXIS)/2, zC, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         lo(IAXIS), hi(IAXIS), &
                         lo(JAXIS), hi(JAXIS)+1, &
                         lo(KAXIS), hi(KAXIS), &
                         out_xMin, out_xMax, out_yMin, out_yMax, 0., 0., &
                         out_flag, out_buffer, out_growthRate, &
                         JAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

#else

   call out_lsDamping(solnData(DFRC_VAR, :, :, :), &
                      solnData(DFUN_VAR, :, :, :), &
                      xC, yC, zC, boundBox, &
                      dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                      lo(IAXIS), hi(IAXIS), &
                      lo(JAXIS), hi(JAXIS), &
                      lo(KAXIS), hi(KAXIS), &
                      out_flag, out_sink, out_buffer, &
                      out_growthRate, &
                      out_xMin, out_xMax, out_yMin, out_yMax, out_zMin, out_zMax)

   call out_velFrcPhased(facexData(VELC_FACE_VAR, :, :, :), &
                         facexData(VFRC_FACE_VAR, :, :, :), &
                         facexData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xC-del(IAXIS)/2, yC, zC, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         lo(IAXIS), hi(IAXIS)+1, &
                         lo(JAXIS), hi(JAXIS), &
                         lo(KAXIS), hi(KAXIS), &
                         out_xMin, out_xMax, out_yMin, out_yMax, out_zMin, out_zMax, &
                         out_flag, out_buffer, out_growthRate, &
                         IAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

   call out_velFrcPhased(faceyData(VELC_FACE_VAR, :, :, :), &
                         faceyData(VFRC_FACE_VAR, :, :, :), &
                         faceyData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xC, yC-del(JAXIS)/2, zC, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         lo(IAXIS), hi(IAXIS), &
                         lo(JAXIS), hi(JAXIS)+1, &
                         lo(KAXIS), hi(KAXIS), &
                         out_xMin, out_xMax, out_yMin, out_yMax, out_zMin, out_zMax, &
                         out_flag, out_buffer, out_growthRate, &
                         JAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

   call out_velFrcPhased(facezData(VELC_FACE_VAR, :, :, :), &
                         facezData(VFRC_FACE_VAR, :, :, :), &
                         facezData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xC, yC, zC-del(KAXIS)/2, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         lo(IAXIS), hi(IAXIS), &
                         lo(JAXIS), hi(JAXIS), &
                         lo(KAXIS), hi(KAXIS)+1, &
                         out_xMin, out_xMax, out_yMin, out_yMax, out_zMin, out_zMax, &
                         out_flag, out_buffer, out_growthRate, &
                         KAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

#endif

   ! Release pointers:

   call Timers_stop("Outlet_setForcing")

   return

end subroutine Outlet_setForcing
