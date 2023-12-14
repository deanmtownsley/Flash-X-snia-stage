!!****if* source/physics/RadTrans/RadTransMain/RadTrans_init
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
!!  NAME 
!!
!!  RadTrans_init
!!
!!  SYNOPSIS
!!
!!  call RadTrans_init()
!!
!!  DESCRIPTION 
!!    Initialize radiative transfer unit
!!
!! ARGUMENTS
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine RadTrans_init()

  use RadTrans_data
  use rt_interface, ONLY : rt_init
  use Driver_interface, ONLY : Driver_getMype, Driver_getComm, &
    Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  implicit none

  call Driver_getMype(MESH_COMM,rt_meshMe)
  call Driver_getNumProcs(MESH_COMM,rt_meshNumProcs)

  call RuntimeParameters_get ("useRadTrans", rt_useRadTrans)
  call RuntimeParameters_get ("gr_useTiling", rt_enableTiling)

  call RuntimeParameters_get("geometry",rt_str_geometry)
  call RuntimeParameters_mapStrToInt(rt_str_geometry, rt_geometry)

  call RuntimeParameters_get("cfl",rt_cfl)

  rt_eosModeGc = MODE_DENS_EI

  rt_gcMask = .FALSE.
  rt_gcMask(DENS_VAR) = .TRUE.
  rt_gcMask(VELX_VAR) = .TRUE.
  rt_gcMask(VELY_VAR) = .TRUE.
  rt_gcMask(VELZ_VAR) = .TRUE.
  rt_gcMask(ENER_VAR) = .TRUE.
#ifdef YE_MSCALAR
  rt_gcMask(YE_MSCALAR) = .TRUE.
#endif
#ifdef SUMY_MSCALAR
  rt_gcMask(SUMY_MSCALAR) = .TRUE.
#endif

  call rt_init

  return

end subroutine RadTrans_init
