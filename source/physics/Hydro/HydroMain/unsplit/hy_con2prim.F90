!!****if* source/physics/Hydro/HydroMain/unsplit/hy_con2prim
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
!! NAME
!!
!!  hy_con2prim
!!
!! SYNOPSIS
!!
!!  hy_con2prim( real(IN)  :: CU(HY_VARINUM),
!!                   real(IN)  :: game,
!!                   real(OUT) :: V(HY_VARINUM))
!!
!! ARGUMENTS
!!
!! CU    - conservative variables
!! game  - gamma for internal energy
!! V     - primitive variables
!!
!! DESCRIPTION
!!
!!  This routine calculates conversion from conservative variables to 
!!  primitive variables.
!!
!!***

Subroutine hy_con2prim(CU,game,V)

  implicit none

#include "Simulation.h"
#include "UHD.h"

  !! Arguments type declaration -----------------
  real ,dimension(HY_VARINUM), intent(IN)  :: CU
  real, intent(IN) :: game
  real ,dimension(HY_VARINUM), intent(OUT) :: V
  !! --------------------------------------------

  real  :: u2,B2

  V(HY_DENS) = CU(HY_DENS)
  V(HY_VELX:HY_VELZ) = CU(HY_XMOM:HY_ZMOM)/CU(HY_DENS)
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  V(HY_MAGX:HY_MAGZ) = CU(HY_MAGX:HY_MAGZ)
#endif

  u2 = dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))
  B2 = 0.
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  B2 = dot_product(V(HY_MAGX:HY_MAGZ),V(HY_MAGX:HY_MAGZ))
#endif
  V(HY_PRES) = (game-1.)*(CU(HY_ENER)-.5 *CU(HY_DENS)*u2-.5 *B2)

End Subroutine hy_con2prim
