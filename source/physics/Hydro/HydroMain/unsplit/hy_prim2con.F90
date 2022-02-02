!!****if* source/physics/Hydro/HydroMain/unsplit/hy_prim2con
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
!!  hy_prim2con
!!
!! SYNOPSIS
!!
!!  hy_prim2con( real(IN)  :: V(HY_VARINUM2),
!!                   real(OUT) :: CU(HY_VARINUM))
!!
!! ARGUMENTS
!!
!! V  - primitive variables + GAMC,GAME
!! CU - conservative variables 
!!
!! DESCRIPTION
!!
!!  This routine calculates conversions from primitive variables to conservative variables.
!!
!!***

Subroutine hy_prim2con(V,CU)

  implicit none

#include "Simulation.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  real ,dimension(HY_VARINUM2), intent(IN)  :: V
  real ,dimension(HY_VARINUM),  intent(OUT) :: CU
  !! --------------------------------------

  real  :: u2,B2

  u2 = dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))
  B2 = 0.
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  B2 = dot_product(V(HY_MAGX:HY_MAGZ),V(HY_MAGX:HY_MAGZ))
#endif

  CU(HY_DENS) = V(HY_DENS)
  CU(HY_XMOM:HY_ZMOM) = V(HY_DENS)*V(HY_VELX:HY_VELZ)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  CU(HY_MAGX:HY_MAGZ) = V(HY_MAGX:HY_MAGZ)
#endif

#ifdef FLASH_UGLM_MHD
  CU(HY_GLMP) = V(HY_GLMP)      ! GLM is currently unsupported
#endif

  CU(HY_ENER) = 0.5*V(HY_DENS)*u2 + V(HY_PRES)/(V(HY_GAME)-1.) + 0.5*B2

End subroutine hy_prim2con
