!!****if* source/physics/Hydro/HydroMain/Spark/prim2con
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
!!  hy_rk_prim2con
!!
!! SYNOPSIS
!!
!!  hy_rk_prim2con( real(IN)  :: V(HY_VARINUM3),
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

subroutine prim2con(V,CU)

  implicit none

#include "Simulation.h"
#include "Spark.h"

  !! Arguments type declaration -----------
  real ,dimension(HY_NUM_VARS), intent(IN)  :: V
  real ,dimension(HY_NUM_FLUX),  intent(OUT) :: CU
  !! --------------------------------------
  real  :: u2, B2
  !$omp declare target
  CU = 0.0

  u2 = dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))
  B2 = 0.
#ifdef SPARK_GLM
  B2 = dot_product(V(HY_MAGX:HY_MAGZ),V(HY_MAGX:HY_MAGZ))
  CU(HY_FMGX:HY_FMGZ) = V(HY_MAGX:HY_MAGZ)
  CU(HY_FPSI) = V(HY_PSIB)
#endif

  CU(HY_MASS) = V(HY_DENS)
  CU(HY_XMOM:HY_ZMOM) = V(HY_DENS)*V(HY_VELX:HY_VELZ)

  CU(HY_ENER) = 0.5*V(HY_DENS)*u2 + V(HY_RHOE) + 0.5*B2

end subroutine prim2con
