!!****if* source/physics/Hydro/HydroMain/unsplit/hy_eigenValue
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
!!  hy_eigenValue
!!
!! SYNOPSIS
!!
!!  hy_eigenValue( real (OUT)          :: EigValue(HY_WAVENUM),
!!                     real (IN)           :: U_normal,
!!                     real (IN)           :: C_fast,
!!                     real (IN), optional :: C_alfn,
!!                     real (IN), optional :: C_slow,
!!                     real (OUT),optional :: C_hyp)
!!                     
!!
!! DESCRIPTION
!!
!!  This routine calculates MHD/Hydro eigenvalues.
!!
!! ARGUMENTS
!!
!!  EigValue   - Eigenvalues (wave speeds)
!!  U_normal   - Fluid velocity in normal direction
!!  C_fast     - Fast magnetoacoustic speed for MHD/Sound speed for Hydro
!!  C_alfn     - Alfven speed (needed for MHD only)
!!  C_slow     - Slow magnetoacoustic speed (needed for MHD only)
!!  C_hyp    - advection wave speed for GLM-MHD
!!
!!***


Subroutine hy_eigenValue(EigValue,U_normal,C_fast,C_alfn,C_slow,C_hyp)
    
  implicit none

#include "Simulation.h"
#include "UHD.h"

  !! Arguments type declaration --------------------------
  real,dimension(HY_WAVENUM), intent(OUT) :: EigValue
  real,intent(IN) :: U_normal,C_fast
  real,intent(IN), optional :: C_alfn,C_slow,C_hyp
  !! -----------------------------------------------------

  EigValue(HY_FASTLEFT) = U_normal-C_fast
  EigValue(HY_SLOWLEFT) = U_normal
  EigValue(HY_ENTROPY)  = U_normal
  EigValue(HY_SLOWRGHT) = U_normal
  EigValue(HY_FASTRGHT) = U_normal+C_fast

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  EigValue(HY_ALFNLEFT) = U_normal-C_alfn
  EigValue(HY_ALFNRGHT) = U_normal+C_alfn
  EigValue(HY_SLOWLEFT) = EigValue(HY_SLOWLEFT)-C_slow
  EigValue(HY_SLOWRGHT) = EigValue(HY_SLOWRGHT)+C_slow
#ifdef FLASH_UGLM_MHD
  EigValue(HY_GLMPLEFT) =-C_hyp
  EigValue(HY_GLMPRGHT) = C_hyp
#endif
#endif


End Subroutine hy_eigenValue
