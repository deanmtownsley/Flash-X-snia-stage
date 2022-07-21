#include "constants.h"
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
#include "Simulation.h"
#include "Spark.h"


module Hydro_data
  implicit none
  save


  real, target, allocatable :: hya_flux(:)
  real, target, allocatable :: hya_flat(:)
  real, target, allocatable :: hya_shck(:)
  real, target, allocatable :: hya_grv(:)

  real, target, allocatable :: hya_rope(:)
  real, target, allocatable :: hya_uPlus(:)
  real, target, allocatable :: hya_uMinus(:)
  real, allocatable, target :: hya_tmpState(:)
  real, allocatable, target :: hya_starState(:)
  real, allocatable, target :: hya_Vc(:)
  real, allocatable, target :: hya_grav(:)  
  real,  allocatable,target :: hya_flat3d(:)
  real, allocatable, target :: hya_flx(:), hya_fly(:), hya_flz(:)
  real, allocatable, dimension(:), target :: hya_fluxBufX, hya_fluxBufY, hya_fluxBufZ
  real, allocatable, dimension(:),target :: hya_farea, hya_cvol
  real, allocatable, dimension(:),target :: hya_xCenter, hya_xLeft, hya_xRight,hya_yCenter, hya_zCenter
  real, allocatable :: hy_mfrac(:), hy_eosData(:)
  !Flux buffers

  real, dimension(MDIM) :: hy_del

  real :: hy_dt, hy_dtmin
  logical :: hy_shockDetectOn
  logical :: hy_useTiling 
  real :: hy_smalldens, hy_smallE, hy_smallpres, hy_smallX, hy_smallu
 


  logical :: hy_fluxCorrect, hy_fluxCorrectPerLevel
  integer, dimension(NFLUXES) :: hy_fluxCorVars
  integer :: hy_geometry

  
  logical :: hy_threadWithinBlock
  logical, dimension(NUNK_VARS) :: hy_gcMask
  ! Additional scratch storage for RK time stepping

  ! Limiter info
  real :: hy_limRad
  real :: hy_cvisc

  real :: hy_tiny=1.e-32
  real :: hy_gravConst, hy_4piGinv

  logical :: hy_hybridRiemann, hy_flattening

  real :: hy_C_hyp, hy_alphaGLM, hy_lChyp
  real :: hy_bref
  ! System of units used

end module Hydro_data


