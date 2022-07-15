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

  
  real, target :: hy_flux(NFLUXES,GRID_IHI_GC+2,GRID_JHI_GC+2*K2D,GRID_KHI_GC+2*K3D)
  real,  target :: hy_rope(NRECON,GRID_IHI_GC+2,GRID_JHI_GC+2*K2D,GRID_KHI_GC+2*K3D)
  real,  target :: hy_uPlus(NRECON,GRID_IHI_GC+2,GRID_JHI_GC+2*K2D,GRID_KHI_GC+2*K3D)
  real,  target :: hy_uMinus(NRECON,GRID_IHI_GC+2,GRID_JHI_GC+2*K2D,GRID_KHI_GC+2*K3D)
  real :: hy_flat(GRID_IHI_GC+2,GRID_JHI_GC+2*K2D,GRID_KHI_GC+2*K3D)
  real :: hy_shck(GRID_IHI_GC+2,GRID_JHI_GC+2*K2D,GRID_KHI_GC+2*K3D)
  real :: hy_grv(GRID_IHI_GC+2,GRID_JHI_GC+2*K2D,GRID_KHI_GC+2*K3D)

  
  real, allocatable, target :: hy_tmpState(:,:,:,:)
  real, allocatable, dimension(:,:,:) :: hy_farea, hy_cvol
  real, allocatable, dimension(:) :: hy_xCenter, hy_xLeft, hy_xRight,hy_yCenter, hy_zCenter
  real, allocatable :: hy_mfrac(:), hy_eosData(:)
  real,  allocatable :: hy_flat3d(:,:,:)
  real, allocatable, target :: hy_flx(:,:,:,:), hy_fly(:,:,:,:), hy_flz(:,:,:,:)
  real, allocatable :: hy_Vc(:,:,:)
  !Flux buffers
  real, allocatable, dimension(:,:,:,:), target :: hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ
  real, allocatable :: hy_grav(:,:,:,:)  
  real, allocatable, target :: hy_starState(:,:,:,:)

  logical :: hydro_GPU_scratch = .False.

  real, dimension(MDIM) :: hy_del
  logical :: scratch_allocated

  integer,  dimension(LOW:HIGH,MDIM) :: hy_dlim, hy_dlimGC
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


