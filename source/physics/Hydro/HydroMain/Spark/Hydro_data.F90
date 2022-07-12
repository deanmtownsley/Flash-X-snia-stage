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

  real,  allocatable, target :: hy_pen(:,:)
  real,  allocatable,target :: hy_pflat(:)
  real,  allocatable,target :: hy_pgrv(:)
  real,  allocatable,target :: hy_pshck(:)
  real, allocatable, target :: hy_flux(:,:,:,:)
  real,  allocatable, target :: hy_snake(:,:,:,:)
  real,  target, allocatable :: hy_uPlus(:,:,:,:)
  real,  target, allocatable :: hy_uMinus(:,:,:,:)
  real,  allocatable :: hy_flat(:,:,:)
  real,  allocatable :: hy_flat3d(:,:,:)
  real,  allocatable :: hy_grv(:,:,:)
  real,  allocatable :: hy_shck(:,:,:)
  real, allocatable, target :: hy_tmpState(:,:,:,:)
  logical :: hydro_GPU_scratch = .False.

  real, dimension(MDIM) :: hy_del
  logical :: scratch_allocated
  real,allocatable :: hy_Vc(:,:,:)

  integer,  dimension(LOW:HIGH,MDIM) :: hy_dlim, hy_dlimGC
  integer :: hy_meshMe
  real :: hy_cfl
  logical :: hy_hydroComputeDtFirstCall
  logical :: hy_updateHydroFluxes
  real :: hy_dt, hy_dtmin
  integer :: hy_gcMaskSize
  integer :: hy_globalComm
  integer :: hy_meshNumProcs
  logical :: hy_restart
  logical :: hy_shockDetectOn
  logical :: hy_useTiling 
  real :: hy_smalldens, hy_smallE, hy_smallpres, hy_smallX, hy_smallu
  real, allocatable, dimension(:,:,:) :: hy_farea, hy_cvol
  real, allocatable, dimension(:) :: hy_xCenter, hy_xLeft, hy_xRight,hy_yCenter, hy_zCenter
  real, allocatable :: hy_mfrac(:), hy_eosData(:)
 
  !One block's worth of fluxes defined generally (not assuming fixed block size mode)
  real, allocatable, target :: hy_flx(:,:,:,:), hy_fly(:,:,:,:), hy_flz(:,:,:,:)
  
  !Flux buffers
  real, allocatable, dimension(:,:,:,:), target :: hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ

  real, allocatable :: hy_grav(:,:,:,:)  

  logical :: hy_useHydro
  logical :: hy_fluxCorrect, hy_fluxCorrectPerLevel
  integer, dimension(NFLUXES) :: hy_fluxCorVars
  integer :: hy_geometry

  
  logical :: hy_threadWithinBlock
  logical, dimension(NUNK_VARS) :: hy_gcMask
  ! Additional scratch storage for RK time stepping
  real, allocatable, target :: hy_starState(:,:,:,:)

  ! Limiter info
  real :: hy_limRad
  real :: hy_cvisc

  real :: hy_tiny=1.e-32
  real :: hy_gravConst, hy_4piGinv

  logical :: hy_hybridRiemann, hy_flattening

  real :: hy_C_hyp, hy_alphaGLM, hy_lChyp
  real :: hy_bref
  ! System of units used
  character(4) :: hy_units

  integer:: hy_sizex,hy_sizey,hy_sizez
  
end module Hydro_data


