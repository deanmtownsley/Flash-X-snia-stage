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
  real, allocatable, dimension(:),target :: hya_xCenter, hya_xLeft, hya_xRight, hya_yCenter, hya_zCenter
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

  !! moved out from static physics routines to make it easier for data packet generation

  integer, dimension(MDIM) :: gCells
  integer, dimension(LOW:HIGH,MDIM,NDIM,MAXSTAGE) :: klim,lim1,lim,limgc
  integer :: dir,stage

  real :: hy_cfl
  logical :: hy_hydroComputeDtFirstCall
  logical :: hy_updateHydroFluxes
  integer :: hy_gcMaskSize
  logical :: hy_restart
  logical :: hy_useHydro, hy_telescoping
  integer :: hy_meshMe, hy_globalComm, hy_meshComm, hy_meshNumProcs

  ! System of units used

  !$omp declare target to &
  !$omp ( hy_cfl, &
  !$omp   hy_hydroComputeDtFirstCall, &
  !$omp   hy_updateHydroFluxes, &
  !$omp   hy_gcMaskSize, &
  !$omp   hy_restart, &
  !$omp   hy_useHydro, hy_telescoping, &
  !$omp   hy_meshMe, hy_globalComm, hy_meshComm, hy_meshNumProcshya_starState, hya_tmpState, &
  !$omp   hya_uPlus, hya_uMinus, hya_Vc, hya_grav, hya_flat3d, hya_flat, hya_grv,&
  !$omp   hya_rope, hya_flux, hya_shck, &
  !$omp   hya_flx, hya_fly, hya_flz, hya_fluxBufX, hya_fluxBufY, hya_fluxBufZ, &
  !$omp   hya_farea, hya_cvol, hya_xCenter, hya_xLeft, hya_xRight, hya_yCenter, hya_zCenter, &
  !$omp   hy_mfrac, hy_eosData, &
  !$omp   hy_del, &
  !$omp   hy_dt, hy_dtmin, &
  !$omp   hy_shockDetectOn, &
  !$omp   hy_useTiling , &
  !$omp   hy_smalldens, hy_smallE, hy_smallpres, hy_smallX, hy_smallu, &
  !$omp   hy_fluxCorrect, hy_fluxCorrectPerLevel, &
  !$omp   hy_fluxCorVars, &
  !$omp   hy_geometry, &
  !$omp   hy_threadWithinBlock, &
  !$omp   hy_gcMask, &
  !$omp   hy_limRad, &
  !$omp   hy_cvisc, &
  !$omp   hy_tiny, &
  !$omp   hy_gravConst, hy_4piGinv, &
  !$omp   hy_hybridRiemann, hy_flattening, &
  !$omp   hy_C_hyp, hy_alphaGLM, hy_lChyp, &
  !$omp   hy_bref, &
  !$omp   klim,lim,limgc,lim1,gCells, dir, stage)

end module Hydro_data


