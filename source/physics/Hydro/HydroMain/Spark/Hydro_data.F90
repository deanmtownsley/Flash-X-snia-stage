!!****if* source/physics/Hydro/HydroMain/Spark/Hydro_data
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
!!    Hydro_data
!!
!!  SYNOPSIS
!!    use Hydro_data
!!
!!  DESCRIPTION
!!    Stores data for Spark Hydro
!!
!!  NOTES
!!
!!***
!!REORDER(4): hy_fluxBuf[XYZ]

#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

module Hydro_data

  implicit none
  save
  !$omp declare target (hy_dlim,hy_starState)
  !$omp declare target (hy_tiny,hy_hybridRiemann,hy_geometry)
  !$omp declare target (hy_C_hyp,hy_smalldens, hy_smallE, hy_smallpres, hy_smallX,hy_cvisc,hy_del)
  !$omp declare target (hy_flx,hy_fly,hy_flz,hy_alphaGLM,hy_grav)

  ! Added for GPU
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

end module Hydro_data

