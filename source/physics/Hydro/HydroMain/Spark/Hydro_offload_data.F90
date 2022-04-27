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


module Hydro_offload_data
  implicit none
  save

  !$omp declare target (hy_flat3d,hy_uPlus,hy_uMinus,hy_snake,hy_tmpState)
  !$omp declare target (hy_flat,hy_grv,hy_shck,hy_flux)
  
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
end module Hydro_offload_data
