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

module Hydro_inhost_data
  implicit none
  save
  real, dimension(NFLUXES) :: hy_flux
  real,  allocatable, target :: hy_pen(:,:)
  real,  dimension(HY_NRECON) :: hy_uPlus
  real,  dimension(HY_NRECON) :: hy_uMinus
  real,  allocatable,target :: hy_pflat(:)
  real,  allocatable,target :: hy_flat3d(:,:,:)
  real,  allocatable,target :: hy_pgrv(:)
  real,  allocatable,target :: hy_pshck(:)
  real, allocatable, target :: hy_tmpState(:,:,:,:)
end module Hydro_inhost_data
