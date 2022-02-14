!!****if* source/physics/HeatAD/HeatADMain/ht_indicators
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
!!
!!
!!*****

subroutine ht_indicators(temp,ix1,ix2,jy1,jy2,kz1,kz2,minaux,maxaux)

#include "Simulation.h"

  implicit none

  !--------Argument List--------!
  real, dimension(:,:,:), intent(in) :: temp
  integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
  real, intent(inout) :: minaux, maxaux

  maxaux  = max(maxaux,maxval(temp (ix1:ix2  , jy1:jy2  , kz1:kz2)))
  minaux  = min(minaux,minval(temp (ix1:ix2  , jy1:jy2  , kz1:kz2)))

  return
end subroutine ht_indicators
