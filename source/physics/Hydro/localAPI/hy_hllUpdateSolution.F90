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

Subroutine hy_hllUpdateSolution(tileLimits, Uin, plo, Uout, flX, flY, flZ, loFl, del, dt)
  implicit none
  
  integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
  integer, intent(IN)  :: plo(*)
  real,    intent(INOUT):: Uin(plo(1):,plo(2):,plo(3):,plo(4):)
  real,    intent(OUT) :: Uout(plo(1):,plo(2):,plo(3):,plo(4):)
  integer, intent(IN)  :: loFl(*)
  real,    intent(IN)  :: flX(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
  real,    intent(IN)  :: flY(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
  real,    intent(IN)  :: flZ(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
  real,    intent(IN)  :: del(1:MDIM)
  real,    intent(IN)  :: dt
  TARGET :: Uin,Uout
    
  Uout(:, :, :, :) = 0.0
End Subroutine hy_hllUpdateSolution
