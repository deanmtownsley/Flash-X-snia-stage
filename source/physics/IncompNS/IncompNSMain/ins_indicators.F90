!!****if* source/physics/IncompNS/IncompNSMain/ins_indicators
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

subroutine ins_indicators(u,v,w,pres,divv,ix1,ix2,jy1,jy2,kz1,kz2,vecminaux,vecmaxaux)

#include "Simulation.h"
#include "constants.h"
  implicit none

  !--------Argument List--------!
  real, dimension(:,:,:), intent(in) :: u,v,w
  real, dimension(:,:,:), intent(in) :: pres,divv
  integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
  real, dimension(5), intent(inout) :: vecminaux, vecmaxaux

  !-------Local Variables--------!
  integer, parameter :: DIV=1, UVEL=2, VVEL=3, WVEL=4, PRS=5

  vecmaxaux(UVEL) = max(vecmaxaux(UVEL),maxval(u    (ix1:ix2+1, jy1:jy2  , kz1:kz2)))
  vecmaxaux(VVEL) = max(vecmaxaux(VVEL),maxval(v    (ix1:ix2  , jy1:jy2+1, kz1:kz2)))
  vecmaxaux(PRS)  = max(vecmaxaux(PRS) ,maxval(pres (ix1:ix2  , jy1:jy2  , kz1:kz2)))
  vecmaxaux(DIV)  = max(vecmaxaux(DIV) ,maxval(divv (ix1:ix2  , jy1:jy2  , kz1:kz2)))

  vecminaux(UVEL) = min(vecminaux(UVEL),minval(u    (ix1:ix2+1, jy1:jy2  , kz1:kz2)))
  vecminaux(VVEL) = min(vecminaux(VVEL),minval(v    (ix1:ix2  , jy1:jy2+1, kz1:kz2)))
  vecminaux(PRS)  = min(vecminaux(PRS) ,minval(pres (ix1:ix2  , jy1:jy2  , kz1:kz2)))
  vecminaux(DIV)  = min(vecminaux(DIV) ,minval(divv (ix1:ix2  , jy1:jy2  , kz1:kz2)))

#if NDIM == MDIM
  vecmaxaux(WVEL) = max(vecmaxaux(WVEL),maxval(w    (ix1:ix2  , jy1:jy2  , kz1:kz2+1)))
  vecminaux(WVEL) = min(vecminaux(WVEL),minval(w    (ix1:ix2  , jy1:jy2  , kz1:kz2+1)))
#else
  vecmaxaux(WVEL) = 0.
  vecminaux(WVEL) = 0.
#endif

  return
end subroutine ins_indicators
