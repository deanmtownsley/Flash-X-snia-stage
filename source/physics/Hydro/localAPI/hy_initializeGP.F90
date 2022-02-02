Subroutine hy_counterGP_init(radiusGP, counterGP, blkLimitsGC)
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

  implicit none

#include "constants.h"

  !!-----Arguments------------------------------------------------
  real,    intent(IN)  :: radiusGP
  integer, intent(OUT) :: counterGP
  integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
  !!--------------------------------------------------------------
  counterGP = 0
end subroutine hy_counterGP_init



! *****************************************************************
Subroutine hy_initGP(RinvGP, WpGP, WmGP, blkLimitsGC)

! * Main init of GP interpolation. Here we set up one-time 
! * initialized arrays and matrices. This is the only place
! * where we invert R_mn(counter,counter). 
! *****************************************************************

  implicit none


  !!-----Arguments------------------------------------------------
  real, intent(INOUT), dimension(:,:) :: RinvGP
  real, intent(INOUT), dimension(:,:) :: WpGP
  real, intent(INOUT), dimension(:,:) :: WmGP
  integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
  !!--------------------------------------------------------------

end subroutine hy_initGP
