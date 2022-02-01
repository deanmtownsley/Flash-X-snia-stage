!!****if* source/physics/Hydro/localAPI/hy_addBiermannBatteryTerms
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
!! NAME
!!
!!  hy_addBiermannBatteryTerms
!!
!!
!! SYNOPSIS
!!
!!  call hy_addBiermannBatteryTerms(integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(LOW:HIGH,MDIM),
!!                            integer(IN) :: ix,
!!                            integer(IN) :: iy,
!!                            integer(IN) :: iz,
!!                            real(IN)    :: Flux,
!!                            integer(IN) :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  Adds Biermann battery terms to total MHD fluxes. This is a stub.
!!
!!  We are solving induction equations
!!
!!  d  / Bx \    d /  0  \    d /  Ez \     d / -Ey \
!! ----| By | + ---| -Ez | + ---|  0  |  + ---|  Ex | = 0, 
!!  dt \ Bz /   dx \  Ey /   dy \ -Ex /    dz \  0  /
!!
!! dEner/dt + grad( E x B ) (mind the sign!)
!!
!!
!! ARGUMENTS
!!
!!  blockID     - a local blockID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section 
!!                of block with the guard cells 
!!  ix,iy,iz    - indices of the line along which the sweep is made
!!  Flux        - array containing MHD fluxes
!!  sweepDir    - direction of sweep
!!
!!***


Subroutine hy_addBiermannBatteryTerms(blockID,blkLimitsGC,ix,iy,iz,Flux,sweepDir)


  implicit none

#include "constants.h"
#include "UHD.h"

  !! Argument List ----------------------------------------------------------
  integer, INTENT(IN) :: blockID,ix,iy,iz
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
  real, dimension(HY_VARINUM), intent(INOUT) :: Flux
  integer, INTENT(IN) :: sweepDir
  !! ----------------------------------------------------------------------

  return
  
End Subroutine hy_addBiermannBatteryTerms

