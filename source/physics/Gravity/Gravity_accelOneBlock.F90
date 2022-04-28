!!****f* source/physics/Gravity/Gravity_accelOneBlock
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
!!  Gravity_accelOneBlock
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneBlock(integer, intent(in) :: tileDesc, 
!!                        integer, intent(in) :: ngcellcomp
!!                        real(:,:,:,:)),intent(out) :: gvec, 
!!                        integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration for this block.  Include ngcell layers outside
!!  block interior.
!!
!!  This routine computes the gravitational acceleration for
!!  zones in a given block. First-order
!!  finite-volume differencing is used everywhere.  It is assumed
!!  here that the requisite number of guard cells have peen appropriately
!!  filled for the variable containting the gravitational potential.
!!
!!  Dean Townsley 2008
!!  Contributed to Flash Center at the University of Chicago 2008
!!  Modified ..., 2017    Flash Center
!!
!! ARGUMENTS
!!
!!  tileDesc           -  metadata for the block, such as its integer bounds, refinement level etc
!!  gvec(:,:,:,:)   -  Array to receive gravitational acceleration
!!                        as as NDIM-dimensional vector.  It is assumed
!!                        the the space provided is the size of the block
!!                        plus all guard cells.  The fist index is the vector
!!                        component and the latter are cell indices.
!!  ngcellcomp         -  Number of layers outside of block interior to
!!                        compute gravity
!!  potentialIndex     -  if specified,  Variable # to take as potential.
!!                        Default is GPOT_VAR for the potential stored in the
!!                        gpot slot of unk, which should correspond to the
!!                        potential at the current timestep.
!!
!!
!!***


subroutine Gravity_accelOneBlock ( tileDesc, ngcellcomp, gvec, potentialIndex)

  use Grid_tile, ONLY : Grid_tile_t
  implicit none

#include "Simulation.h"
#include "constants.h"

  type(Grid_tile_t), intent(in)   :: tileDesc
  integer, intent(in)                  :: ngcellcomp
  real, dimension(:,:,:,:),intent(out) :: gvec
  integer, intent(in),optional         :: potentialIndex

  gvec(:,:,:,:) = 0.0

  return
   
end subroutine Gravity_accelOneBlock
