!!****if* source/Particles/localAPI/pt_dpdNonBondedForces
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
!!    pt_DPDNonBondedForces
!!
!! SYNOPSIS
!!    pt_dpdNonBondedForces( real,   INTENT(IN)     ::pos
!!                         real,   INTENT(IN)     ::v
!!                         real,   INTENT(IN)     ::btypes
!!                         real,   INTENT(IN)     ::parent 
!!                         real,   INTENT(IN)     ::parentType
!!                         real  , INTENT(OUT)    ::fvec
!!                         integer,INTENT(IN)     ::p_count)
!!
!! DESCRIPTION
!!
!!    Calculates the bonded and nonbonded forces between particles on each processor
!!    The Aij array holds the repulsion parameters between all beads types
!!    The property BDT_PART_PROP for each particle contains its bead type, 
!!    therefore, for calculating the force between the two particles i and j  
!!    The value stored in Aij(particle(BDT_PART_PROP,i),particle(BDT_PART_PROP,j)) will be used  
!!    to evaluate the forces between these two particles
!!
!!
!!
!! ARGUMENTS
!!
!!  particles :A real type array holding the local particles (real and virtual) on this processor
!! 
!!            for this version of the routine
!!
!! updateRefine : is true if the routine wished to retain the already
!!                initialized particles instead of reinitializing them
!!                as the grid refine.
!!
!!***

subroutine pt_dpdNonBondedForces(pos,v,btypes,parents,parentType, &
     internalIndex,fvec,p_count)
  
  implicit none
#include "Simulation.h"

  integer,INTENT(IN) :: p_count
  real,dimension(NDIM,p_count),INTENT(IN) :: pos,v
  real,dimension(p_count),INTENT(IN) :: btypes,parents,parentType,internalIndex
  real,dimension(NDIM,p_count),INTENT(OUT) :: fvec

  fvec = 0.0

  return 

end subroutine pt_dpdNonBondedForces
