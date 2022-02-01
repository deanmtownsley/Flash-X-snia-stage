!!****if* source/Particles/localAPI/pt_prepareEsti
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
!!  pt_prepareEsti
!!
!! SYNOPSIS
!!
!!  pt_prepareEsti(real(in) :: dtOld,
!!                    real(in) :: dtNew,
!!                    real(inout):: particles(:,p_count),
!!                    integer(in):: p_count,
!!                    integer(in):: ind)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the passive particle module.
!!  This portion cleans up and prepares for the next time step.
!!  
!!
!! ARGUMENTS
!!
!!   dtOld -- previous time interval
!!   dtNew -- current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind  -- index of this particle type into pt_typeInfo
!!  
!! PARAMETERS
!!
!!
!!***

!===============================================================================

subroutine pt_prepareEsti (dtOld,dtNew,particles,p_count, ind)
  
  implicit none
  
#include "Simulation.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles
  
  return
  
end subroutine pt_prepareEsti

