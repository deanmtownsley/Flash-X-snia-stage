!!****f* source/Driver/Driver_superTimeStep
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
!!  Driver_superTimeStep
!!
!! SYNOPSIS
!!
!!  Driver_superTimeStep()
!!
!! DESCRIPTION
!!
!! This routine implements the super time steppping advancement algorithm
!! to overcome small diffusive time scales in explicit formulation
!!
!!  
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_superTimeStep(dt,nuSTS,nstepSTS,nstepTotalSTS,dt_subSTS)

  implicit none
  real, intent(IN)    :: dt,nuSTS
  integer, intent(IN) :: nstepSTS,nstepTotalSTS
  real, intent(OUT)   :: dt_subSTS

  return
end subroutine Driver_superTimeStep
