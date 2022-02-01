!!****if* source/Particles/localAPI/pt_utComputeError
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
!!  pt_utComputeError
!!
!! SYNOPSIS
!!
!!  pt_utComputeError(real(in) :: dtOld,
!!                           real(in) :: dtNew,
!!                           real(in) :: t)
!!
!! DESCRIPTION
!!
!!  Compute error, that is, compare actual to analytic solution.
!!
!!  Compares particles' POS{X,Y,Z}_PART_PROP and POSANAL{X,Y,Z}_PART_PROP
!!  properties.
!!
!! ARGUMENTS
!!
!!   dtOld -- previous time increment
!!   dtNew -- current time increment
!!   t     -- time at which solutions are compared
!!  
!!***

!===============================================================================

subroutine pt_utComputeError (dtOld,dtNew,t)
    
  
  implicit none

  
  real, INTENT(in)  :: dtOld, dtNew, t
     
  return
!!------------------------------------------------------------------------------
  
end subroutine pt_utComputeError


