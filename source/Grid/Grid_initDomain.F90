!!****f* source/Grid/Grid_initDomain
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
!!  Grid_initDomain
!!
!! SYNOPSIS
!!
!!  call Grid_initDomain(logical(IN)  :: restart,
!!                       logical(INOUT) :: particlesInitialized)
!!
!!
!! DESCRIPTION
!!  Initialize the Grid's data stuctures for the discretized mesh and
!!  apply the initial conditions.
!!
!! ARGUMENTS
!!  restart  : whether starting from initial conditions or from checkpoint
!!  particlesInitialized : is true if particle positions were initialized before returning
!!                         from this routine
!!
!!***

subroutine Grid_initDomain( restart,particlesInitialized)

implicit none
  logical, intent(IN) :: restart
  logical, intent(INOUT) :: particlesInitialized
  
end subroutine Grid_initDomain
