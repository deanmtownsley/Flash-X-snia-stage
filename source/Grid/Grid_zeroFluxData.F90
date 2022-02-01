!!****f* source/Grid/Grid_zeroFluxData
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
!!  Grid_zeroFluxData
!!
!! SYNOPSIS
!!  call Grid_zeroFluxData()
!!
!! DESCRIPTION
!!  Request that the Grid unit zero all flux data managed by the unit.
!!
!!  This is only implemented, and should only be used, when the following
!!  two conditions apply:
!!  1. The Grid implementation is Amrex;
!!  2. level-wide fluxes are in use.
!!  Otherwise, a stub version that does nothing is in effect.
!!***

subroutine Grid_zeroFluxData
  implicit none

  RETURN
end subroutine Grid_zeroFluxData

