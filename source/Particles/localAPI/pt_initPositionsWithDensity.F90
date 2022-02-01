!!****if* source/Particles/localAPI/pt_initPositionsWithDensity
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
!!  pt_initPositionsWithDensity
!!
!! SYNOPSIS
!!
!!  call pt_initPositionsWithDensity(integer(in)  :: blockID,
!!                                   logical(OUT) :: success)
!!
!! DESCRIPTION
!!
!!    Initialize particle locations.  This version sets up passive tracer
!!      particles that are distributed randomly according to the gas density     
!!
!! ARGUMENTS
!!
!!   blockID : ID of block in current processor
!!  success:   returns .TRUE. if positions for all particles
!!             that should be assigned to this block have been
!!             successfully initialized.
!!
!! PARAMETERS
!!
!!  pt_numParticlesWanted:   integer  Approximate number of tracer particles to use
!!                                throughout domain ??
!!
!!***
!========================================================================

subroutine pt_initPositionsWithDensity (tileDesc,success)

  use Grid_tile, ONLY : Grid_tile_t

  implicit none

  type(Grid_tile_t), INTENT(in) :: tileDesc
  logical, INTENT(out) :: success

  success = .false.
  return

end subroutine pt_initPositionsWithDensity


