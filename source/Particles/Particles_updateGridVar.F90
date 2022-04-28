!!****f* source/Particles/Particles_updateGridVar
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
!!    Particles_updateGridVar
!!
!! SYNOPSIS
!!    Particles_updateGridVar(integer(IN)    :: partProp,
!!                            integer(IN)    :: varGrid,
!!                            integer(IN),optional    :: mode)
!!
!! DESCRIPTION
!!
!!    Updates the given grid variable with data from the given particle
!!    property.
!!
!! ARGUMENTS
!!               partProp:  Index of particle attribute to interpolate onto 
!!                          the mesh
!!               varGrid:   Index of gridded variable to receive interpolated
!!                          quantity
!!               mode:      (Optional) If zero (default), zero varGrid first;
!!                          if present and nonzero, do not zero varGrid first
!!                          but add data from particles to existing grid data.
!!
!! PARAMETERS
!! 
!!***
  
subroutine Particles_updateGridVar(partProp, varGrid, mode)

  implicit none

  integer, INTENT(in) :: partProp, varGrid
  integer, INTENT(in), optional :: mode

  return

end subroutine Particles_updateGridVar
