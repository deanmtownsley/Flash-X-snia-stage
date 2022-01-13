!!****f* source/Grid/Grid_conserveFluxes
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!! NAME
!!  Grid_conserveFluxes
!!
!! SYNOPSIS
!!
!!  call Grid_conserveFluxes(integer(IN) :: axis,
!!                           integer(IN) :: coarse_level)
!!  
!! DESCRIPTION 
!!  
!!  Flux conservation is necessary when 2 blocks of differing
!!  levels (meaning having different grid spacings) border 
!!  one another. 
!!  
!!  This routine can perform flux conservation on the finest
!!  blocks, the most typical usage for the Paramesh Grid or on
!!  blocks of a certain level.
!!  
!!  The routine overwrites the flux arrays maintained by the Grid
!!  
!! ARGUMENTS 
!!
!!
!!  axis - conserve fluxes in just one direction if 
!!         IAXIS, JAXIS, KAXIS, or in all directions if ALLDIR.
!!         These constants are defined in constants.h.
!!
!!  coarse_level - refinement level. Selects the level (coarse level) for
!!          which fluxes are updated.
!!          Can be UNSPEC_LEVEL for all levels (except, as an
!!          optimizing shortcut, the highest possible one).
!!
!! NOTES
!!
!!  Some implementations might ignore the level arguments, and always
!!  act as if it were UNSPEC_LEVEL.
!!
!!***

subroutine Grid_conserveFluxes(axis, coarse_level)
  implicit none
  
  integer, intent(IN)                   :: axis
  integer, intent(IN)                   :: coarse_level
end subroutine Grid_conserveFluxes

