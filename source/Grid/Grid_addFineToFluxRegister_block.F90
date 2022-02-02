!!****f* source/Grid/Grid_addFineToFluxRegister_block
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
!!  Grid_addFineToFluxRegister_block
!!
!! SYNOPSIS
!!  call Grid_addFineToFluxRegister(Grid_tile_t(IN) :: blockDesc,
!!                        optional, logical(IN)     :: isDensity(:),
!!                        optional, real(IN)        :: coefficient,
!!                        optional, logical(IN)     :: zeroFullRegister)
!!
!! DESCRIPTION 
!!  Each flux register is associated with a fine and a coarse level.  In normal
!!  use, client code could add flux data from both levels into the flux register
!!  for use with adjusting flux data on the coarse level.
!!
!!  This routine allows client code to request that the Grid unit add fine data
!!  from a single block in the Grid unit's flux data structures to the
!!  contents of the associated flux registers that are associated with that same
!!  block.  This routine is clearly intended for use with AMR.  Note that the
!!  flux registers may choose to only store flux data that exists at 
!!  fine/coarse boundaries.
!!
!!  All data stored in the Grid unit's flux data structures as flux densities
!!  will automatically be transformed to flux before applying to the flux
!!  register.
!!
!!  Additionally, a multiplicative scale factor may be applied to all flux data before
!!  passing the data to the flux register.
!!
!!  It is assumed that before calling this routine, the client code has already
!!  written flux data for the given block to the Grid's data structures using
!!  the Grid_tile_t's getDataPtr routine.
!!
!!  IMPORTANT: This routine should *only* be used on blocks and not 
!!             when tiling is enabled.
!!
!! ARGUMENTS
!!  blockDesc - the block whose flux data should be stored
!!  isDensity - a mask that identifies which physical flux quantities are
!!              actually stored in the Grid unit's flux data structures as
!!              flux densities.  If no mask is given, it is assumed that data
!!              is stored as flux.
!!              !DEV: CURRENTLY THE isDensity ARGUMENT IS NOT SUPPORTED.
!!  coefficient - a scaling parameter to apply to all flux data before applying
!!                the data to the flux register.
!!              CURRENTLY THE coefficient ARGUMENT IS IGNORED BY THE IMPLEMENTATION.
!!  zeroFullRegister - zero the current fine and coarse data for the block
!!                     in the register before adding the indicated flux data to
!!                     the register.  If this parameter is not given, then the
!!                     current data is not zeroed.
!!
!! NOTES
!!
!!   DEV: This interface may be obsolete.
!!   DEV: Currently only implemented for the Amrex Grid.
!!
!! SEE ALSO
!!   Grid_addCoarseToFluxRegister
!!   Grid_overwriteFluxes
!!
!!***

subroutine Grid_addFineToFluxRegister_block(blockDesc, isDensity, coefficient, &
                                           zeroFullRegister)
  use Driver_interface, ONLY : Driver_abort
  use Grid_tile,        ONLY : Grid_tile_t

  implicit none

  type(Grid_tile_t), intent(IN)           :: blockDesc
  logical,           intent(IN), optional :: isDensity(:)
  real,              intent(IN), optional :: coefficient
  logical,           intent(IN), optional :: zeroFullRegister

  call Driver_abort("[Grid_addFineToFluxRegister_block] Prototype stub.  Do NOT use!")
end subroutine Grid_addFineToFluxRegister_block

