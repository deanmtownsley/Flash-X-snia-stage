!!****f* source/Grid/Grid_addFineToFluxRegister
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
!!  Grid_addFineToFluxRegister
!!
!! SYNOPSIS
!!  call Grid_addFineToFluxRegister(integer(IN) :: fine_level,
!!                        optional, logical(IN) :: isDensity(:),
!!                        optional, real(IN)    :: coefficient,
!!                        optional, logical(IN) :: zeroFullRegister)
!!
!! DESCRIPTION 
!!  Each flux register is associated with a fine and a coarse level.  In normal
!!  use, client code could add flux data from both levels into the flux register
!!  for use with adjusting flux data on the coarse level.
!!
!!  This routine allows client code to request that the Grid unit add fine data
!!  from the Grid unit's flux data structures to the contents of the associated
!!  flux registers.  This routine is clearly intended for use with AMR.  Note
!!  that the flux registers may choose to only store flux data that exists at 
!!  fine/coarse boundaries.
!!
!!  All data stored in the Grid unit's flux data structures as flux densities
!!  will automatically be transformed to flux before applying to the flux
!!  register.
!!
!!  Additionally, a multiple scale factor may be applied to all flux data before
!!  passing the data to the flux register.
!!
!!  It is assumed that before calling this routine, the client code has already
!!  written flux data to Grid's data structures using the Grid_getFluxPtr
!!  interface.
!!
!! ARGUMENTS
!!  fine_level - the 1-based level index (1 is the coarsest level) indicating
!!               which level's data should be added to the flux register as
!!               fine data.
!!  isDensity - a mask that identifies which physical flux quantities are
!!              actually stored in the Grid unit's flux data structures as
!!              flux densities.  If no mask is given, it is assumed that data
!!              is stored as flux.
!!              !DEV: CURRENTLY THE isDensity ARGUMENT IS NOT SUPPORTED.
!!  coefficient - a scaling parameter to apply to all flux data before applying
!!                the data to the flux register.
!!  zeroFullRegister - zero the current fine and coarse data in the register
!!                     before adding the indicated flux data to the register.
!!                     If this parameter is not given, then the current data is
!!                     not zeroed.
!!
!! NOTES
!!
!!   DEV: This interface may be obsolete.
!!   DEV: Currently only implemented for the Amrex Grid.
!!
!!   This interface is called to implement Grid_putFluxData, which is only
!!   used in configurations with per-level fluxes, for the Amrex Grid.
!!
!! SEE ALSO
!!   Grid_getFluxPtr/Grid_releaseFluxPtr
!!   Grid_zeroFluxRegister
!!   Grid_addCoarseToFluxRegister
!!   Grid_overwriteFluxes
!!
!!***


subroutine Grid_addFineToFluxRegister(fine_level, isDensity, coefficient, &
                                      zeroFullRegister)
  use Driver_interface, ONLY : Driver_abort

  
  implicit none

  integer, intent(IN)           :: fine_level
  logical, intent(IN), optional :: isDensity(:)
  real,    intent(IN), optional :: coefficient
  logical, intent(IN), optional :: zeroFullRegister

  call Driver_abort("[Grid_addFineToFluxRegister] Prototype stub.  Do NOT use!")
end subroutine Grid_addFineToFluxRegister

