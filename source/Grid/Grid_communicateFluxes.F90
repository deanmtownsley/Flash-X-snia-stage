!!****f* source/Grid/Grid_communicateFluxes
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
!!  Grid_communicateFluxes
!!
!! SYNOPSIS
!!  call Grid_communicateFluxes(integer(IN) :: axis,
!!                              integer(IN) :: coarse_level)
!!
!! DESCRIPTION
!!
!!  This call launches the communication phase for flux correction.
!!
!!  The communication happens entirely within SPFS, and does not
!!  involve any data structures that should be directly visible to
!!  client code. Instead, client code uses subroutines like
!!  Grid_putFluxData_block to transfer flux data to SPFS before a
!!  communication phase, and Grid_correctFluxData to get data back
!!  from SPFS after the communication.
!!
!!  The coarse_level argument selects a pair of refinement levels
!!  whose blocks may participate in the flux correction.  The
!!  coarse_level argument identifies the refinement level on the
!!  coarse size of the fine/coarse boundaries for which data is
!!  exchanged.
!!
!!  It is assumed that before calling this routine, the code has
!!  already stored the correct fluxes for all relevant blocks at the
!!  finer level into the semi-permanent flux storage, using
!!  Grid_putFluxData or a variant thereof.
!!
!!  This routine will comunicate the data necessary, including
!!  communication between MPI ranks, so that the flux data on the
!!  coarse side of fine/coarse boundaries can be overwritten by
!!  presumably more accurate values from computations on the fine
!!  side.
!!
!! ARGUMENTS
!!  axis - direction: IAXIS, JAXIS, KAXIS, or ALLDIR.
!!  coarse_level - the 1-based level index of the coarse blocks,
!!                 or the special value UNSPEC_LEVEL.
!!
!! NOTES
!!
!!   SPFS means semi-permanent flux storage. When using a Grid
!!   implementation based on AMReX, SPFS is implemented by an AMReX
!!   flux register class, such as FlashFluxRegister.
!!
!!   The symbol UNSPEC_LEVEL is defined in constants.h.
!!
!!   The special value UNSPEC_LEVEL should be taken to mean
!!   "communicate between all levels". It is not supported
!!   by the AMReX implementation.
!!
!!  This routine should only get called when Flash-X is configured to
!!  NOT use level-wide flux arrays, i.e., when levelFlux source directories
!!  are not included and the NoLevelwideFluxes setup variable is False.
!!
!! SEE ALSO
!!
!!  Grid_putFluxData
!!  Grid_putFluxData_block
!!  Grid_correctFluxData
!!
!!***

recursive subroutine Grid_communicateFluxes(axis, coarse_level)

    implicit none

    integer, intent(IN)                   :: axis
    integer, intent(IN)                   :: coarse_level

end subroutine Grid_communicateFluxes

