!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_finishPotential
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
!!     Gravity_finishPotential
!!
!! SYNOPSIS
!!
!!  call Gravity_finishPotential(   optional,integer(IN) :: potentialIndex)
!!
!! DESCRIPTION
!!
!!      This routine computes the gravitational potential on all
!!      blocks specified in the list, for the gravity implementations
!!      (i.e., various Poisson implementations), which make use of it
!!      in computing the gravitational acceleration.
!!
!!      The resulting potential can be considered as zone-averaged.
!!
!!      Supported boundary conditions are isolated (0) and
!!      periodic (1).  The same boundary conditions are applied
!!      in all directions.  For some implementations of Gravity,
!!      additional combinations of boundary conditions may be supported.
!!
!! ARGUMENTS
!!
!!   potentialIndex : If present, determines which variable in UNK to use
!!                    for storing the updated potential.  If not present,
!!                    GPOT_VAR is assumed.
!!                    Presence or absense of this optional dummy argument
!!                    also determines whether some side effects are enabled
!!                    or disabled; see discussion of two modes under NOTES
!!                    below.
!!
!! NOTES
!!
!!  Gravity_finishPotential can operate in one of two modes:
!!  * automatic mode  - when called without the optional potentialIndex.
!!    Such a call will usually be made once per time step, usually
!!    from the main time advancement loop in Driver_evolveAll.
!!    Various side effects are enabled in this mode, see SIDE EFFECT below.
!!
!!  * explicit mode  - when called with the optional potentialIndex.
!!    The potential is stored in the variable explicitly given, and
!!    side effects like saving the previous potential in GPOL_VAR
!!    and updating some sink particle state and properties are
!!    suppressed.
!!
!! SIDE EFFECTS
!!
!!  Updates certain variables in permanent UNK storage to contain the
!!  gravitational potential.  Invokes a solver (of the Poisson equation)
!!  if necessary. On return, if potentialIndex is not present,
!!     GPOT_VAR:  contains potential for the current simulation time.
!!     GPOL_VAR (if defined): contains potential at the previous simulation time.
!!  On return, if potentialIndex is present, the UNK variable given by
!!  potentialIndex contains the newly computed potential.
!!
!!  May affect other variables related to particle properties if particles
!!  are included in the simulation.  In particular,
!!     PDEN_VAR (if defined): may get updated to the current density from
!!                particles if particles have mass.
!!
!!  There are additional side effects if sink particles are used.
!!  These effects happen by the call to Particles_sinkAccelGasOnSinksAndSinksOnGas,
!!  which may update sink particle properties and additional UNK variables that store
!!  accelerations. The calls are only made in automatic mode.
!!
!!  May modify certain variables used for intermediate results by the solvers
!!  invoked. The list of variables depends on the Gravity implementation.
!!  The following information is subject to change without notice.
!!  For the Multigrid implementation:
!!     ISLS_VAR (residual)
!!     ICOR_VAR (correction)
!!     IMGM_VAR (image mass)
!!     IMGP_VAR (image potential)
!!  For the Multipole implementation:
!!     (none)
!!
!!***

!!REORDER(4): solnVec

subroutine Gravity_finishPotential( potentialIndex)


  use Gravity_data, ONLY : grav_poisfact, grav_temporal_extrp, grav_boundary, &
       grav_unjunkPden, &
       useGravity, updateGravity, grv_meshComm
  use Driver_interface, ONLY : Driver_abort
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Particles_interface, ONLY: Particles_updateGridVar
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
       Grid_getTileIterator, Grid_releaseTileIterator, &
       Grid_notifySolnDataUpdate, &
       Grid_finalizePoisson
  use Grid_tile,     ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
  
#include "Flashx_mpi_implicitNone.fh"
#include "Simulation.h"
#include "constants.h"

  integer, intent(IN), optional :: potentialIndex


  real, POINTER, dimension(:,:,:,:) :: solnVec

  integer       :: ierr

  real          :: redshift=0, oldRedshift=0
  real          :: scaleFactor, oldScaleFactor
  real          :: invscale
  integer       :: bcTypes(6)
  real          :: bcValues(2,6) = 0.
  integer       :: density
  integer       :: newPotVar
  logical       :: saveLastPot
  type(Grid_tile_t)     :: tileDesc
  type(Grid_iterator_t) :: itor

  integer :: i, j, k

  nullify(solnVec)

  saveLastPot = (.NOT. present(potentialIndex))
  if (present(potentialIndex)) then
     newPotVar = potentialIndex
  else
     newPotVar = GPOT_VAR
  end if

!!$  call Cosmology_getRedshift(redshift)
!!$  call Cosmology_getOldRedshift(oldRedshift)
!!$
  scaleFactor = 1./(1.+redshift)
  oldScaleFactor = 1./(1.+oldRedshift)

  invscale = 1./scaleFactor**3

!=========================================================================

  if(.not.useGravity) return

  if(.not.updateGravity) return

  call Timers_start("gravity")

  invscale=grav_poisfact*invscale
  call Grid_notifySolnDataUpdate( (/newPotVar/) )
  call Grid_finalizePoisson (newPotVar, density, bcTypes, bcValues, &
       invscale)


! Un-junk PDEN if it exists and if requested.

#ifdef PDEN_VAR
  if (grav_unjunkPden) then
     density = PDEN_VAR
#ifdef DENS_VAR           
     call Driver_abort("[Gravity_finishPotential] Not tested I guess!")
     ! I (JO) added this line to create the iterator.  I don't know if it is
     ! correct or if this code has ever been called.
     call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnVec, CENTER)
        ! If tiling is used here, we probably need to write this as an
        ! explicit loop nest over the tile's indices
        solnVec(density,:,:,:) = solnVec(density,:,:,:) + &
             solnVec(DENS_VAR,:,:,:)
        call tileDesc%releaseDataPtr(solnVec, CENTER)
        call itor%next()
     enddo
     call Grid_releaseTileIterator(itor)
     if (density .NE. PDEN_VAR) call Grid_notifySolnDataUpdate( (/density/) )
#endif
  end if
#endif

!!$  if (.NOT. present(potentialIndex)) then
!!$     ! Compute acceleration of the sink particles caused by gas and vice versa
!!$     call Particles_sinkAccelGasOnSinksAndSinksOnGas()
!!$  end if
  
  
#ifdef USEBARS
  call MPI_Barrier (grv_meshComm, ierr)
#endif  
  call Timers_stop ("gravity")
  
  return
end subroutine Gravity_finishPotential
