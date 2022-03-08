!!****if* source/Grid/GridMain/AMR/Amrex/Grid_initDomain
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
!!  Grid_initDomain(logical(IN)    :: restart,
!!                  logical(INOUT) :: particlesInitialized)
!!
!! DESCRIPTION
!!  Create the coarsest mesh, initialize all the mesh data structures,
!!  apply initial conditions, and run EoS on interior and guardcells to make
!!  them thermodynamically consistent.
!!
!!  User-defined refinement critera is applied to determine the 
!!  blocks that require refinement.  All new child blocks are filled
!!  with the initial conditions, guardells are filled for cell-centered
!!  variables, and EoS is run on interiors and guardcells.
!!
!!  In simulations with particles, under certain conditions particle
!!  positions will also be initialized.  Currently this is the case
!!  if and only if the runtime parameter refine_on_particle_count is
!!  true.
!!
!! ARGUMENTS
!!  restart : is true if the execution is starting from a checkpoint
!!            file, otherwise false.
!!  particlesInitialized : is true if particle positions were initialized before returning
!!                         from this routine
!!
!! NOTES
!!  When restarting from a checkpoint file, block interiors are assumed to
!!  have been filled when this interface is called. The EOS is not called on
!!  the block interiors in this implementation for use with Paramesh. It is
!!  assumed that data is already thermodynamically consistent, because
!!  that is how data are written to checkpoint files.
!!
!!***

#include "Simulation.h"
#include "constants.h"

subroutine Grid_initDomain(restart,particlesInitialized)
  use amrex_fort_module,    ONLY : wp => amrex_real
  use amrex_amr_module,     ONLY : amrex_init_from_scratch, &
                                   amrex_max_level, amrex_geom, amrex_ref_ratio
  use gr_amrexInterface, ONLY : gr_fillPhysicalBC, gr_preinterpolationWork, &
                              gr_postinterpolationWork
  use amrex_fillpatch_module,    ONLY : amrex_fillpatch

  use Eos_interface, ONLY : Eos_wrapped
  use Grid_data,            ONLY : gr_doFluxCorrection
  use gr_physicalMultifabs, ONLY : unk, &
                                   gr_scratchCtr, &
                                   facevarx, facevary, facevarz, &
                                   fluxes, &
                                   flux_registers
  use Grid_data, ONLY : gr_globalNumBlocks, gr_globalDomain, gr_interpolator, &
                      lo_bc_amrex, hi_bc_amrex, gr_eosModeInit, &
                      gr_maxRefine, gr_doFluxCorrection

  use Driver_interface,     ONLY : Driver_abort
  use Grid_iterator,             ONLY : Grid_iterator_t
  use Grid_tile,                 ONLY : Grid_tile_t
  use Grid_interface,            ONLY : Grid_getTileIterator, &
                                        Grid_releaseTileIterator
  implicit none
  
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
  real(wp), contiguous, pointer :: dp(:,:,:,:)

  logical, intent(IN)    :: restart
  logical, intent(INOUT) :: particlesInitialized
  real :: time
  integer :: i
  real(wp), parameter :: T_INIT = 0.0_wp
  !!!!!----- ALLOCATE DATA STRUCTURES
  ! multifabs 
  !
  ! NOTE: We implement these with the 0-based level indexing scheme native to
  ! AMReX instead of the 1-based level indexing scheme of FLASH.
  !   => all code dealing with multifabs arrays must consider the need for 
  !      index translation

  ! DEV: TODO Implement parameters
  if (.NOT. restart) then
    allocate(gr_scratchCtr(0:amrex_max_level))

#if NFLUXES > 0
# ifdef USE_LEVELWIDE_FLUXES
    allocate(fluxes(0:amrex_max_level, 1:NDIM))
# endif

    ! Flux registers
    !
    ! By definition, each flux register is associated with two levels 
    ! (a coarse level and the next finest level).  Following the AMReX 
    ! interface and the AMReX tutorials, the 0-based level index used here
    ! is associated with the level index of the fine level 
    !    (e.g. flux_register(1) is the flux register for the coarse level 0
    !          and the fine level 1).
    if (gr_doFluxCorrection) then
      allocate(flux_registers(1:amrex_max_level))
    end if
#endif
    !  This creates all refinement levels needed based on ICs,
    !  runs EoS on interiors, fills GCs, and runs EoS on GCs.
    !  All this is done through the callback functions.
    allocate(unk     (0:amrex_max_level))

#if NFACE_VARS > 0
    !call Driver_abort("[Grid_initDomain] Face-centered variables not yet supported with AMReX Grid implementation")
    allocate(facevarx(0:amrex_max_level))
#if NDIM >= 2
    allocate(facevary(0:amrex_max_level))
#endif
#if NDIM == 3
    allocate(facevarz(0:amrex_max_level))
#endif
#endif

    ! DEV: T_INIT is set to zero now as AMReX does not care about this value.
    ! It might be useful to pass non-zero T_INIT for time-dependent refinement criteria;
    ! this might be time a part of \FlashOfTheFuture with simulated time passed as T_INIT
    ! that is read from checkpoint file. Useful with user-defined callbacks
    call amrex_init_from_scratch(T_INIT)
  else ! this is a restart 
    call Simulation_initRestart()
  end if

  !Initialize for Grid sovlers
    call gr_solversInit()

  particlesInitialized = .FALSE.
end subroutine Grid_initDomain

