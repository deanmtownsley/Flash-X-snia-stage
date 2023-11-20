!!****if* source/Grid/GridMain/AMR/Paramesh4/Incomp/Grid_updateRefinement
!!
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  Grid_updateRefinement
!!
!!
!! SYNOPSIS
!!
!!  call Grid_updateRefinement(integer(IN) :: nstep,
!!                             real(IN)    :: time,
!!                    OPTIONAL,logical(OUT):: gridChanged)
!!
!!
!! DESCRIPTION
!!
!!  Apply the user-defined refinment critera to determine which blocks need
!!  to be refined and derefined.  Once the blocks are marked, call
!!  gr_updateGridRefinement to carry out the rest of the routine.
!!  The internal routine does the refinements (amr_refine_derefine)
!!  During this
!!  stage, the blocks are redistributed across processors (if needed).
!!  After the refinement, the newly created child blocks are filled via
!!  prolongation from the coarse parents.  This prolongation step can use
!!  prolongation routines from paramesh, or defined by the user
!!  Once the prolongation is done, the guardcells are filled.  Finally, the
!!  EOS is called on the block interiors to make them thermodynamically
!!  consistent. The internal routine also calls Particles_updateRefinement to
!!  move the particles to the correct block after the grid refines.
!!
!!
!!
!! ARGUMENTS
!!
!!  nstep : current step number
!!  time  : current evolution time
!!  gridChanged : returns TRUE if grid may actually have changed.
!!
!!***

subroutine Grid_updateRefinement(nstep, time, gridChanged)

   use Grid_data, ONLY: gr_nrefs, gr_meshMe
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Grid_interface, ONLY: Grid_markRefineDerefine
   use gr_interface, ONLY: gr_updateRefinement
   use paramesh_interfaces, ONLY: amr_restrict
   use Logfile_interface, ONLY: Logfile_open, Logfile_close
   use Simulation_interface, ONLY: Simulation_wantsRebalance
   use tree, ONLY: grid_changed
   use physicaldata, only: interp_mask_facex, interp_mask_facey, interp_mask_facez, &
                           interp_mask_facex_res, interp_mask_facey_res, interp_mask_facez_res

!!$  use tree, ONLY : newchild, lnblocks
!!$  use paramesh_interfaces, ONLY : amr_refine_derefine, &
!!$                                  amr_prolong

   implicit none

#include "constants.h"
#include "Simulation.h"

   integer, intent(in) :: nstep
   real, intent(in) :: time
   logical, intent(out), OPTIONAL :: gridChanged

!!$  integer :: i

!!$  integer :: count

   integer :: iopt, iempty, logUnit
   logical :: logUnitLocal = .true.
   logical :: force_rebalance

   !=============================================================================

   !Before we do any refinement, store a globalID made up of the block's cornerID,
   !refinement level and nodetype.  This is necessary to refine and derefine the
   !particles.  If particles are not included, this routine will be a stub
   call gr_ptFillBlkParticleInfo()

   !Do this to initialize the value of grid_changed before any refinement is done.
   !If refinement is applied grid_changed becomes 1.
   !grid_changed is used in other parts of the code as a reference and therefore this intialization
   !is neccessary to make sure correct logic is propagated.
   grid_changed = 0

   ! We only consider refinements every nrefs timesteps, or if Simulation
   ! is requesting a rebalance.
   force_rebalance = Simulation_wantsRebalance(nstep, time)
   if ((mod(nstep, gr_nrefs) == 0) .or. force_rebalance) then

      call Timers_start("tree")  !1 of 2

      ! Make restriction of velocities and rhs linear (conservative)
      interp_mask_facex(VELC_FACE_VAR) = 2
      interp_mask_facex_res(VELC_FACE_VAR) = 1
      interp_mask_facey(VELC_FACE_VAR) = 2
      interp_mask_facey_res(VELC_FACE_VAR) = 1
      interp_mask_facez(VELC_FACE_VAR) = 2
      interp_mask_facez_res(VELC_FACE_VAR) = 1

      interp_mask_facex(HVN1_FACE_VAR) = 2
      interp_mask_facex_res(HVN1_FACE_VAR) = 1
      interp_mask_facey(HVN1_FACE_VAR) = 2
      interp_mask_facey_res(HVN1_FACE_VAR) = 1
      interp_mask_facez(HVN1_FACE_VAR) = 2
      interp_mask_facez_res(HVN1_FACE_VAR) = 1

      iopt = 1; iempty = 1
      call amr_restrict(gr_meshMe, iopt, iempty)

      call Timers_start("markRefineDerefine")
      call Grid_markRefineDerefine()
      call Timers_stop("markRefineDerefine")

      call Timers_stop("tree")  !1 of 2 (We restart in gr_updateRefinement)
      !internal routine that does the actual amr refinement and
      !other housekeeping

      call gr_updateRefinement(gridChanged, force_rebalance)

   else
      if (present(gridChanged)) gridChanged = .FALSE.
   end if

   return
end subroutine Grid_updateRefinement
