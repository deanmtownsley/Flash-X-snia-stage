!!****if* source/Grid/GridMain/AMR/Paramesh4/Incomp/Grid_markRefineDerefine
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
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!
!! DESCRIPTION
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub. The AMReX-based Grid
!!  currently does not understand this way of implementing
!!  refinement criteria either, and uses callbacks instead,
!!  which are implemented as private functions of the Grid unit.
!!
!!  With the PARAMESH-based Grid implementation,
!!  this routine is normally called by the implementation of
!!  Grid_updateRefinement. It may also get called repeatedly
!!  during the initial construction of the Grid from
!!  Grid_initDomain.
!!
!! ARGUMENTS
!!
!!  none
!!
!! SEE ALSO
!!
!!  Grid_updateRefinement
!!  Grid_initDomain
!!  gr_expandDomain
!!
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

   use Grid_data, ONLY: gr_refine_cutoff, gr_derefine_cutoff, &
                        gr_refine_filter, &
                        gr_numRefineVars, gr_refine_var, gr_refineOnParticleCount, &
                        gr_enforceMaxRefinement, gr_maxRefine, &
                        gr_lrefineMaxByTime, &
                        gr_lrefineMaxRedDoByTime, &
                        gr_lrefineMaxRedDoByLogR, &
                        gr_lrefineCenterI, gr_lrefineCenterJ, gr_lrefineCenterK, &
                        gr_eosModeNow
   use tree, ONLY: newchild, refine, derefine, stay, nodetype, lrefine_max
   use Logfile_interface, ONLY: Logfile_stampVarMask
   use Grid_interface, ONLY: Grid_fillGuardCells
   use gr_interface, ONLY: gr_markRefineDerefine
   implicit none

#include "constants.h"
#include "Simulation.h"

   real :: ref_cut, deref_cut, ref_filter
   integer       :: l, i, iref
   logical, save :: gcMaskArgsLogged = .FALSE.
   integer, save :: eosModeLast = 0
   logical :: doEos = .true.
   integer, parameter :: maskSize = NUNK_VARS
   logical, dimension(maskSize) :: gcMask
   real, dimension(MAXBLOCKS) :: err

   if (gr_lrefineMaxRedDoByTime) then
      call gr_markDerefineByTime()
   end if

   if (gr_lrefineMaxByTime) then
      call gr_setMaxRefineByTime()
   end if

   if (gr_eosModeNow .NE. eosModeLast) then
      gcMaskArgsLogged = .FALSE.
      eosModeLast = gr_eosModeNow
   end if

   ! that are implemented in this file need values in guardcells

   gcMask = .false.
   do i = 1, gr_numRefineVars
      iref = gr_refine_var(i)
      if (iref > 0) gcMask(iref) = .TRUE.
   end do

   if (.NOT. gcMaskArgsLogged) then
      call Logfile_stampVarMask(gcMask, .true., '[Grid_markRefineDerefine]', 'gcArgs')
   end if

!!$  force_consistency = .FALSE.
   call Grid_fillGuardCells(CENTER, ALLDIR, doEos=.true., &
                            maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true., doLogMask=.NOT. gcMaskArgsLogged, &
                            selectBlockType=ACTIVE_BLKS)
   gcMaskArgsLogged = .TRUE.
!!$  force_consistency = .TRUE.

   newchild(:) = .FALSE.
   refine(:) = .FALSE.
   derefine(:) = .FALSE.
   stay(:) = .FALSE.

   do l = 1, gr_numRefineVars
      iref = gr_refine_var(l)
      ref_cut = gr_refine_cutoff(l)
      deref_cut = gr_derefine_cutoff(l)
      ref_filter = gr_refine_filter(l)
      err(:) = 0.0

#ifdef DFUN_VAR
      if (iref == DFUN_VAR) then
         call gr_markVarBounds(iref, min(ref_cut, deref_cut), max(ref_cut, deref_cut), lrefine_max)

      else
#endif
         call gr_estimateError(err, iref, ref_filter)
         call gr_markRefineDerefine(err, ref_cut, deref_cut)

#ifdef DFUN_VAR
      end if
#endif

   end do

   if (gr_refineOnParticleCount) call gr_ptMarkRefineDerefine()

   if (gr_enforceMaxRefinement) call gr_enforceMaxRefine(gr_maxRefine)

   if (gr_lrefineMaxRedDoByLogR) &
      call gr_unmarkRefineByLogRadius(gr_lrefineCenterI, &
                                      gr_lrefineCenterJ, gr_lrefineCenterK)

   ! When the flag arrays are passed to Paramesh for processing, only leaf
   ! blocks should be marked. - KW
   where (nodetype(:) .NE. LEAF)
      refine(:) = .false.
      derefine(:) = .false.
   end where

   return
end subroutine Grid_markRefineDerefine

