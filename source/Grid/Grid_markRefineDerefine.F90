!!****f* source/Grid/Grid_markRefineDerefine
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
!! SIDE EFFECTS
!!  
!!  This routine works by modifying the global (processor-local) flag
!!  arrays
!!      newchild, refine, derefine, stay
!!  imported from the AMR implementation.
!!
!! NOTES
!!  
!!  The default implementation of this routine uses the Flash runtime
!!  parameters refine_var_{1,2,3,4}, refine_cutoff_{1,2,3,4},
!!  derefine_cutoff_{1,2,3,4}, and refine_filter_{1,2,3,4} to determine
!!  refinement.
!!  
!!  A non-directional guardcell fill for all CENTER variables that are
!!  referenced by any of those refine_var_# runtime parameters
!!  (including EOS calls for guardcells, if any refinement variables
!!  refine_var_# require them to be current) is performed in the
!!  default implementation of this routine.  Since the checking of
!!  refinement criteria may reference the values of some variables in
!!  inactive blocks, the caller of Grid_markRefineDerefine should
!!  ensure that those inactive blocks have values that are up to date
!!  (at least for the variables to be used for refinement criteria);
!!  this is done by an explicit restriction call in the default APR
!!  implementation of Grid_updateRefinement.
!!  
!!  Users creating their own implementation of this interface or of
!!  Grid_updateRefinement should make sure that the above remains
!!  true; that is, they should include the appropriate call(s) to
!!  Grid_fillGuardCells (and to a restriction routine if necessary) in
!!  their implementation.
!!
!! SEE ALSO
!!
!!  Grid_updateRefinement
!!  Grid_initDomain
!!  Grid_fillGuardCells
!!***



subroutine Grid_markRefineDerefine()

implicit none
   
end subroutine Grid_markRefineDerefine














