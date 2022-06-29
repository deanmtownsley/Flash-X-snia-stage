!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! This subroutine is a callback that is registered with Milhoja's grid backend at
!! initialization.  The backend may call this subroutine many times during the
!! process of grid refinement so that Flash-X may communicate to the backend which
!! blocks in the given level require refinement.  The final refinement decisions
!! are made by the backend based on the information gathered with this callback.
!!
!! This routine iterates across all blocks in the given level and determines if the
!! current block needs refinement.  A block is marked for refinement if the block's
!! error estimate for any refinement variable is greater than the variable's
!! associated refinement cutoff value.
!!
!! This subroutine is presently a noop and as a result, it is one factor that
!! limits the Milhoja Grid unit implementation to pseudo-UG.
!!
!! @todo Should this routine be adapted so that it uses the generic Milhoja
!!       interface?  That seems unnecessary as this routine will likely only 
!!       ever be used specifically with AMReX.  I wonder if this should be ported
!!       to C so that it can be included in the AMReX-based implementation of the
!!       Milhoja Grid class.  It could be provided as a general-use, default
!!       error checking routine.  Calling code could simply select which error
!!       routine to use.  Calling code would still need to provide their own
!!       version of gr_estimateBlkError.  This should be studied along with 
!!       the need for some Flash-X users to supply their own routine.  The types
!!       have been set to MILHOJA.  However, this is fundamentally an AMReX 
!!       callback and it makes more conceptual sense for these to be AMReX types.
!!
!! @param lev  The 0-based level index
!! @param tags  C-pointer to an AMReX tagbox array.  The elements of this are tag
!!          boxes.  The cells of these tagboxes are set to communicate a need
!!          to refine the associated block.
!! @param time    Not used in this default implementation of Flash-X
!! @param tagval  For full, rich AMReX tagging, this values should be assigned to
!!            each cell that has insufficient resolution.
!! @param clearval  For full, rich AMReX tagging, this values should be assigned to
!!              each cell that has sufficient resolution.
subroutine gr_markRefineDerefineCallback(lev, tags, time, tagval, clearval) bind(c)
   use iso_c_binding, ONLY : C_CHAR, &
                             C_PTR

   use milhoja_types_mod, ONLY : MILHOJA_INT, &
                                 MILHOJA_REAL

   implicit none

   integer(MILHOJA_INT), intent(IN), value :: lev
   type(C_PTR),          intent(IN), value :: tags
   real(MILHOJA_REAL),   intent(IN), value :: time
   character(C_CHAR),    intent(IN), value :: tagval
   character(C_CHAR),    intent(IN), value :: clearval

    RETURN
!   write(*,'(A)') "[gr_markRefineDerefineCallback] AMReX wants marking"
end subroutine gr_markRefineDerefineCallback

