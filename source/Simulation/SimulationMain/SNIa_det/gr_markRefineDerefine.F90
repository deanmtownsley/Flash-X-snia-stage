!! Replacement of single-variable refinement marking routine from
!!  source/Grid/GridMain/paramesh/paramesh4/gr_markRefineDerefine.F90
!! see that file for main functional description
!!
!! This version is modified so that the parent-child reconciliation is
!! performed outside this subroutine (see notes at end).  Thus each
!! block is checked and marked if refinement is required based on the
!! input parameters (see main documentation in above subroutine for their
!! meaning).  The refinement test, based on a combination of derivatives
!! is the same as the default routine.
!!
!! Originial code: FLASH X
!! Modifications: Dean M. Townsley 2009-2026
!!

subroutine gr_markRefineDerefine(error, refine_cutoff,derefine_cutoff)

  use tree, ONLY: lnblocks, parent, nchild, child, newchild, nodetype, &
                  lrefine, refine, derefine, stay
  implicit none

#include "Simulation.h"
#include "constants.h"  

  real, intent(IN) :: error(MAXBLOCKS)
  real, intent(IN) :: refine_cutoff, derefine_cutoff

  integer lb
  
! MARK FOR REFINEMENT OR DEREFINEMENT

! Note that max refine/derefine and reconciliation between
!  parents and children must be done separately after all
!  refinement tests have been applied

  do lb = 1,lnblocks

     if (nodetype(lb).eq.1 .or. nodetype(lb).eq.2) then

        ! default state is derefine
        ! we must indicate if derefinement is forbidden
        ! or if refinement is required
        ! this is effectively or'ed over all criteria
        ! any one can prevent derefinement or request refinement

        ! test for refinement
        if (error(lb) .gt. refine_cutoff) then
           derefine(lb) = .false.
           refine(lb) = .true.
        else if ( error(lb).gt.derefine_cutoff ) then
           derefine(lb) = .false.
        endif
        
     end if
     
  end do

  !=========================================================================
  return
end subroutine gr_markRefineDerefine














