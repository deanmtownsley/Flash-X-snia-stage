!
! Dean Townsley 2008
! Klaus Weide 2013
!

#include "FortranLangFeatures.fh"

module fl_effInterface
  interface fl_effInit
     subroutine fl_effInit()
        implicit none
     end subroutine fl_effInit
  end interface

  interface fl_effFinalize
     subroutine fl_effFinalize()
        implicit none
     end subroutine
  end interface

  interface fl_effects
     subroutine fl_effects( solnData, flamdot, dt, tileDesc)
        use Grid_tile, ONLY : Grid_tile_t
        implicit none
        real, dimension(:,:,:,:), POINTER_INTENT_IN :: solnData
        real, dimension(:,:,:), intent(in)  :: flamdot
        real, intent(in) :: dt
        type(Grid_tile_t), intent(in) :: tileDesc
                ! Applies ancillary effects of the flame (e.g. energy release,
                ! composition change) by updating variables in the block passed
                ! via solnData
                ! Depending on the implementation, energy may be deposited
                ! in a unit called later (like Burn)
     end subroutine
  end interface

end module fl_effInterface
