#include "Simulation.h"
subroutine eos_externalComputeAbarZbar( solnScalars, abarData, zbarData)

  use fl_effData, only : fl_eff_ye_u, fl_eff_ye_b, fl_eff_sumy_u, fl_eff_sumy_b

  real, dimension(:,:), intent(in) :: solnScalars
  real, dimension(:), intent(out)  :: abarData, zbarData

  ! index of flam variable in scalars array
  integer :: flami = FLAM_MSCALAR-SPECIES_BEGIN+1

  real :: flam, ye, yi
  integer i

  do i = 1, ubound(solnScalars, 2)
     flam = solnScalars(flami,i)
     abarData(i) = 1.0 / ( (1.0-flam)*fl_eff_sumy_u + flam*fl_eff_sumy_b )
     zbarData(i) = abarData(i)* ( (1.0-flam)*fl_eff_ye_u + flam*fl_eff_ye_b )
  enddo

end subroutine
