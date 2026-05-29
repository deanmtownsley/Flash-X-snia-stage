module sim_local_interface

  interface sim_interpolate1dWd
     subroutine sim_interpolate1dWd( radius, dr, dens, temp, xhe4, xc12, xo16, xne20, fluff)
     implicit none
     real, intent(in) :: radius, dr
     real, intent(out):: dens, temp, xhe4, xc12, xo16, xne20
     logical, intent(out) :: fluff
     end subroutine
  end interface

  interface sim_find_ign_state_keeping_pres
     subroutine sim_find_ign_state_keeping_pres(ign_temp, thermstate, abund)
     implicit none
#include "Eos.h"
#include "Flash.h"
     real, intent(in) :: ign_temp
     real, intent(inout), dimension(EOS_NUM) :: thermstate
     real, intent(in), dimension(SPECIES_BEGIN:SPECIES_END) :: abund
     end subroutine
  end interface

end module

