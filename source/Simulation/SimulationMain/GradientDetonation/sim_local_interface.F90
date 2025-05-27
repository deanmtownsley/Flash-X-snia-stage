module sim_local_interface

  interface sim_find_ign_state_keeping_pres
     subroutine sim_find_ign_state_keeping_pres(ign_temp, thermstate, abund)
     implicit none
#include "Eos.h"
#include "Simulation.h"
     real, intent(in) :: ign_temp
     real, intent(inout), dimension(EOS_NUM) :: thermstate
     real, intent(in), dimension(SPECIES_BEGIN:SPECIES_END) :: abund
     end subroutine
  end interface

end module

