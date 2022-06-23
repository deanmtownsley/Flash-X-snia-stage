subroutine Simulation_init()
    use Simulation_data
    use RuntimeParameters_interface, only: RuntimeParameters_get
    use MoL_interface, only: MoL_registerVariable

#include "Simulation.h"
#include "constants.h"

    implicit none

    call RuntimeParameters_get("sim_center", sim_center)
    call RuntimeParameters_get("sim_width",  sim_width)
    call RuntimeParameters_get("sim_height", sim_height)

    call RuntimeParameters_get("sim_speed",  sim_speed)
    call RuntimeParameters_get("sim_beta",   sim_beta)

    call RuntimeParameters_get("sim_cfl",    sim_cfl)

    call MoL_registerVariable("phi0", PHI0_VAR, PHI0_RHS)
    call MoL_registerVariable("phi1", PHI1_VAR, PHI1_RHS)

end subroutine Simulation_init
