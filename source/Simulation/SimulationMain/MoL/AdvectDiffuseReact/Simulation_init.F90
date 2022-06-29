subroutine Simulation_init()
    use Simulation_data
    use RuntimeParameters_interface, only: RuntimeParameters_get
    use MoL_interface, only: MoL_registerVariable

#include "Simulation.h"
#include "constants.h"

    implicit none

    call RuntimeParameters_get("sim_a", sim_a)
    call RuntimeParameters_get("sim_b", sim_b)

    call RuntimeParameters_get("sim_alpha",   sim_alpha)
    call RuntimeParameters_get("sim_epsilon", sim_epsilon)
    call RuntimeParameters_get("sim_rho",     sim_rho)

    call RuntimeParameters_get("sim_k", sim_k)

    call MoL_registerVariable("u", U_VAR, U_RHS)
    call MoL_registerVariable("v", V_VAR, V_RHS)
    call MoL_registerVariable("w", W_VAR, W_RHS)
    
end subroutine Simulation_init
