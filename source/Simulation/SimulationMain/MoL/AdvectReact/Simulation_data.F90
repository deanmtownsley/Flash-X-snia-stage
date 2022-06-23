module Simulation_data

    implicit none

    real, save :: sim_cfl

    real, save :: sim_center, sim_width, sim_height
    real, save :: sim_speed, sim_beta

    integer, save :: PHI0_RHS, PHI1_RHS

end module Simulation_data
