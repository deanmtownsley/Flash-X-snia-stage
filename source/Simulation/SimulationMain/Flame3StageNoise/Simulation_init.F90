!!****if* source/Simulation/SimulationMain/Flame3StageNoise/Simulation_init

!!

!! NAME

!!

!!  Simulation_init

!!

!! SYNOPSIS

!!

!!  Simulation_init()

!!

!! DESCRIPTION

!!  Initialize private data for the 3-stage flame test setup

!!

!! ARGUMENTS

!!

!!

!!

!!***

subroutine Simulation_init()

    use Simulation_data
    use Flame_interface, ONLY : Flame_rhJumpReactive, Flame_getWidth, Flame_laminarSpeed
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use bn_paraInterface, ONLY : bn_paraFuelAshProperties
  
    implicit none
  
#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
  
    real :: laminarWidth
    real :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a
    real, dimension(2) :: info
  
    !--------------------------------------------------------
    !  initialize runtime parameters and some other constants
    !--------------------------------------------------------
    call RuntimeParameters_get( 'rho_ambient', sim_rhoAmbient)
    call RuntimeParameters_get( 't_ambient', sim_tAmbient)
    
    call RuntimeParameters_get( 'ignite', sim_ignite)
    call RuntimeParameters_get( 'frac_perturb', sim_fracPerturb)
    
    call RuntimeParameters_get( 'pseudo_1d', sim_pseudo1d)
    call RuntimeParameters_get( 'xctr_perturb', sim_xctrPerturb)
    call RuntimeParameters_get( 'yctr_perturb', sim_yctrPerturb)
    call RuntimeParameters_get( 'zctr_perturb', sim_zctrPerturb)
    call RuntimeParameters_get( 'theta', sim_theta)
    
    call RuntimeParameters_get( 'c_frac', sim_cFrac)
    call RuntimeParameters_get( 'ne_frac', sim_neFrac)
    
    !  this is grid info, no accessor functions available
    call RuntimeParameters_get( 'xmin', sim_xmin)
    call RuntimeParameters_get( 'xmax', sim_xmax)
    
    call RuntimeParameters_get( 'ymin', sim_ymin)
    call RuntimeParameters_get( 'ymax', sim_ymax)
    
    ! only need to get width of artificial flame once
    call Flame_getWidth(sim_laminarWidth)
  
    !--------------------------------------------------------
    !  find unburned and burned states
    !--------------------------------------------------------
    !  retrieve properties of fuel and 1st stage ash
    call bn_paraFuelAshProperties(sim_cFrac, sim_neFrac, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)
    ! put this information in an eos datastructure and save qbar
    sim_eosData_u(EOS_DENS) = sim_rhoAmbient
    sim_eosData_u(EOS_TEMP) = sim_tAmbient
    sim_eosData_u(EOS_ABAR) = 1.e0 / yi_f
    sim_eosData_u(EOS_ZBAR) = ye_f * sim_eosData_u(EOS_ABAR)
    sim_qbar_u = qbar_f
  
    ! now determine praperties of final NSE burned state
    call Flame_rhJumpReactive(sim_eosData_u, sim_qbar_u, sim_eosData_nse, sim_qbar_nse, MODE_DENS_TEMP)
    sim_deltae_nse = (sim_qbar_nse - sim_qbar_u) * 9.6485e17
  
    info(1) = sim_cFrac
    info(2) = sim_neFrac
    call Flame_laminarSpeed(sim_eosData_u(EOS_DENS), sim_flamespeed, info=info)
  
  end subroutine Simulation_init
  