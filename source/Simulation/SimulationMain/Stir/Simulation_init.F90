!!****if* source/Simulation/SimulationMain/Stir/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!! SYNOPSIS
!!  Simulation_init()
!!
!! DESCRIPTION
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Driver_data, ONLY: dr_globalMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Simulation.h"

  call RuntimeParameters_get('sim_dens', sim_dens)
  call RuntimeParameters_get('sim_cs', sim_cs)
#ifdef MAGZ_VAR
  call RuntimeParameters_get('sim_magz', sim_magz)
#endif

  if (dr_globalMe .eq. MASTER_PE) then
    write(*,'(A)') 'Initializing the Stir turbulence driving setup...'
#ifdef MAGZ_VAR
    print *, 'Constant magnetic field was set to ', sim_magz, ' in the z-direction.'
#endif
  endif

end subroutine Simulation_init
