!!****if* source/Simulation/SimulationMain/Sedov/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init()
  use Particles_data, ONLY : pt_meshMe
  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_meshMe
  use reader, ONLY : read_array_from_file
  USE Grid_interface, ONLY : Grid_getDomainBoundBox 

  implicit none
#include "constants.h"
#include "Flash.h"

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_tempAmbient', sim_tempAmbient)
  call RuntimeParameters_get('sim_YeAmbient', sim_YeAmbient)

end subroutine Simulation_init
