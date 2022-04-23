!!****if* source/Simulation/SimulationMain/incompFlow/FlowBoiling/Simulation_init
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! ARGUMENTS
!!
!!   none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for INS-isotropic turbulence problem.
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine Simulation_init()

   use Simulation_data
   use Driver_interface, ONLY: Driver_getMype, Driver_abort
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use sim_heaterInterface, ONLY: sim_heaterInit
   use sim_outletInterface, ONLY: sim_outletInit
   use sim_inletInterface, ONLY: sim_inletInit

   implicit none

   call Driver_getMype(MESH_COMM, sim_meshMe)

   call RuntimeParameters_get('xmin', sim_xMin)
   call RuntimeParameters_get('ymin', sim_yMin)
   call RuntimeParameters_get('xmax', sim_xMax)
   call RuntimeParameters_get('ymax', sim_yMax)
   call RuntimeParameters_get('zmin', sim_zmin)
   call RuntimeParameters_get('zmax', sim_zmax)
   call RuntimeParameters_get('ins_gravX', sim_gravX)
   call RuntimeParameters_get('ins_gravY', sim_gravY)
   call RuntimeParameters_get('ins_gravZ', sim_gravZ)

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'sim_gravX=', sim_gravX
      write (*, *) 'sim_gravY=', sim_gravY
      write (*, *) 'sim_gravZ=', sim_gravZ
   end if

   call sim_heaterInit()
   call sim_inletInit()
   call sim_outletInit()

end subroutine Simulation_init
