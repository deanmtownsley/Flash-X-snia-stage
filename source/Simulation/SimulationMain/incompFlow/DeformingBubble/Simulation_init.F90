!!****if* source/Simulation/SimulationMain/incompFlow/DeformingBubble/Simulation_init
!!
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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

   use Driver_interface, ONLY: Driver_getMype
   use Simulation_data, ONLY: sim_xMin, sim_yMin, &
                              sim_xMax, sim_yMax, &
                              sim_zMin, sim_zMax, &
                              sim_meshMe, sim_reInitFlow, &
                              sim_numBubbles, sim_bubbleLoc, &
                              sim_refineMax

   use RuntimeParameters_interface, ONLY: RuntimeParameters_get

   implicit none
   integer :: ib, jb, kb, ibubble

   call Driver_getMype(MESH_COMM, sim_meshMe)

   call RuntimeParameters_get('xmin', sim_xMin)
   call RuntimeParameters_get('ymin', sim_yMin)
   call RuntimeParameters_get('xmax', sim_xMax)
   call RuntimeParameters_get('ymax', sim_yMax)
   call RuntimeParameters_get('zmin', sim_zMin)
   call RuntimeParameters_get('zmax', sim_zMax)

   call RuntimeParameters_get('sim_reInitFlow', sim_reInitFlow)
   call RuntimeParameters_get('lrefine_max', sim_refineMax)

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'sim_reInitFlow =', sim_reInitFlow
      write (*, *) 'sim_refineMax =', sim_refineMax
   end if

   ! Initialize dimensional scales
#if NDIM < MDIM
   sim_numBubbles(IAXIS) = int(sim_xMax-sim_xMin)
   sim_numBubbles(JAXIS) = int(sim_yMax-sim_yMin)
   sim_numBubbles(KAXIS) = 1

#else
   sim_numBubbles(IAXIS) = int(sim_xMax-sim_xMin)
   sim_numBubbles(JAXIS) = int(sim_yMax-sim_yMin)
   sim_numBubbles(KAXIS) = int(sim_zMax-sim_zMin)

#endif

   allocate (sim_bubbleLoc(MDIM, product(sim_numBubbles)))

   ibubble = 0

   do kb = 1, sim_numBubbles(KAXIS)
      do jb = 1, sim_numBubbles(JAXIS)
         do ib = 1, sim_numBubbles(IAXIS)

            ibubble = ibubble+1

#if NDIM < MDIM
            sim_bubbleLoc(:, ibubble) = (/(ib-1)+0.75, &
                                          (jb-1)+0.75, &
                                          0./)
#else
            sim_bubbleLoc(:, ibubble) = (/(ib-1)+0.75, &
                                          (jb-1)+0.75, &
                                          (kb-1)+0.50/)
#endif

         end do
      end do
   end do

end subroutine Simulation_init
