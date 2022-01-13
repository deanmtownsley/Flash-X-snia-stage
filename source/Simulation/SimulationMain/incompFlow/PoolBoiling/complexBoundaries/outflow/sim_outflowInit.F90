!!***if* source/Simulation/SimulationMain/incompFlow/PoolBoiling/localAPI/outflow/sim_outflowInit
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!
!!***

#include "constants.h"

subroutine sim_outflowInit()

  use Simulation_data,             ONLY : sim_meshMe
  use sim_outflowData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  call RuntimeParameters_get('sim_lsSink',sim_lsSink)
  call RuntimeParameters_get('sim_outflowRegion',sim_outflowRegion)
  call RuntimeParameters_get('ins_gravX',sim_gravX)
  call RuntimeParameters_get('ins_gravY',sim_gravY)
  call RuntimeParameters_get('ins_gravZ',sim_gravZ)

  if (sim_meshMe .eq. MASTER_PE) then
     write(*,*) 'sim_lsSink=',sim_lsSink
     write(*,*) 'sim_outflowRegion=',sim_outflowRegion
     write(*,*) 'sim_gravX=',sim_gravX
     write(*,*) 'sim_gravY=',sim_gravY
     write(*,*) 'sim_gravZ=',sim_gravZ
  end if

  sim_outflowVel = 0.

end subroutine sim_outflowInit
