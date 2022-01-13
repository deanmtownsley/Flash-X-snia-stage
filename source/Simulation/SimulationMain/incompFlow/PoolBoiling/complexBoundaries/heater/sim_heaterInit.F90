!!***if* source/Simulation/SimulationMain/incompFlow/PoolBoiling/localAPI/heater/sim_heaterInit
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

subroutine sim_heaterInit()

  use Simulation_data,             ONLY : sim_meshMe
  use sim_heaterData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use sim_heaterInterface,         ONLY : sim_heaterRead

  implicit none
  character(len=20) :: heaterFile
  integer           :: htr

  call RuntimeParameters_get('sim_numHeaters',sim_numHeaters)
  call RuntimeParameters_get('sim_nucSeedRadius', sim_nucSeedRadius)

  if (sim_meshMe .eq. MASTER_PE) then
     write(*,*) 'sim_numHeaters=',sim_numHeaters
     write(*,*) 'sim_nucSeedRadius=',sim_nucSeedRadius
  end if

  allocate(sim_heaterInfo(sim_numHeaters))

  do htr=1,sim_numHeaters
     write(heaterFile,"(A,I4.4,A)") 'sim_heater.',htr,'.h5'
     call sim_heaterRead(htr,heaterFile)
  end do

  return

end subroutine sim_heaterInit
