!!****f* source/physics/Hydro/Hydro
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
!!
!! NAME
!!
!!  Hydro
!!
!! SYNOPSIS
!!
!!  Hydro( real(IN)    :: timeEndAdv, 
!!         real(IN)    :: dt, 
!!         real(IN)    :: dtOld, 
!!         integer(IN) :: sweepOrder )
!!
!! DESCRIPTION
!!
!!  Computer one timestep advance for Hydro   
!!
!! ARGUMENTS
!!
!!  timeEndAdv -  end time
!!  dt -          timestep
!!  dtOld -       old timestep
!!  sweepOrder -  direction of hydro sweep, can be: SWEEP_XYZ or SWEEP_ZYX
!!                as defined in  constants.h if directionally split implementation is used
!!
!!
!!***


subroutine Hydro(simTime, dt, dtOld, sweeporder)
implicit none
#include "Simulation.h"
#include "constants.h"
  
  real,    INTENT(IN) :: simTime, dt, dtOld
  integer, optional, intent(IN) :: sweeporder
end subroutine Hydro
