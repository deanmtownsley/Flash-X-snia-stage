!!****f* source/Simulation/Simulation_wantsRebalance
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
!!  Simulation_wantsRebalance
!!
!! SYNOPSIS
!!
!!  logical force_rebalance = Simulation_wantsRebalance(integer(IN)  :: nstep
!!                                                      real(IN)     :: time )
!!
!! DESCRIPTION
!!
!!  The Simulation_wantsRebalance interface provides a way to
!!  customize how often Grid_updateRefinement requests a rebalance.
!!  This provides a way to redistribute load even if no blocks are
!!  being refined or derefined.
!!
!! ARGUMENTS
!!   nstep : current step
!!   time : current time
!!
!!***

logical function Simulation_wantsRebalance(nstep,time)

  implicit none
  integer, intent(in) :: nstep
  real, intent(in) :: time

  Simulation_wantsRebalance = .FALSE.

  ! User can customize requirement for requesting a rebalance.
  ! For example, request rebalance every 10 timesteps:
  ! if(mod(nstep, 10) == 0) Simulation_wantsRebalance = .TRUE.

end function Simulation_wantsRebalance
