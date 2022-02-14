!!****if* source/Simulation/SimulationMain/unitTest/Roots/x3Polynomials/Driver_initAll
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
!!  Driver_initAll
!!
!! SYNOPSIS
!!
!!  Driver_initAll ()
!!
!! DESCRIPTION
!!
!!  Stripped down version for testing single units.
!!
!! NOTES
!!
!!***

subroutine Driver_initAll ()
  
  use Driver_data,                 ONLY : dr_elapsedWCTime,        &
                                          dr_globalComm,           &
                                          dr_globalMe,             &
                                          dr_globalNumProcs,       &
                                          dr_initialWCTime,        &
                                          dr_particlesInitialized, &
                                          dr_restart
  use Driver_interface,            ONLY : Driver_init,               &
                                          Driver_initNumericalTools, &
                                          Driver_setupParallelEnv
  use RuntimeParameters_interface, ONLY : RuntimeParameters_init
  use Simulation_interface,        ONLY : Simulation_init
  use Logfile_interface,           ONLY : Logfile_init

  implicit none       
  
  call dr_set_rlimits            (dr_globalMe)
  call RuntimeParameters_init    (dr_restart)
  call Driver_setupParallelEnv   ()
  call Timers_init               (dr_initialWCTime)
  call Logfile_init              ()
  call Driver_init               ()
  call Driver_initNumericalTools ()
  call Simulation_init           ()

  return
end subroutine Driver_initAll
