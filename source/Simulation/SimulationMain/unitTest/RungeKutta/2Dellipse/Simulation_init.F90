!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/Simulation_init
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
!! SYNOPSIS
!!
!!  Simulation_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the parameters for the Runge Kutta 2D ellipse unit test.
!!
!!***

subroutine Simulation_init ()

  use  Simulation_data
  use  RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use  Logfile_interface,           ONLY: Logfile_stamp
  use  Driver_interface,            ONLY: Driver_abort 
  use  sim_interface,               ONLY: sim_calculateInitialData

  implicit none

#include "constants.h"
#include "Simulation.h"

  call RuntimeParameters_get ('sim_errorFraction',            sim_errorFraction          )
  call RuntimeParameters_get ('sim_RungeKuttaMethod',         sim_RungeKuttaMethod       )
  call RuntimeParameters_get ('sim_x0',                       sim_x0                     )
  call RuntimeParameters_get ('sim_y0',                       sim_y0                     )
  call RuntimeParameters_get ('sim_ellipseAspectRatio',       sim_ellipseAspectRatio     )
  call RuntimeParameters_get ('sim_stepSize',                 sim_stepSize               )
  call RuntimeParameters_get ('sim_numberOfEllipses',         sim_numberOfEllipses       )
!
!
!    ...Check parameters.
!
!
  if (sim_x0 /= 1.0 .or. sim_y0 /= 1.0) then
      call Driver_abort ('[Simulation_init] ERROR: Particle position not at (1,1)!')
  end if

  if (sim_stepSize <= 0.0) then
      call Driver_abort ('[Simulation_init] ERROR: Step size =< 0 !')
  end if
!
!
!    ...Set conversion factor(s).
!
!
  sim_deg2rad = acos (-1.0) / 180.0
!
!
!    ...Calculates the initital data.
!
!
  call sim_calculateInitialData ()
!
!
!    ...Ready!
!
!
  return
end subroutine Simulation_init
