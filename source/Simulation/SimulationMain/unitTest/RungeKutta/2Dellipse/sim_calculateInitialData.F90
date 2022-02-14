!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/sim_calculateInitialData
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
!!  sim_calculateInitialData
!! 
!! SYNOPSIS
!!
!!  call sim_calculateInitialData
!!
!! DESCRIPTION
!!
!!  Calculates the positive 'k' value for the elliptical 2D ODE equation from the given aspect
!!  ratio A:
!!
!!                 k = (A^2 + 1) / (A^2 - 1)
!!
!!  The routine also gives the parameters characterizing the elliptical path that will be traced.
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_calculateInitialData ()

  use  Simulation_data, ONLY: sim_ellipseMajorSemiAxis, &
                              sim_ellipseMinorSemiAxis, &
                              sim_ellipseRotationAngle, &
                              sim_ellipseCenterX,       &
                              sim_ellipseCenterY,       &
                              sim_ellipseAspectRatio,   &
                              sim_k

  implicit none

  real :: AxA
  real :: sqrt2
!
!
!   ...Calculate the 'k' value.
!
!
  AxA = sim_ellipseAspectRatio * sim_ellipseAspectRatio

  sim_k = (AxA + 1.0) / (AxA - 1.0)
!
!
!   ...Calculate the parameters that will characterize the elliptical path.
!
!
  sqrt2 = sqrt (2.0)

  sim_ellipseMinorSemiAxis = sqrt2
  sim_ellipseMajorSemiAxis = sqrt2 * sim_ellipseAspectRatio
  sim_ellipseRotationAngle = 45.0                             ! in degrees
  sim_ellipseCenterX       = 0.0
  sim_ellipseCenterY       = 0.0
!
!
!   ...Ready!
!
!
  return
end subroutine sim_calculateInitialData
