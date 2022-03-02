!!  Multiphase_interface
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
!! SYNOPSIS
!!
!!  use Multiphase_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Multiphase module
!! module that defines its public interfaces.
!!
!!***

Module Multiphase_interface

  implicit none

#include "constants.h"
#include "Simulation.h"


  interface
  subroutine Multiphase_init(restart)
  implicit none
  logical, intent(in) :: restart
  end subroutine Multiphase_init
  end interface

  interface
  subroutine Multiphase_advection()
  implicit none
  end subroutine Multiphase_advection
  end interface

  interface
  subroutine Multiphase_solve(dt)
  implicit none
  real,    INTENT(IN) :: dt
  end subroutine Multiphase_solve
  end interface

  interface
  subroutine Multiphase_redistance(iteration)
  implicit none
  integer, INTENT(IN) :: iteration
  end subroutine Multiphase_redistance
  end interface

  interface
  subroutine Multiphase_indicators()
  implicit none
  end subroutine Multiphase_indicators
  end interface

  interface
  subroutine Multiphase_finalize()
  implicit none
  end subroutine Multiphase_finalize
  end interface

  interface
  subroutine Multiphase_setFluidProps()
  implicit none
  end subroutine Multiphase_setFluidProps
  end interface

  interface
  subroutine Multiphase_setThermalProps()
  implicit none
  end subroutine Multiphase_setThermalProps
  end interface

  interface
  subroutine Multiphase_setPressureJumps()
  implicit none
  end subroutine Multiphase_setPressureJumps
  end interface

  interface
  subroutine Multiphase_getScalarProp(name, value)
  implicit none
  character(len=*), intent(in)  :: name
  real, intent(out)             :: value
  end subroutine Multiphase_getScalarProp
  end interface

  interface
  subroutine Multiphase_thermalForcing()
  implicit none
  end subroutine Multiphase_thermalForcing
  end interface

  interface
  subroutine Multiphase_divergence()
  implicit none
  end subroutine Multiphase_divergence
  end interface

  interface
  subroutine Multiphase_extrapFluxes(iteration)
  implicit none
  integer, INTENT(IN) :: iteration
  end subroutine Multiphase_extrapFluxes
  end interface

  interface
  subroutine Multiphase_setMassFlux()
  implicit none
  end subroutine Multiphase_setMassFlux
  end interface

  interface
  subroutine Multiphase_velForcing(dt)
  implicit none
  real, intent(in) :: dt
  end subroutine Multiphase_velForcing
  end interface

  interface
  subroutine Multiphase_reInitGridVars()
  implicit none
  end subroutine Multiphase_reInitGridVars
  end interface

  interface
  subroutine Multiphase_getGridVar(name, value)
  implicit none
  character(len=*), intent(in)  :: name
  integer, intent(out)          :: value
  end subroutine Multiphase_getGridVar
  end interface

end module Multiphase_interface
