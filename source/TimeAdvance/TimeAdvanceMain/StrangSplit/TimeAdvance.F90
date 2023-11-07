!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief TimeAdvance implementation

!> @ingroup TimeAdvanceMain/StraongSplit
!!
!! @brief Implements TimeAdvance
!!
!! @details
!! This is the Strang Split implementation of time integration
!!
!! @stubref{TimeAdvance}

!!#define DEBUG_ADVANCE
#ifdef DEBUG_ALL
#define DEBUG_ADVANCE
#endif

#include "Simulation.h"
#include "constants.h"

subroutine TimeAdvance(dt, dtold, time)

   use Hydro_interface, ONLY: Hydro, Hydro_gravPotIsAlreadyUpdated
   use Gravity_interface, ONLY: Gravity_potential
   use RadTrans_interface, ONLY: RadTrans
   use Particles_interface, ONLY: Particles_advance, Particles_dump
   use Burn_interface, ONLY: Burn
   use Deleptonize_interface, ONLY: Deleptonize
   use Timers_interface, ONLY: Timers_start, Timers_stop
   implicit none

   real, intent(IN) :: dt, dtold, time

   call Hydro(time, dt, dtOld)
#ifdef DEBUG_ADVANCE
   print *, 'returned from hydro '
#endif

#ifndef DRIVER_DIFFULAST
   ! 3. Diffusive processes:
   call RadTrans(dt)
#ifdef DEBUG_ADVANCE
   print *, 'returned from RadTrans'
#endif
#endif

   ! 4. Add source terms:
   call Timers_start("sourceTerms")
   call Burn(dt)
   call Deleptonize(.false., dt, time)

   call Timers_stop("sourceTerms")
#ifdef DEBUG_ADVANCE
   print *, 'returned from sourceTerms'
#endif

#ifdef DRIVER_DIFFULAST
   ! 3. Diffusive processes: *** CHANGED ORDER !!! ***
   !    Radiation, viscosity, conduction, & magnetic registivity
   call RadTrans(dt)
#endif

   ! #. Advance Particles
   call Timers_start("Particles_advance")
   call Particles_advance(dtOld, dt)
   call Timers_stop("Particles_advance")
#ifdef DEBUG_ADVANCE
   print *, 'return from Particles_advance '  ! DEBUG
#endif

   !Allows evolution of gravitational potential for Spark Hydro
   ! #. Calculate gravitational potentials
   if (.NOT. Hydro_gravPotIsAlreadyUpdated()) then
      call Timers_start("Gravity potential")
      call Gravity_potential()
      call Timers_stop("Gravity potential")
#ifdef DEBUG_ADVANCE
      print *, 'return from Gravity_potential '  ! DEBUG
#endif
   end if

end subroutine TimeAdvance
