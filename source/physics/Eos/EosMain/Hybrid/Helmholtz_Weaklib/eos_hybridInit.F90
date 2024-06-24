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
!> @ingroup physics_Eos
!!
!! @brief  Initialize the Hybrid (Helmholtz+WeakLib) EOS.
!!


subroutine eos_hybridInit()

   use Eos_data, ONLY: eos_type, eos_doYe
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use IO_interface, ONLY: IO_getPrevScalar
   use Driver_interface, ONLY: Driver_abort

   use eos_hybridData, ONLY: eos_hybTransitionDensLo, eos_hybTransitionDensHi, &
                             eos_hybDeltaE_WL, eos_hybBoverA
   ! use Eos_wlInterface, only: Eos_wlEnerShift
   use Multispecies_interface, only: Multispecies_getPropertyVector

#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

#ifdef FLASH_MULTISPECIES
#include "Multispecies.h"
#endif

   implicit none

   logical :: restart
   integer :: ierr, toterr

   real, dimension(NSPECIES) :: As, Bs

   ! set eos type
   eos_type = EOS_HYB
   eos_doYe = .true.
   !get the runtime parameters for the transition point
   call RuntimeParameters_get('restart', restart)
   if (restart) then
      toterr = 0
      call IO_getPrevScalar("hyb_transDensHi", eos_hybTransitionDensHi, ierr)
      toterr = toterr + ierr
      call IO_getPrevScalar("hyb_transDensLo", eos_hybTransitionDensLo, ierr)
      toterr = toterr + ierr

      if (toterr /= 0) then
         call Driver_abort("[eos_hybridInit] Could not initialize Hybrid EoS from checkpoint!")
      end if
   else
      call RuntimeParameters_get('eos_hybTransitionDensHi', eos_hybTransitionDensHi)
      call RuntimeParameters_get('eos_hybTransitionDensLo', eos_hybTransitionDensLo)
   end if

   ! Binding energy / mass number for each evolved species
   ! These will be multiplied against each species mass fraction
   ! when calculating the energy shift for Helmholtz
#ifdef FLASH_MULTISPECIES
   if (NSPECIES .gt. 0) then
      call Multispecies_getPropertyVector(A, As)
      call Multispecies_getPropertyVector(EB, Bs)

      eos_hybBoverA = Bs/As
   end if
#endif

   ! Setting this directly here to avoid confusion as to which
   ! energy shift in WeakLib this refers to that the `eos_wlEnerShift`
   ! subroutine returns
   eos_hybDeltaE_WL = 8.9 ! MeV/nucleon
end subroutine eos_hybridInit
