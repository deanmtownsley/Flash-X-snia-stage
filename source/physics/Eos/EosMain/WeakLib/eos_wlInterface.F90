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
!! @brief  header file for the WeakLib  public interfaces.
!!

Module eos_wlInterface

   implicit none

   interface
      subroutine eos_wlOneZone(xDens, xTemp, xYe, xEner, xPres, xEntr, xdedt, xCs2, xXp, xXn, xXa, xXh, xAbar, xVar, varID, mode)
         implicit none
         real, intent(INOUT) :: xDens, xYe
         real, intent(INOUT) :: xTemp, xEner, xEntr, xPres
         integer, intent(IN) :: mode, varID
         real, intent(OUT) :: xXp, xXn, xXa, xXh, xdedt, xCs2, xVar, xAbar
         real :: xZbar, xMu_e, xMu_n, xMu_p, xMuhat
      end subroutine eos_wlOneZone
   end interface

   interface
      subroutine eos_wlPotentials(xDens, xTemp, xYe, xMu_n, xMu_p, xMu_e)
         implicit none
         real, intent(in) :: xDens, xTemp, xYe
         real, intent(out) :: xMu_n, xMu_p
         real, intent(out), optional :: xMu_e
      end subroutine eos_wlPotentials
   end interface

   interface
      subroutine eos_wlDetectBounce(postBounce, bounceTime, centralDens, centralEntr)
         implicit none
         logical, intent(OUT) :: postBounce
         real, optional, intent(OUT) :: bounceTime, centralDens, centralEntr
      end subroutine eos_wlDetectBounce
   end interface

   interface
      subroutine eos_wlEnerShift(energyShift)
         implicit none
         real, intent(OUT) :: energyShift
      end subroutine eos_wlEnerShift
   end interface

end Module eos_wlInterface

