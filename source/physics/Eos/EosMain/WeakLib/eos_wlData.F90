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
!! @brief  General parameters for EOS WeakLib
!!

module eos_wlData

 use wlEquationOfStateTableModule

 implicit none

 character(len=80),save :: eos_file

 type(EquationOfStateTableType), target, save :: EosNewTable
 type(EquationOfStateTableType), pointer :: eos_pointer

 logical, save :: eos_wl_muShift

 integer, save :: nVariables

 logical, save :: eos_postBounce = .FALSE.
 real, save :: eos_bounceTime = 0.0
 real, save :: eos_centralDens, eos_centralEntr
 integer, save :: eos_nstep

 ! The entropy within radius eos_shockEntrRad at which bounce will be flagged:
 real, parameter :: eos_shockEntr = 3.0
 real, parameter :: eos_shockEntrRad = 3.0e6
 ! The minimum central density at which bounce will be flagged:
 real, parameter :: eos_bounceDens = 2.0e14

end module eos_wlData
