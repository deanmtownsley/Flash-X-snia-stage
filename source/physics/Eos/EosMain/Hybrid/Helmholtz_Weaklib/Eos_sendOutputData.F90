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
!! @brief sends the scalar variables 
!!         to the IO unit, to be written to a checkpoint file
!! 
!! @subref{Eos_sendOutputData}
!!

subroutine Eos_sendOutputData()

   use IO_interface, ONLY: IO_setScalar
   use eos_hybridData, ONLY: eos_hybTransitionDensLo, eos_hybTransitionDensHi

   implicit none

   call IO_setScalar("hyb_transDensHi", eos_hybTransitionDensHi)
   call IO_setScalar("hyb_transDensLo", eos_hybTransitionDensLo)

end subroutine Eos_sendOutputData
