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
!! @brief ml_warn implementation

!> @ingroup MoLMain
!!
!! @brief Implements ml_warn
!!
!! @stubref{ml_warn}
subroutine ml_warn(msg)
   use MoL_data, only: ml_abortOnWarn, ml_verbosity, ml_mpiRank

   use Driver_interface, only: Driver_abort

#include "MoL.h"
#include "constants.h"

   implicit none

   character(len=*), intent(in) :: msg

   if (ml_abortOnWarn) then
      ! Always print warning if abort-on-warning is turned on
      call Driver_abort("[MoL] WARNING: "//msg)
   end if

   if ((ml_verbosity .ge. MOL_VERBOSE_WARN) .and. (ml_mpiRank .eq. MASTER_PE)) then
      print *, "[MoL] WARNING: "//msg
   end if
end subroutine ml_warn
