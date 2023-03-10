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
!! @brief MoL_getOrder implementation

!> @ingroup MoLMR
!!
!! @brief Implements MoL_getOrder
!!
!! @stubref{MoL_getOrder}
function MoL_getOrder() result(order)
   use ml_mrData, only: ml_slowOrder, ml_fastOrder

   implicit none

   integer :: order

   order = min(ml_slowOrder, ml_fastOrder)
end function MoL_getOrder
