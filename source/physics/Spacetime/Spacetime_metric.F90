!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief Stub file for `Spacetime_metric_t` derived-type

!> @ingroup Spacetime
!! @brief Contains a derived-type stub for managing metric-related quantities
module Spacetime_metric

   implicit none

   !> @brief Encapsulate commonly used metric quantities and operations
   !!
   !! @details
   !! Metric quantities:
   !! - Lapse function (scalar)
   !! - Shift vector (spatial, contravariant components)
   !! - Spatial metric (spatial, covariant, symmetric components)
   !! - Inverse spatial metric (spatial, contravariant, symmetric components)
   !! - Extrinsic curvature (spatial, covariant, symmetric components)
   !! - Square-root of the spatial metric determinant (scalar)
   !!
   !! Metric operations:
   !! - Raising/lowering rank-1 tensor components
   !! - Rank-1 tensor contractions
   type :: Spacetime_metric_t
   end type Spacetime_metric_t

end module Spacetime_metric
