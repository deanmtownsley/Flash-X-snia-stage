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
!! @brief Spacetime_getMetric implementation

!> @ingroup Minkowksi
!!
!! @brief Obtain the metric quantities at the provided location
!!
!! @details
!!
!! @stubref{Spacetime_getMetric}
subroutine Spacetime_getMetric(metric, tileDesc, solnData, loc)
   use Spacetime_interface, only: Spacetime_getSpatialMetric
   use Spacetime_metric, only: Spacetime_metric_t
   use Grid_tile, only: Grid_tile_t

#include "constants.h"

   implicit none

   type(Spacetime_metric_t), intent(out) :: metric
   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   metric%alp = 1d0

   metric%betax = 0d0
   metric%betay = 0d0
   metric%betaz = 0d0

   metric%kxx = 0d0
   metric%kxy = 0d0
   metric%kxz = 0d0
   metric%kyy = 0d0
   metric%kyz = 0d0
   metric%kzz = 0d0

   call Spacetime_getSpatialMetric(metric%gxx, metric%gxy, metric%gxz, &
                                   metric%gyy, metric%gyz, metric%gzz, &
                                   tileDesc, solnData, loc)

   metric%guxx = 1d0/metric%gxx
   metric%guxy = 0d0
   metric%guxz = 0d0
   metric%guyy = 1d0/metric%gyy
   metric%guyz = 0d0
   metric%guzz = 1d0/metric%gzz

   metric%sdetg = sqrt(this%gxx*this%gyy*this%gzz + 2d0*this%gxy*this%gxz*this%gyz &
                       - this%gxx*this%gyz**2 - this%gyy*this%gxz**2 - this%gzz*this%gxy**2)
end subroutine Spacetime_getMetric
