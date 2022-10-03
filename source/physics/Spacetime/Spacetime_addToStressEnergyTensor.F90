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
!! @brief Spacetime_addToStressEnergyTensor stub

!> @ingroup Spacetime
!!
!! @brief Add contributions to the total stress energy tensor
!!
!! @details
!! @anchor Spacetime_addToStressEnergyTensor_stub
!!
!! This procedure will add the given contribution to the overall stress-energy
!! tensor.  Each physics unit will be responsible for making calls to this
!! procedure to add its own contribution.  This procedure will take the
!! zeroth-, first-, and second-rank projections of the (densitiezed) stress
!! energy tensor in the computational frame of an Eulerian observer, i.e.
!!
!! @f[
!!    \tilde{T}_{\mu\nu} = \tilde{E} n_\mu n_\nu + 2 \tilde{S}_{(\mu}n_{\nu)} + \tilde{S}_{\mu\nu}
!! @f]
!!
!! where @f$ n^\mu @f$  is the four-velocity of an Eulerian observer, and
!!
!! @f{align}{
!!    \tilde{E} &= \sqrt{\gamma}E\\
!!    \tilde{S}_{\mu} &= \sqrt{\gamma}S_{\mu}\\
!!    \tilde{S}_{\mu\nu} &= \sqrt{\gamma}S_{\mu\nu}
!! @f}
!!
!! are the "densitized" projections projections, with the square-root of the
!! determinant of the spatial metric @f$ \sqrt{gamma} @f$
!!
!! @param  E                        Rank-0 projection @f$ E @f$
!! @param  Sxx,Sy,Sz                Covariant spatial components of the
!!                                  rank-1 projection @f$ S_{i} @f$
!! @param  Sxx,Sxy,Sxz,Syy,Syz,Szz  Symmetric covariant spatial components of the
!!                                  rank-2 projection @f$ S_{ij} @f$
!! @param  tileDesc                 Current tile descriptor
!! @param  solnData                 Pointer to variables in UNK for the current tile
!! @param  loc                      Location (i,j,k) in the current tile
subroutine Spacetime_addToStressEnergyTensor(E, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                                             tileDesc, solnData, loc)
   use Grid_tile, only: Grid_tile_t

#include "constants.h"

   implicit none

   real, intent(in) :: E
   real, intent(in) :: Sx, Sy, Sz
   real, intent(in) :: Sxx, Sxy, Sxz, Syy, Syz, Szz
   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   return
end subroutine Spacetime_addToStressEnergyTensor
