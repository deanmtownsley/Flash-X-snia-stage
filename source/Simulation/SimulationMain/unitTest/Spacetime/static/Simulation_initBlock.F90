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
!! @brief `Simulation_initBlock` implementation for unitTest/Spacetime/static

!> @brief Fill in the initial data for the static Spacetime unit test
!!
!! @details
!!
!! This procedure will fill in the ADM variables contained in @ref SpacetimeStatic
!! with a Schwarschild black hole solution using Kerr-Schild coordinates in a
!! Cartesian basis.  For simplicity, a black hole of unitary mass (in geometric
!! units) will be assumed.
!!
!! @note The provided setup line and `flash.par` file will always be origin
!!       avoid (i.e. r = 0 will never occur at a cell center)
!!
!! @param vars      Pointer to the solution data
!! @param tileDesc  Current tile descriptor
subroutine Simulation_initBlock(vars, tileDesc)
   use Grid_tile, only: Grid_tile_t
   use Grid_interface, only: Grid_getCellCoords

#include "ADM.h"
#include "constants.h"

   implicit none

   real, dimension(:, :, :, :), pointer :: vars
   type(Grid_tile_t), intent(in) :: tileDesc

   integer :: i, j, k

   real, dimension(:), allocatable :: x, y, z
   real :: H, r, alp, lx, ly, lz

   real, parameter :: M = 1.0

   allocate (x(tileDesc%limits(LOW, IAXIS):tileDesc%limits(HIGH, IAXIS)), &
             y(tileDesc%limits(LOW, JAXIS):tileDesc%limits(HIGH, JAXIS)), &
             z(tileDesc%limits(LOW, KAXIS):tileDesc%limits(HIGH, KAXIS)))
   x = 0.0; y = 0.0; z = 0.0

   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                           tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), x)
#if (NDIM > 1)
   call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                           tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), y)
#endif
#if (NDIM > 2)
   call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                           tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), z)
#endif

   do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
      do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
         do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
            ! Radius
            r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)

            ! mass and coordinate ratios
            H = M/r

            lx = x(i)/r
            ly = y(j)/r
            lz = z(k)/r

            ! Keeping a copy of the lapse outside of vars(...) for later use
            alp = 1.0/sqrt(1.0 + 2.0*H)

            ! Lapse
            vars(ALP_VAR, i, j, k) = alp

            ! Shift
            vars(BETAX_VAR, i, j, k) = 2.0*H*alp**2*lx
            vars(BETAY_VAR, i, j, k) = 2.0*H*alp**2*ly
            vars(BETAZ_VAR, i, j, k) = 2.0*H*alp**2*lz

            ! Spatial metric (including 0 explicitly for ease of reading)
            vars(GXX_VAR, i, j, k) = 1.0 + 2.0*H*lx*lx
            vars(GXY_VAR, i, j, k) = 0.0 + 2.0*H*lx*ly
            vars(GXZ_VAR, i, j, k) = 0.0 + 2.0*H*lx*lz
            vars(GYY_VAR, i, j, k) = 1.0 + 2.0*H*ly*ly
            vars(GYZ_VAR, i, j, k) = 0.0 + 2.0*H*ly*lz
            vars(GZZ_VAR, i, j, k) = 1.0 + 2.0*H*lz*lz

            ! Extrinsic curvature (including 0 explicitly for ease of reading)
            vars(KXX_VAR, i, j, k) = 2.0*H*alp/r*(1.0 - (2.0 + H)*lx*lx)
            vars(KXY_VAR, i, j, k) = 2.0*H*alp/r*(0.0 - (2.0 + H)*lx*ly)
            vars(KXZ_VAR, i, j, k) = 2.0*H*alp/r*(0.0 - (2.0 + H)*lx*lz)
            vars(KYY_VAR, i, j, k) = 2.0*H*alp/r*(1.0 - (2.0 + H)*ly*ly)
            vars(KYZ_VAR, i, j, k) = 2.0*H*alp/r*(0.0 - (2.0 + H)*ly*lz)
            vars(KZZ_VAR, i, j, k) = 2.0*H*alp/r*(1.0 - (2.0 + H)*lz*lz)
         end do ! i
      end do ! j
   end do ! k

   deallocate (x, y, z)

end subroutine Simulation_initBlock
