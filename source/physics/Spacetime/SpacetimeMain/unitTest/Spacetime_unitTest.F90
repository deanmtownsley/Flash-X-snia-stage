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
!! @brief Spacetime_unitTest implementation

!> @ingroup SpacetimeMain
!!
!! @details
!!
!! This implementation will verify that the Schwarzschild initial data
!! in Cartesion Kerr-Schild coordinates for a black hole of unitary
!! mass in geometric units is properly set and stored by comparing
!! the values obtained by the ADM accessors matches the expected
!! solution
!!
!! @stubref{Spacetime_unitTest}
subroutine Spacetime_unitTest(fileUnit, perfect)
   use Grid_tile, only: Grid_tile_t
   use Grid_iterator, only: Grid_iterator_t
   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator, &
                             Grid_getCellCoords
   use ut_testDriverMod, only: start_test_run, finish_test_run, assertEqual

#include "ADM.h"
#include "constants.h"

   implicit none

   integer, intent(in) :: fileUnit
   logical, intent(inout) :: perfect

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc
   real, dimension(:, :, :, :), pointer :: vars

   integer :: i, j, k

   real, dimension(:), allocatable :: x, y, z
   real :: H, r, lx, ly, lz

   real :: alpActual
   real :: betaxActual, betayActual, betazActual
   real :: gxxActual, gxyActual, gxzActual, gyyActual, gyzActual, gzzActual
   real :: KxxActual, KxyActual, KxzActual, KyyActual, KyzActual, KzzActual

   real, parameter :: M = 1d0

   nullify (vars)

   ! We should probably pass fileUnit here to make the output consistent
   ! with what may be expected by a unitTest simulation
   call start_test_run()

   call Grid_getTileIterator(itor, LEAF)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      allocate (x(tileDesc%limits(LOW, IAXIS):tileDesc%limits(HIGH, IAXIS)), &
                y(tileDesc%limits(LOW, JAXIS):tileDesc%limits(HIGH, JAXIS)), &
                z(tileDesc%limits(LOW, KAXIS):tileDesc%limits(HIGH, KAXIS)))
      x = 0d0; y = 0d0; z = 0d0

      call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                              tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), x)
#if NDIM > 1
      call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                              tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), y)
#endif
#if NDIM > 2
      call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                              tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), z)
#endif

      call tileDesc%getDataPtr(vars, CENTER)

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

               ! Lapse
               alpActual = 1d0/sqrt(1d0 + 2d0*H)

               ! Shift
               betaxActual = 2d0*H*alpActual**2*lx
               betayActual = 2d0*H*alpActual**2*ly
               betazActual = 2d0*H*alpActual**2*lz

               ! Spatial metric (including 0 explicitly for ease of reading)
               gxxActual = 1d0 + 2d0*H*lx*lx
               gxyActual = 0d0 + 2d0*H*lx*ly
               gxzActual = 0d0 + 2d0*H*lx*lz
               gyyActual = 1d0 + 2d0*H*ly*ly
               gyzActual = 0d0 + 2d0*H*ly*lz
               gzzActual = 1d0 + 2d0*H*lz*lz

               ! Extrinsic curvature (including 0 explicitly for ease of reading)
               KxxActual = 2d0*H*alpActual/r*(1d0 - (2d0 + H)*lx*lx)
               KxyActual = 2d0*H*alpActual/r*(0d0 - (2d0 + H)*lx*ly)
               KxzActual = 2d0*H*alpActual/r*(0d0 - (2d0 + H)*lx*lz)
               KyyActual = 2d0*H*alpActual/r*(1d0 - (2d0 + H)*ly*ly)
               KyzActual = 2d0*H*alpActual/r*(0d0 - (2d0 + H)*ly*lz)
               KzzActual = 2d0*H*alpActual/r*(1d0 - (2d0 + H)*lz*lz)

               ! Verify that these match
               call assertEqual(vars(ALP_VAR, i, j, k), alpActual, "Lapse is incorrect")

               call assertEqual(vars(BETAX_VAR, i, j, k), betaxActual, "Shift vector x-component is incorrect")
               call assertEqual(vars(BETAY_VAR, i, j, k), betayActual, "Shift vector y-component is incorrect")
               call assertEqual(vars(BETAZ_VAR, i, j, k), betazActual, "Shift vector z-component is incorrect")

               call assertEqual(vars(GXX_VAR, i, j, k), gxxActual, "Spatial metric xx-component is incorrect")
               call assertEqual(vars(GXY_VAR, i, j, k), gxyActual, "Spatial metric xy-component is incorrect")
               call assertEqual(vars(GXZ_VAR, i, j, k), gxzActual, "Spatial metric xz-component is incorrect")
               call assertEqual(vars(GYY_VAR, i, j, k), gyyActual, "Spatial metric yy-component is incorrect")
               call assertEqual(vars(GYZ_VAR, i, j, k), gyzActual, "Spatial metric yz-component is incorrect")
               call assertEqual(vars(GZZ_VAR, i, j, k), gzzActual, "Spatial metric zz-component is incorrect")

               call assertEqual(vars(KXX_VAR, i, j, k), KxxActual, "Extrinsic curvature xx-component is incorrect")
               call assertEqual(vars(KXY_VAR, i, j, k), KxyActual, "Extrinsic curvature xy-component is incorrect")
               call assertEqual(vars(KXZ_VAR, i, j, k), KxzActual, "Extrinsic curvature xz-component is incorrect")
               call assertEqual(vars(KYY_VAR, i, j, k), KyyActual, "Extrinsic curvature yy-component is incorrect")
               call assertEqual(vars(KYZ_VAR, i, j, k), KyzActual, "Extrinsic curvature yz-component is incorrect")
               call assertEqual(vars(KZZ_VAR, i, j, k), KzzActual, "Extrinsic curvature zz-component is incorrect")
            end do ! i
         end do ! j
      end do ! k

      call tileDesc%releaseDataPtr(vars, CENTER)

      deallocate (x, y, z)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)

   perfect = perfect .and. finish_test_run()

end subroutine Spacetime_unitTest
