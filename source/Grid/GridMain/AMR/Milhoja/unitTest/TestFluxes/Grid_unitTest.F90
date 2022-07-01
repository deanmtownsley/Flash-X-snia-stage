#include "Simulation.h"
#include "constants.h"

!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! Refer to the test's design doc for more information.
!!
!! @param fileUnit   Ignored.  All output is written to stdout.
!! @param perfect    True if no errors occurred; False, otherwise.
subroutine Grid_unitTest(fileUnit, perfect)
    use Grid_interface,  ONLY : Grid_getTileIterator, &
                                Grid_releaseTileIterator
    use Grid_iterator,   ONLY : Grid_iterator_t
    use Grid_tile,       ONLY : Grid_tile_t
    use ut_testDriverMod

    implicit none

    integer, intent(in)    :: fileUnit
    logical, intent(inout) :: perfect

    type(Grid_iterator_t)         :: itor
    type(Grid_tile_t)             :: tileDesc
    real,                 pointer :: fluxData(:, :, :, :)
    integer                       :: lBdd(MDIM+1)
    integer                       :: uBdd(MDIM+1)

    integer :: i, j, k, var

    NULLIFY(fluxData)

    CALL start_test_run

    CALL assertEqual(NFLUXES, 3, "Invalid N flux variables")

    !!!!! CONFIRM WRITING OF DATA
    CALL Grid_getTileIterator(itor, ALL_BLKS, tiling=.TRUE.)
    do while (itor%isValid())
        CALL itor%currentTile(tileDesc)

        associate(lo => tileDesc%limits(LOW,  :), &
                  hi => tileDesc%limits(HIGH, :))
            CALL tileDesc%getDataPtr(fluxData, FLUXX)
            CALL assertTrue(ASSOCIATED(fluxData), 'fluxData not acquired for X')
            lBdd(:) = lBound(fluxData)
            CALL assertEqual(lBdd(IAXIS),  lo(IAXIS), "Invalid fluxX lo(x)")
            CALL assertEqual(lBdd(JAXIS),  lo(JAXIS), "Invalid fluxX lo(y)")
            CALL assertEqual(lBdd(KAXIS),  lo(KAXIS), "Invalid fluxX lo(z)")
            CALL assertEqual(lBdd(MDIM+1), 1,         "Invalid fluxX lo(var)")

            uBdd(:) = uBound(fluxData)
            CALL assertEqual(uBdd(IAXIS),  hi(IAXIS)+1, "Invalid fluxX hi(x)")
            CALL assertEqual(uBdd(JAXIS),  hi(JAXIS),   "Invalid fluxX hi(y)")
            CALL assertEqual(uBdd(KAXIS),  hi(KAXIS),   "Invalid fluxX hi(z)")
            CALL assertEqual(uBdd(MDIM+1), NFLUXES,     "Invalid N flux vars")

            do           var = 1, NFLUXES
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)+1
                            fluxData(i, j, k, var) = -1.1 * var
                        end do
                    end do
                end do
            end do
            CALL tileDesc%releaseDataPtr(fluxData, FLUXX)
            CALL assertFalse(ASSOCIATED(fluxData), 'fluxData not released')

            CALL tileDesc%getDataPtr(fluxData, FLUXY)
#if NDIM >= 2
            CALL assertTrue(ASSOCIATED(fluxData), 'fluxData not acquired for Y')
            lBdd(:) = lBound(fluxData)
            CALL assertEqual(lBdd(IAXIS),  lo(IAXIS), "Invalid fluxY lo(x)")
            CALL assertEqual(lBdd(JAXIS),  lo(JAXIS), "Invalid fluxY lo(y)")
            CALL assertEqual(lBdd(KAXIS),  lo(KAXIS), "Invalid fluxY lo(z)")
            CALL assertEqual(lBdd(MDIM+1), 1,         "Invalid fluxY lo(var)")

            uBdd(:) = uBound(fluxData)
            CALL assertEqual(uBdd(IAXIS),  hi(IAXIS),   "Invalid fluxY hi(x)")
            CALL assertEqual(uBdd(JAXIS),  hi(JAXIS)+1, "Invalid fluxY hi(y)")
            CALL assertEqual(uBdd(KAXIS),  hi(KAXIS),   "Invalid fluxY hi(z)")
            CALL assertEqual(uBdd(MDIM+1), NFLUXES,     "Invalid N flux vars")

            do           var = 1, NFLUXES
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)+1
                        do i = lo(IAXIS), hi(IAXIS)
                            fluxData(i, j, k, var) = -2.1 * var
                        end do
                    end do
                end do
            end do
#else
            CALL assertFalse(ASSOCIATED(fluxData), 'fluxData not null for Y')
#endif
            CALL tileDesc%releaseDataPtr(fluxData, FLUXY)
            CALL assertFalse(ASSOCIATED(fluxData), 'fluxData not released')

            CALL tileDesc%getDataPtr(fluxData, FLUXZ)
#if NDIM == 3
            CALL assertTrue(ASSOCIATED(fluxData), 'fluxData not acquired for Z')
            lBdd(:) = lBound(fluxData)
            CALL assertEqual(lBdd(IAXIS),  lo(IAXIS), "Invalid fluxZ lo(x)")
            CALL assertEqual(lBdd(JAXIS),  lo(JAXIS), "Invalid fluxZ lo(y)")
            CALL assertEqual(lBdd(KAXIS),  lo(KAXIS), "Invalid fluxZ lo(z)")
            CALL assertEqual(lBdd(MDIM+1), 1,         "Invalid fluxZ lo(var)")

            uBdd(:) = uBound(fluxData)
            CALL assertEqual(uBdd(IAXIS),  hi(IAXIS),   "Invalid fluxZ hi(x)")
            CALL assertEqual(uBdd(JAXIS),  hi(JAXIS),   "Invalid fluxZ hi(y)")
            CALL assertEqual(uBdd(KAXIS),  hi(KAXIS)+1, "Invalid fluxZ hi(z)")
            CALL assertEqual(uBdd(MDIM+1), NFLUXES,     "Invalid N flux vars")

            do           var = 1, NFLUXES
                do         k = lo(KAXIS), hi(KAXIS)+1
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)
                            fluxData(i, j, k, var) = -3.1 * var
                        end do
                    end do
                end do
            end do
#else
            CALL assertFalse(ASSOCIATED(fluxData), 'fluxData not null for Z')
#endif
            CALL tileDesc%releaseDataPtr(fluxData, FLUXZ)
            CALL assertFalse(ASSOCIATED(fluxData), 'fluxData not released')
        end associate

        CALL itor%next()
    end do
    CALL Grid_releaseTileIterator(itor)

    CALL Grid_getTileIterator(itor, ALL_BLKS, tiling=.TRUE.)
    do while (itor%isValid())
        CALL itor%currentTile(tileDesc)

        associate(lo => tileDesc%limits(LOW,  :), &
                  hi => tileDesc%limits(HIGH, :))
            CALL tileDesc%getDataPtr(fluxData, FLUXX)
            do           var = 1, NFLUXES
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)+1
                            CALL assertEqual(fluxData(i, j, k, var), &
                                             -1.1 * var, &
                                             "Incorrect data in fluxX")
                        end do
                    end do
                end do
            end do
            CALL tileDesc%releaseDataPtr(fluxData, FLUXX)

#if NDIM >= 2
            CALL tileDesc%getDataPtr(fluxData, FLUXY)
            do           var = 1, NFLUXES
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)+1
                        do i = lo(IAXIS), hi(IAXIS)
                            CALL assertEqual(fluxData(i, j, k, var), &
                                             -2.1 * var, &
                                             "Incorrect data in fluxY")
                        end do
                    end do
                end do
            end do
            CALL tileDesc%releaseDataPtr(fluxData, FLUXY)
#endif
#if NDIM == 3
            CALL tileDesc%getDataPtr(fluxData, FLUXZ)
            do           var = 1, NFLUXES
                do         k = lo(KAXIS), hi(KAXIS)+1
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)
                            CALL assertEqual(fluxData(i, j, k, var), &
                                             -3.1 * var, &
                                             "Incorrect data in fluxZ")
                        end do
                    end do
                end do
            end do
            CALL tileDesc%releaseDataPtr(fluxData, FLUXZ)
#endif
        end associate

        CALL itor%next()
    end do
    CALL Grid_releaseTileIterator(itor)

    perfect = finish_test_run()
end subroutine Grid_unitTest

