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
!! @brief Burn_makeAbundanceYeConsistent implementation

!> @ingroup physics_sourceTerms_Burn
!!
!! @brief Implements Burn_makeAbundanceYeConsistent
!!
!! @details
!! @stubref{Burn_makeAbundanceYeConsistent}
!!
!! This assumes all changes in Ye are represented by changes to free nucleons
!!

!!REORDER(4): Uin

subroutine Burn_makeAbundanceYeConsistent()

#include "Simulation.h"
#include "constants.h"

  use Burn_data, ONLY : bn_nuclearDensMax, bn_smallx
  use Driver_interface, ONLY : Driver_abort
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

  real, pointer, dimension(:,:,:,:) :: Uin

  type(Grid_iterator_t)  :: itor
  type(Grid_tile_t) :: tileDesc

  integer, dimension(MDIM) :: lo, hi
  integer :: i, j, k

  real :: xn, xp, dye

#ifdef DYE_VAR
#ifdef N_SPEC
#ifdef P_SPEC

  call Grid_getTileIterator(itor, LEAF, tiling=.false. )
  do while(itor%isValid())
     call itor%CurrentTile(tileDesc)

     lo = tileDesc % limits(LOW ,1:MDIM)
     hi = tileDesc % limits(HIGH,1:MDIM)

     ! Get a pointer to solution data
     call tileDesc%getDataPtr(Uin, CENTER)

     do k = lo(KAXIS), hi(KAXIS)
        do j = lo(JAXIS), hi(JAXIS)
           do i = lo(IAXIS), hi(IAXIS)

              ! Check if composition needs to be adjusted
              if ( Uin(DENS_VAR,i,j,k) <= bn_nuclearDensMax ) then

                 ! Get old n and p mass fractions
                 xn = Uin(N_SPEC,i,j,k)
                 xp = Uin(P_SPEC,i,j,k)

                 ! Assume Ye change captured by changing free nucleon mass fractions
                 dye = Uin(DYE_VAR,i,j,k)
                 xn = MAX( xn - dye, bn_smallx )
                 xp = MAX( xp + dye, 1.0 )

                 ! Apply the update
                 Uin(N_SPEC,i,j,k) = xn
                 Uin(P_SPEC,i,j,k) = xp

              end if

           end do
        end do
     end do

     call tileDesc%releaseDataPtr(Uin, CENTER)
     nullify( Uin )

     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

#endif
#endif
#endif

end subroutine Burn_makeAbundanceYeConsistent
