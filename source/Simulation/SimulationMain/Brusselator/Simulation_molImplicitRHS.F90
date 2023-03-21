!!****if* source/Simulation/SimulationMain/Brusselator/Simulation_molImplicitRHS
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!  NAME
!!
!!      Simulation_molImplicitRHS
!!
!!  SYNOPSIS
!!
!!      call Simulation_molImplicitRHS(real,    intent(in) :: t,
!                                      integer, intent(in) :: activeRHS
!!                                     real,    intent(in) :: dtWeight)
!!
!!  DESCRIPTION
!!
!!      Calculate implicit RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      t         : Current time
!!      activeRHS : RHS data struct to fill
!!      dtWeight  : Weighted timestep (e.g. for flux corrections)
!!
!!***

!!REORDER(4): vars,rhs

subroutine Simulation_molImplicitRHS(t, activeRHS, dtWeight)
   use Simulation_data, only: sim_alpha, U_RHS, V_RHS, W_RHS

   use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_iterator, only: Grid_iterator_t
   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "MoL.h"
#include "constants.h"

   implicit none

   real, intent(in) :: t
   integer, intent(in) :: activeRHS
   real, intent(in) :: dtWeight

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   real, dimension(:, :, :, :), pointer :: rhs, vars
   integer, dimension(LOW:HIGH, MDIM) :: lim, bcs
   integer :: i, j, k
   real :: del(MDIM), idx2

   nullify (rhs); nullify (vars)

   call Grid_getTileIterator(itor, LEAF)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      call tileDesc%faceBCs(bcs)

      lim = tileDesc%limits

      if (bcs(LOW, IAXIS) .ne. NOT_BOUNDARY) lim(LOW, IAXIS) = lim(LOW, IAXIS) + 1
      if (bcs(HIGH, IAXIS) .ne. NOT_BOUNDARY) lim(HIGH, IAXIS) = lim(HIGH, IAXIS) - 1

      call tileDesc%deltas(del)
      idx2 = 1.0/(del(IAXIS)**2)

      ! Note: In the following, the request for MOL_EVOLVED will
      !       always obtain a pointer to the variables in UNK; this
      !       call simply forwards to the tile descriptors `getDataPtr`.
      call MoL_getDataPtr(tileDesc, vars, MOL_EVOLVED)
      call MoL_getDataPtr(tileDesc, rhs, activeRHS)

      do k = lim(LOW, KAXIS), lim(HIGH, KAXIS)
         do j = lim(LOW, JAXIS), lim(HIGH, JAXIS)
            do i = lim(LOW, IAXIS), lim(HIGH, IAXIS)
               rhs(U_RHS, i, j, k) = rhs(U_RHS, i, j, k) &
                                     + sim_alpha*(vars(U_VAR, i + 1, j, k) &
                                                  - 2.0*vars(U_VAR, i, j, k) &
                                                  + vars(U_VAR, i - 1, j, k))*idx2

               rhs(V_RHS, i, j, k) = rhs(V_RHS, i, j, k) &
                                     + sim_alpha*(vars(V_VAR, i + 1, j, k) &
                                                  - 2.0*vars(V_VAR, i, j, k) &
                                                  + vars(V_VAR, i - 1, j, k))*idx2

               rhs(W_RHS, i, j, k) = rhs(W_RHS, i, j, k) &
                                     + sim_alpha*(vars(W_VAR, i + 1, j, k) &
                                                  - 2.0*vars(W_VAR, i, j, k) &
                                                  + vars(W_VAR, i - 1, j, k))*idx2
            end do ! i
         end do ! j
      end do ! k

      call MoL_releaseDataPtr(tileDesc, rhs, activeRHS)
      call MoL_releaseDataPtr(tileDesc, vars, MOL_EVOLVED)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)
end subroutine Simulation_molImplicitRHS
