!!****if* source/Simulation/SimulationMain/Brusselator/Simulation_molExplicitRHS
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
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
!!      Simulation_molExplicitRHS
!!
!!  SYNOPSIS
!!
!!      call Simulation_molExplicitRHS(real, intent(in) :: t)
!!
!!  DESCRIPTION
!!
!!      Calculate explicit RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      t : Current time
!!
!!***
subroutine Simulation_molExplicitRHS(t)
   use Simulation_data, only: sim_rho, U_RHS, V_RHS, W_RHS

   use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator, &
                             Grid_fillGuardCells
   use Grid_iterator, only: Grid_iterator_t
   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "MoL.h"
#include "constants.h"

   implicit none

   real, intent(in) :: t

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   real, dimension(:, :, :, :), pointer :: rhs, vars
   integer :: i, j, k, ip, im

   integer, dimension(LOW:HIGH, MDIM) :: lim, bcs
   real :: del(MDIM), idx

   nullify (rhs); nullify (vars)

   if (sim_rho .lt. 0d0) then
      ip = 0
      im = -1
   else
      ip = 1
      im = 0
   end if

   ! Required for applying the upwind stencil
   call Grid_fillGuardCells(CENTER, ALLDIR)

   call Grid_getTileIterator(itor, LEAF)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      call tileDesc%faceBCs(bcs)

      lim = tileDesc%limits

      if (bcs(LOW, IAXIS) .ne. NOT_BOUNDARY) lim(LOW, IAXIS) = lim(LOW, IAXIS) + 1
      if (bcs(HIGH, IAXIS) .ne. NOT_BOUNDARY) lim(HIGH, IAXIS) = lim(HIGH, IAXIS) - 1

      call tileDesc%deltas(del)
      idx = 1d0/del(IAXIS)

      ! Note: In the following, the request for MOL_EVOLVED will
      !       always obtain a pointer to the variables in UNK; this
      !       call simply forwards to the tile descriptors `getDataPtr`.
      !       The request for MOL_RHS will behave in one of two ways,
      !       depending on the requirements of the selected integrator:
      !         - `rhs` will point to the current integration stage
      !            RHS memory structure as determined internally in MoL,
      !            and this will be typically be to stage-specific and
      !            type-specific (explicit, implicit, etc.)
      !         - `rhs` will point to the same (always the first and
      !           provided by default in MoL) RHS memory structure, and
      !           if the integrator requires saving this state, it will
      !           make a copy of the state into another block of memory
      !           that is not directly accessible to the user via
      !           requests for MOL_RHS
      call MoL_getDataPtr(tileDesc, vars, MOL_EVOLVED)
      call MoL_getDataPtr(tileDesc, rhs, MOL_RHS)

      do k = lim(LOW, KAXIS), lim(HIGH, KAXIS)
         do j = lim(LOW, JAXIS), lim(HIGH, JAXIS)
            do i = lim(LOW, IAXIS), lim(HIGH, IAXIS)
               rhs(U_RHS, i, j, k) = rhs(U_RHS, i, j, k) &
                                     + sim_rho*(vars(U_VAR, i + ip, j, k) - vars(U_VAR, i + im, j, k))*idx

               rhs(V_RHS, i, j, k) = rhs(V_RHS, i, j, k) &
                                     + sim_rho*(vars(V_VAR, i + ip, j, k) - vars(V_VAR, i + im, j, k))*idx

               rhs(W_RHS, i, j, k) = rhs(W_RHS, i, j, k) &
                                     + sim_rho*(vars(W_VAR, i + ip, j, k) - vars(W_VAR, i + im, j, k))*idx
            end do ! i
         end do ! j
      end do ! k

      call MoL_releaseDataPtr(tileDesc, rhs, MOL_RHS)
      call MoL_releaseDataPtr(tileDesc, vars, MOL_EVOLVED)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)

end subroutine Simulation_molExplicitRHS
