!!****if* source/Simulation/SimulationMain/Brusselator/Simulation_molExplicitRHS
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
!!      Simulation_molExplicitRHS
!!
!!  SYNOPSIS
!!
!!      call Simulation_molExplicitRHS(real,    intent(in) :: t,
!                                      integer, intent(in) :: activeRHS
!!                                     real,    intent(in) :: dtWeight)
!!
!!  DESCRIPTION
!!
!!      Calculate explicit RHS terms
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

subroutine Simulation_molExplicitRHS(t, activeRHS, dtWeight)
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
   integer, intent(in) :: activeRHS
   real, intent(in) :: dtWeight

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   real, dimension(:, :, :, :), pointer :: rhs, vars
   integer :: i, j, k, ip, im

   integer, dimension(LOW:HIGH, MDIM) :: lim, bcs
   real :: del(MDIM), idx

   nullify (rhs); nullify (vars)

   if (sim_rho .lt. 0.0) then
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
      idx = 1.0/del(IAXIS)

      ! Note: In the following, the request for MOL_EVOLVED will
      !       always obtain a pointer to the variables in UNK; this
      !       call simply forwards to the tile descriptors `getDataPtr`.
      call MoL_getDataPtr(tileDesc, vars, MOL_EVOLVED)
      call MoL_getDataPtr(tileDesc, rhs, activeRHS)

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

      call MoL_releaseDataPtr(tileDesc, rhs, activeRHS)
      call MoL_releaseDataPtr(tileDesc, vars, MOL_EVOLVED)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)

end subroutine Simulation_molExplicitRHS
