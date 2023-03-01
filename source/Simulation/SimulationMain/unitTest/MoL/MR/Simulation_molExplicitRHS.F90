!!****if* source/Simulation/SimulationMain/unitTest/MoL/MR/Simulation_molExplicitRHS
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
subroutine Simulation_molExplicitRHS(t, activeRHS, dtWeight)
   use Simulation_data, only: V_RHS

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
   integer :: i, j, k

   integer, dimension(LOW:HIGH, MDIM) :: lim

   nullify (rhs); nullify (vars)

   ! No guard-cell filling necessary - just a bunch of local equations to solve

   call Grid_getTileIterator(itor, LEAF)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      lim = tileDesc%limits

      call MoL_getDataPtr(tileDesc, vars, MOL_EVOLVED)
      call MoL_getDataPtr(tileDesc, rhs, activeRHS)

      do k = lim(LOW, KAXIS), lim(HIGH, KAXIS)
         do j = lim(LOW, JAXIS), lim(HIGH, JAXIS)
            do i = lim(LOW, IAXIS), lim(HIGH, IAXIS)
               ! Only V has a slow explicit term
               rhs(V_RHS, i, j, k) = rhs(V_RHS, i, j, k) - 0.5d0*sin(t)/vars(V_VAR, i, j, k)
            end do ! i
         end do ! j
      end do ! k

      call MoL_releaseDataPtr(tileDesc, rhs, activeRHS)
      call MoL_releaseDataPtr(tileDesc, vars, MOL_EVOLVED)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)

end subroutine Simulation_molExplicitRHS
