!!****if* source/Simulation/SimulationMain/Brusselator/Simulation_molFastRHS
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
!!      Simulation_molFastRHS
!!
!!  SYNOPSIS
!!
!!      call Simulation_molFastRHS(real,    intent(in) :: t,
!                                  integer, intent(in) :: activeRHS
!!                                 real,    intent(in) :: dtWeight)
!!
!!  DESCRIPTION
!!
!!      Calculate fast RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      t         : Current time
!!      activeRHS : RHS data struct to fill
!!      dtWeight  : Weighted timestep (e.g. for flux corrections)
!!
!!***
subroutine Simulation_molFastRHS(t, activeRHS, dtWeight)
   use Simulation_data, only: sim_alpha, sim_beta, sim_epsilon, sim_lambdaF, sim_lambdaS, U_RHS

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
   integer, dimension(LOW:HIGH, MDIM) :: lim
   integer :: i, j, k
   real :: u, v, A, B, cost, cosbt

   A = 0.5d0*sim_lambdaF
   B = 0.5d0*(1d0 - sim_epsilon)*(sim_lambdaF - sim_lambdaS)/sim_alpha

   cost = cos(t)
   cosbt = cos(sim_beta*t)

   nullify (rhs); nullify (vars)

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
               ! Only V has a slow, implicit term
               u = vars(U_VAR, i, j, k)
               v = vars(V_VAR, i, j, k)

               rhs(U_RHS, i, j, k) = rhs(U_RHS, i, j, k) &
                                     + (A*(u**2 - cosbt - 3d0) - 0.5d0*sim_beta*sin(sim_beta*t))/u &
                                     + B*(v**2 - cost - 2d0)/v
            end do ! i
         end do ! j
      end do ! k

      call MoL_releaseDataPtr(tileDesc, rhs, activeRHS)
      call MoL_releaseDataPtr(tileDesc, vars, MOL_EVOLVED)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)
end subroutine Simulation_molFastRHS
