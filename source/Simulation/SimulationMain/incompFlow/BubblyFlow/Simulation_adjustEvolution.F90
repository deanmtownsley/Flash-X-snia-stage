!!****f* source/Simulation/incompFlow/BubblyFlow/Simulation_adjustEvolution
!! NOTICE
!!  Copyright 2024 UChicago Argonne, LLC and contributors
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
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!
!! ARGUMENTS
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
!!REORDER(4): solnData

#include "constants.h"
#include "Simulation.h"
#include "FortranLangFeatures.fh"

subroutine Simulation_adjustEvolution(nstep, dt, stime)

   use Grid_interface, ONLY: Grid_getTileIterator, Grid_releaseTileIterator
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Grid_iterator, ONLY: Grid_iterator_t
   use Grid_tile, ONLY: Grid_tile_t

   implicit none

   integer, intent(in) :: nstep
   real, intent(in) :: dt
   real, intent(in) :: stime

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc
   real, pointer, dimension(:, :, :, :) :: solnData
   integer :: i, j, k
   real :: del(MDIM)
   real, parameter :: lsBuffer = 1.0, lsForce = 0.2

   nullify (solnData)

   call Timers_start("Simulation_adjustEvolution")

   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
      call tileDesc%deltas(del)
      call tileDesc%getDataPtr(solnData, CENTER)

      do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
         do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
            do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)

               ! TODO: Add a repulsive force algorithm
               !------------------------------------------------------------
               !if ((solnData(DFUN_VAR, i, j, k) .lt. 0.0) .and. &
               !    (abs(solnData(DFUN_VAR, i, j, k)) .lt. lsBuffer)) then
               !   solnData(DFRC_VAR, i, j, k) = solnData(DFRC_VAR, i, j, k)- &
               !                                 lsForce*(solnData(NRMX_VAR, i, j, k)+solnData(NRMY_VAR, i, j, k))
               !
               !end if
               !------------------------------------------------------------

            end do
         end do
      end do

      call tileDesc%releaseDataPtr(solnData, CENTER)
      call itor%next()

   end do
   call Grid_releaseTileIterator(itor)

   call Timers_stop("Simulation_adjustEvolution")

end subroutine Simulation_adjustEvolution
