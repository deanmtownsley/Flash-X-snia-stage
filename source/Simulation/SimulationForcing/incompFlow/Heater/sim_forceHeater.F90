!!****f* source/Simulation/SimulationForcing/incompFlow/Heater/sim_forceHeater
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
!! NAME
!!  sim_forceHeater
!!
!!
!! SYNOPSIS
!!  sim_forceHeater( integer(IN) :: nstep,
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
#include "constants.h"
#include "Simulation.h"

subroutine sim_forceHeater(nstep, dt, stime)

   use Grid_interface, ONLY: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_iterator, ONLY: Grid_iterator_t
   use Grid_tile, ONLY: Grid_tile_t

   use sim_heaterInterface, ONLY: sim_heaterCheckSites, sim_heaterLSReInit, &
                                  sim_heaterTagSites
   implicit none
   include "Flashx_mpi.h"
   integer, intent(in) :: nstep
   real, intent(in) :: dt
   real, intent(in) :: stime

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer :: ierr

   ! Check Nucleation Sites
   !------------------------------------------------------------
   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
      !---------------------------------------------------------
      call sim_heaterCheckSites(tileDesc)
      !---------------------------------------------------------
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)

   ! Tag Nucleation Sites for renucleation
   !-------------------------------------------------------------
   call sim_heaterTagSites(stime)

   ! Re-initialize Level-Set Function
   !-------------------------------------------------------------
   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
      !---------------------------------------------------------
      call sim_heaterLSReInit(tileDesc, stime)
      !---------------------------------------------------------
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)

end subroutine sim_forceHeater
