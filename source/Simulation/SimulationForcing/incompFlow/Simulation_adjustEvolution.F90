!!****f* source/Simulation/SimulationForcing/incompFlow/Simulation_adjustEvolution
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
#include "constants.h"
#include "Simulation.h"

subroutine Simulation_adjustEvolution(nstep, dt, stime)

   use Grid_interface, ONLY: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_iterator, ONLY: Grid_iterator_t
   use Grid_tile, ONLY: Grid_tile_t

   use IncompNS_interface, ONLY: IncompNS_setVectorProp
   use Simulation_data, ONLY: sim_meshMe

#ifdef SIMULATION_FORCE_HEATER
   use sim_heaterInterface, ONLY: sim_heaterCheckSites, sim_heaterLSReInit, &
                                  sim_heaterTagSites
#endif

#ifdef SIMULATION_FORCE_OUTFLOW
   use sim_outflowInterface, ONLY: sim_outflowSetForcing
   use sim_outflowData, ONLY: sim_outflowVel
#endif

   implicit none
   include "Flashx_mpi.h"
   integer, intent(in) :: nstep
   real, intent(in) :: dt
   real, intent(in) :: stime

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer :: ierr
   integer :: htr, isite
   real :: velOutAux(LOW:HIGH, MDIM)

#ifdef SIMULATION_FORCE_HEATER

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
#endif

#ifdef SIMULATION_FORCE_OUTFLOW
   velOutAux = 0.

   ! Set Outflow Forcing
   !-------------------------------------------------------------
   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
      !---------------------------------------------------------
      call sim_outflowSetForcing(tileDesc, velOutAux, dt)
      !---------------------------------------------------------
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)

   ! Consolidate data
   !-------------------------------------------------------------
   sim_outflowVel = 0.

   call MPI_Allreduce(velOutAux, sim_outflowVel, (HIGH - LOW + 1)*MDIM, FLASH_REAL, &
                      MPI_MAX, MPI_COMM_WORLD, ierr)

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'Outflow Velocity Low  =', sim_outflowVel(LOW, :)
      write (*, *) 'Outflow Velocity High =', sim_outflowVel(HIGH, :)
   end if

   call IncompNS_setVectorProp("Outflow_Vel_Low", sim_outflowVel(LOW, :))
   call IncompNS_setVectorProp("Outflow_Vel_High", sim_outflowVel(HIGH, :))
#endif

end subroutine Simulation_adjustEvolution
