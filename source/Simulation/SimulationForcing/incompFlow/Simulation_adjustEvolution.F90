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

#ifdef SIMULATION_FORCE_OUTLET
   use sim_outletInterface, ONLY: sim_outletSetForcing
   use sim_outletData, ONLY: sim_QMean, sim_QAux, sim_volMean, sim_volAux, &
                             sim_QMeanLiq, sim_QMeanGas, sim_QAuxLiq, sim_QAuxGas, &
                             sim_volMeanLiq, sim_volMeanGas, sim_volAuxLiq, sim_volAuxGas
#endif

#ifdef SIMULATION_FORCE_INLET
   use sim_inletInterface, ONLY: sim_inletSetForcing
#endif

   implicit none
   include "Flashx_mpi.h"
   integer, intent(in) :: nstep
   real, intent(in) :: dt
   real, intent(in) :: stime

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer :: ierr

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

#ifdef SIMULATION_FORCE_INLET
   ! Set Inlet Forcing
   !-------------------------------------------------------------
   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
      !---------------------------------------------------------
      call sim_inletSetForcing(tileDesc, dt)
      !---------------------------------------------------------
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)
#endif

#ifdef SIMULATION_FORCE_OUTLET
   sim_QAux = 0.
   sim_QAuxLiq = 0.
   sim_QAuxGas = 0.

   sim_volAux = 0.
   sim_volAuxLiq = 0.
   sim_volAuxGas = 0.

   ! Set Outlet Forcing
   !-------------------------------------------------------------
   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
      !---------------------------------------------------------
      call sim_outletSetForcing(tileDesc, dt)
      !---------------------------------------------------------
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)

   ! Consolidate data
   !-------------------------------------------------------------
   sim_QMean = 0.
   sim_QMeanLiq = 0.
   sim_QMeanGas = 0.

   sim_volMean = 0.
   sim_volMeanLiq = 0.
   sim_volMeanGas = 0.

#ifdef SIMULATION_OUTLET_PHASED

   call MPI_Allreduce(sim_QAuxLiq, sim_QMeanLiq, MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(sim_QAuxGas, sim_QMeanGas, MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(sim_volAuxLiq, sim_volMeanLiq, MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(sim_volAuxGas, sim_volMeanGas, MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   sim_QMeanLiq = sim_QMeanLiq/(sim_volMeanLiq + 1e-13)
   sim_QMeanGas = sim_QMeanGas/(sim_volMeanGas + 1e-13)

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'Mean Liq Velocity =', sim_QMeanLiq
      write (*, *) '--------------------------------------------------------'
      write (*, *) 'Mean Gas Velocity =', sim_QMeanGas
   end if

#else

   call MPI_Allreduce(sim_QAux, sim_QMean, MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(sim_volAux, sim_volMean, MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   sim_QMean = sim_QMean/(sim_volMean + 1e-13)

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'Mean Velocity =', sim_QMean
   end if

#endif

#endif

end subroutine Simulation_adjustEvolution
