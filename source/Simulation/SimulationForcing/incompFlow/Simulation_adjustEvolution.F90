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
   use sim_outletData, ONLY: sim_QOut, sim_QAux, sim_volOut, sim_volAux, &
                             sim_QOutLiq, sim_QOutGas, sim_QAuxLiq, sim_QAuxGas, &
                             sim_volOutLiq, sim_volOutGas, sim_volAuxLiq, sim_volAuxGas
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
   sim_QOut = 0.
   sim_QOutLiq = 0.
   sim_QOutGas = 0.

   sim_volOut = 0.
   sim_volOutLiq = 0.
   sim_volOutGas = 0.

#ifdef SIMULATION_OUTLET_PHASED

   call MPI_Allreduce(sim_QAuxLiq, sim_QOutLiq, (HIGH - LOW + 1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(sim_QAuxGas, sim_QOutGas, (HIGH - LOW + 1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(sim_volAuxLiq, sim_volOutLiq, (HIGH - LOW + 1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(sim_volAuxGas, sim_volOutGas, (HIGH - LOW + 1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   sim_QOutLiq = sim_QOutLiq/(sim_volOutLiq + 1e-13)
   sim_QOutGas = sim_QOutGas/(sim_volOutGas + 1e-13)

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'Outlet Liq Velocity LOW,', sim_QOutLiq(LOW, :)
      write (*, *) 'Outlet Liq Velocity HIGH,', sim_QOutLiq(HIGH, :)
      write (*, *) '--------------------------------------------------------'
      write (*, *) 'Outlet Gas Velocity LOW,', sim_QOutGas(LOW, :)
      write (*, *) 'Outlet Gas Velocity HIGH,', sim_QOutGas(HIGH, :)
   end if

#else

   call MPI_Allreduce(sim_QAux, sim_QOut, (HIGH - LOW + 1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(sim_volAux, sim_volOut, (HIGH - LOW + 1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   sim_QOut = sim_QOut/(sim_volOut + 1e-13)

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'Outlet Velocity LOW,', sim_QOut(LOW, :)
      write (*, *) 'Outlet Velocity HIGH,', sim_QOut(HIGH, :)
   end if

#endif

#endif

end subroutine Simulation_adjustEvolution
