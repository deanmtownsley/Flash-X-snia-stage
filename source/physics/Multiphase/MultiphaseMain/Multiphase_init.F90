!!****if* source/physics/Multiphase/MultiphaseMain/Multiphase_init
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
!!
!!  Multiphase_init
!!
!!
!! SYNOPSIS
!!
!!  call Multiphase_init()
!!
!!
!! DESCRIPTION
!!
!!
!!***

subroutine Multiphase_init(restart)

   use Multiphase_data
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_getMype, Driver_getNumProcs, &
                               Driver_getComm
   use mph_interface, ONLY: mph_init
   use mph_evapInterface, ONLY: mph_evapInit
   use IncompNS_interface, ONLY: IncompNS_getGridVar, IncompNS_getScalarProp
   use HeatAD_interface, ONLY: HeatAD_getGridVar, HeatAD_getScalarProp

   implicit none
   include 'Flashx_mpi.h'

#include "constants.h"
#include "Simulation.h"
#include "Multiphase.h"

   logical, intent(in) :: restart

   logical :: useIncompNS, useHeatAD

   call RuntimeParameters_get("useMultiphase", mph_useMultiphase)

   if (.NOT. mph_useMultiphase) RETURN

   call Driver_getMype(MESH_COMM, mph_meshMe)
   call Driver_getNumProcs(MESH_COMM, mph_meshNumProcs)
   call Driver_getComm(MESH_COMM, mph_meshComm)

   call RuntimeParameters_get("mph_rhoGas", mph_rhoGas)
   call RuntimeParameters_get("mph_muGas", mph_muGas)
   call RuntimeParameters_get("mph_invWeber", mph_invWeber)
   call RuntimeParameters_get("mph_Tsat", mph_Tsat)
   call RuntimeParameters_get("mph_thcoGas", mph_thcoGas)
   call RuntimeParameters_get("mph_CpGas", mph_CpGas)
   call RuntimeParameters_get("mph_Stefan", mph_Stefan)
   call RuntimeParameters_get("mph_lsIt", mph_lsIt)
   call RuntimeParameters_get("mph_extpIt", mph_extpIt)
   call RuntimeParameters_get("mph_iPropSmear", mph_iPropSmear)
   call RuntimeParameters_get("ins_invReynolds", mph_invReynolds)

   mph_Prandtl = 1.

   call RuntimeParameters_get("useHeatAD", useHeatAD)
   if (useHeatAD) call RuntimeParameters_get("ht_Prandtl", mph_Prandtl)

   if (mph_meshMe .eq. MASTER_PE) then
      write (*, *) 'mph_rhoGas=', mph_rhoGas
      write (*, *) 'mph_muGas=', mph_muGas
      write (*, *) 'mph_invWeber=', mph_invWeber
      write (*, *) 'mph_Tsat=', mph_Tsat
      write (*, *) 'mph_thcoGas=', mph_thcoGas
      write (*, *) 'mph_CpGas=', mph_CpGas
      write (*, *) 'mph_iPropSmear=', mph_iPropSmear
      write (*, *) 'mph_Stefan=', mph_Stefan
      write (*, *) 'mph_lsIt=', mph_lsIt
      write (*, *) 'mph_extpIt=', mph_extpIt
      write (*, *) 'mph_invReynolds=', mph_invReynolds
      write (*, *) 'mph_Prandtl=', mph_Prandtl
   end if

   call IncompNS_getGridVar('Face_Velocity', mph_iVelFVar)
   call IncompNS_getGridVar('Face_Density', mph_iRhoFVar)
   call IncompNS_getGridVar('Face_Pressure_Jump', mph_iJumpVar)
   call IncompNS_getGridVar('Center_Viscosity', mph_iMuCVar)
   call IncompNS_getGridVar('Center_Divergence', mph_iDivCVar)

   call mph_evapInit()
   call mph_init()

end subroutine Multiphase_init
