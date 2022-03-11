!!****if* source/physics/IncompNS/IncompNSMain/IncompNS_init
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
!!  IncompNS_init
!!
!!
!! SYNOPSIS
!!
!!  call IncompNS_init()
!!
!!
!! DESCRIPTION
!!
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initAll.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine IncompNS_init(restart)

   use IncompNS_data
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_getMype, Driver_getNumProcs, &
                               Driver_getComm, Driver_getNstep, Driver_abort
   use ins_interface, ONLY: ins_init

   implicit none

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

   logical, intent(IN) :: restart

   logical :: useMultiphase

   call Driver_getMype(MESH_COMM, ins_meshMe)
   call Driver_getNumProcs(MESH_COMM, ins_meshNumProcs)
   call Driver_getComm(MESH_COMM, ins_meshComm)

   call RuntimeParameters_get("useIncompNS", ins_useIncompNS)
   if (.NOT. ins_useIncompNS) RETURN

   call RuntimeParameters_get("ins_cflFlg", ins_cflflg)
   call RuntimeParameters_get("ins_cfl", ins_cfl)
   call RuntimeParameters_get("ins_isgs", ins_isgs)
   call RuntimeParameters_get("ins_invReynolds", ins_invReynolds)
   call RuntimeParameters_get("ins_sigma", ins_sigma)
   call RuntimeParameters_get("ins_dtSpec", ins_dtspec)
   call RuntimeParameters_get("ins_velProlongMethod", ins_prol_method)

   ins_rhoGas = 1.
   ins_muGas = 1.

   call RuntimeParameters_get("useMultiphase", useMultiphase)
   if (useMultiphase) then
      call RuntimeParameters_get("mph_rhoGas", ins_rhoGas)
      call RuntimeParameters_get("mph_muGas", ins_muGas)
   end if

   if (ins_meshMe .eq. MASTER_PE) then
      write (*, *) 'ins_cfl   =', ins_cfl
      write (*, *) 'ins_isgs  =', ins_isgs
      write (*, *) 'ins_invReynolds =', ins_invReynolds
      write (*, *) 'ins_sigma =', ins_sigma
      write (*, *) 'ins_dtSpec=', ins_dtspec
      write (*, *) 'ins_velProlongMethod=', ins_prol_method
      write (*, *) 'ins_rhoGas=', ins_rhoGas
      write (*, *) 'ins_muGas=', ins_muGas
   end if

   call RuntimeParameters_get("ins_pressureCorrect", ins_prescorr)
   ins_prescoeff = 0.
   if (ins_prescorr) ins_prescoeff = 1.

   ! Read gravity acceleration components:
   call RuntimeParameters_get("ins_gravX", ins_gravX)
   call RuntimeParameters_get("ins_gravY", ins_gravY)
   call RuntimeParameters_get("ins_gravZ", ins_gravZ)

   ! Read pressure gradients if necessary, constant mass simulation data:
   call RuntimeParameters_get("ins_dpdx", ins_dpdx)
   call RuntimeParameters_get("ins_dpdy", ins_dpdy)
   call RuntimeParameters_get("ins_dpdz", ins_dpdz)

   call ins_init()

end subroutine IncompNS_init
