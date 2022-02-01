!!****if* source/physics/ImBound/ImBoundMain/ImBound_init
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
!!  ImBound_init
!!
!!
!! SYNOPSIS
!!
!!  call ImBound_init()
!!  
!!
!! DESCRIPTION
!! 
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine ImBound_init(restart)

  use ImBound_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface,            ONLY : Driver_getMype, Driver_getNumProcs, &
                                          Driver_getComm

  implicit none
  include 'Flashx_mpi.h'
  logical, intent(in) :: restart

  call RuntimeParameters_get("useImBound", ib_useImBound)

  if (.NOT. ib_useImBound) RETURN

  call Driver_getMype(MESH_COMM, ib_meshMe)
  call Driver_getNumProcs(MESH_COMM, ib_meshNumProcs)
  call Driver_getComm(MESH_COMM, ib_meshComm)

  call RuntimeParameters_get("ib_rhoSolid",ib_rhoSolid)
  call RuntimeParameters_get("ib_muSolid",ib_muSolid)
  call RuntimeParameters_get("ib_thcoSolid",ib_thcoSolid)
  call RuntimeParameters_get("ib_CpSolid",ib_CpSolid)
  call RuntimeParameters_get("ib_lsIt",ib_lsIt)

 if (ib_meshMe .eq. MASTER_PE) then
     write(*,*) 'ib_rhoSolid=',ib_rhoSolid
     write(*,*) 'ib_muSolid=',ib_muSolid
     write(*,*) 'ib_thcoSolid=',ib_thcoSolid
     write(*,*) 'ib_CpSolid=',ib_CpSolid
     write(*,*) 'ib_lsIt=',ib_lsIt
  endif

end subroutine ImBound_init
