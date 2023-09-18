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
   use ib_interface, ONLY: ib_readBody, ib_annBuildTree
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_getMype, Driver_getNumProcs, &
                               Driver_getComm

   implicit none
   include 'Flashx_mpi.h'
   logical, intent(in) :: restart

   character(len=30) :: bodyFile
   integer :: ibd
   call RuntimeParameters_get("useImBound", ib_useImBound)

   if (.NOT. ib_useImBound) RETURN

   call Driver_getMype(MESH_COMM, ib_meshMe)
   call Driver_getNumProcs(MESH_COMM, ib_meshNumProcs)
   call Driver_getComm(MESH_COMM, ib_meshComm)

   call RuntimeParameters_get("ib_lsIt", ib_lsIt)
   call RuntimeParameters_get("ib_numBodies", ib_numBodies)
   call RuntimeParameters_get("ib_bodyName", ib_bodyName)
   call RuntimeParameters_get("ib_enableSelectiveMapping", ib_enableSelectiveMapping)
   call RuntimeParameters_get("ib_bruteForceMapping", ib_bruteForceMapping)
   call RuntimeParameters_get("ib_annQueries", ib_annQueries)

   if (ib_meshMe .eq. MASTER_PE) then
      write (*, *) 'ib_lsIt=', ib_lsIt
      write (*, *) 'ib_numBodies=', ib_numBodies
      write (*, *) 'ib_bodyName=', ib_bodyName
      write (*, *) 'ib_enableSelectiveMapping', ib_enableSelectiveMapping
      write (*, *) 'ib_bruteForceMapping', ib_bruteForceMapping
      write (*, *) 'ib_annQueries', ib_annQueries
   end if

   allocate (ib_bodyInfo(ib_numBodies))
   allocate (ib_annIdx(ib_annQueries))

   do ibd = 1, ib_numBodies
      write (bodyFile, "(A,A,I4.4)") trim(ib_bodyName), '_hdf5_ibd_', ibd
      call ib_readBody(ib_bodyInfo(ibd), bodyFile)
      call ib_annBuildTree(ib_bodyInfo(ibd))
   end do

end subroutine ImBound_init
