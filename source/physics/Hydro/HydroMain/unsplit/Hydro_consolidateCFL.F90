!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_consolidateCFL
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!! NAME
!!
!!  Hydro_consolidateCFL
!!
!! SYNOPSIS
!!
!!  call Hydro_consolidateCFL()
!!
!! DESCRIPTION
!!
!!  Consolidate the CFL factors of different MPI tasks by a global
!!  reduction operation.
!!
!!  Different MPI tasks may have computed different values for an updated Hydro
!!  CFL factor. This routine should be called before Hydro_computeDt is called
!!  to compute an updated Hydro time step.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! SIDE EFFECTS
!!
!!  The Hydro_data variable hy_cfl is updated.
!!
!! NOTES
!!
!!  This routine must be called collectively by all MPI tasks in the global
!!  communicator (hy_globalComm).
!!
!!  For split Hydro implementations, there is only a stub implementation of
!!  this interface which does not do any computation or communication.
!!  That is okay since unsplit Hydro implementations are not expected to
!!  modify the CFL factor that is given as a runtime parameter.
!!***

subroutine Hydro_consolidateCFL()
  use Grid_interface, ONLY : Grid_fillGuardCells
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Driver_interface,  ONLY : Driver_getSimTime
  use Hydro_data, ONLY : hy_globalComm, hy_cfl, hy_useVaryingCFL
  use Hydro_data, ONLY : hy_cflStencil, hy_reduceGcellFills
  use Hydro_data, ONLY : hy_simTime, hy_simGeneration, hy_dtminValid
#include "Flash_mpi_implicitNone.fh"
#include "Simulation.h"
#include "constants.h"
  integer :: error
  real :: newCfl
  real :: nowTime
  integer :: nowGeneration
#ifdef CFL_VAR
  logical       :: gcMask(NUNK_VARS)
  logical, save :: gcMaskLogged = .FALSE.
#endif

#ifdef CFL_VAR
  if (hy_cflStencil > 0) then
     if (hy_reduceGcellFills) then
        call Grid_fillGuardCells(CENTER_FACES,ALLDIR,minLayers=hy_cflStencil,&
             unitReadsMeshDataOnly=.true.)
     else
        gcMask(:)       = .FALSE.
        gcMask(CFL_VAR) = .TRUE.
        if (.NOT.gcMaskLogged) then
           call Logfile_stampVarMask(gcMask, .FALSE., '[Hydro_consolidateCFL]', 'gcMask')
        end if
        call Grid_fillGuardCells(CENTER,ALLDIR,minLayers=hy_cflStencil,&
             maskSize=NUNK_VARS,mask=gcMask,&
             doLogMask=.NOT.gcMaskLogged)

        gcMaskLogged = .TRUE.
     end if
  end if
#endif

  if (hy_useVaryingCFL) then
     call MPI_AllReduce (hy_cfl, newCfl, 1, & 
          FLASH_REAL, MPI_MIN, hy_globalComm, error)

     hy_cfl = newCfl
  end if

  ! Check whether the grid has changed since the last time Hydro was called
  call Driver_getSimTime(nowTime, nowGeneration)
  if (nowTime == hy_simTime .AND. nowGeneration > hy_simGeneration) then
     hy_dtminValid = .FALSE.
  end if

end subroutine Hydro_consolidateCFL
